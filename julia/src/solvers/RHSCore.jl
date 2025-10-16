"""
RHSCore.jl

calRHS!関数のコア処理（初期化 + 6面一様境界条件）を提供

3つのソルバー（DHCP, Sensitivity, Adjoint）で共通する
RHSベクトル計算のコア部分を一元管理する。

共通化範囲:
- wk.bの初期化（ゼロクリア）
- 6面の一様熱流束境界条件（X±, Y±, Z±）

各ソルバー固有処理:
- 分布境界条件（calRHS!内で個別実装）
- 最終RHS計算（内部熱源項の有無）
"""

module RHSCore

using FLoops
using ThreadsX

import ..Commons
using ..Commons: WorkBuffers, get_backend, HEAT_FLUX, CONVECTION, myfill!

import ..BoundaryConditions
using ..BoundaryConditions: BoundaryCondition, BoundaryConditionSet

export calRHS_core!

"""
    calRHS_core!(wk, HF, dx, dy, Δt, dz, par)
      -> (SZ, dx1, dy1, z_st, z_ed, ddt)

RHSベクトルのコア計算（3ソルバー共通部分）

初期化とX/Y/Z方向の一様熱流束境界条件を処理する。
分布境界条件と最終RHS計算は呼び出し側で実装する。

# 引数
- `b, θ
- `dx::Float64`: x方向格子幅 [m]
- `dy::Float64`: y方向格子幅 [m]
- `Δt::Float64`: 時間刻み [s]
- `dz::Vector{Float64}`: z方向格子幅配列 [m]
- `par::String`: 並列化バックエンド（"sequential" or "thread"）

# 処理内容
1. wk.bをゼロクリア（全要素を0.0に初期化）
2. 6面の一様熱流束境界条件を適用:

# 注意
- ガイドセル（境界外）は処理対象外
"""
function calRHS_core!(
            b::AbstractArray{T,3},
            θ::AbstractArray{T,3},
            dx::T,
            dy::T,
            dz::AbstractVector{T},
            bc_set::BoundaryConditionSet,
            par::String)::Nothing where {T <: AbstractFloat}

  # bをゼロクリア
  myfill!(b, zero(T), par)

  # X軸負方向面 (i=1)
  apply_face_bc!(b, θ, dx, dy, dz, bc_set.x_minus, :x_minus, par)

  # X軸正方向面 (i=SZ[1])
  apply_face_bc!(b, θ, dx, dy, dz, bc_set.x_plus, :x_plus, par)

  # Y軸負方向面 (j=1)
  apply_face_bc!(b, θ, dx, dy, dz, bc_set.y_minus, :y_minus, par)

  # Y軸正方向面 (j=SZ[2])
   apply_face_bc!(b, θ, dx, dy, dz, bc_set.y_plus, :y_plus, par)

  # Z軸負方向面 (k=1)
  apply_face_bc!(b, θ, dx, dy, dz, bc_set.z_minus, :z_minus, par)

  # Z軸正方向面 (k=SZ[3])
  apply_face_bc!(b, θ, dx, dy, dz, bc_set.z_plus, :z_plus, par)

  return nothing
end


"""
個別の境界面に境界条件を適用
@param θ 温度
@param λ 熱伝導率
@param ρ 密度
@param cp 比熱
@param mask マスク
@param bc 境界条件
@param face_type 面のタイプ (:x_minus, :x_plus, :y_minus, :y_plus, :z_minus, :z_plus)
"""
function apply_face_bc!(b::AbstractArray{T,3},
                        θ::AbstractArray{T,3},
                        dx::T,
                        dy::T,
                        dz::AbstractVector{T},
                        bc::BoundaryCondition,
                        face_type::Symbol,
                        par::String)::Nothing where {T <: AbstractFloat}

    if bc.type == HEAT_FLUX
        RHS_heat_flux!(b, dx, dy, dz, bc, face_type, par)

    elseif bc.type == CONVECTION
        RHS_convection!(b, θ, dx, dy, dz, bc, face_type, par)
    end

    return nothing
end


"""
熱流束境界条件のRHS項への適用
"""
function RHS_heat_flux!(b::AbstractArray{T,3},
                        dx::T,
                        dy::T,
                        dz::AbstractVector{T},
                        bc::BoundaryCondition,
                        face_type::Symbol,
                        par::String)::Nothing where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(b)
    dx1 = inv(dx)
    dy1 = inv(dy)

    if face_type == :x_minus
      let i=2, a = bc.heat_flux * dx1
        @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
          b[i,j,k] += a
        end
      end
    elseif face_type == :x_plus
      let i = SZ[1]-1, a = bc.heat_flux * dx1
        @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
          b[i,j,k] -= a
        end
      end
    elseif face_type == :y_minus
      let j = 2, a = bc.heat_flux * dy1
        @floop backend for k in 2:SZ[3]-1, i in 2:SZ[1]-1
          b[i,j,k] += a
        end
      end
    elseif face_type == :y_plus
      let j = SZ[2]-1, a = bc.heat_flux * dy1
        @floop backend for k in 2:SZ[3]-1, i in 2:SZ[1]-1
          b[i,j,k] -= a
        end
      end
    elseif face_type == :z_minus
      let k = 2,  a = bc.heat_flux / dz[k]
        @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
          b[i,j,k] += a
        end
      end
    elseif face_type == :z_plus
      let k = SZ[3]-1, a = bc.heat_flux / dz[k]
        @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
          b[i,j,k] -= a
        end
      end
    end

    return nothing
end


"""
熱伝達境界条件のRHS項への適用
"""
function RHS_convection!(b::AbstractArray{T,3},
                        θ::AbstractArray{T,3},
                        dx::T,
                        dy::T,
                        dz::AbstractVector{T},
                        bc::BoundaryCondition,
                        face_type::Symbol,
                        par::String)::Nothing where {T <: AbstractFloat}
  backend = get_backend(par)
  SZ = size(b)
  dx1 = inv(dx)
  dy1 = inv(dy)

  if face_type == :x_minus
    let i=2, a = bc.heat_transfer_coefficient * dx1, ta = bc.ambient_temperature
      @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        b[i,j,k] -= a * (ta - θ[i,j,k])
      end
    end
  elseif face_type == :x_plus
    let i = SZ[1]-1, a = bc.heat_transfer_coefficient * dx1, ta = bc.ambient_temperature
      @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        b[i,j,k] += a * (ta - θ[i,j,k])
      end
    end
  elseif face_type == :y_minus
    let j = 2, a = bc.heat_transfer_coefficient * dy1, ta = bc.ambient_temperature
      @floop backend for k in 2:SZ[3]-1, i in 2:SZ[1]-1
        b[i,j,k] -= a * (ta - θ[i,j,k])
      end
    end
  elseif face_type == :y_plus
    let j = SZ[2]-1, a = bc.heat_transfer_coefficient * dy1, ta = bc.ambient_temperature
      @floop backend for k in 2:SZ[3]-1, i in 2:SZ[1]-1
        b[i,j,k] += a * (ta - θ[i,j,k])
      end
    end
  elseif face_type == :z_minus
    let k = 2, a = bc.heat_transfer_coefficient / dz[k], ta = bc.ambient_temperature
      @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        b[i,j,k] -= a * (ta - θ[i,j,k])
      end
    end
  elseif face_type == :z_plus
    let k = SZ[3]-1, a = bc.heat_transfer_coefficient / dz[k], ta = bc.ambient_temperature
      @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
        b[i,j,k] += a * (ta - θ[i,j,k])
      end
    end
  end

  return nothing
end




end # module RHSCore
