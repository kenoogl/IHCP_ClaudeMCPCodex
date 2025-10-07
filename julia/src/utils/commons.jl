"""
commons.jl

共通ユーティリティ

"""
module Commons

export initialize_cells!, compute_z_range, λf, WorkBuffers, get_backend,
       BoundaryType, ISOTHERMAL, HEAT_FLUX, CONVECTION, FloatMin

# ゼロ除算回避用の小さな数値（BiCGSTAB法などで使用）
const FloatMin = 1.0e-37
      
"""
Harmonic mean with mask correction (Heat3ds互換)
@param a left value
@param b right value
@param ma mask for left (1.0=interior, 0.0=boundary)
@param mb mask for right (1.0=interior, 0.0=boundary)
"""
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))


# 境界条件タイプの列挙型
@enum BoundaryType begin
    ISOTHERMAL   # 等温条件 (Dirichlet)
    HEAT_FLUX    # 熱流束条件 (Neumann)
    CONVECTION   # 熱伝達条件 (Robin)
end


"""
並列動作バックエンドを返す
"""
function get_backend(par::String)
    return (par == "thread") ? ThreadedEx() : SequentialEx()
end


"""
 Heat3D 配列
"""
struct WorkBuffers
    θ      ::Array{Float64,3}
    b      ::Array{Float64,3}
    mask   ::Array{Float64,3}
    cp     ::Array{Float64,3}
    λ      ::Array{Float64,3}
    α      ::Array{Float64,3}
    pcg_p  ::Array{Float64,3}
    pcg_p_ ::Array{Float64,3}
    pcg_r  ::Array{Float64,3}
    pcg_r0 ::Array{Float64,3}
    pcg_q  ::Array{Float64,3}
    pcg_s  ::Array{Float64,3}
    pcg_s_ ::Array{Float64,3}
    pcg_t_ ::Array{Float64,3}
    hsrc   ::Array{Float64,3}
end

"""
 Heat3D 配列確保
"""
function WorkBuffers(mx::Int64, my::Int64, mz::Int64)
  WorkBuffers(
    zeros(Float64, mx, my, mz), # θ
    zeros(Float64, mx, my, mz), # b
     ones(Float64, mx, my, mz), # mask
    zeros(Float64, mx, my, mz), # cp
    zeros(Float64, mx, my, mz), # λ
    zeros(Float64, mx, my, mz), # α
    zeros(Float64, mx, my, mz), # pcg_p
    zeros(Float64, mx, my, mz), # pcg_p_
    zeros(Float64, mx, my, mz), # pcg_r
    zeros(Float64, mx, my, mz), # pcg_r0
    zeros(Float64, mx, my, mz), # pcg_q
    zeros(Float64, mx, my, mz), # pcg_s
    zeros(Float64, mx, my, mz), # pcg_s_
    zeros(Float64, mx, my, mz), # pcg_t_ 
    zeros(Float64, mx, my, mz)  # hsrc
  )
end


"""
  initialize_guard_cells!(θ, λ, mask, θ_init, λ_init)

ガイドセルを含む配列を初期化（内点のデータをコピー）

# Arguments
- θ: 温度配列 (ni+2, nj+2, nk+2)
- λ: 熱伝導率配列 (ni+2, nj+2, nk+2)
- mask: マスク配列 (ni+2, nj+2, nk+2)
- θ_init: 初期温度（内点） (ni, nj, nk)
- λ_init: 初期熱伝導率（内点） (ni, nj, nk)

# Side Effects
- θ, λ, mask配列を初期化
- 全領域のマスクを1.0（内点）に設定
- 内点データをガイドセル配列の内部にコピー
- ガイドセル部分の値は未初期化（境界条件設定で上書き）

# Notes
- 境界条件の詳細設定は BoundaryConditions.apply_boundary_conditions! を使用
- この関数は配列の基本的な初期化のみを行う
"""
function initialize_cells!(
  θ::Array{Float64,3},
  λ::Array{Float64,3},
  mask::Array{Float64,3},
  θ_init::Array{Float64,3},
  λ_init::Array{Float64,3}
)
  SZ = size(θ)
  @assert size(λ) == SZ "λ size mismatch"
  @assert size(mask) == SZ "mask size mismatch"

  ni, nj, nk = size(θ_init)
  @assert SZ == (ni+2, nj+2, nk+2) "Array size mismatch: expected ($(ni+2), $(nj+2), $(nk+2)), got $SZ"

  # 全領域を内点として初期化
  mask .= 1.0

  # 内点の温度と熱伝導率を設定（インデックス2:ni+1, 2:nj+1, 2:nk+1）
  for k in 1:nk, j in 1:nj, i in 1:ni
    θ[i+1, j+1, k+1] = θ_init[i, j, k]
    λ[i+1, j+1, k+1] = λ_init[i, j, k]
  end

  # ガイドセル部分（境界）は未初期化
  # BoundaryConditions.apply_boundary_conditions! で設定

  return nothing
end


"""
  compute_z_range(nk, z_minus_bc, z_plus_bc) -> Vector{Int}

境界条件に基づいてZ方向の計算範囲を決定

# Arguments
- nk: 計算内点数（配列サイズはnk+2）
- z_minus_bc: Z-方向（下側）の境界条件タイプ
- z_plus_bc: Z+方向（上側）の境界条件タイプ

# Returns
- [z_start, z_end]: 計算範囲のインデックス（ガイドセル配列基準）

# 配列インデックス構造 (MZ = nk+2)
- k=1: 下側ガイドセル
- k=2: 下側境界セル
- k=3～nk+1: 計算内点
- k=nk+2: 上側ガイドセル

# 計算範囲の決定ルール
- Z-方向が等温条件 → z_start = 3 (k=2は温度固定、計算対象外)
- Z-方向が熱流束/熱伝達 → z_start = 2 (k=2も計算対象)
- Z+方向が等温条件 → z_end = nk (k=nk+1は温度固定、計算対象外)
- Z+方向が熱流束/熱伝達 → z_end = nk+1 (k=nk+1も計算対象)

# Examples
```julia
compute_z_range(10, ISOTHERMAL, HEAT_FLUX)  # [3, 11]
compute_z_range(10, HEAT_FLUX, ISOTHERMAL)  # [2, 10]
compute_z_range(10, HEAT_FLUX, CONVECTION)  # [2, 11]
compute_z_range(10, ISOTHERMAL, ISOTHERMAL) # [3, 10]
```
"""
function compute_z_range(nk::Int, z_minus_bc::BoundaryType, z_plus_bc::BoundaryType)
  # Z-方向: 等温条件なら内点から開始（k=3）
  z_start = (z_minus_bc == ISOTHERMAL) ? 3 : 2

  # Z+方向: 等温条件なら境界セルを除外（k=nkまで）
  z_end = (z_plus_bc == ISOTHERMAL) ? nk : nk+1

  return [z_start, z_end]
end


end # module GridTransform
