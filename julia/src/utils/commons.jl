"""
commons.jl

共通ユーティリティ

"""
module Commons

using FLoops

export initialize_cells!, λf, WorkBuffers, get_backend,
       BoundaryType, ISOTHERMAL, HEAT_FLUX, CONVECTION, FloatMin, reset_work_buffers!

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
 WorkBuffersのリセット（CGM反復間でクリーンな状態を保証）
"""
function reset_work_buffers!(wk::WorkBuffers)
  # 反復ソルバー用配列をゼロクリア
  wk.pcg_p  .= 0.0
  wk.pcg_p_ .= 0.0
  wk.pcg_r  .= 0.0
  wk.pcg_r0 .= 0.0
  wk.pcg_q  .= 0.0
  wk.pcg_s  .= 0.0
  wk.pcg_s_ .= 0.0
  wk.pcg_t_ .= 0.0

  # θとbもリセット（solve内で上書きされるが、念のため）
  wk.θ    .= 0.0
  wk.b    .= 0.0
  wk.hsrc .= 0.0

  # maskは1.0に保持（境界条件で上書きされる）
  wk.mask .= 1.0

  # cp、λは物性値計算で毎回設定されるのでリセット不要
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


end # module GridTransform
