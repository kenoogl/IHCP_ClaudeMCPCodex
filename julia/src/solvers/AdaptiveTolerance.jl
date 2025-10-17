"""
AdaptiveTolerance.jl

適応的収束判定モジュール

時間発展ソルバー（DHCP, Adjoint, Sensitivity）向けの適応的収束基準。
解の相対変化量に基づいてPBICGSTAB!のtolを動的に調整する。

理論的根拠:
- Eisenstat-Walker inexact Newton (1996)
- 時間離散化誤差との整合 (Knoll & Keyes 2004)
- 温度場変化量と収束性 (Elman et al. 2014)

参考文献: docs/adaptive_tolerance.md
"""

module AdaptiveTolerance

using LinearAlgebra

export AdaptiveToleranceParams, compute_adaptive_tol, compute_relative_change

"""
適応的収束判定のパラメータ

# フィールド
- `tol_min::T`: 最小許容誤差（精度保証の下限）
- `tol_max::T`: 最大許容誤差（緩和の上限）
- `κ_θ::T`: 相対変化量のスケール係数（推奨: 0.2）
- `warmup_steps::Int`: 固定tolを使用する初期ステップ数（推奨: 3）
- `ε::T`: ゼロ除算防止用の小さい値（推奨: 1e-12）

# 推奨パラメータ
- DHCP/Sensitivity: tol_min=1e-7, tol_max=1e-5
- Adjoint: tol_min=1e-9, tol_max=1e-5
"""
struct AdaptiveToleranceParams{T <: AbstractFloat}
  tol_min::T
  tol_max::T
  κ_θ::T
  warmup_steps::Int
  ε::T

  function AdaptiveToleranceParams{T}(;
    tol_min::T = T(1e-7),
    tol_max::T = T(1e-5),
    κ_θ::T = T(0.2),
    warmup_steps::Int = 3,
    ε::T = T(1e-12)
  ) where {T <: AbstractFloat}
    new{T}(tol_min, tol_max, κ_θ, warmup_steps, ε)
  end
end

# Float64コンストラクタ（デフォルト）
AdaptiveToleranceParams(; kwargs...) = AdaptiveToleranceParams{Float64}(; kwargs...)

"""
    compute_adaptive_tol(θ_current, θ_previous, t, tol_default, params)

適応的収束基準の計算（解の相対変化ベース）

# アルゴリズム
1. 初期ステップ（t ≤ warmup_steps）: 固定tol使用
2. 相対変化量の計算: δ = ||θⁿ - θⁿ⁻¹||₂ / (||θⁿ||₂ + ε)
3. 適応的tol: tol_adaptive = κ_θ · δ
4. クランプ: clamp(tol_adaptive, tol_min, tol_max)

# 引数
- `θ_current::AbstractArray{T,3}`: 現在の温度場（内点のみ、ni×nj×nk）
- `θ_previous::AbstractArray{T,3}`: 前ステップの温度場（内点のみ、ni×nj×nk）
- `t::Int`: 現在の時間ステップ（1始まり）
- `tol_default::T`: デフォルト収束基準（warmup期間中に使用）
- `params::AdaptiveToleranceParams{T}`: パラメータ

# 戻り値
- `tol::T`: 適応的収束基準

# 使用例
```julia
params = AdaptiveToleranceParams{Float64}(tol_min=1e-7, tol_max=1e-5, κ_θ=0.2)
tol = compute_adaptive_tol(T_all[:,:,:,t], T_all[:,:,:,t-1], t, 1e-6, params)
```
"""
function compute_adaptive_tol(
  θ_current::AbstractArray{T,3},
  θ_previous::AbstractArray{T,3},
  t::Int,
  tol_default::T,
  params::AdaptiveToleranceParams{T}
) where {T <: AbstractFloat}
  # 初期ステップは固定tol使用（移動平均が安定しない期間）
  if t <= params.warmup_steps
    return tol_default
  end

  # 相対変化量の計算
  δ = compute_relative_change(θ_current, θ_previous, params.ε)

  # 適応的tol計算（δが小さい→収束基準を緩和）
  tol_adaptive = params.κ_θ * δ

  # クランプ（精度保証と過剰緩和防止）
  return clamp(tol_adaptive, params.tol_min, params.tol_max)
end

"""
    compute_relative_change(θ_current, θ_previous, ε)

解の相対変化量を計算

δ = ||θⁿ - θⁿ⁻¹||₂ / (||θⁿ||₂ + ε)

# 引数
- `θ_current::AbstractArray{T,3}`: 現在の温度場（ni×nj×nk）
- `θ_previous::AbstractArray{T,3}`: 前ステップの温度場（ni×nj×nk）
- `ε::T`: ゼロ除算防止用の小さい値

# 戻り値
- `δ::T`: 相対変化量（0から1の範囲、通常は0.001〜0.1程度）

# 注意
- L2ノルムを使用（全格子点の二乗和の平方根）
- εによるゼロ除算防止（定常状態でも安定動作）
"""
function compute_relative_change(
  θ_current::AbstractArray{T,3},
  θ_previous::AbstractArray{T,3},
  ε::T
) where {T <: AbstractFloat}
  # ノルムの計算（LinearAlgebra.normは高速だが、ここでは明示的に記述）
  diff_norm = zero(T)
  current_norm = zero(T)

  @inbounds for i in eachindex(θ_current, θ_previous)
    diff = θ_current[i] - θ_previous[i]
    diff_norm += diff * diff
    current_norm += θ_current[i] * θ_current[i]
  end

  diff_norm = sqrt(diff_norm)
  current_norm = sqrt(current_norm)

  # 相対変化量（ε でゼロ除算防止）
  δ = diff_norm / (current_norm + ε)

  return δ
end

end # module AdaptiveTolerance
