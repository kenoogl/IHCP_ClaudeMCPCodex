"""
StoppingCriteria.jl

Phase 4: CGM最適化の停止判定モジュール

CGM反復の停止条件を判定する機能を提供：
1. Discrepancy停止判定（成功収束）
2. プラトー検出（計算ボトルネック）
3. 反復数制限

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- global_CGM_time() (1436-1454行)

主要関数:
- check_discrepancy: Discrepancy条件満足判定
- check_plateau: プラトー（目的関数の停滞）検出
- check_stopping_criteria: 統合停止判定
"""

module StoppingCriteria

using Printf

export check_discrepancy, check_plateau, check_stopping_criteria, StoppingStatus

"""
停止判定の結果を格納する構造体

Fields:
  should_stop: 停止すべきか (true/false)
  reason: 停止理由（文字列）
  criterion: 判定基準名（:discrepancy, :plateau, :max_iter, :none）
"""
struct StoppingStatus
  should_stop::Bool
  reason::String
  criterion::Symbol
end


"""
    check_discrepancy(J, residual, epsilon, sigma) -> Bool

Discrepancy停止判定

判定条件:
  J < ε  かつ  max|residual| ≤ σ

Args:
  J: 目的関数値 (残差二乗和)
  residual: 残差場 (nt-1, ni, nj) [K]
  epsilon: Discrepancy基準値 (通常 M*σ²*(nt-1))
  sigma: 測定誤差標準偏差 [K]

Returns:
  satisfied: 条件満足の可否
"""
function check_discrepancy(J::Float64, residual::Array{Float64,3}, epsilon::Float64, sigma::Float64)
  max_abs_residual = maximum(abs.(residual))
  satisfied = (J < epsilon) && (max_abs_residual <= sigma)
  return satisfied
end


"""
    check_plateau(J_hist, P, eta, eps) -> (Bool, Float64)

プラトー（停滞）検出

近P回の平均相対下降率が閾値η未満の場合にプラトー判定

判定条件:
  mean(rel_drops) < η
  ここで rel_drops[i] = max(0, (J[i-1] - J[i]) / |J[i-1]|)

Args:
  J_hist: 目的関数履歴 (vector)
  P: 平均化ウィンドウサイズ（通常10）
  eta: 相対下降閾値（通常1e-4）
  eps: ゼロ除算防止の微小値（通常1e-12）

Returns:
  is_plateau: プラトー検出の可否
  rel_drop_avg: 平均相対下降率（計算できない場合は NaN）
"""
function check_plateau(J_hist::Vector{Float64}, P::Int, eta::Float64, eps::Float64=1e-12)
  n_hist = length(J_hist)

  # 履歴が不十分
  if n_hist < P + 1
    return false, NaN
  end

  # 近P回の相対下降率計算
  drops = Float64[]
  for i in (n_hist - P + 1):n_hist
    prev_J = J_hist[i - 1]
    curr_J = J_hist[i]
    rel_drop = max(0.0, (prev_J - curr_J) / (abs(prev_J) + eps))
    push!(drops, rel_drop)
  end

  rel_drop_avg = sum(drops) / P
  is_plateau = (rel_drop_avg < eta)

  return is_plateau, rel_drop_avg
end


"""
    check_stopping_criteria(J, J_hist, residual, it, params) -> StoppingStatus

CGM統合停止判定

優先順位:
  1. Discrepancy条件満足（成功収束）
  2. プラトー検出（計算ボトルネック）
  3. 最大反復数到達

Args:
  J: 現在の目的関数値
  J_hist: 目的関数履歴 (vector)
  residual: 現在の残差場 (nt-1, ni, nj) [K]
  it: 現在の反復数（0-indexed）
  params: 停止判定パラメータ（NamedTuple）
    - epsilon: Discrepancy基準値
    - sigma: 測定誤差標準偏差 [K]
    - min_iter: 最小反復数（早期停止防止）
    - P: プラトー検出ウィンドウサイズ
    - eta: プラトー相対下降閾値
    - max_iter: 最大反復数
    - eps: ゼロ除算防止値（オプション、デフォルト1e-12）

Returns:
  status: StoppingStatus構造体
"""
function check_stopping_criteria(
  J::Float64,
  J_hist::Vector{Float64},
  residual::Array{Float64,3},
  it::Int,
  params::NamedTuple
)
  epsilon = params.epsilon
  sigma = params.sigma
  min_iter = params.min_iter
  P = params.P
  eta = params.eta
  max_iter = params.max_iter
  eps = haskey(params, :eps) ? params.eps : 1e-12

  # 1. Discrepancy判定（最小反復数以上）
  if it >= min_iter
    if check_discrepancy(J, residual, epsilon, sigma)
      max_residual = maximum(abs.(residual))
      reason = @sprintf("[停止] Discrepancy条件満足: J=%.4e < %.4e かつ max|ΔT|=%.3e ≤ σ=%.1f",
                        J, epsilon, max_residual, sigma)
      return StoppingStatus(true, reason, :discrepancy)
    end
  end

  # 2. プラトー判定（最小反復数以上）
  if it >= min_iter
    is_plateau, rel_drop_avg = check_plateau(J_hist, P, eta, eps)
    if is_plateau
      reason = @sprintf("[停止] プラトー検出: rel_drop_avg=%.3e < eta=%.1e（直近%dステップの平均進展が微小）",
                        rel_drop_avg, eta, P)
      return StoppingStatus(true, reason, :plateau)
    end
  end

  # 3. 最大反復数判定（更新完了後にチェック）
  if it + 1 >= max_iter
    reason = @sprintf("[停止] 最大反復数到達: iter=%d >= max_iter=%d", it + 1, max_iter)
    return StoppingStatus(true, reason, :max_iter)
  end

  # 継続
  return StoppingStatus(false, "[継続] 停止条件未満", :none)
end

end # module StoppingCriteria
