"""
SolverDiagnostics.jl

ソルバー診断情報を格納する構造体

Python-Julia比較検証用に、各ソルバーの詳細情報を記録する。
"""

module SolverDiagnostics

export SolverDiag, TimeStepDiag, CGMIterDiag
export add_timestep!, summary

"""
1時間ステップの診断情報
"""
mutable struct TimeStepDiag
  timestep::Int              # 時間ステップ番号
  iterations::Int            # 反復回数
  initial_residual::Float64  # 初期残差ノルム
  final_residual::Float64    # 最終残差ノルム
  converged::Bool            # 収束フラグ
end

"""
ソルバー診断情報（複数時間ステップ分）
"""
mutable struct SolverDiag
  solver_name::String                    # ソルバー名（"DHCP", "Adjoint", "Sensitivity"）
  timesteps::Vector{TimeStepDiag}        # 各時間ステップの診断情報
  elapsed_time::Float64                   # 総実行時間

  function SolverDiag(name::String)
    new(name, TimeStepDiag[], 0.0)
  end
end

"""
1時間ステップの診断情報を追加
"""
function add_timestep!(diag::SolverDiag, timestep::Int, iters::Int,
                       init_resid::Float64, final_resid::Float64, converged::Bool)
  push!(diag.timesteps, TimeStepDiag(timestep, iters, init_resid, final_resid, converged))
end

"""
診断情報のサマリーを返す

Returns:
  Dict with keys:
    - total_iters: 総反復回数
    - avg_iters: 平均反復回数
    - max_iters: 最大反復回数
    - min_iters: 最小反復回数
    - elapsed_time: 総実行時間
    - convergence_rate: 収束率（収束した時間ステップの割合）
    - timesteps: 時間ステップ数
"""
function summary(diag::SolverDiag)
  if isempty(diag.timesteps)
    return Dict(
      "total_iters" => 0,
      "avg_iters" => 0.0,
      "max_iters" => 0,
      "min_iters" => 0,
      "elapsed_time" => diag.elapsed_time,
      "convergence_rate" => 0.0,
      "timesteps" => 0
    )
  end

  iters_list = [ts.iterations for ts in diag.timesteps]
  total_iters = sum(iters_list)
  avg_iters = total_iters / length(iters_list)
  max_iters = maximum(iters_list)
  min_iters = minimum(iters_list)

  converged_count = count(ts -> ts.converged, diag.timesteps)
  convergence_rate = converged_count / length(diag.timesteps)

  return Dict(
    "total_iters" => total_iters,
    "avg_iters" => avg_iters,
    "max_iters" => max_iters,
    "min_iters" => min_iters,
    "elapsed_time" => diag.elapsed_time,
    "convergence_rate" => convergence_rate,
    "timesteps" => length(diag.timesteps)
  )
end

"""
CGM1反復の診断情報
"""
mutable struct CGMIterDiag
  iteration::Int                         # CGM反復番号
  J::Float64                             # 目的関数値
  beta::Float64                          # ステップサイズ
  gamma::Float64                         # 共役勾配係数
  wall_time::Float64                     # 壁時間
  dhcp_summary::Union{Dict,Nothing}      # DHCP診断サマリー
  adjoint_summary::Union{Dict,Nothing}   # Adjoint診断サマリー
  sensitivity_summary::Union{Dict,Nothing}  # Sensitivity診断サマリー

  function CGMIterDiag(iter::Int)
    new(iter, 0.0, 0.0, 0.0, 0.0, nothing, nothing, nothing)
  end
end

end  # module SolverDiagnostics
