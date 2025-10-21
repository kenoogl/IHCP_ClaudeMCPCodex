#!/usr/bin/env julia
"""
ソルバー・前処理組み合わせ比較スクリプト

summary.md（shared/results/solver_comparison/summary.md）の条件で
全ソルバー・前処理の組み合わせをテストします。

テスト条件:
- 格子サイズ: 80×100×20
- 時間ステップ: 10
- CGM反復: 1回

ソルバー:
- PCG (共役勾配法)
- PBICGSTAB (前処理付きBiCGSTAB!)

前処理:
- none (前処理なし)
- diagonal (対角前処理、Python版互換)
- gs (Gauss-Seidel前処理)

実行方法:
  julia --project=julia julia/benchmarks/solver_preconditioner_comparison.jl
"""

using Printf
using Statistics
using Dates
using IHCP_CGM

# テスト条件（summary.mdと同一）
const NI = 80
const NJ = 100
const NK = 20
const NT = 10
const CGM_ITERATIONS = 1

# ソルバー・前処理の組み合わせ
const SOLVERS = [:PCG, :PBICGSTAB]
const PRECONDITIONERS = [:none, :diagonal, :gs]

"""
単一の組み合わせをテスト

Args:
  solver: :PCG または :PBICGSTAB
  precond: :none, :diagonal, :gs

Returns:
  結果辞書（cgm_time, dhcp_iters, adjoint_iters, sensitivity_iters）
"""
function test_single_combination(solver::Symbol, precond::Symbol)
  println("\n" * "="^70)
  println("テスト: $solver + $precond")
  println("="^70)

  # パラメータ設定
  dx = 0.12e-3
  dy = 0.12e-3
  dt = 1e-3
  rho = 7900.0

  # z方向格子（非均等）
  dz = vcat(
    fill(5.0e-5, 5),
    fill(1.0e-4, 5),
    fill(2.0e-4, NK - 10)
  )

  # 境界条件係数（HC配列）
  hc_array = fill(0.0, 6)  # [hc_xm, hc_xp, hc_ym, hc_yp, hc_zm, hc_zp]

  # 熱物性値係数（SUS304）
  cp_coeffs = [420.0, 0.15, 0.0, 0.0]
  k_coeffs = [14.0, 0.01, 0.0, 0.0]

  # 初期温度場（300K）
  T_init = fill(300.0, NI, NJ, NK)

  # 測定温度データ（ダミー、300K）
  T_measure = fill(300.0, NT, NI, NJ)

  # 初期熱流束推定値（0.0）
  q_init = zeros(NT - 1, NI, NJ)

  # CGMソルバーオプション
  cgm_options = Dict(
    "max_iterations" => CGM_ITERATIONS,
    "tolerance" => 1e-6,
    "solver" => string(solver),
    "smoother" => precond,
    "parallel" => "thread"
  )

  try
    # CGM実行（時間計測）
    start_time = time()
    q_estimated, _, convergence_info = solve_sliding_window_cgm(
      T_measure,
      T_init,
      q_init,
      cp_coeffs,
      k_coeffs,
      rho,
      dx,
      dy,
      dz,
      dt,
      hc_array;
      cgm_options=cgm_options,
      window_size=NT,
      overlap_size=0,
      verbose=false
    )
    cgm_time = time() - start_time

    # 反復回数の集計
    n_windows = length(convergence_info)
    dhcp_iters = []
    adjoint_iters = []
    sensitivity_iters = []

    for window_info in convergence_info
      for iter_info in window_info
        push!(dhcp_iters, iter_info["dhcp_iterations"])
        push!(adjoint_iters, iter_info["adjoint_iterations"])
        push!(sensitivity_iters, iter_info["sensitivity_iterations"])
      end
    end

    # 平均反復回数
    avg_dhcp = mean(dhcp_iters)
    avg_adjoint = mean(adjoint_iters)
    avg_sensitivity = mean(sensitivity_iters)

    println("\n結果:")
    println("  CGM実行時間: $(round(cgm_time, digits=2))秒")
    println("  DHCP平均反復: $(round(avg_dhcp, digits=1))")
    println("  Adjoint平均反復: $(round(avg_adjoint, digits=1))")
    println("  Sensitivity平均反復: $(round(avg_sensitivity, digits=1))")

    return Dict(
      "success" => true,
      "cgm_time" => cgm_time,
      "avg_dhcp_iters" => avg_dhcp,
      "avg_adjoint_iters" => avg_adjoint,
      "avg_sensitivity_iters" => avg_sensitivity
    )

  catch e
    println("\n✗ エラー発生: $e")
    return Dict(
      "success" => false,
      "error" => string(e)
    )
  end
end

"""
全組み合わせをテストして結果を表示
"""
function run_all_tests()
  println("\n" * "="^70)
  println("ソルバー・前処理組み合わせ比較テスト")
  println("="^70)
  println("条件:")
  println("  - 格子サイズ: $(NI)×$(NJ)×$(NK)")
  println("  - 時間ステップ: $(NT)")
  println("  - CGM反復: $(CGM_ITERATIONS)回")
  println()

  results = Dict()

  for solver in SOLVERS
    for precond in PRECONDITIONERS
      key = "$(solver)_$(precond)"
      results[key] = test_single_combination(solver, precond)
    end
  end

  # 結果サマリー表示
  println("\n" * "="^70)
  println("結果サマリー")
  println("="^70)
  println()
  println("| ソルバー | 前処理 | CGM時間 | DHCP反復 | Adjoint反復 | Sensitivity反復 |")
  println("|---------|--------|---------|----------|-------------|----------------|")

  for solver in SOLVERS
    for precond in PRECONDITIONERS
      key = "$(solver)_$(precond)"
      res = results[key]

      if res["success"]
        cgm_time_str = @sprintf("%.2f秒", res["cgm_time"])
        dhcp_str = @sprintf("%.1f", res["avg_dhcp_iters"])
        adjoint_str = @sprintf("%.1f", res["avg_adjoint_iters"])
        sensitivity_str = @sprintf("%.1f", res["avg_sensitivity_iters"])
        println("| $solver | $precond | $cgm_time_str | $dhcp_str | $adjoint_str | $sensitivity_str |")
      else
        println("| $solver | $precond | エラー | - | - | - |")
      end
    end
  end

  println()
  return results
end

# メイン実行
if abspath(PROGRAM_FILE) == @__FILE__
  results = run_all_tests()
  println("\n✓ テスト完了")
end
