"""
適応的収束判定の動作確認テスト

小規模格子（10×10×5、5ステップ）で適応的tolの動作を確認：
1. 固定tol vs 適応的tolの反復回数比較
2. verbose出力でtol値の変化を確認
3. 温度場の精度を確認
"""

using Pkg
Pkg.activate("julia")

using IHCP_CGM
using Printf
using Statistics

println("="^70)
println("適応的収束判定の動作確認テスト")
println("="^70)
println()

# 小規模格子設定
ni, nj, nk = 10, 10, 5
nt = 6  # 初期+5ステップ
dx = 0.00012
dy = 0.00012
dt = 0.001

# Z方向格子
Z = Float64[0.0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0007]
ΔZ = diff(Z)

# 熱物性値
ρ = 7900.0
cp_coeffs = [439.0, 0.18483, 0.0, 0.0]
k_coeffs = [9.248, 0.01571, 0.0, 0.0]

# 初期温度場（一定温度）
T_initial = fill(300.0, ni, nj, nk)

# 熱流束（一定値）
q_surface = fill(1e6, ni, nj, nt-1)

# ワークバッファ
wk = WorkBuffers(ni+2, nj+2, nk+2)

println("格子設定:")
println("  ni×nj×nk = $(ni)×$(nj)×$(nk)")
println("  時間ステップ数: $(nt)")
println("  dx = $(dx) m, dy = $(dy) m, dt = $(dt) s")
println()

# テスト1: 固定tol（デフォルト）
println("="^70)
println("テスト1: 固定tol（デフォルト、adaptive_tol=false）")
println("="^70)
println()

T_buffer_fixed = zeros(Float64, ni, nj, nk, nt)
iter_buffer_fixed = zeros(Int, nt)

t_start = time()
T_fixed, iters_fixed = solve_dhcp!(
  T_initial,
  q_surface,
  wk,
  nt,
  ρ,
  cp_coeffs,
  k_coeffs,
  dx,
  dy,
  Z,
  ΔZ,
  dt;
  rtol=1e-6,
  maxiter=20000,
  verbose=true,
  par="sequential",
  T_buffer=T_buffer_fixed,
  iter_buffer=iter_buffer_fixed,
  use_previous_solution=true,
  extrapolation=:none,
  smoother=:gs,
  adaptive_tol=false
)
t_fixed = time() - t_start

println()
println("固定tol結果:")
println("  実行時間: $(round(t_fixed, digits=3)) s")
println("  総反復回数: $(sum(iters_fixed))")
println("  平均反復回数: $(round(mean(iters_fixed), digits=1))")
println("  反復回数詳細: $(iters_fixed)")
println("  最終温度範囲: [$(round(minimum(T_fixed), digits=2)), $(round(maximum(T_fixed), digits=2))] K")
println()

# テスト2: 適応的tol
println("="^70)
println("テスト2: 適応的tol（adaptive_tol=true）")
println("="^70)
println()

# 適応的tolパラメータ設定
adaptive_params = AdaptiveToleranceParams{Float64}(
  tol_min=1e-7,
  tol_max=1e-5,
  κ_θ=0.2,
  warmup_steps=3,
  ε=1e-12
)

println("適応的tolパラメータ:")
println("  tol_min: $(adaptive_params.tol_min)")
println("  tol_max: $(adaptive_params.tol_max)")
println("  κ_θ: $(adaptive_params.κ_θ)")
println("  warmup_steps: $(adaptive_params.warmup_steps)")
println()

T_buffer_adaptive = zeros(Float64, ni, nj, nk, nt)
iter_buffer_adaptive = zeros(Int, nt)

t_start = time()
T_adaptive, iters_adaptive = solve_dhcp!(
  T_initial,
  q_surface,
  wk,
  nt,
  ρ,
  cp_coeffs,
  k_coeffs,
  dx,
  dy,
  Z,
  ΔZ,
  dt;
  rtol=1e-6,
  maxiter=20000,
  verbose=true,
  par="sequential",
  T_buffer=T_buffer_adaptive,
  iter_buffer=iter_buffer_adaptive,
  use_previous_solution=true,
  extrapolation=:none,
  smoother=:gs,
  adaptive_tol=true,
  adaptive_tol_params=adaptive_params
)
t_adaptive = time() - t_start

println()
println("適応的tol結果:")
println("  実行時間: $(round(t_adaptive, digits=3)) s")
println("  総反復回数: $(sum(iters_adaptive))")
println("  平均反復回数: $(round(mean(iters_adaptive), digits=1))")
println("  反復回数詳細: $(iters_adaptive)")
println("  最終温度範囲: [$(round(minimum(T_adaptive), digits=2)), $(round(maximum(T_adaptive), digits=2))] K")
println()

# 比較
println("="^70)
println("比較結果")
println("="^70)
println()

iter_reduction = (sum(iters_fixed) - sum(iters_adaptive)) / sum(iters_fixed) * 100
time_reduction = (t_fixed - t_adaptive) / t_fixed * 100

println("反復回数削減:")
println("  固定tol: $(sum(iters_fixed)) 回")
println("  適応的tol: $(sum(iters_adaptive)) 回")
println("  削減率: $(round(iter_reduction, digits=2))%")
println()

println("実行時間:")
println("  固定tol: $(round(t_fixed, digits=3)) s")
println("  適応的tol: $(round(t_adaptive, digits=3)) s")
println("  削減率: $(round(time_reduction, digits=2))%")
println()

# 精度確認
temp_diff = maximum(abs.(T_fixed .- T_adaptive))
temp_rel_error = temp_diff / maximum(abs.(T_fixed)) * 100

println("精度:")
println("  最大絶対誤差: $(round(temp_diff, digits=6)) K")
println("  最大相対誤差: $(round(temp_rel_error, digits=4))%")
println()

if temp_rel_error < 0.1
  println("✅ 精度は十分保たれています（< 0.1%）")
else
  println("⚠️  精度が低下している可能性があります（> 0.1%）")
end

println()
println("="^70)
println("テスト完了")
println("="^70)
