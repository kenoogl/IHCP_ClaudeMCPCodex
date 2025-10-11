#!/usr/bin/env julia
"""
Performance benchmark for tuning4 branch (allocation reduction)
Baseline comparison with commit 22fde2d

Test conditions:
- Grid: 80×100×20 (N=160,000)
- Time steps: nt=10, dt=0.001s
- CGM iterations: 1 (for performance measurement)
- CG tolerance: rtol_dhcp=1e-6, rtol_adjoint=1e-8
"""

using Pkg
Pkg.activate(".")

using Printf

# ソースコード読み込み
include("src/IHCP_CGM.jl")
using .IHCP_CGM

# テスト条件設定
ni, nj, nk = 80, 100, 20
nt = 10
dt = 0.001  # 1ms

dx = 0.12e-3  # 0.12mm
dy = 0.12e-3  # 0.12mm

# Z方向格子（非均等、表面側に集中）
Z = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
     0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0] .* 0.7e-3
ΔZ = diff(Z)

println("="^80)
println("Performance Benchmark - tuning4 Branch (Allocation Reduction)")
println("="^80)
println("Grid: $ni×$nj×$nk (N=$(ni*nj*nk))")
println("Time steps: nt=$nt, dt=$(dt)s")
println("CGM iterations: 1 (performance test)")
println("="^80)
println()

# 初期条件
T_init = fill(300.0, ni, nj, nk)  # 300K

# 測定温度（裏面S2、ダミーデータ）
Y_obs = fill(300.0, ni, nj, nt)
for t in 2:nt
  Y_obs[:, :, t] .= 300.0 .+ (t-1) * 0.1 .+ randn(ni, nj) * 0.01
end

# 初期熱流束（表面S1）
q = fill(500.0, ni, nj, nt-1)  # 500 W/m²

# 熱物性値（SUS304、600K基準）
T_ref = 600.0
rho = thermal_properties_calculator(T_ref, "density")
cp_coeffs = thermal_properties_calculator(T_ref, "specific_heat", "coeffs")
k_coeffs = thermal_properties_calculator(T_ref, "thermal_conductivity", "coeffs")

println("Starting CGM solve...")
println()

# 性能測定実行
GC.gc()  # GC実行してクリーンな状態で開始
start_time = time()

result = solve_cgm!(
  T_init, q, Y_obs,
  maxiter_cgm=1,
  sigma=0.01,
  rtol_dhcp=1e-6,
  rtol_adjoint=1e-8,
  maxiter_cg=20000,
  verbose=true,
  par="sequential"
)

total_time = time() - start_time

println()
println("="^80)
println("Performance Summary")
println("="^80)
@printf("Total runtime: %.3f seconds (%.2f minutes)\n", total_time, total_time/60)
println()

# 結果取得
q_final = result.q_optimal
T_cal = result.T_cal
J_hist = result.J_history

@printf("Final objective function J: %.6e\n", J_hist[end])
@printf("q range: [%.2f, %.2f] W/m²\n", minimum(q_final), maximum(q_final))
@printf("T range: [%.2f, %.2f] K\n", minimum(T_cal), maximum(T_cal))

println()
println("="^80)
println("Benchmark completed successfully")
println("="^80)
