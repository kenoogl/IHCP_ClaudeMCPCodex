"""
Phase 2.2 メモリレイアウト最適化 性能測定
完全一致検証サイズ: 80×100×20格子、357ステップ
"""

using LinearAlgebra
using Statistics
using Printf

# BLAS逐次実行設定（ベースライン条件）
BLAS.set_num_threads(1)

# モジュール読み込み
using IHCP_CGM

# ============================================================
# 完全一致検証パラメータ（FINAL_VERIFICATION_SUMMARY.mdより）
# ============================================================
ni, nj, nk = 80, 100, 20
nt = 357
dx, dy, dt = 0.12e-3, 0.12e-3, 1e-3

# dz配列（表面集中、幾何級数）
z_total = 0.7e-3
base = 1.6
dz = zeros(nk)
r_sum = sum(base^(i-1) for i in 1:nk)
dz[1] = z_total / r_sum
for i in 2:nk
  dz[i] = dz[i-1] * base
end

dz_b = copy(dz)
dz_b[1] = Inf
dz_t = copy(dz)
dz_t[end] = Inf

rho = 7900.0
cp_coeffs = [0.0, 0.0, 0.0, 500.0]  # 定数500 J/(kg·K)
k_coeffs = [0.0, 0.0, 0.0, 15.0]    # 定数15 W/(m·K)

# 初期条件
T_initial = fill(300.0, ni, nj, nk)

# 境界条件（熱流束一定）
q_surface = fill(1000.0, ni, nj, nt-1)

println("="^70)
println("Phase 2.2 メモリレイアウト最適化 性能測定")
println("="^70)
println("問題サイズ: $(ni)×$(nj)×$(nk)  (N=$(ni*nj*nk)格子点)")
println("時間ステップ: $(nt)")
println("総時間: $(nt*dt) s")
println("BLAS threads: $(BLAS.get_num_threads())")
println("="^70)

# ============================================================
# 1. コンパイル実行（初回JIT）
# ============================================================
println("\n[1] コンパイル実行（初回JIT）...")
@time T_all_warmup = solve_dhcp!(
  T_initial, q_surface[:, :, 1:10], 11, rho, cp_coeffs, k_coeffs,
  dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000, verbose=false
)
println("  コンパイル完了")

# ============================================================
# 2. DHCP単一ステップ性能測定（10回平均）
# ============================================================
println("\n[2] DHCP単一ステップ性能測定（10時間ステップ、5回繰り返し）")
dhcp_times = Float64[]
for trial in 1:5
  t = @elapsed begin
    T_all = solve_dhcp!(
      T_initial, q_surface[:, :, 1:10], 11, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000, verbose=false
    )
  end
  push!(dhcp_times, t * 1000)  # msに変換
  @printf("  試行%d: %.1f ms\n", trial, t * 1000)
end

dhcp_mean = mean(dhcp_times)
dhcp_std = std(dhcp_times)
dhcp_min = minimum(dhcp_times)
dhcp_max = maximum(dhcp_times)

println("\nDHCP統計（10ステップ）:")
@printf("  平均: %.1f ms\n", dhcp_mean)
@printf("  標準偏差: %.1f ms\n", dhcp_std)
@printf("  最小: %.1f ms\n", dhcp_min)
@printf("  最大: %.1f ms\n", dhcp_max)
@printf("  単一ステップあたり: %.1f ms\n", dhcp_mean / 10)

# ============================================================
# 3. フルスケール計算性能測定（357ステップ）
# ============================================================
println("\n[3] フルスケール計算性能測定（357ステップ、1回実行）")
println("  DHCP検証計算開始...")

GC.gc()  # GC実行してクリーンな状態で計測
full_time = @elapsed begin
  T_all_full = solve_dhcp!(
    T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000, verbose=false
  )
end

println("  DHCP検証計算完了")
@printf("  総実行時間: %.1f s\n", full_time)
@printf("  単一ステップあたり: %.1f ms\n", full_time * 1000 / (nt - 1))

# ============================================================
# 4. 結果サマリー
# ============================================================
println("\n" * "="^70)
println("性能測定結果サマリー")
println("="^70)
@printf("DHCP単一ステップ平均: %.1f ms\n", dhcp_mean / 10)
@printf("DHCP 357ステップ総時間: %.1f s\n", full_time)
@printf("メモリ使用量推定: %.1f MB (T_all配列)\n",
  ni * nj * nk * nt * 8 / 1024 / 1024)

println("\n[参考] Phase 2.2以前のベースライン（phase2_ilu_failure_analysis.mdより）:")
println("  DHCP単一ステップ: 325 ms（対角前処理、N=160,000）")
println("  期待改善率: 1.5-2倍（メモリレイアウト最適化）")

if dhcp_mean / 10 < 325
  improvement = 325 / (dhcp_mean / 10)
  @printf("\n✅ 改善達成: %.2f倍高速化\n", improvement)
else
  @printf("\n⚠️  改善未達: %.2f倍（ベースラインより遅い）\n", (dhcp_mean / 10) / 325)
end

println("="^70)
