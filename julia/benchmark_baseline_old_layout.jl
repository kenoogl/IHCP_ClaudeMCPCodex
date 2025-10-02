"""
Phase 2.2 ベースライン測定: 旧メモリレイアウト

完全一致検証サイズ（80×100×20、N=160,000）での性能測定
旧レイアウト: (nt, ni, nj, nk) および (nt-1, ni, nj)

注意: このスクリプトは古いコミット (06e5c56以前) で実行する必要があります
"""

using LinearAlgebra
using Statistics
using Printf

# シングルスレッド実行（公平な比較のため）
BLAS.set_num_threads(1)

# モジュールのパスを追加
push!(LOAD_PATH, joinpath(@__DIR__, "src"))

# モジュール読み込み
include("src/IHCP_CGM.jl")
using .IHCP_CGM

# ===========================
# 問題設定（完全一致検証サイズ）
# ===========================
ni, nj, nk = 80, 100, 20
nt = 357  # 356時間ステップ + 初期時刻

# 格子幅
dx = 0.12e-3  # m
dy = 0.12e-3  # m
dt = 1e-3     # s

# z方向格子（幾何級数）
depth = 2.4e-3  # m
ratio = 1.2
dz = zeros(Float64, nk)
dz[1] = depth * (1 - ratio) / (1 - ratio^nk)
for k in 2:nk
  dz[k] = dz[k-1] * ratio
end

# 界面距離
dz_b = zeros(Float64, nk)
dz_t = zeros(Float64, nk)
dz_b[1] = dz[1] / 2
for k in 2:nk
  dz_b[k] = dz_b[k-1] + dz[k-1]/2 + dz[k]/2
end
for k in 1:nk-1
  dz_t[k] = dz_b[k+1]
end
dz_t[nk] = dz_b[nk] + dz[nk]/2

# 材料物性（一定値）
rho = 7900.0  # kg/m³
cp_coeffs = [0.0, 0.0, 0.0, 500.0]  # J/(kg·K)
k_coeffs = [0.0, 0.0, 0.0, 16.0]    # W/(m·K)

# 初期条件・境界条件
T_initial = fill(300.0, ni, nj, nk)
q_surface = fill(1000.0, nt-1, ni, nj)  # 旧レイアウト: (nt-1, ni, nj)

println("="^70)
println("Phase 2.2 ベースライン測定: 旧メモリレイアウト")
println("="^70)
println("問題サイズ: $(ni)×$(nj)×$(nk)  (N=$(ni*nj*nk))")
println("時間ステップ数: $(nt-1)")
println("配列レイアウト: T_all=(nt,ni,nj,nk), q_surface=(nt-1,ni,nj)")
println("="^70)

# ===========================
# コンパイル実行（10ステップ）
# ===========================
println("\n[1/3] コンパイル実行（10ステップ）...")
q_warmup = q_surface[1:10, :, :]
T_all_warmup = solve_dhcp!(
  T_initial, q_warmup, 11,
  dx, dy, dz, dz_b, dz_t, dt,
  rho, cp_coeffs, k_coeffs;
  rtol=1e-6, maxiter=20000
)
println("  コンパイル完了")

# ===========================
# 単一ステップ性能測定（10ステップ×5試行）
# ===========================
println("\n[2/3] 単一ステップ性能測定（10ステップ×5試行）...")
dhcp_times = Float64[]

for trial in 1:5
  q_test = q_surface[1:10, :, :]
  T_test = copy(T_initial)

  t_start = time()
  T_all_test = solve_dhcp!(
    T_test, q_test, 11,
    dx, dy, dz, dz_b, dz_t, dt,
    rho, cp_coeffs, k_coeffs;
    rtol=1e-6, maxiter=20000
  )
  t_elapsed = (time() - t_start) * 1000  # ms

  push!(dhcp_times, t_elapsed)
  println("  試行 $trial: $(round(t_elapsed, digits=1)) ms")
end

avg_time = mean(dhcp_times)
std_time = std(dhcp_times)
per_step_time = avg_time / 10

println("\n  平均時間（10ステップ）: $(round(avg_time, digits=1)) ± $(round(std_time, digits=1)) ms")
println("  ステップあたり: $(round(per_step_time, digits=1)) ms/step")

# ===========================
# フルスケール実行（357ステップ）
# ===========================
println("\n[3/3] フルスケール実行（357ステップ）...")
println("  予測時間: $(round(per_step_time * (nt-1) / 1000, digits=1)) 秒")

T_full = copy(T_initial)
t_start = time()
T_all_full = solve_dhcp!(
  T_full, q_surface, nt,
  dx, dy, dz, dz_b, dz_t, dt,
  rho, cp_coeffs, k_coeffs;
  rtol=1e-6, maxiter=20000
)
full_time = (time() - t_start)

full_per_step = full_time * 1000 / (nt-1)

println("  実測時間: $(round(full_time, digits=2)) 秒")
println("  ステップあたり: $(round(full_per_step, digits=1)) ms/step")

# ===========================
# 結果サマリー
# ===========================
println("\n" * "="^70)
println("旧メモリレイアウト ベースライン測定結果")
println("="^70)
println("問題サイズ: $(ni)×$(nj)×$(nk)  (N=$(ni*nj*nk))")
println("時間ステップ数: $(nt-1)")
println("")
println("単一ステップ性能:")
println("  平均: $(round(per_step_time, digits=1)) ms/step")
println("  標準偏差: $(round(std_time/10, digits=1)) ms/step")
println("")
println("フルスケール実行:")
println("  総時間: $(round(full_time, digits=2)) 秒")
println("  ステップあたり: $(round(full_per_step, digits=1)) ms/step")
println("  スループット: $(round((nt-1)/full_time, digits=2)) steps/sec")
println("="^70)

# 結果をファイルに保存
open("benchmark_baseline_results.txt", "w") do io
  println(io, "Phase 2.2 旧メモリレイアウト ベースライン測定")
  println(io, "実行日時: $(now())")
  println(io, "")
  println(io, "問題サイズ: $(ni)×$(nj)×$(nk)  (N=$(ni*nj*nk))")
  println(io, "時間ステップ数: $(nt-1)")
  println(io, "配列レイアウト: T_all=(nt,ni,nj,nk), q_surface=(nt-1,ni,nj)")
  println(io, "")
  println(io, "単一ステップ性能: $(round(per_step_time, digits=1)) ± $(round(std_time/10, digits=1)) ms/step")
  println(io, "フルスケール性能: $(round(full_per_step, digits=1)) ms/step")
  println(io, "総実行時間: $(round(full_time, digits=2)) 秒")
end

println("\n結果を benchmark_baseline_results.txt に保存しました")
