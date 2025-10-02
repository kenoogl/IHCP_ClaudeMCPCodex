"""
benchmark_large.jl

大規模問題でのベンチマーク
配列形状: 新レイアウト (ni, nj, nt)
"""

using Printf

# 中規模問題パラメータ（計算時間を考慮）
ni, nj, nk = 40, 50, 10  # 空間を小規模の2倍に設定
nt = 50  # 時間ステップ数

# 空間格子
dx = 0.00012  # m
dy = 0.0001671  # m

# z方向非均等格子（簡略化10層）
dz = Float64[
  5.0e-6, 5.0e-6, 5.0e-6,
  1.0e-5, 1.0e-5, 1.0e-5,
  2.0e-5, 2.0e-5,
  5.0e-5, 5.0e-5
]

dz_b = zeros(Float64, nk)
dz_t = zeros(Float64, nk)

dz_b[1] = dz[1] / 2.0
for k in 2:nk
  dz_b[k] = (dz[k-1] + dz[k]) / 2.0
end

for k in 1:nk-1
  dz_t[k] = (dz[k] + dz[k+1]) / 2.0
end
dz_t[nk] = dz[nk] / 2.0

# 時間刻み
dt = 0.001  # s

# 物性値（SUS304）
rho = 7823.494  # kg/m³
cp_coeffs = Float64[0.0, 0.0, 0.0, 500.0]
k_coeffs = Float64[0.0, 0.0, 0.0, 16.0]

# 初期温度
T_init = fill(300.0, ni, nj, nk)

# スライディングウィンドウパラメータ（計算時間短縮）
window_size = 30
overlap = 10
q_init_value = 1.0e5  # W/m²
cgm_iteration = 3

println("="^70)
println("中規模問題ベンチマーク（新レイアウト版）")
println("="^70)
println("\n問題サイズ:")
println("  ni×nj×nk×nt = $ni×$nj×$nk×$nt")
println("  配列要素数: 温度場 = $(ni*nj*nk*nt) ($(ni*nj*nk*nt*8÷1024÷1024) MB)")
println("  ウィンドウサイズ: $window_size, オーバーラップ: $overlap")
println("  CGM反復数: $cgm_iteration")

# モジュール読み込み
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
include("src/IHCP_CGM.jl")
using .IHCP_CGM

println("\n" * "="^70)
println("ベンチマーク実行")
println("="^70)

# ウォームアップ実行（コンパイル時間除外）
println("\nウォームアップ実行中...")
Y_obs_warmup = zeros(Float64, ni, nj, nt)  # 新レイアウト
for t in 1:nt
  Y_obs_warmup[:, :, t] .= 300.0 + 0.1 * t
end
_, _ = solve_sliding_window_cgm(
  Y_obs_warmup, T_init,
  dx, dy, dz, dz_b, dz_t, dt,
  rho, cp_coeffs, k_coeffs,
  window_size, overlap, q_init_value, cgm_iteration
)
println("ウォームアップ完了")

# 本測定（3回実行して平均）
n_runs = 3
times = Float64[]

println("\n本測定開始（$n_runs 回実行）...")

for run in 1:n_runs
  # 観測データ作成（新レイアウト）
  Y_obs = zeros(Float64, ni, nj, nt)
  for t in 1:nt
    Y_obs[:, :, t] .= 300.0 + 0.1 * t
  end

  # 計算時間測定
  println("\nRun $run 開始...")
  elapsed = @elapsed begin
    q_global, windows_info = solve_sliding_window_cgm(
      Y_obs, T_init,
      dx, dy, dz, dz_b, dz_t, dt,
      rho, cp_coeffs, k_coeffs,
      window_size, overlap, q_init_value, cgm_iteration
    )
  end

  push!(times, elapsed)
  @printf("  Run %d: %.2f 秒\n", run, elapsed)
end

# 統計
mean_time = sum(times) / length(times)
std_time = sqrt(sum((times .- mean_time).^2) / length(times))
min_time = minimum(times)
max_time = maximum(times)

println("\n" * "="^70)
println("ベンチマーク結果（中規模問題）")
println("="^70)
@printf("\n平均時間: %.2f 秒\n", mean_time)
@printf("標準偏差: %.2f 秒\n", std_time)
@printf("最小時間: %.2f 秒\n", min_time)
@printf("最大時間: %.2f 秒\n", max_time)

println("\n結果を benchmark_large_new.txt に保存...")
open("benchmark_large_new.txt", "w") do f
  println(f, "中規模問題ベンチマーク（新レイアウト版）")
  println(f, "配列形状: (ni, nj, nt) = ($ni, $nj, $nt)")
  println(f, "問題サイズ: ni=$ni, nj=$nj, nk=$nk, nt=$nt")
  println(f, "ウィンドウサイズ: $window_size, オーバーラップ: $overlap")
  println(f, "CGM反復数: $cgm_iteration")
  println(f, "")
  println(f, "実行時間:")
  for (i, t) in enumerate(times)
    @printf(f, "  Run %d: %.2f 秒\n", i, t)
  end
  println(f, "")
  @printf(f, "平均時間: %.2f 秒\n", mean_time)
  @printf(f, "標準偏差: %.2f 秒\n", std_time)
  @printf(f, "最小時間: %.2f 秒\n", min_time)
  @printf(f, "最大時間: %.2f 秒\n", max_time)
end

println("="^70)
println("ベンチマーク完了")
println("="^70)
