"""
verify_layout_old.jl

旧レイアウト版(main)の検証スクリプト

配列形状:
- Y_obs: (nt, ni, nj)  # 旧レイアウト
- T_init: (ni, nj, nk)
- q_init: (nt-1, ni, nj)  # 旧レイアウト
"""

using NPZ

# 検証用小規模パラメータ（計算時間短縮）
ni, nj, nk = 20, 20, 10
nt = 50

# 空間格子
dx = 0.00012  # m
dy = 0.0001671  # m

# z方向非均等格子（10層に削減）
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
cp_coeffs = Float64[0.0, 0.0, 0.0, 500.0]  # 定数500 J/(kg·K)
k_coeffs = Float64[0.0, 0.0, 0.0, 16.0]    # 定数16 W/(m·K)

# 初期温度
T_init = fill(300.0, ni, nj, nk)

# 観測データ（ダミー、旧レイアウト）
Y_obs = zeros(Float64, nt, ni, nj)  # (nt, ni, nj) 旧レイアウト
for t in 1:nt
  Y_obs[t, :, :] .= 300.0 + 0.1 * t
end

# スライディングウィンドウパラメータ（小規模問題用）
window_size = 30
overlap = 10
q_init_value = 1.0e5  # W/m²
cgm_iteration = 3

println("=== 旧レイアウト版検証開始 ===")
println("配列形状:")
println("  Y_obs: $(size(Y_obs))  # 期待: (nt=$nt, ni=$ni, nj=$nj)")
println("  T_init: $(size(T_init))  # 期待: (ni=$ni, nj=$nj, nk=$nk)")

# モジュール読み込み
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
include("src/IHCP_CGM.jl")
using .IHCP_CGM

# スライディングウィンドウCGM実行
println("\nスライディングウィンドウCGM開始...")
q_global, windows_info = solve_sliding_window_cgm(
  Y_obs, T_init,
  dx, dy, dz, dz_b, dz_t, dt,
  rho, cp_coeffs, k_coeffs,
  window_size, overlap, q_init_value, cgm_iteration
)

println("\n計算完了")
println("  q_global形状: $(size(q_global))")
println("  ウィンドウ数: $(length(windows_info))")

# 結果保存
output_file = "layout_old_result.npz"
npzwrite(output_file, Dict(
  "q_global" => q_global,
  "ni" => ni,
  "nj" => nj,
  "nk" => nk,
  "nt" => nt
))

println("\n結果保存: $output_file")
println("=== 旧レイアウト版検証完了 ===")
