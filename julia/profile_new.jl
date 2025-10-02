"""
profile_new.jl

新レイアウト版（tuning2）のプロファイリング
"""

using Profile
using Printf

# 中規模問題パラメータ
ni, nj, nk = 40, 50, 10
nt = 50

# 空間格子
dx = 0.00012
dy = 0.0001671

# z方向非均等格子
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

dt = 0.001
rho = 7823.494
cp_coeffs = Float64[0.0, 0.0, 0.0, 500.0]
k_coeffs = Float64[0.0, 0.0, 0.0, 16.0]

T_init = fill(300.0, ni, nj, nk)

window_size = 30
overlap = 10
q_init_value = 1.0e5
cgm_iteration = 3

println("="^70)
println("新レイアウト版プロファイリング")
println("="^70)
println("\n問題サイズ: ni=$ni, nj=$nj, nk=$nk, nt=$nt")

# モジュール読み込み
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
include("src/IHCP_CGM.jl")
using .IHCP_CGM

# ウォームアップ実行
println("\nウォームアップ実行中...")
Y_obs_warmup = zeros(Float64, ni, nj, nt)
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

# 本測定（プロファイリング）
Y_obs = zeros(Float64, ni, nj, nt)
for t in 1:nt
  Y_obs[:, :, t] .= 300.0 + 0.1 * t
end

println("\nプロファイリング開始...")
Profile.clear()
@profile begin
  q_global, windows_info = solve_sliding_window_cgm(
    Y_obs, T_init,
    dx, dy, dz, dz_b, dz_t, dt,
    rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, cgm_iteration
  )
end
println("プロファイリング完了")

# プロファイル結果を保存
println("\nプロファイル結果を profile_new.txt に保存中...")
open("profile_new.txt", "w") do io
  Profile.print(io, format=:flat, sortedby=:count)
end

# トップ20関数の詳細を出力
println("\n" * "="^70)
println("トップ20関数（呼び出し回数順）")
println("="^70)
Profile.print(format=:flat, sortedby=:count, maxdepth=20)

println("\n" * "="^70)
println("プロファイル結果保存完了: profile_new.txt")
println("="^70)
