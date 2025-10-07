"""
profile_new_mf.jl

マトリクスフリー版PBICGSTAB!のプロファイリング
"""

using Profile
using Printf

# 中規模問題パラメータ
ni, nj, nk = 40, 50, 10
nt = 50

# 空間格子
dx = 0.00012
dy = 0.0001671

# z方向非均等格子（IHCP形式）
dz = Float64[
  5.0e-6, 5.0e-6, 5.0e-6,
  1.0e-5, 1.0e-5, 1.0e-5,
  2.0e-5, 2.0e-5,
  5.0e-5, 5.0e-5
]
dz_b = copy(dz)  # 簡易版: セル中心幅と同じ
dz_t = copy(dz)

dt = 0.001
ρ = 7823.494
cp_coeffs = Float64[0.0, 0.0, 0.0, 500.0]
k_coeffs = Float64[0.0, 0.0, 0.0, 16.0]

T_init = fill(300.0, ni, nj, nk)

window_size = 30
overlap = 10
q_init_value = 1.0e5
cgm_iteration = 3

println("="^70)
println("マトリクスフリー版PBICGSTAB!プロファイリング")
println("="^70)
println("\n問題サイズ: ni=$ni, nj=$nj, nk=$nk, nt=$nt")

# モジュール読み込み
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
include("src/IHCP_CGM.jl")
using .IHCP_CGM
using .IHCP_CGM.GridTransform: convert_to_guard_cell_grid

# Heat3ds形式へ変換（ガイドセル含む）
Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

# WorkBuffers作成（ガイドセル含む）
wk = WorkBuffers(ni+2, nj+2, nk+2)

# 表面熱流束（簡易版）
q_surface = zeros(Float64, ni, nj, nt-1)
for t in 1:(nt-1)
  q_surface[:, :, t] .= q_init_value
end

# ウォームアップ実行
println("\nウォームアップ実行中...")
T_warmup = solve_dhcp!(
  T_init, q_surface, wk,
  nt, ρ, cp_coeffs, k_coeffs,
  dx, dy, Z, ΔZ, dt,
  rtol=1e-6, maxiter=20000, verbose=false, par="sequential"
)
println("ウォームアップ完了")

# 本測定（プロファイリング）
println("\nプロファイリング開始...")
Profile.clear()

# WorkBuffers初期化
wk = WorkBuffers(ni+2, nj+2, nk+2)

@profile begin
  T_result = solve_dhcp!(
    T_init, q_surface, wk,
    nt, ρ, cp_coeffs, k_coeffs,
    dx, dy, Z, ΔZ, dt,
    rtol=1e-6, maxiter=20000, verbose=false, par="sequential"
  )
end
println("プロファイリング完了")

# プロファイル結果を保存
println("\nプロファイル結果を profile_new_mf.txt に保存中...")
open("profile_new_mf.txt", "w") do io
  Profile.print(io, format=:flat, sortedby=:count)
end

# トップ20関数の詳細を出力
println("\n" * "="^70)
println("トップ20関数（呼び出し回数順）")
println("="^70)
Profile.print(format=:flat, sortedby=:count, maxdepth=20)

println("\n" * "="^70)
println("プロファイル結果保存完了: profile_new_mf.txt")
println("="^70)
