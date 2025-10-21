#!/usr/bin/env julia
"""
テストデータ読み込み例

Python版ベンチマークで生成された約1分実行用テストデータを読み込む
"""

using NPZ

# プロジェクトルート
PROJECT_ROOT = dirname(dirname(@__DIR__))

# テストデータ読み込み
data_path = joinpath(PROJECT_ROOT, "shared/data/T_measure_test_1min.npy")

println("=" ^ 80)
println("テストデータ読み込み")
println("=" ^ 80)
println("\nデータファイル: ", data_path)

if !isfile(data_path)
    error("データファイルが見つかりません: $data_path")
end

# データ読み込み
Y_obs = npzread(data_path)

println("\nデータ情報:")
println("  形状: ", size(Y_obs))
println("  データ型: ", eltype(Y_obs))
println("  サイズ: ", round(sizeof(Y_obs) / 1e6, digits=2), " MB")
println("  温度範囲: ", round(minimum(Y_obs), digits=2), "K ~ ", round(maximum(Y_obs), digits=2), "K")

# データ形状の確認
nt, ni, nj = size(Y_obs)
println("\n格子情報:")
println("  時間ステップ数: nt = ", nt)
println("  空間格子サイズ: ni × nj = ", ni, " × ", nj)
println("  総測定点数: ", ni * nj)

# 推奨パラメータ
window_size = 71
overlap = 17

println("\n推奨スライディングウィンドウパラメータ:")
println("  ウィンドウサイズ: window_size = ", window_size)
println("  オーバーラップ: overlap = ", overlap)
println("  予想ウィンドウ数: ", ceil(Int, (nt - overlap) / (window_size - overlap)))

# 時間刻み
dt = 0.001  # 1ms

println("\n時間情報:")
println("  時間刻み: dt = ", dt * 1000, " ms")
println("  総時間: ", nt * dt, " s = ", nt * dt * 1000, " ms")

println("\n" * "=" ^ 80)
println("Python版との性能比較の準備完了")
println("Python版予想実行時間: 約54秒（1CGM反復）")
println("=" ^ 80)

# 初期温度場の設定（実際の計算では必要）
nz = 20  # z方向格子数
T_measure_init = Y_obs[1, :, :]  # 最初の時間ステップ
T0 = repeat(T_measure_init, 1, 1, nz)  # (ni, nj, nz)

println("\n初期温度場:")
println("  形状: ", size(T0))
println("  データ型: ", eltype(T0))
println("  温度範囲: ", round(minimum(T0), digits=2), "K ~ ", round(maximum(T0), digits=2), "K")

# 初期熱流束（ゼロ初期化）
q_init = zeros(nt - 1, ni, nj)

println("\n初期熱流束:")
println("  形状: ", size(q_init))
println("  データ型: ", eltype(q_init))

println("\n" * "=" ^ 80)
println("データ読み込み完了")
println("=" ^ 80)
