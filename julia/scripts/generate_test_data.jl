"""
generate_test_data.jl

テスト用データ生成スクリプト

小規模な合成データを生成して、システム全体の動作確認を行う
"""

using NPZ
using Printf
using Random
using Statistics

"""
    generate_test_data()

テスト用の合成データを生成
"""
function generate_test_data()
    println("="^60)
    println("テストデータ生成")
    println("="^60)

    # 乱数シード固定（再現性のため）
    Random.seed!(42)

    # 問題サイズ
    ni, nj, nk = 3, 3, 5
    nt = 50

    println("\n問題サイズ:")
    println("  空間格子: ($ni, $nj, $nk)")
    println("  時間ステップ: $nt")

    # 初期温度場（一様温度 + 微小ノイズ）
    T_init = 300.0 .+ 0.1 .* randn(ni, nj, nk)

    println("\n初期温度場:")
    println(@sprintf("  形状: %s", size(T_init)))
    println(@sprintf("  範囲: [%.2f, %.2f] K", minimum(T_init), maximum(T_init)))
    println(@sprintf("  平均: %.2f K", mean(T_init)))

    # 観測温度（底面）: 線形加熱 + ノイズ
    # 時間とともに温度が上昇するシナリオ
    Y_obs = zeros(ni, nj, nt)

    for t in 1:nt
        # 線形加熱: 0.1 K/ステップ
        heating_rate = 0.1
        base_temp = 300.0 + (t - 1) * heating_rate

        # 空間的なばらつき（中心が最も高温）
        for i in 1:ni, j in 1:nj
            # 中心からの距離
            cx, cy = (ni + 1) / 2, (nj + 1) / 2
            dist = sqrt((i - cx)^2 + (j - cy)^2)
            spatial_factor = exp(-0.5 * dist)

            # ノイズ追加（測定誤差を模擬）
            noise = 0.05 * randn()

            Y_obs[i, j, t] = base_temp * (0.8 + 0.2 * spatial_factor) + noise
        end
    end

    println("\n観測温度:")
    println(@sprintf("  形状: %s", size(Y_obs)))
    println(@sprintf("  範囲: [%.2f, %.2f] K", minimum(Y_obs), maximum(Y_obs)))
    println(@sprintf("  初期平均: %.2f K", mean(Y_obs[:, :, 1])))
    println(@sprintf("  最終平均: %.2f K", mean(Y_obs[:, :, end])))

    # データ保存
    output_dir = joinpath(@__DIR__, "..", "data")
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    T_init_file = joinpath(output_dir, "test_T_init.npy")
    Y_obs_file = joinpath(output_dir, "test_Y_obs.npy")

    println("\nデータ保存中...")
    npzwrite(T_init_file, T_init)
    println(@sprintf("  ✓ 初期温度場: %s (%.2f KB)", T_init_file, filesize(T_init_file) / 1024))

    npzwrite(Y_obs_file, Y_obs)
    println(@sprintf("  ✓ 観測温度: %s (%.2f KB)", Y_obs_file, filesize(Y_obs_file) / 1024))

    println("\n" * "="^60)
    println("✓ テストデータ生成完了")
    println("="^60)

    return T_init_file, Y_obs_file
end

# スクリプトとして実行された場合
if abspath(PROGRAM_FILE) == @__FILE__
    generate_test_data()
end
