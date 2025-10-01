#!/usr/bin/env julia
"""
Julia版縮小テストデータ実行スクリプト（1分想定）

Python版との比較検証用にJulia版を実行し、結果を保存する

実行方法:
  cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia
  julia --project=. examples/run_test_1min.jl

出力:
  - shared/results/julia_test_1min.npz: 計算結果
    - q_result: 熱流束結果 (nt-1, ni, nj)
    - elapsed_time: 実行時間 [秒]
    - J_final: 最終目的関数値
    - windows_count: ウィンドウ数
"""

using NPZ
using Printf
using Statistics

# プロジェクトルート
PROJECT_ROOT = dirname(dirname(@__DIR__))

# IHCP_CGMモジュールのインクルード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

function main()
    println("=" ^ 80)
    println("Julia版縮小テストデータ実行（1分想定）")
    println("=" ^ 80)

    # テストデータ読み込み
    data_path = joinpath(PROJECT_ROOT, "shared/data/T_measure_test_1min.npy")
    println("\nデータファイル: ", data_path)

    if !isfile(data_path)
        error("データファイルが見つかりません: $data_path")
    end

    Y_obs = npzread(data_path)
    println("データ形状: ", size(Y_obs))
    @printf("データサイズ: %.2f MB\n", sizeof(Y_obs) / 1e6)

    # パラメータ設定（Python版と同じ）
    nt, ni, nj = size(Y_obs)
    window_size = 71
    overlap = 17
    cgm_iteration = 1  # 最小限の反復数
    q_init_value = 0.0  # 初期熱流束
    dt = 0.001  # 1ms

    println("\nパラメータ: nt=$nt, window=$window_size, overlap=$overlap, CGM=$cgm_iteration")

    # 空間格子設定（Python版と同じ）
    dx = 0.00012  # m
    dy = 0.00016712741767680456  # m
    nz = 20

    # z方向格子（Python版と同じ）
    dz = [
        7.32951631e-05, 6.30857315e-05, 5.42983923e-05, 4.67350594e-05,
        4.02252384e-05, 3.46221835e-05, 2.97995895e-05, 2.56487444e-05,
        2.20760789e-05, 1.90010572e-05, 1.63543615e-05, 1.40763294e-05,
        1.21156090e-05, 1.04280013e-05, 8.97546388e-06, 7.72525335e-06,
        6.64918718e-06, 5.72300844e-06, 4.92583902e-06, 4.23970893e-06
    ]

    dz_t = [
        1.04838029e-04, 5.86920619e-05, 5.05167258e-05, 4.34801489e-05,
        3.74237109e-05, 3.22108865e-05, 2.77241670e-05, 2.38624116e-05,
        2.05385680e-05, 1.76777093e-05, 1.52153454e-05, 1.30959692e-05,
        1.12718051e-05, 9.70173258e-06, 8.35035862e-06, 7.18722027e-06,
        6.18609781e-06, 5.32442373e-06, 6.70262844e-06, Inf
    ]

    dz_b = [
        Inf, 1.04838029e-04, 5.86920619e-05, 5.05167258e-05,
        4.34801489e-05, 3.74237109e-05, 3.22108865e-05, 2.77241670e-05,
        2.38624116e-05, 2.05385680e-05, 1.76777093e-05, 1.52153454e-05,
        1.30959692e-05, 1.12718051e-05, 9.70173258e-06, 8.35035862e-06,
        7.18722027e-06, 6.18609781e-06, 5.32442373e-06, 6.70262844e-06
    ]

    # 熱物性値（Python版と同じ）
    rho = 7823.493962874829  # kg/m³
    cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
    k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

    # 初期温度場の設定
    T_measure_init = Y_obs[1, :, :]  # Julia: 1始まり
    T0 = repeat(T_measure_init, outer=(1, 1, nz))  # (ni, nj, nz)

    @printf("初期温度場: %s, T_range=[%.2f, %.2f]K\n",
            size(T0), minimum(T0), maximum(T0))

    # 計算実行
    println("\n" * "=" ^ 80)
    println("Julia版計算開始")
    println("=" ^ 80)

    start_time = time()

    q_result, windows_info = solve_sliding_window_cgm(
        Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
        window_size, overlap, q_init_value, cgm_iteration
    )

    elapsed = time() - start_time

    println("=" ^ 80)
    @printf("総実行時間: %.2f秒\n", elapsed)
    @printf("1時間ステップあたり: %.1fms\n", elapsed / nt * 1000)

    # 結果保存
    output_path = joinpath(PROJECT_ROOT, "shared/results/julia_test_1min.npz")
    mkpath(dirname(output_path))

    windows_count = length(windows_info)
    J_final = windows_count > 0 ? windows_info[end].J_final : 0.0

    # NPZ形式で保存（Python版と互換）
    npzwrite(output_path, Dict(
        "q_result" => q_result,
        "elapsed_time" => elapsed,
        "windows_count" => windows_count,
        "J_final" => J_final,
        # パラメータ
        "nt" => nt,
        "window_size" => window_size,
        "overlap" => overlap,
        "cgm_iteration" => cgm_iteration,
        "dt" => dt,
        "dx" => dx,
        "dy" => dy,
        "dz" => dz,
        "dz_b" => dz_b,
        "dz_t" => dz_t,
        "rho" => rho,
        "cp_coeffs" => cp_coeffs,
        "k_coeffs" => k_coeffs
    ))

    println("\n結果保存: ", output_path)
    println("  q_result: ", size(q_result))
    @printf("  elapsed: %.2fs\n", elapsed)
    println("  windows: ", windows_count)
    @printf("  J_final: %.6e\n", J_final)
    @printf("  q_range: [%.3e, %.3e] W/m²\n", minimum(q_result), maximum(q_result))
    @printf("  q_mean: %.3e W/m²\n", mean(q_result))
    @printf("  q_std: %.3e W/m²\n", std(q_result))

    println("\n" * "=" ^ 80)
    println("完了")
    println("=" ^ 80)

    # Python版との比較情報
    python_result_path = joinpath(PROJECT_ROOT, "shared/results/python_test_1min.npz")
    if isfile(python_result_path)
        println("\n" * "=" ^ 80)
        println("Python版との比較")
        println("=" ^ 80)

        python_data = npzread(python_result_path)
        python_elapsed = python_data["elapsed_time"]
        python_q = python_data["q_result"]

        @printf("実行時間:\n")
        @printf("  Python: %.2fs\n", python_elapsed)
        @printf("  Julia:  %.2fs\n", elapsed)
        @printf("  比率:   %.2fx (Julia/Python)\n", elapsed / python_elapsed)
        @printf("  改善:   %.1f%% %s\n",
                abs(1 - elapsed / python_elapsed) * 100,
                elapsed < python_elapsed ? "高速化" : "低速化")

        # 数値比較
        diff = q_result - python_q
        @printf("\n数値比較:\n")
        @printf("  最大絶対誤差: %.6e\n", maximum(abs.(diff)))
        @printf("  平均絶対誤差: %.6e\n", mean(abs.(diff)))
        @printf("  相対誤差:     %.6e\n", maximum(abs.(diff ./ (python_q .+ 1e-10))))
        @printf("  相関係数:     %.6f\n", cor(vec(q_result), vec(python_q)))

        println("=" ^ 80)
    end

    return nothing
end

# メイン実行
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
