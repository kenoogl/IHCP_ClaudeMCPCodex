#!/usr/bin/env julia
"""
Julia版超小規模テスト（デバッグ用）

nt=20程度の超小規模データで実行時間を測定
"""

using NPZ
using Printf

PROJECT_ROOT = dirname(dirname(@__DIR__))

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

function main()
    println("=" ^ 80)
    println("Julia版超小規模テスト")
    println("=" ^ 80)

    # テストデータ読み込み（最初の20ステップのみ）
    data_path = joinpath(PROJECT_ROOT, "shared/data/T_measure_test_1min.npy")
    Y_obs_full = npzread(data_path)
    Y_obs = Y_obs_full[1:20, :, :]  # 超小規模: 20ステップのみ

    println("\nデータ形状: ", size(Y_obs))

    # パラメータ設定
    nt, ni, nj = size(Y_obs)
    window_size = 10  # 小さなウィンドウ
    overlap = 3
    cgm_iteration = 1
    q_init_value = 0.0
    dt = 0.001

    println("パラメータ: nt=$nt, window=$window_size, overlap=$overlap, CGM=$cgm_iteration")

    # 空間格子設定
    dx = 0.00012
    dy = 0.00016712741767680456
    nz = 20

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

    rho = 7823.493962874829
    cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
    k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

    # 初期温度場
    T_measure_init = Y_obs[1, :, :]
    T0 = repeat(T_measure_init, outer=(1, 1, nz))

    # 計算実行
    println("\n計算開始...")
    start_time = time()

    q_result, windows_info = solve_sliding_window_cgm(
        Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
        window_size, overlap, q_init_value, cgm_iteration
    )

    elapsed = time() - start_time

    @printf("\n実行時間: %.2f秒\n", elapsed)
    @printf("1時間ステップあたり: %.1fms\n", elapsed / nt * 1000)
    @printf("推定357ステップ実行時間: %.1f秒\n", elapsed / nt * 357)

    println("\n完了")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
