#!/usr/bin/env julia
"""
Julia版完全一致検証用実行スクリプト

このスクリプトは、TOMLパラメータファイルから完全に同じパラメータを読み込んで
スライディングウィンドウCGM計算を実行し、Python版との完全一致検証用データを生成します。

実行方法:
  julia examples/run_exact_match.jl
"""

using NPZ
using TOML
using Printf
using LinearAlgebra

# プロジェクトのソースコードをロード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM


function load_verification_params()
    """共通パラメータファイル（TOML）を読み込み"""
    config_path = joinpath(@__DIR__, "../../shared/config/verification_params.toml")

    if !isfile(config_path)
        error("パラメータファイルが見つかりません: $config_path")
    end

    params = TOML.parsefile(config_path)
    return params
end


function prepare_grid_arrays(params)
    """格子配列を準備（infを処理）"""
    dz = Float64.(params["grid_z"]["dz"])

    # dz_tとdz_bの'inf'文字列を処理
    dz_t_raw = params["grid_z"]["dz_t"]
    dz_t = [x == "inf" ? Inf : Float64(x) for x in dz_t_raw]

    dz_b_raw = params["grid_z"]["dz_b"]
    dz_b = [x == "inf" ? Inf : Float64(x) for x in dz_b_raw]

    return dz, dz_t, dz_b
end


function main()
    """メイン実行関数"""
    println("="^60)
    println("Julia版 完全一致検証実行")
    println("="^60)

    start_time_total = time()

    # ========================================
    # 1. パラメータ読み込み
    # ========================================
    println("\n[1/6] パラメータ読み込み...")
    params = load_verification_params()

    ni = params["problem"]["ni"]
    nj = params["problem"]["nj"]
    nk = params["problem"]["nk"]
    nt = params["problem"]["nt"]
    dt = params["problem"]["dt"]
    dx = params["problem"]["dx"]
    dy = params["problem"]["dy"]

    rho = params["material"]["rho_reference_value"]
    cp_coeffs = Float64.(params["material"]["cp_coeffs"])
    k_coeffs = Float64.(params["material"]["k_coeffs"])

    dz, dz_t, dz_b = prepare_grid_arrays(params)

    window_size = params["cgm"]["window_size"]
    overlap = params["cgm"]["overlap"]
    cgm_iteration = params["cgm"]["cgm_iteration"]
    q_init_value = params["cgm"]["q_init_value"]

    T_init_value = params["initial"]["T_init_value"]

    rtol_dhcp = params["numerical"]["rtol_dhcp"]
    maxiter_dhcp = params["numerical"]["maxiter_dhcp"]
    rtol_adjoint = params["numerical"]["rtol_adjoint"]
    maxiter_adjoint = params["numerical"]["maxiter_adjoint"]

    @printf("  問題サイズ: ni=%d, nj=%d, nk=%d, nt=%d\n", ni, nj, nk, nt)
    @printf("  時間刻み: dt=%.6e s\n", dt)
    @printf("  空間刻み: dx=%.6e m, dy=%.6e m\n", dx, dy)
    @printf("  材料密度: rho=%.6f kg/m³\n", rho)
    @printf("  CGMパラメータ: window_size=%d, overlap=%d, iteration=%d\n",
            window_size, overlap, cgm_iteration)

    # ========================================
    # 2. データ読み込み
    # ========================================
    println("\n[2/6] 測定データ読み込み...")
    data_path = joinpath(@__DIR__, "../../shared/data/T_measure_test_1min.npy")

    if !isfile(data_path)
        error("測定データが見つかりません: $data_path")
    end

    Y_obs = npzread(data_path)
    println("  データ形状: $(size(Y_obs))")
    println("  データ型: $(eltype(Y_obs))")
    @printf("  温度範囲: %.2f ~ %.2f K\n", minimum(Y_obs), maximum(Y_obs))

    # データサイズの一致確認
    if size(Y_obs) != (nt, ni, nj)
        error("データ形状が不一致: $(size(Y_obs)) != ($nt, $ni, $nj)")
    end

    # ========================================
    # 3. 初期温度場の作成
    # ========================================
    println("\n[3/6] 初期温度場作成...")

    # 厳密に定義された初期温度（全格子点で一定値）
    T_init = fill(T_init_value, ni, nj, nk)

    @printf("  初期温度: %.1f K\n", T_init_value)
    println("  T_init形状: $(size(T_init))")
    println("  T_init型: $(eltype(T_init))")

    # ========================================
    # 4. スライディングウィンドウCGM計算
    # ========================================
    println("\n[4/6] スライディングウィンドウCGM計算開始...")
    expected_windows = div(nt - 1, window_size - overlap) + 1
    @printf("  予想ウィンドウ数: 約%d個\n", expected_windows)

    start_time_cgm = time()

    # スライディングウィンドウ計算
    q_result, windows_info = solve_sliding_window_cgm(
        Y_obs,
        T_init,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        dt,
        rho,
        cp_coeffs,
        k_coeffs,
        window_size,
        overlap,
        q_init_value,
        cgm_iteration;
        rtol_dhcp=rtol_dhcp,
        maxiter_dhcp=maxiter_dhcp,
        rtol_adjoint=rtol_adjoint,
        maxiter_adjoint=maxiter_adjoint
    )

    elapsed_cgm = time() - start_time_cgm

    println("\n  CGM計算完了！")
    @printf("  計算時間: %.2f秒\n", elapsed_cgm)
    println("  結果形状: $(size(q_result))")
    @printf("  熱流束範囲: %.2e ~ %.2e W/m²\n", minimum(q_result), maximum(q_result))
    @printf("  実際のウィンドウ数: %d個\n", length(windows_info))

    # ========================================
    # 5. 結果検証（DHCP forward計算）
    # ========================================
    println("\n[5/6] 結果検証（DHCP順解析）...")
    start_time_verify = time()

    T_verify = solve_dhcp_multiple_timesteps(
        T_init,
        q_result,
        nt,
        rho,
        cp_coeffs,
        k_coeffs,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        dt;
        rtol=rtol_dhcp,
        maxiter=maxiter_dhcp
    )

    elapsed_verify = time() - start_time_verify

    # 検証誤差計算（表面温度）
    T_calc_surface = T_verify[:, :, :, 1]  # Julia: 底面は1番目
    residual = T_calc_surface .- Y_obs
    rms_error = sqrt(mean(residual.^2))
    max_error = maximum(abs.(residual))

    @printf("  検証計算時間: %.2f秒\n", elapsed_verify)
    @printf("  温度誤差（RMS）: %.4e K\n", rms_error)
    @printf("  温度誤差（最大）: %.4e K\n", max_error)

    # ========================================
    # 6. 結果保存
    # ========================================
    println("\n[6/6] 結果保存...")
    output_dir = joinpath(@__DIR__, "../../shared/results/validation")
    mkpath(output_dir)

    output_path = joinpath(output_dir, "julia_exact.npz")

    npzwrite(output_path, Dict(
        # 入力データ
        "Y_obs" => Y_obs,
        "T_init" => T_init,
        # パラメータ
        "dt" => dt,
        "dx" => dx,
        "dy" => dy,
        "dz" => dz,
        "dz_t" => dz_t,
        "dz_b" => dz_b,
        "rho" => rho,
        "cp_coeffs" => cp_coeffs,
        "k_coeffs" => k_coeffs,
        # 結果
        "q_result" => q_result,
        "T_verify" => T_verify,
        # 統計情報
        "elapsed_time_cgm" => elapsed_cgm,
        "elapsed_time_verify" => elapsed_verify,
        "rms_error" => rms_error,
        "max_error" => max_error,
        "n_windows" => length(windows_info)
    ))

    elapsed_total = time() - start_time_total

    file_size_mb = filesize(output_path) / 1024 / 1024
    @printf("  保存完了: %s\n", output_path)
    @printf("  ファイルサイズ: %.2f MB\n", file_size_mb)

    # ========================================
    # サマリー表示
    # ========================================
    println("\n" * "="^60)
    println("実行完了サマリー")
    println("="^60)
    @printf("総実行時間: %.2f秒\n", elapsed_total)
    @printf("  CGM計算:  %.2f秒 (%.1f%%)\n", elapsed_cgm, elapsed_cgm/elapsed_total*100)
    @printf("  検証計算: %.2f秒 (%.1f%%)\n", elapsed_verify, elapsed_verify/elapsed_total*100)
    println("\n結果統計:")
    @printf("  熱流束範囲: %.4e ~ %.4e W/m²\n", minimum(q_result), maximum(q_result))
    @printf("  温度誤差（RMS）: %.4e K\n", rms_error)
    @printf("  温度誤差（最大）: %.4e K\n", max_error)
    println("\n次のステップ:")
    println("  比較実行: python python/validation/compare_exact_match.py")
    println("="^60)

    return q_result, T_verify
end


# スクリプトとして実行された場合
if abspath(PROGRAM_FILE) == @__FILE__
    try
        q_result, T_verify = main()
    catch e
        println("\n❌ エラーが発生しました: ", e)
        Base.showerror(stdout, e, catch_backtrace())
        println()
        exit(1)
    end
end
