#!/usr/bin/env julia
"""
Julia版 300ステップ一致検証テストスクリプト

このスクリプトは、Python版との計算結果一致確認用に、T_measure_700um_1ms.npyの
最初の300ステップを使用してスライディングウィンドウCGM計算を実行します。

実行方法:
  JULIA_NUM_THREADS=1 julia julia/examples/test_python_julia_match_300steps.jl
"""

using NPZ
using Printf
using LinearAlgebra
using Statistics

# プロジェクトのソースコードをロード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM


function main()
    """メイン実行関数"""
    println("="^80)
    println("Julia版 300ステップ一致検証テスト")
    println("="^80)
    println("Juliaスレッド数: ", Threads.nthreads())

    start_time_total = time()

    # ========================================
    # 1. パラメータ設定
    # ========================================
    println("\n[1/6] パラメータ設定...")

    # 時間パラメータ
    nt = 300
    dt = 0.001  # 1 ms

    # 空間パラメータ（python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.pyから）
    ni = 80
    nj = 100
    nk = 20
    dx = 0.00012  # m
    dy = 0.00012 * sin(deg2rad(80)) / sin(deg2rad(45))  # m
    Lz = 0.0005  # m
    stretch_factor = 3

    # Z方向格子（非均等）
    z_faces_normalized = range(1, 0, length=nk+1)
    z_faces = Lz .- (Lz / (exp(stretch_factor) - 1)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1)

    z_centers = zeros(Float64, nk)
    z_centers[1] = z_faces[1]
    z_centers[end] = z_faces[end]
    z_centers[2:end-1] = (z_faces[2:end-2] .+ z_faces[3:end-1]) ./ 2

    dz = diff(z_faces)

    dz_t = zeros(Float64, nk)
    dz_t[end] = Inf
    dz_t[1:end-1] = z_centers[2:end] .- z_centers[1:end-1]

    dz_b = zeros(Float64, nk)
    dz_b[1] = Inf
    dz_b[2:end] = z_centers[2:end] .- z_centers[1:end-1]

    # 材料物性値（cp_coeffs, k_coeffsはpython/original/から取得）
    rho = 7823.493962874829
    cp_coeffs = Float64[2.0092965857857302e-10, -3.426055709046039e-7, 0.13492793577599277, 469.85285976023886]
    k_coeffs = Float64[4.799122446183406e-12, -8.182993477116119e-9, 0.016176544533897493, 8.117517482517492]

    # CGMパラメータ
    window_size = 71
    overlap = 17
    cgm_iteration = 2
    q_init_value = 0.0

    # 数値ソルバーパラメータ
    rtol_dhcp = 1.0e-6
    maxiter_dhcp = 20000
    rtol_adjoint = 1.0e-8
    maxiter_adjoint = 20000

    @printf("  時間ステップ数: nt=%d\n", nt)
    @printf("  時間刻み: dt=%.6f s\n", dt)
    @printf("  空間刻み: dx=%.6e m, dy=%.6e m\n", dx, dy)
    @printf("  材料密度: rho=%.6f kg/m³\n", rho)
    @printf("  CGMパラメータ: window_size=%d, overlap=%d, iteration=%d\n",
            window_size, overlap, cgm_iteration)

    # ========================================
    # 2. データ読み込み
    # ========================================
    println("\n[2/6] 測定データ読み込み...")
    data_path = joinpath(@__DIR__, "../../shared/data/T_measure_700um_1ms.npy")

    if !isfile(data_path)
        error("測定データが見つかりません: $data_path")
    end

    # 全データ読み込み
    T_measure_full = npzread(data_path)
    println("  全データ形状: $(size(T_measure_full))")

    # 300ステップに制限
    # Python形式 (nt, ni, nj) → Julia形式 (ni, nj, nt) に変換
    Y_obs_python = T_measure_full[1:nt, :, :]
    Y_obs = permutedims(Y_obs_python, (2, 3, 1))
    println("  使用データ形状（Julia形式）: $(size(Y_obs))")
    println("  データ型: $(eltype(Y_obs))")
    @printf("  温度範囲: %.2f ~ %.2f K\n", minimum(Y_obs), maximum(Y_obs))

    # ========================================
    # 3. 初期温度場の作成
    # ========================================
    println("\n[3/6] 初期温度場作成...")

    # 最初のフレームから初期温度を設定
    T_measure_init = Y_obs[:, :, 1]
    T_init = repeat(T_measure_init, outer=(1, 1, 1))
    T_init = repeat(reshape(T_measure_init, ni, nj, 1), outer=(1, 1, nk))
    T_init = convert(Array{Float64,3}, T_init)

    println("  T_init形状: $(size(T_init))")
    println("  T_init型: $(eltype(T_init))")
    @printf("  初期温度範囲: %.2f ~ %.2f K\n", minimum(T_init), maximum(T_init))

    # ========================================
    # 4. WorkBuffers と Z配列の準備
    # ========================================
    println("\n[4/6] WorkBuffers と Z配列の準備...")

    # convert_to_guard_cell_gridを使用してZ, ΔZを計算
    Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

    # WorkBuffers作成（ガイドセル含むサイズ: ni+2, nj+2, nk+2）
    wk = IHCP_CGM.WorkBuffers(ni+2, nj+2, nk+2)

    println("  Z配列形状: $(size(Z))")
    println("  ΔZ配列形状: $(size(ΔZ))")
    println("  WorkBuffers作成完了")

    # ========================================
    # 5. スライディングウィンドウCGM計算
    # ========================================
    println("\n[5/6] スライディングウィンドウCGM計算開始...")
    expected_windows = div(nt - 1, window_size - overlap) + 1
    @printf("  予想ウィンドウ数: 約%d個\n", expected_windows)

    start_time_cgm = time()

    # スライディングウィンドウ計算
    q_result, windows_info = solve_sliding_window_cgm(
        Y_obs,
        T_init,
        wk,
        dx,
        dy,
        Z,
        ΔZ,
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
    # 6. 結果検証（DHCP forward計算）
    # ========================================
    println("\n[6/6] 結果検証（DHCP順解析）...")
    start_time_verify = time()

    # WorkBuffers作成（検証用、ガイドセル含むサイズ）
    wk_verify = IHCP_CGM.WorkBuffers(ni+2, nj+2, nk+2)

    T_verify = IHCP_CGM.solve_dhcp!(
        T_init,
        q_result,
        wk_verify,
        nt,
        rho,
        cp_coeffs,
        k_coeffs,
        dx,
        dy,
        Z,
        ΔZ,
        dt;
        rtol=rtol_dhcp,
        maxiter=maxiter_dhcp
    )

    elapsed_verify = time() - start_time_verify

    # 検証誤差計算（表面温度）
    # Y_obsは表面（Z上面、k=nk）の温度、形状: (ni, nj, nt)
    # T_verifyは全格子点の温度、形状: (ni, nj, nk, nt)
    T_calc_surface = T_verify[:, :, nk, :]  # 表面（k=nk）の温度を取得
    residual = T_calc_surface .- Y_obs
    rms_error = sqrt(mean(residual.^2))
    max_error = maximum(abs.(residual))

    @printf("  検証計算時間: %.2f秒\n", elapsed_verify)
    @printf("  温度誤差（RMS）: %.4e K\n", rms_error)
    @printf("  温度誤差（最大）: %.4e K\n", max_error)

    # ========================================
    # 7. 結果保存
    # ========================================
    println("\n[7/7] 結果保存...")
    output_dir = joinpath(@__DIR__, "../../shared/results")
    mkpath(output_dir)

    # 個別に保存（.npy形式、非圧縮）
    q_path = joinpath(output_dir, "julia_300steps_q.npy")
    T_path = joinpath(output_dir, "julia_300steps_T.npy")

    npzwrite(q_path, q_result)
    npzwrite(T_path, T_verify)

    # メタデータ付き圧縮版も保存
    npz_path = joinpath(output_dir, "julia_300steps.npz")

    npzwrite(npz_path, Dict(
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
        # CGMパラメータ
        "window_size" => window_size,
        "overlap" => overlap,
        "cgm_iteration" => cgm_iteration,
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

    file_size_q = filesize(q_path) / 1024 / 1024
    file_size_T = filesize(T_path) / 1024 / 1024
    file_size_npz = filesize(npz_path) / 1024 / 1024

    println("  保存完了:")
    @printf("    - %s (%.2f MB)\n", q_path, file_size_q)
    @printf("    - %s (%.2f MB)\n", T_path, file_size_T)
    @printf("    - %s (%.2f MB)\n", npz_path, file_size_npz)

    # ========================================
    # サマリー表示
    # ========================================
    println("\n" * "="^80)
    println("実行完了サマリー")
    println("="^80)
    @printf("総実行時間: %.2f秒\n", elapsed_total)
    @printf("  CGM計算:  %.2f秒 (%.1f%%)\n", elapsed_cgm, elapsed_cgm/elapsed_total*100)
    @printf("  検証計算: %.2f秒 (%.1f%%)\n", elapsed_verify, elapsed_verify/elapsed_total*100)
    println("\n結果統計:")
    @printf("  熱流束範囲: %.4e ~ %.4e W/m²\n", minimum(q_result), maximum(q_result))
    @printf("  温度誤差（RMS）: %.4e K\n", rms_error)
    @printf("  温度誤差（最大）: %.4e K\n", max_error)
    println("\n次のステップ:")
    println("  結果比較: julia julia/examples/compare_python_julia_300steps.jl")
    println("="^80)

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
