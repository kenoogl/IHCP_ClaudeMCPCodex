#!/usr/bin/env julia
"""
ベースライン性能測定スクリプト

Phase 1.1 並列化導入前のベースライン性能を測定します。

測定項目:
1. build_dhcp_system! - 係数構築（DHCP）
2. assemble_dhcp_matrix - 疎行列組み立て（DHCP）
3. solve_dhcp! - DHCP順解析（単一ステップ）
4. build_adjoint_system! - 係数構築（Adjoint）
5. assemble_adjoint_matrix - 疎行列組み立て（Adjoint）
6. solve_adjoint! - Adjoint随伴解析（単一ステップ）
7. 完全一致検証（小規模問題）

実行方法:
  julia --project=. benchmarks/baseline_benchmark.jl
"""

using BenchmarkTools
using Printf
using Statistics
using Dates
using LinearAlgebra

# プロジェクトのソースコードをロード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

# 結果保存先
const RESULTS_DIR = joinpath(@__DIR__, "../benchmarks/results")

"""
    setup_test_problem(size::Symbol) -> Dict

テスト問題のセットアップ

# 引数
- size: :small, :medium, :large
"""
function setup_test_problem(size::Symbol)
    if size == :small
        ni, nj, nk = 10, 10, 10
        nt = 50
    elseif size == :medium
        ni, nj, nk = 40, 50, 15
        nt = 100
    elseif size == :large
        ni, nj, nk = 80, 100, 20
        nt = 357
    else
        error("不明な問題サイズ: $size")
    end

    # パラメータ
    dx = 0.12e-3
    dy = 0.12e-3
    dt = 1e-3
    rho = 7900.0

    # z方向格子（非均等）
    dz = vcat(
        fill(5.0e-5, 5),
        fill(1.0e-4, 5),
        fill(2.0e-4, nk - 10)
    )
    dz_t = vcat(fill(1.0e-4, nk - 1), [Inf])
    dz_b = vcat([Inf], fill(1.0e-4, nk - 1))

    # 熱物性値係数（SUS304）
    cp_coeffs = [420.0, 0.15, 0.0, 0.0]
    k_coeffs = [14.0, 0.01, 0.0, 0.0]

    # 初期温度場（300K）
    T_init = fill(300.0, ni, nj, nk)

    # 熱流束（一定値）
    q_surface = fill(1.0e5, nt - 1, ni, nj)

    # 観測温度（ダミーデータ）
    Y_obs = 300.0 .+ 50.0 .* rand(nt, ni, nj)

    return Dict(
        "ni" => ni, "nj" => nj, "nk" => nk, "nt" => nt,
        "dx" => dx, "dy" => dy, "dz" => dz, "dz_b" => dz_b, "dz_t" => dz_t,
        "dt" => dt, "rho" => rho,
        "cp_coeffs" => cp_coeffs, "k_coeffs" => k_coeffs,
        "T_init" => T_init, "q_surface" => q_surface, "Y_obs" => Y_obs
    )
end


"""
    benchmark_dhcp_coefficient_building(prob::Dict) -> BenchmarkTools.Trial

DHCP係数構築のベンチマーク
"""
function benchmark_dhcp_coefficient_building(prob::Dict)
    ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
    T_init = prob["T_init"]
    q_surf = prob["q_surface"][1, :, :]
    dx, dy, dz = prob["dx"], prob["dy"], prob["dz"]
    dz_b, dz_t = prob["dz_b"], prob["dz_t"]
    dt, rho = prob["dt"], prob["rho"]
    cp_coeffs, k_coeffs = prob["cp_coeffs"], prob["k_coeffs"]

    # 熱物性値計算
    cp, k = IHCP_CGM.thermal_properties_calculator(T_init, cp_coeffs, k_coeffs)

    println("  測定中: build_dhcp_system! (ni=$ni, nj=$nj, nk=$nk, N=$(ni*nj*nk))...")

    trial = @benchmark IHCP_CGM.build_dhcp_system!(
        $T_init, $q_surf, $rho, $cp, $k, $dx, $dy, $dz, $dz_b, $dz_t, $dt
    ) samples=100 seconds=60

    return trial
end


"""
    benchmark_dhcp_matrix_assembly(prob::Dict) -> BenchmarkTools.Trial

DHCP疎行列組み立てのベンチマーク
"""
function benchmark_dhcp_matrix_assembly(prob::Dict)
    ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
    T_init = prob["T_init"]
    q_surf = prob["q_surface"][1, :, :]
    dx, dy, dz = prob["dx"], prob["dy"], prob["dz"]
    dz_b, dz_t = prob["dz_b"], prob["dz_t"]
    dt, rho = prob["dt"], prob["rho"]
    cp_coeffs, k_coeffs = prob["cp_coeffs"], prob["k_coeffs"]

    # 熱物性値計算
    cp, k = IHCP_CGM.thermal_properties_calculator(T_init, cp_coeffs, k_coeffs)

    # 係数構築
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = IHCP_CGM.build_dhcp_system!(
        T_init, q_surf, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    println("  測定中: assemble_dhcp_matrix (N=$(ni*nj*nk))...")

    trial = @benchmark IHCP_CGM.assemble_dhcp_matrix(
        $ni, $nj, $nk, $a_w, $a_e, $a_s, $a_n, $a_b, $a_t, $a_p
    ) samples=100 seconds=60

    return trial
end


"""
    benchmark_dhcp_single_step(prob::Dict) -> BenchmarkTools.Trial

DHCP単一ステップ求解のベンチマーク
"""
function benchmark_dhcp_single_step(prob::Dict)
    ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
    T_init = prob["T_init"]
    q_surface = prob["q_surface"][1:2, :, :]  # 2ステップ分
    dx, dy, dz = prob["dx"], prob["dy"], prob["dz"]
    dz_b, dz_t = prob["dz_b"], prob["dz_t"]
    dt, rho = prob["dt"], prob["rho"]
    cp_coeffs, k_coeffs = prob["cp_coeffs"], prob["k_coeffs"]

    println("  測定中: solve_dhcp! (単一ステップ, N=$(ni*nj*nk))...")

    trial = @benchmark IHCP_CGM.solve_dhcp!(
        $T_init, $q_surface, 2, $rho, $cp_coeffs, $k_coeffs,
        $dx, $dy, $dz, $dz_b, $dz_t, $dt
    ) samples=50 seconds=60

    return trial
end


"""
    benchmark_adjoint_coefficient_building(prob::Dict) -> BenchmarkTools.Trial

Adjoint係数構築のベンチマーク
"""
function benchmark_adjoint_coefficient_building(prob::Dict)
    ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
    T_init = prob["T_init"]
    T_cal_bottom = zeros(ni, nj)
    Y_obs = zeros(ni, nj)
    dx, dy, dz = prob["dx"], prob["dy"], prob["dz"]
    dz_b, dz_t = prob["dz_b"], prob["dz_t"]
    dt, rho = prob["dt"], prob["rho"]
    cp_coeffs, k_coeffs = prob["cp_coeffs"], prob["k_coeffs"]

    # 熱物性値計算
    cp, k = IHCP_CGM.thermal_properties_calculator(T_init, cp_coeffs, k_coeffs)

    println("  測定中: build_adjoint_system! (ni=$ni, nj=$nj, nk=$nk, N=$(ni*nj*nk))...")

    trial = @benchmark IHCP_CGM.build_adjoint_system!(
        $T_init, $T_cal_bottom, $Y_obs, $rho, $cp, $k, $dx, $dy, $dz, $dz_b, $dz_t, $dt
    ) samples=100 seconds=60

    return trial
end


"""
    benchmark_adjoint_single_step(prob::Dict) -> BenchmarkTools.Trial

Adjoint単一ステップ求解のベンチマーク
"""
function benchmark_adjoint_single_step(prob::Dict)
    ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
    T_all = zeros(2, ni, nj, nk)
    T_all[1, :, :, :] = prob["T_init"]
    T_all[2, :, :, :] = prob["T_init"] .+ 10.0
    residual = zeros(2, ni, nj)
    dx, dy, dz = prob["dx"], prob["dy"], prob["dz"]
    dz_b, dz_t = prob["dz_b"], prob["dz_t"]
    dt, rho = prob["dt"], prob["rho"]
    cp_coeffs, k_coeffs = prob["cp_coeffs"], prob["k_coeffs"]

    println("  測定中: solve_adjoint! (単一ステップ, N=$(ni*nj*nk))...")

    trial = @benchmark IHCP_CGM.solve_adjoint!(
        $T_all, $residual, 2, $rho, $cp_coeffs, $k_coeffs,
        $dx, $dy, $dz, $dz_b, $dz_t, $dt
    ) samples=50 seconds=60

    return trial
end


"""
    run_small_exact_match_benchmark() -> Dict

小規模問題での完全一致検証ベンチマーク
"""
function run_small_exact_match_benchmark()
    println("\n" * "="^60)
    println("小規模問題での完全一致検証ベンチマーク")
    println("="^60)

    # 小規模問題セットアップ
    ni, nj, nk = 20, 25, 15
    nt = 71
    window_size = 71
    overlap = 0

    dx = 0.12e-3
    dy = 0.12e-3
    dt = 1e-3
    rho = 7900.0

    dz = vcat(fill(5.0e-5, 5), fill(1.0e-4, 5), fill(2.0e-4, 5))
    dz_t = vcat(fill(1.0e-4, nk - 1), [Inf])
    dz_b = vcat([Inf], fill(1.0e-4, nk - 1))

    cp_coeffs = [420.0, 0.15, 0.0, 0.0]
    k_coeffs = [14.0, 0.01, 0.0, 0.0]

    T_init = fill(300.0, ni, nj, nk)
    Y_obs = 300.0 .+ 50.0 .* rand(nt, ni, nj)
    q_init_value = 0.0
    cgm_iteration = 15

    @printf("  問題サイズ: ni=%d, nj=%d, nk=%d, nt=%d\n", ni, nj, nk, nt)
    @printf("  CGMパラメータ: window_size=%d, iteration=%d\n", window_size, cgm_iteration)

    println("\n測定開始...")
    start_time = time()

    # CGM計算
    println("  [1/2] スライディングウィンドウCGM計算...")
    time_cgm_start = time()
    q_result, windows_info = IHCP_CGM.solve_sliding_window_cgm(
        Y_obs, T_init,
        dx, dy, dz, dz_b, dz_t, dt,
        rho, cp_coeffs, k_coeffs,
        window_size, overlap, q_init_value, cgm_iteration
    )
    time_cgm = time() - time_cgm_start

    # DHCP検証
    println("  [2/2] DHCP検証計算...")
    time_verify_start = time()
    T_verify = IHCP_CGM.solve_dhcp_multiple_timesteps(
        T_init, q_result, nt, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, dt
    )
    time_verify = time() - time_verify_start

    elapsed_total = time() - start_time

    # 結果サマリー
    println("\n結果:")
    @printf("  CGM計算時間:   %.2f秒\n", time_cgm)
    @printf("  DHCP検証時間:  %.2f秒\n", time_verify)
    @printf("  総実行時間:    %.2f秒\n", elapsed_total)
    @printf("  熱流束範囲:    %.2e ~ %.2e W/m²\n", minimum(q_result), maximum(q_result))

    return Dict(
        "time_cgm" => time_cgm,
        "time_verify" => time_verify,
        "time_total" => elapsed_total,
        "q_min" => minimum(q_result),
        "q_max" => maximum(q_result),
        "problem_size" => (ni, nj, nk, nt)
    )
end


"""
    print_trial_summary(name::String, trial::BenchmarkTools.Trial, N::Int)

ベンチマーク結果サマリーを表示
"""
function print_trial_summary(name::String, trial::BenchmarkTools.Trial, N::Int)
    min_time = minimum(trial.times) / 1e9  # ns → s
    median_time = median(trial.times) / 1e9
    mean_time = mean(trial.times) / 1e9
    max_time = maximum(trial.times) / 1e9

    println("\n  結果: $name")
    @printf("    最小時間: %.6f 秒\n", min_time)
    @printf("    中央値:   %.6f 秒\n", median_time)
    @printf("    平均値:   %.6f 秒\n", mean_time)
    @printf("    最大時間: %.6f 秒\n", max_time)
    @printf("    サンプル数: %d\n", length(trial.times))
    if N > 0
        @printf("    スループット: %.2e 格子点/秒\n", N / median_time)
    end
end


"""
    save_benchmark_results(results::Dict, filename::String)

ベンチマーク結果をファイルに保存
"""
function save_benchmark_results(results::Dict, filename::String)
    mkpath(RESULTS_DIR)
    filepath = joinpath(RESULTS_DIR, filename)

    open(filepath, "w") do io
        println(io, "="^60)
        println(io, "ベースライン性能測定結果")
        println(io, "="^60)
        println(io, "測定日時: $(now())")
        println(io, "Julia バージョン: $(VERSION)")
        println(io, "スレッド数: $(Threads.nthreads())")
        println(io, "BLASスレッド数: $(BLAS.get_num_threads())")
        println(io, "="^60)

        for (key, value) in results
            if value isa BenchmarkTools.Trial
                min_time = minimum(value.times) / 1e9
                median_time = median(value.times) / 1e9
                mean_time = mean(value.times) / 1e9

                println(io, "\n[$key]")
                @printf(io, "  最小時間: %.6f 秒\n", min_time)
                @printf(io, "  中央値:   %.6f 秒\n", median_time)
                @printf(io, "  平均値:   %.6f 秒\n", mean_time)
                @printf(io, "  サンプル数: %d\n", length(value.times))
            elseif value isa Dict
                println(io, "\n[$key]")
                for (k, v) in value
                    println(io, "  $k: $v")
                end
            else
                println(io, "\n[$key]: $value")
            end
        end

        println(io, "\n" * "="^60)
    end

    @info "結果を保存しました: $filepath"
end


"""
    main()

メイン実行関数
"""
function main()
    println("="^60)
    println("ベースライン性能測定")
    println("="^60)
    println("測定日時: $(now())")
    println("Julia バージョン: $(VERSION)")
    println("スレッド数: $(Threads.nthreads())")
    println("BLASスレッド数: $(BLAS.get_num_threads())")
    println("="^60)

    results = Dict{String, Any}()

    # 問題サイズ選択
    println("\n測定する問題サイズ:")
    println("  1. Small  (10×10×10)")
    println("  2. Medium (40×50×15)")
    println("  3. Large  (80×100×20, 完全一致検証サイズ)")
    println("  4. すべて")

    print("\n選択 [1-4, デフォルト: 1]: ")
    choice_input = readline()
    choice = isempty(choice_input) ? 1 : parse(Int, choice_input)

    sizes = if choice == 1
        [:small]
    elseif choice == 2
        [:medium]
    elseif choice == 3
        [:large]
    elseif choice == 4
        [:small, :medium, :large]
    else
        error("不正な選択: $choice")
    end

    # 各問題サイズでベンチマーク実行
    for size in sizes
        println("\n" * "="^60)
        println("問題サイズ: $size")
        println("="^60)

        prob = setup_test_problem(size)
        ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
        N = ni * nj * nk

        @printf("格子点数: %d×%d×%d (N=%d)\n", ni, nj, nk, N)

        # 1. DHCP係数構築
        println("\n[1/6] DHCP係数構築ベンチマーク")
        trial_dhcp_coeff = benchmark_dhcp_coefficient_building(prob)
        print_trial_summary("build_dhcp_system!", trial_dhcp_coeff, N)
        results["$(size)_dhcp_coefficient"] = trial_dhcp_coeff

        # 2. DHCP疎行列組み立て
        println("\n[2/6] DHCP疎行列組み立てベンチマーク")
        trial_dhcp_matrix = benchmark_dhcp_matrix_assembly(prob)
        print_trial_summary("assemble_dhcp_matrix", trial_dhcp_matrix, 0)
        results["$(size)_dhcp_matrix"] = trial_dhcp_matrix

        # 3. DHCP単一ステップ
        println("\n[3/6] DHCP単一ステップベンチマーク")
        trial_dhcp_step = benchmark_dhcp_single_step(prob)
        print_trial_summary("solve_dhcp! (1ステップ)", trial_dhcp_step, 0)
        results["$(size)_dhcp_step"] = trial_dhcp_step

        # 4. Adjoint係数構築
        println("\n[4/6] Adjoint係数構築ベンチマーク")
        trial_adj_coeff = benchmark_adjoint_coefficient_building(prob)
        print_trial_summary("build_adjoint_system!", trial_adj_coeff, N)
        results["$(size)_adjoint_coefficient"] = trial_adj_coeff

        # 5. Adjoint単一ステップ
        println("\n[5/6] Adjoint単一ステップベンチマーク")
        trial_adj_step = benchmark_adjoint_single_step(prob)
        print_trial_summary("solve_adjoint! (1ステップ)", trial_adj_step, 0)
        results["$(size)_adjoint_step"] = trial_adj_step

        println()
    end

    # 小規模完全一致検証ベンチマーク
    if :small in sizes || choice == 4
        exact_match_results = run_small_exact_match_benchmark()
        results["exact_match_small"] = exact_match_results
    end

    # 結果保存
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    filename = "baseline_benchmark_$(timestamp).txt"
    save_benchmark_results(results, filename)

    println("\n" * "="^60)
    println("✓ ベースライン性能測定完了")
    println("="^60)

    return results
end


# スクリプトとして実行された場合
if abspath(PROGRAM_FILE) == @__FILE__
    try
        results = main()
    catch e
        println("\n❌ エラーが発生しました: ", e)
        Base.showerror(stdout, e, catch_backtrace())
        println()
        exit(1)
    end
end
