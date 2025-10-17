#!/usr/bin/env julia
"""
簡易性能測定スクリプト（現在のAPI対応版）

現在のAPIに合わせたベンチマーク測定

測定項目:
1. solve_cgm! - CGMソルバー（小・中・大規模問題）
2. 実行時間と反復回数
3. 以前の結果との比較

実行方法:
  julia --project=. benchmarks/simple_benchmark.jl
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
const RESULTS_DIR = joinpath(@__DIR__, "results")

"""
    build_z_grid_simple(nk::Int, Lz::Float64, stretch_factor::Float64)

z方向格子生成（run_10steps_fullsize_test.jlと同じ方式）
"""
function build_z_grid_simple(nk::Int, Lz::Float64, stretch_factor::Float64)
  # Step 1: z_faces生成
  z_faces_normalized = range(1.0, stop = 0.0, length = nk + 1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  # Step 2: dzの計算
  dz = diff(z_faces)

  # Step 3: ΔZ配列（ガイドセル込み、nk+2個）
  ΔZ = zeros(nk + 2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # Step 4: z_centers（ガイドセル込み、nk+2個）
  z_centers = zeros(nk + 2)
  z_centers[1] = z_faces[1]
  z_centers[nk+2] = z_faces[nk+1]
  z_centers[2] = z_faces[1]
  z_centers[nk+1] = z_faces[nk+1]

  if nk > 2
    for k in 3:nk
      z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
    end
  end

  return z_centers, ΔZ, z_faces, dz
end

"""
    setup_test_problem(size::Symbol) -> Dict

テスト問題のセットアップ

# 引数
- size: :small, :medium, :large
"""
function setup_test_problem(size::Symbol)
    if size == :small
        ni, nj, nk = 10, 10, 10
        nt = 10
    elseif size == :medium
        ni, nj, nk = 40, 50, 15
        nt = 10
    elseif size == :large
        ni, nj, nk = 80, 100, 20
        nt = 10
    else
        error("不明な問題サイズ: $size")
    end

    # パラメータ
    dx = 0.12e-3
    dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
    dt = 1e-3
    rho = 7823.493962874829

    # z方向格子（非均等、run_10steps方式）
    Lz = 0.5e-3
    stretch_factor = 3.0
    z_centers, ΔZ, z_faces, dz = build_z_grid_simple(nk, Lz, stretch_factor)

    # 熱物性値係数（SUS304、多項式）
    cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
    k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

    # 初期温度場（550K）
    T_init = fill(550.0, ni, nj, nk)

    # 観測温度（ダミーデータ、(ni, nj, nt)形式）
    Y_obs = 550.0 .+ 20.0 .* rand(ni, nj, nt)

    # 初期熱流束（(ni, nj, nt-1)形式）
    q_init = zeros(Float64, ni, nj, nt - 1)

    # WorkBuffers
    work = WorkBuffers(ni + 2, nj + 2, nk + 2)

    return Dict(
        "ni" => ni, "nj" => nj, "nk" => nk, "nt" => nt,
        "dx" => dx, "dy" => dy, "z_centers" => z_centers, "ΔZ" => ΔZ,
        "dt" => dt, "rho" => rho,
        "cp_coeffs" => cp_coeffs, "k_coeffs" => k_coeffs,
        "T_init" => T_init, "Y_obs" => Y_obs, "q_init" => q_init,
        "work" => work
    )
end


"""
    benchmark_cgm_solver(prob::Dict) -> Dict

CGMソルバーのベンチマーク
"""
function benchmark_cgm_solver(prob::Dict)
    ni, nj, nk, nt = prob["ni"], prob["nj"], prob["nk"], prob["nt"]
    T_init = prob["T_init"]
    Y_obs = prob["Y_obs"]
    q_init = prob["q_init"]
    work = prob["work"]
    dx, dy = prob["dx"], prob["dy"]
    z_centers, ΔZ = prob["z_centers"], prob["ΔZ"]
    dt, rho = prob["dt"], prob["rho"]
    cp_coeffs, k_coeffs = prob["cp_coeffs"], prob["k_coeffs"]

    println("  測定中: solve_cgm! (ni=$ni, nj=$nj, nk=$nk, nt=$nt)...")

    # CGMパラメータ（Phase 1-B最適化適用）
    cgm_params = (
        max_iter = 1,  # 性能測定のため1回
        sigma = 1.8,
        rtol_dhcp = 1.0e-6,
        rtol_adjoint = 1.0e-8,
        maxiter_cg = 20000,
        dire_reset_every = 5,
        eps = 1.0e-12,
        min_iter = 10,
        P = 10,
        eta = 1.0e-4,
        beta_max = 1.0e8,
        verbose = false,
        dhcp_extrapolation = :quadratic,
        adjoint_residual_scale = 0.5
    )

    # ウォームアップ実行（1回）
    println("    ウォームアップ実行...")
    q_result, T_result, J_history = solve_cgm!(
        copy(T_init), copy(Y_obs), copy(q_init), work,
        dx, dy, z_centers, ΔZ, dt, rho, cp_coeffs, k_coeffs;
        params = cgm_params
    )

    # 本測定（3回）
    println("    本測定中（3回）...")
    times = Float64[]

    for i in 1:3
        start_time = time()
        q_result, T_result, J_history = solve_cgm!(
            copy(T_init), copy(Y_obs), copy(q_init), work,
            dx, dy, z_centers, ΔZ, dt, rho, cp_coeffs, k_coeffs;
            params = cgm_params
        )
        elapsed = time() - start_time
        push!(times, elapsed)
        println(@sprintf("      試行%d: %.6f 秒", i, elapsed))
    end

    return Dict(
        "times" => times,
        "min" => minimum(times),
        "median" => median(times),
        "mean" => mean(times),
        "max" => maximum(times),
        "J_history_length" => length(J_history)
    )
end


"""
    print_benchmark_summary(name::String, result::Dict, N::Int)

ベンチマーク結果サマリーを表示
"""
function print_benchmark_summary(name::String, result::Dict, N::Int)
    println("\n  結果: $name")
    @printf("    最小時間: %.6f 秒\n", result["min"])
    @printf("    中央値:   %.6f 秒\n", result["median"])
    @printf("    平均値:   %.6f 秒\n", result["mean"])
    @printf("    最大時間: %.6f 秒\n", result["max"])
    @printf("    サンプル数: %d\n", length(result["times"]))
    if N > 0
        @printf("    スループット: %.2e 格子点/秒\n", N / result["median"])
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
        println(io, "簡易性能測定結果（現在のAPI対応版）")
        println(io, "="^60)
        println(io, "測定日時: $(now())")
        println(io, "Julia バージョン: $(VERSION)")
        println(io, "スレッド数: $(Threads.nthreads())")
        println(io, "BLASスレッド数: $(BLAS.get_num_threads())")
        println(io, "="^60)

        for (key, value) in results
            if value isa Dict
                println(io, "\n[$key]")
                for (k, v) in value
                    if k == "times"
                        println(io, "  $k: $(join([@sprintf("%.6f", t) for t in v], ", "))")
                    else
                        println(io, "  $k: $v")
                    end
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
    compare_with_previous(current_results::Dict)

以前の結果と比較
"""
function compare_with_previous(current_results::Dict)
    # 2025年10月2日のベースライン結果（Large問題）
    previous_results = Dict(
        "large_dhcp_coefficient" => 0.003746,  # 中央値（秒）
        "large_dhcp_step" => 0.324853,
        "large_adjoint_coefficient" => 0.003746,
        "large_adjoint_step" => 0.070252
    )

    # 22fde2dの10ステップフルサイズ結果
    previous_fullsize = Dict(
        "total_time" => 1072.43,  # 秒
        "cgm_time" => 1069.53,
        "dhcp_time" => 321.406,
        "adjoint_time" => 447.191,
        "sensitivity_time" => 300.695
    )

    println("\n" * "="^60)
    println("以前の結果との比較")
    println("="^60)

    println("\n【2025年10月2日 ベースライン】")
    println("  DHCP係数構築（Large）: 3.746 ms")
    println("  DHCP単一ステップ（Large）: 324.853 ms")
    println("  Adjoint係数構築（Large）: 3.746 ms")
    println("  Adjoint単一ステップ（Large）: 70.252 ms")

    println("\n【22fde2d 10ステップフルサイズ（80×100×20, 10ステップ）】")
    println("  総実行時間: 1072.43 秒（約17.9分）")
    println("  CGM時間: 1069.53 秒（99.7%）")
    println("  DHCP時間: 321.406 秒（30.0%）")
    println("  Adjoint時間: 447.191 秒（41.7%）")
    println("  Sensitivity時間: 300.695 秒（28.0%）")

    if haskey(current_results, "large_cgm")
        current_large = current_results["large_cgm"]
        println("\n【今回の測定（Large, 10ステップ）】")
        @printf("  CGM実行時間（1反復）: %.3f 秒\n", current_large["median"])
        @printf("  推定10反復時間: %.3f 秒（%.2f分）\n",
                current_large["median"] * 10,
                current_large["median"] * 10 / 60)

        # 22fde2dとの比較
        speedup = previous_fullsize["cgm_time"] / current_large["median"]
        println(@sprintf("  22fde2dとの速度比: %.2fx", speedup))

        improvement_pct = (1 - current_large["median"] / previous_fullsize["cgm_time"]) * 100
        println(@sprintf("  改善率: %.2f%%", improvement_pct))
    end
end


"""
    main()

メイン実行関数
"""
function main()
    # BLAS逐次実行に設定（ベースライン測定）
    BLAS.set_num_threads(1)

    println("="^60)
    println("簡易性能測定（現在のAPI対応版）")
    println("="^60)
    println("測定日時: $(now())")
    println("Julia バージョン: $(VERSION)")
    println("スレッド数: $(Threads.nthreads())")
    println("BLASスレッド数: $(BLAS.get_num_threads()) [逐次実行]")
    println("="^60)

    results = Dict{String, Any}()

    # 問題サイズ選択
    println("\n測定する問題サイズ:")
    println("  1. Small  (10×10×10, 10ステップ)")
    println("  2. Medium (40×50×15, 10ステップ)")
    println("  3. Large  (80×100×20, 10ステップ) [推奨]")
    println("  4. すべて")

    print("\n選択 [1-4, デフォルト: 3]: ")
    choice_input = readline()
    choice = isempty(choice_input) ? 3 : parse(Int, choice_input)

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
        ni, nj, nk, nt = prob["ni"], prob["nj"], prob["nk"], prob["nt"]
        N = ni * nj * nk

        @printf("格子点数: %d×%d×%d (N=%d), 時間ステップ: %d\n", ni, nj, nk, N, nt)

        # CGMソルバーベンチマーク
        println("\nCGMソルバーベンチマーク")
        result_cgm = benchmark_cgm_solver(prob)
        print_benchmark_summary("solve_cgm! ($(nt)ステップ, 1反復)", result_cgm, N)
        results["$(size)_cgm"] = result_cgm

        println()
    end

    # 結果保存
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    filename = "simple_benchmark_$(timestamp).txt"
    save_benchmark_results(results, filename)

    # 以前の結果と比較
    compare_with_previous(results)

    println("\n" * "="^60)
    println("✓ 簡易性能測定完了")
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
