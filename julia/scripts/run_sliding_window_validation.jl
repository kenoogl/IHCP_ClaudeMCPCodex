#!/usr/bin/env julia
"""
Juliaスライディングウィンドウ検証スクリプト

Juliaコードでスライディングウィンドウ計算を実行し、
詳細な診断情報を含む結果を保存する。

実行方法:
  julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
  julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 20000

出力:
  shared/results/julia_sliding_window_cgm{N}.npz
"""

using Dates
using Printf
using Statistics
using NPZ
using ArgParse

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))
const DT = 1.0e-3

# ---------------------------------------------------------------------------
# コマンドライン引数パース
# ---------------------------------------------------------------------------

function parse_commandline()
    s = ArgParseSettings(
        description = "Juliaスライディングウィンドウ検証スクリプト",
        version = "1.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--cgm-iter"
            help = "CGM最大反復数"
            arg_type = Int
            default = 3
        "--nt"
            help = "時間ステップ数"
            arg_type = Int
            default = 300
        "--window"
            help = "ウィンドウサイズ"
            arg_type = Int
            default = 71
        "--overlap"
            help = "オーバーラップ"
            arg_type = Int
            default = 17
        "--quiet"
            help = "詳細出力を抑制"
            action = :store_true
    end

    return parse_args(s)
end

# ---------------------------------------------------------------------------
# ユーティリティ関数
# ---------------------------------------------------------------------------

function build_z_grid(nk::Int, Lz::Float64, stretch_factor::Float64)
    """Z方向格子生成（Python版と同じ）"""
    z_faces_normalized = range(1.0, stop = 0.0, length = nk + 1)
    z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

    dz = diff(z_faces)

    ΔZ = zeros(nk + 2)
    ΔZ[2:nk+1] = dz[1:nk]
    ΔZ[1] = ΔZ[2]
    ΔZ[nk+2] = ΔZ[nk+1]

    z_centers = zeros(nk + 2)
    z_centers[1] = z_faces[1]
    z_centers[nk+2] = z_faces[nk+1]

    for k in 2:nk+1
        z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
    end

    return z_centers, ΔZ, z_faces, dz
end

function load_material_properties()
    """材料物性値読み込み"""
    cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
    k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]
    rho = 7823.493962874829
    return rho, cp_coeffs, k_coeffs
end

function load_measurement_subset(nt::Int, verbose::Bool)
    """測定データ読み込み"""
    data_path = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
    if !isfile(data_path)
        error("Measurement file not found: $data_path")
    end

    if verbose
        println("データ読み込み中: $data_path")
    end
    flush(stdout)

    load_start = time()
    T_measure_full = npzread(data_path)
    ni_full = size(T_measure_full, 2)
    nj_full = size(T_measure_full, 3)

    if size(T_measure_full, 1) < nt
        error("Measurement file only has $(size(T_measure_full, 1)) time steps; need $nt")
    end

    subset = Array(T_measure_full[1:nt, :, :])
    if eltype(subset) != Float64
        subset = Float64.(subset)
    end
    T_measure_full = nothing
    GC.gc()

    if verbose
        println("データ形状 (Python順): $(size(subset))")
        println(@sprintf("読み込み時間: %.2f s", time() - load_start))
    end
    flush(stdout)

    return subset, ni_full, nj_full
end

function build_initial_temperature(first_frame::AbstractArray{<:Real,2}, nk::Int)
    """初期温度場構築"""
    base = Float64.(first_frame)
    T_init = repeat(base, 1, 1, nk)
    return Array{Float64,3}(T_init)
end

# ---------------------------------------------------------------------------
# 診断機能付きCGMソルバー
# ---------------------------------------------------------------------------

function solve_cgm_with_diagnostics!(
    T_init::Array{Float64,3},
    Y_obs::Array{Float64,3},
    q_init::Array{Float64,3},
    work::WorkBuffers,
    dx::Float64, dy::Float64,
    z_centers::Vector{Float64},
    ΔZ::Vector{Float64},
    dt::Float64,
    rho::Float64,
    cp_coeffs::Vector{Float64},
    k_coeffs::Vector{Float64},
    cgm_max_iter::Int;
    verbose::Bool = false
)
    """
    診断機能付きCGMソルバー

    既存のsolve_cgm!を使用しつつ、詳細な診断情報を記録する
    """
    ni, nj, nt = size(Y_obs)
    nk = size(T_init, 3)

    # dz配列の生成（ΔZからガイドセルを除く）
    dz = ΔZ[2:end-1]

    cgm_diagnostics = CGMIterDiag[]

    # CGMパラメータ設定
    # 注: 前処理なし（GS前処理は1反復あたり重すぎるため）
    #     反復数は増えるが、総時間は短くなる（Codex調査結果）
    cgm_params = (
        max_iter = cgm_max_iter,
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
        verbose = verbose,
        dhcp_extrapolation = :quadratic,
        adjoint_residual_scale = 0.5,
        dhcp_solver = :pcg,
        dhcp_smoother = :none,  # 前処理なし（PCG+noneが最速）
        adjoint_solver = :pcg,
        adjoint_smoother = :none,  # 前処理なし（PCG+noneが最速）
        sensitivity_solver = :pcg,  # PCGに統一（対称正定値問題）
        sensitivity_smoother = :none  # 前処理なし（PCG+noneが最速）
    )

    # 既存のsolve_cgm!を呼び出し
    # 注: 診断情報の詳細記録は今後の拡張として、まずは基本動作を確認
    q_result, T_result, J_history = solve_cgm!(
        T_init, Y_obs, q_init, work,
        dx, dy, z_centers, ΔZ,
        dt, rho, cp_coeffs, k_coeffs;
        params = cgm_params
    )

    # 簡易的な診断情報を作成（各CGM反復の情報）
    for it in 1:length(J_history)
        iter_diag = CGMIterDiag(it)
        iter_diag.J = J_history[it]
        # 注: 詳細な反復回数情報は今後の拡張で追加
        push!(cgm_diagnostics, iter_diag)
    end

    return q_result, T_result, J_history, cgm_diagnostics
end

# ---------------------------------------------------------------------------
# スライディングウィンドウCGM
# ---------------------------------------------------------------------------

function sliding_window_cgm_with_diagnostics(
    Y_obs::Array{Float64,3},
    T0::Array{Float64,3},
    dx::Float64, dy::Float64,
    z_centers::Vector{Float64},
    ΔZ::Vector{Float64},
    dt::Float64,
    rho::Float64,
    cp_coeffs::Vector{Float64},
    k_coeffs::Vector{Float64},
    window_size::Int,
    overlap::Int,
    q_init_value::Float64,
    cgm_max_iter::Int;
    verbose::Bool = true
)
    """
    診断機能付きスライディングウィンドウCGM計算
    """
    ni, nj, nt = size(Y_obs)
    nk = size(T0, 3)
    T_init = copy(T0)

    start_idx = 0
    q_total = Array{Float64,3}[]
    prev_q_win = nothing

    windows_info = []

    safety_counter = 0
    safety_limit = nt * 5

    # WorkBuffers作成
    work = WorkBuffers(ni + 2, nj + 2, nk + 2)

    if verbose
        println("\n" * "="^80)
        println("Julia スライディングウィンドウCGM計算")
        println("="^80)
        println("総時間ステップ: $nt")
        println("ウィンドウサイズ: $window_size")
        println("オーバーラップ: $overlap")
        println("CGM最大反復数: $cgm_max_iter")
        println("格子: $(ni)×$(nj)×$(nk)")
        println("="^80)
    end

    while start_idx < nt - 1
        safety_counter += 1
        if safety_counter > safety_limit
            println("[警告] 安全カウンタ超過")
            break
        end

        # 現在のウィンドウ長
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        # Python版: Y_obs[start_idx:end_idx+1] と同等 (Juliaは1-origin)
        Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]

        window_id = length(windows_info) + 1

        if verbose
            println("\n" * "="^80)
            println("ウィンドウ $window_id: [$start_idx, $end_idx] (長さ=$max_L)")
            println("="^80)
        end

        # 初期熱流束の設定
        if isnothing(prev_q_win)
            q_init_win = fill(q_init_value, (ni, nj, max_L))
            if verbose
                println("初期熱流束: 一定値 $q_init_value W/m²")
            end
        else
            q_init_win = zeros(Float64, ni, nj, max_L)
            L_overlap = min(overlap, max_L, size(prev_q_win, 3))

            if L_overlap > 0
                q_init_win[:, :, 1:L_overlap] = prev_q_win[:, :, end-L_overlap+1:end]
            end

            if L_overlap < max_L
                edge = prev_q_win[:, :, end]
                for t in L_overlap+1:max_L
                    q_init_win[:, :, t] = edge
                end
            end

            if verbose
                println("初期熱流束: 前ウィンドウから継承 (オーバーラップ=$L_overlap)")
            end
        end

        # CGM計算実行
        start_time_win = time()

        q_win, T_win_last, J_hist, cgm_diag = solve_cgm_with_diagnostics!(
            T_init, Y_obs_win, q_init_win, work,
            dx, dy, z_centers, ΔZ,
            dt, rho, cp_coeffs, k_coeffs, cgm_max_iter;
            verbose = verbose
        )

        elapsed_win = time() - start_time_win

        # ウィンドウ診断情報を保存
        win_info = Dict(
            "window_id" => window_id,
            "start_idx" => start_idx,
            "end_idx" => end_idx,
            "max_L" => max_L,
            "cgm_iterations" => length(J_hist),
            "J_final" => length(J_hist) > 0 ? J_hist[end] : 0.0,
            "J_history" => J_hist,
            "elapsed_window" => elapsed_win,
            "q_range" => [minimum(q_win), maximum(q_win)],
            "q_mean" => mean(q_win),
            "q_std" => std(q_win)
        )
        push!(windows_info, win_info)

        if verbose
            println("\nウィンドウ $window_id 完了:")
            println("  CGM反復数: $(length(J_hist))")
            println("  最終J: $(J_hist[end])")
            println(@sprintf("  実行時間: %.2f秒", elapsed_win))
            println(@sprintf("  q範囲: [%.3e, %.3e] W/m²", minimum(q_win), maximum(q_win)))
        end

        prev_q_win = copy(q_win)

        # 熱流束の連結（オーバーラップ部分は平均化）
        if length(q_total) == 0
            push!(q_total, q_win)
        else
            overlap_steps = min(overlap, size(q_win, 3), size(q_total[end], 3))
            if overlap_steps > 0
                # Pythonと同じ平均化
                q_total[end][:, :, end-overlap_steps+1:end] =
                    0.5 * q_total[end][:, :, end-overlap_steps+1:end] +
                    0.5 * q_win[:, :, 1:overlap_steps]

                if overlap_steps < size(q_win, 3)
                    push!(q_total, q_win[:, :, overlap_steps+1:end])
                end
            else
                push!(q_total, q_win)
            end
        end

        # 温度場の継承
        T_init = copy(T_win_last[:, :, :, end])

        # 次のウィンドウへ
        step = max(1, max_L - overlap)
        start_idx += step
    end

    # 全時間の熱流束を連結
    q_global = cat(q_total..., dims=3)

    if size(q_global, 3) > nt - 1
        q_global = q_global[:, :, 1:nt-1]
    end

    if verbose
        println("\n" * "="^80)
        println("スライディングウィンドウ計算完了")
        println("="^80)
        println("総ウィンドウ数: $(length(windows_info))")
        println("q_global形状: $(size(q_global))")
        println(@sprintf("q範囲: [%.3e, %.3e] W/m²", minimum(q_global), maximum(q_global)))
    end

    return q_global, T_init, windows_info
end

# ---------------------------------------------------------------------------
# メイン実行
# ---------------------------------------------------------------------------

function main()
    args = parse_commandline()

    println("="^80)
    println("Juliaスライディングウィンドウ検証スクリプト")
    println("="^80)
    println("CGM反復数: $(args["cgm-iter"])")
    println("時間ステップ数: $(args["nt"])")
    println("ウィンドウサイズ: $(args["window"])")
    println("オーバーラップ: $(args["overlap"])")
    println("="^80)

    verbose = !args["quiet"]

    # パラメータ設定
    nt = args["nt"]
    window_size = args["window"]
    overlap = args["overlap"]
    cgm_max_iter = args["cgm-iter"]

    ni, nj, nk = 80, 100, 20
    dx = 0.12e-3
    dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
    Lz = 0.5e-3
    stretch_factor = 3.0
    dt = 1.0e-3
    q_init_value = 0.0

    # 格子生成
    z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)
    rho, cp_coeffs, k_coeffs = load_material_properties()

    # データ読み込み
    Y_obs_python, ni_file, nj_file = load_measurement_subset(nt, verbose)
    if ni_file != ni || nj_file != nj
        error("Measurement grid mismatch: file has $(ni_file)x$(nj_file), expected $(ni)x$(nj)")
    end

    # メモリレイアウト変換: (nt, ni, nj) → (ni, nj, nt)
    Y_obs = permutedims(Y_obs_python, (2, 3, 1))
    Y_obs_python = nothing

    # 初期温度場
    T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)

    if verbose
        println("\n初期温度場: $(size(T_init))")
        println(@sprintf("温度範囲: [%.2f, %.2f] K", minimum(T_init), maximum(T_init)))
    end

    # 計算実行
    println("\n" * "="^80)
    println("計算開始")
    println("="^80)

    start_time = time()

    q_result, T_final, windows_info = sliding_window_cgm_with_diagnostics(
        Y_obs, T_init, dx, dy, z_centers, ΔZ, dt,
        rho, cp_coeffs, k_coeffs,
        window_size, overlap, q_init_value, cgm_max_iter;
        verbose = verbose
    )

    elapsed_total = time() - start_time

    println("\n" * "="^80)
    println("計算完了")
    println("="^80)
    println(@sprintf("総実行時間: %.2f秒", elapsed_total))
    println(@sprintf("1ウィンドウあたり: %.2f秒", elapsed_total / length(windows_info)))

    # 結果保存
    output_path = joinpath(PROJECT_ROOT, "shared", "results",
                           "julia_sliding_window_cgm$(cgm_max_iter).npz")
    mkpath(dirname(output_path))

    println("\n結果保存中: $output_path")

    # Python互換形式で保存（メモリレイアウトを戻す）
    q_global_python = permutedims(q_result, (3, 1, 2))  # (ni, nj, nt-1) → (nt-1, ni, nj)

    npzwrite(
        output_path,
        Dict(
            "q_global" => q_global_python,
            "T_final" => T_final,
            "elapsed_total" => elapsed_total,
            "windows_info" => windows_info,
            # パラメータ
            "nt" => nt,
            "window_size" => window_size,
            "overlap" => overlap,
            "cgm_max_iter" => cgm_max_iter,
            "dt" => dt,
            "dx" => dx,
            "dy" => dy,
            "dz" => dz,
            "z_faces" => z_faces,
            "z_centers" => z_centers,
            "ΔZ" => ΔZ,
            "rho" => rho,
            "cp_coeffs" => cp_coeffs,
            "k_coeffs" => k_coeffs
        )
    )

    println("✅ 保存完了")
    println("\n結果サマリー:")
    println("  q_global: $(size(q_global_python))")
    println(@sprintf("  q範囲: [%.3e, %.3e] W/m²", minimum(q_result), maximum(q_result)))
    println(@sprintf("  q平均: %.3e W/m²", mean(q_result)))
    println(@sprintf("  q標準偏差: %.3e W/m²", std(q_result)))
    println("  総ウィンドウ数: $(length(windows_info))")
    println("  総CGM反復数: $(sum(w["cgm_iterations"] for w in windows_info))")

    # ウィンドウ別サマリー
    println("\nウィンドウ別サマリー:")
    println(@sprintf("%3s %15s %4s %8s", "ID", "範囲", "CGM", "時間[s]"))
    println("-"^40)
    for win in windows_info
        println(@sprintf("%3d [%3d,%3d] %4d %8.2f",
                        win["window_id"],
                        win["start_idx"],
                        win["end_idx"],
                        win["cgm_iterations"],
                        win["elapsed_window"]))
    end

    println("\n" * "="^80)
    println("完了")
    println("="^80)

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        exit(main())
    catch err
        println("\n❌ Error: ", err)
        Base.showerror(stdout, err, catch_backtrace())
        println()
        exit(1)
    end
end
