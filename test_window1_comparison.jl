#!/usr/bin/env julia
"""
ウィンドウ1のJulia計算テスト

Python版と全く同じ入力でウィンドウ1を計算し、結果を比較する。
"""

using NPZ
using Printf
using Statistics

include("julia/src/IHCP_CGM.jl")
using .IHCP_CGM

const PROJECT_ROOT = @__DIR__
const DT = 1.0e-3

function build_z_grid(nk::Int, Lz::Float64, stretch_factor::Float64)
    """Z方向格子生成"""
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

function main()
    println("="^80)
    println("ウィンドウ1のJulia計算テスト")
    println("="^80)

    # パラメータ設定
    ni, nj, nk = 80, 100, 20
    dx = 0.12e-3
    dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
    Lz = 0.5e-3
    stretch_factor = 3.0
    dt = 1.0e-3

    # 格子生成
    z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)

    # 熱物性値
    cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
    k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]
    rho = 7823.493962874829

    # データ読み込み
    data_path = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
    Y_obs_full = npzread(data_path)
    println("データ形状: $(size(Y_obs_full))")

    # ウィンドウ1の設定
    start_idx = 0
    end_idx = 70
    max_L = end_idx - start_idx  # 70ステップの熱流束

    # メモリレイアウト変換: (nt, ni, nj) → (ni, nj, nt)
    Y_obs_full_julia = permutedims(Y_obs_full, (2, 3, 1))

    # Y_obsのスライス（Julia: 1-indexedなので+1）
    # Python: Y_obs[0:71] = Y_obs[0], Y_obs[1], ..., Y_obs[70] = 71ステップ
    # Julia:  Y_obs[:, :, 1:71] = 同じく71ステップ
    Y_obs_win = Y_obs_full_julia[:, :, start_idx+1:end_idx+1]  # 71ステップ (1:71)
    println("Y_obs_win形状: $(size(Y_obs_win))")
    println(@sprintf("Y_obs_win範囲: [%.2f, %.2f] K", minimum(Y_obs_win), maximum(Y_obs_win)))

    # 初期温度場
    T_init = repeat(Y_obs_win[:, :, 1], 1, 1, nk)
    println("\nT_init形状: $(size(T_init))")
    println(@sprintf("T_init範囲: [%.2f, %.2f] K", minimum(T_init), maximum(T_init)))

    # 初期熱流束（全てゼロ）
    q_init = zeros(Float64, ni, nj, max_L)
    println("\nq_init形状: $(size(q_init))")
    println("q_init値: 全て 0.0 W/m²")

    # WorkBuffers作成
    work = WorkBuffers(ni + 2, nj + 2, nk + 2)

    # CGMパラメータは直接指定

    println("\n" * "="^80)
    println("CGM計算開始（3反復）")
    println("="^80 * "\n")

    # CGM計算（max_iter=3を指定）
    q_result, T_result, J_hist = solve_cgm!(
        T_init, Y_obs_win, q_init, work,
        dx, dy, z_centers, ΔZ,
        dt, rho, cp_coeffs, k_coeffs;
        params = (max_iter = 3, verbose = true)
    )

    println("\n" * "="^80)
    println("CGM計算完了")
    println("="^80)
    println("q_result形状: $(size(q_result))")
    println(@sprintf("q_result範囲: [%.6e, %.6e] W/m²", minimum(q_result), maximum(q_result)))
    println(@sprintf("q_result平均: %.6e W/m²", mean(q_result)))
    println(@sprintf("q_result標準偏差: %.6e W/m²", std(q_result)))
    println("\nT_result形状: $(size(T_result))")
    println(@sprintf("T_result範囲: [%.2f, %.2f] K", minimum(T_result), maximum(T_result)))
    println("\nJ履歴: $J_hist")

    # 結果保存（Python互換形式）
    output_path = joinpath(PROJECT_ROOT, "shared", "results", "window1_julia_test.npz")

    # メモリレイアウトを戻す
    Y_obs_win_python = permutedims(Y_obs_win, (3, 1, 2))
    q_result_python = permutedims(q_result, (3, 1, 2))
    T_result_python = permutedims(T_result, (4, 1, 2, 3))

    npzwrite(
        output_path,
        Dict(
            "Y_obs_win" => Y_obs_win_python,
            "T_init" => T_init,
            "q_init" => permutedims(q_init, (3, 1, 2)),
            "q_result" => q_result_python,
            "T_result" => T_result_python,
            "J_hist" => J_hist
        )
    )
    println("\n結果保存: $output_path")
end

main()
