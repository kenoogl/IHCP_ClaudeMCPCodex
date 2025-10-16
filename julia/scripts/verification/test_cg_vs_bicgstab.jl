#!/usr/bin/env julia
"""
CG法とBICGSTAB法の反復回数を比較
"""

using IterativeSolvers
using LinearAlgebra
using Printf

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

function test_cg_vs_bicgstab()
    println("="^80)
    println("CG vs BICGSTAB 反復回数比較")
    println("="^80)

    # 小さなテスト問題
    ni, nj, nk = 20, 25, 10
    N = ni * nj * nk

    println("格子サイズ: $(ni)×$(nj)×$(nk) = $(N)")

    # グリッド設定
    dx = 0.12e-3
    dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
    Lz = 0.5e-3
    stretch_factor = 3.0
    dt = 1.0e-3

    # Z方向格子
    z_faces_normalized = range(1.0, stop=0.0, length=nk+1)
    z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)
    dz = Base.diff(z_faces)
    z_centers = similar(dz)
    z_centers[1] = z_faces[1]
    z_centers[end] = z_faces[end]
    if nk > 2
        z_centers[2:end-1] = (z_faces[2:end-2] .+ z_faces[3:end-1]) ./ 2
    end

    dz_top = zeros(Float64, nk)
    dz_top[end] = Inf
    if nk > 1
        dz_top[1:end-1] = z_centers[2:end] .- z_centers[1:end-1]
    end

    dz_bottom = zeros(Float64, nk)
    dz_bottom[1] = Inf
    if nk > 1
        dz_bottom[2:end] = z_centers[2:end] .- z_centers[1:end-1]
    end

    z_centers_grid, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

    # 熱物性値
    rho = 7823.493962874829
    cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
    k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

    # 初期温度場
    T_init = ones(Float64, ni, nj, nk) .* 550.0
    q_surface = zeros(Float64, ni, nj)

    # WorkBuffers
    work = WorkBuffers(ni+2, nj+2, nk+2)
    work.θ[2:end-1, 2:end-1, 2:end-1] .= T_init

    # 熱物性値計算
    calculate_thermal_properties!(work.cp, work.λ, work.θ, cp_coeffs, k_coeffs)

    # RHSベクトル構築（簡略版）
    work.b[2:end-1, 2:end-1, 2:end-1] .= T_init .* (rho .* work.cp[2:end-1, 2:end-1, 2:end-1] ./ dt)

    Δh = (dx, dy, 0.0)
    z_range = [1, nk]
    HT = zeros(Float64, 6)

    # マスク初期化
    work.mask .= 1.0
    work.mask[1, :, :] .= 0.0
    work.mask[end, :, :] .= 0.0
    work.mask[:, 1, :] .= 0.0
    work.mask[:, end, :] .= 0.0
    work.mask[:, :, 1] .= 0.0
    work.mask[:, :, end] .= 0.0

    println("\n1. BICGSTAB法（現在のJulia実装、前処理: Jacobi）")
    work_copy = deepcopy(work)
    t_start = time()
    success, iters, res0 = PBiCGSTAB!(work_copy, Δh, dt, z_centers_grid, dz_grid, z_range, HT, rho,
                                      tol=1e-6, maxItr=5000, smoother=:jacobi, par="sequential")
    t_bicgstab = time() - t_start
    @printf("  反復回数: %d\n", iters)
    @printf("  収束: %s\n", success ? "成功" : "失敗")
    @printf("  時間: %.3f秒\n", t_bicgstab)
    @printf("  初期残差: %.6e\n", res0)

    println("\n2. IterativeSolvers.jl の cg!（対称行列用）")
    println("  注意: 行列が対称でない場合、CG法は失敗する可能性があります")

    # TODO: CG法のテストには明示的な行列が必要
    # マトリクスフリー実装からは直接比較が難しい
    println("  → マトリクスフリー実装のため、直接比較は困難")
    println("  → 係数行列の対称性をまず確認する必要があります")

    println("\n" * "="^80)
    println("結論:")
    println("  BICGSTAB反復回数: $(iters)回")
    println("  係数行列が対称なら、CG法で約2-3倍高速化が期待できます")
    println("="^80)
end

test_cg_vs_bicgstab()
