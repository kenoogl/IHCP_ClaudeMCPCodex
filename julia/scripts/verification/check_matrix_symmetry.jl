#!/usr/bin/env julia
"""Julia版の係数行列が対称かどうかを確認"""

include("../src/IHCP_CGM.jl")
using .IHCP_CGM
using LinearAlgebra
using Printf

function check_matrix_symmetry()
    println("="^80)
    println("Julia版係数行列の対称性チェック")
    println("="^80)

    # 小さなテスト問題
    ni, nj, nk = 10, 10, 5
    N = ni * nj * nk

    println("行列サイズ: $(N)x$(N) (ni=$(ni), nj=$(nj), nk=$(nk))")

    # グリッド設定
    dx = 0.12e-3
    dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
    Lz = 0.5e-3
    stretch_factor = 3.0

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
    T_init = ones(Float64, ni, nj, nk) .* 500.0

    # 係数行列を明示的に構築（マトリクスフリーから変換）
    dt = 1.0e-3

    # WorkBuffersを準備
    work = WorkBuffers(ni+2, nj+2, nk+2)

    # 初期化
    work.θ[2:end-1, 2:end-1, 2:end-1] .= T_init
    work.b .= 0.0

    # 熱物性値計算
    T_with_guard = work.θ
    calculate_thermal_properties!(work.cp, work.λ, T_with_guard, cp_coeffs, k_coeffs)

    # HT（境界条件、ここではゼロ）
    HT = zeros(Float64, 6)

    # 係数行列を明示的に構築
    println("\n係数行列を明示的に構築中...")

    A = zeros(Float64, N, N)

    Δh = (dx, dy, 0.0)  # dzは使わない
    z_range = [1, nk]

    # 各格子点について行列要素を計算
    for linear_idx in 1:N
        # 線形インデックスを3Dインデックスに変換
        k = div(linear_idx - 1, ni * nj) + 1
        remainder = mod(linear_idx - 1, ni * nj)
        j = div(remainder, ni) + 1
        i = mod(remainder, ni) + 1

        # ガイドセル込みのインデックス
        ig = i + 1
        jg = j + 1
        kg = k + 1

        λ0 = work.λ[ig, jg, kg]
        r_ρc = 1.0 / (rho * work.cp[ig, jg, kg])
        m0 = work.mask[ig, jg, kg]

        mw = 1.0 - work.mask[ig-1, jg, kg]
        me = 1.0 - work.mask[ig+1, jg, kg]
        ms = 1.0 - work.mask[ig, jg-1, kg]
        mn = 1.0 - work.mask[ig, jg+1, kg]
        mb = work.mask[ig, jg, kg-1]
        mt = work.mask[ig, jg, kg+1]

        dx2 = 1.0 / (dx * dx)
        dy2 = 1.0 / (dy * dy)
        dx1 = 1.0 / dx
        dy1 = 1.0 / dy
        ddt = 1.0 / dt

        # 界面での熱伝導率（調和平均）
        λf = (a, b, ma, mb) -> begin
            if ma ≈ 0.0 || mb ≈ 0.0
                return 0.0
            else
                return 2.0 * a * b / (a + b + 1e-30)
            end
        end

        axm = λf(work.λ[ig-1,jg,kg], λ0, work.mask[ig-1,jg,kg], m0) * dx2 + mw*dx1*HT[1]*r_ρc
        axp = λf(work.λ[ig+1,jg,kg], λ0, work.mask[ig+1,jg,kg], m0) * dx2 + me*dx1*HT[2]*r_ρc
        aym = λf(work.λ[ig,jg-1,kg], λ0, work.mask[ig,jg-1,kg], m0) * dy2 + ms*dy1*HT[3]*r_ρc
        ayp = λf(work.λ[ig,jg+1,kg], λ0, work.mask[ig,jg+1,kg], m0) * dy2 + mn*dy1*HT[4]*r_ρc

        zb = (z_centers_grid[kg]-z_centers_grid[kg-1])*mb + (1.0-mb)*dz_grid[kg]
        zt = (z_centers_grid[kg+1]-z_centers_grid[kg])*m0 + (1.0-m0)*dz_grid[kg]

        azm = (work.λ[ig,jg,kg-1]/ (dz_grid[kg]*zb))*mb + (1.0-mb)/dz_grid[kg]*HT[5]*r_ρc
        azp = (λ0        / (dz_grid[kg]*zt))*mt + (1.0-mt)/dz_grid[kg]*HT[6]*r_ρc

        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm + ddt)*m0

        # 対角要素
        A[linear_idx, linear_idx] = -dd * m0

        # 隣接要素
        if i > 1
            neighbor_idx = linear_idx - 1
            A[linear_idx, neighbor_idx] = axm * m0
        end

        if i < ni
            neighbor_idx = linear_idx + 1
            A[linear_idx, neighbor_idx] = axp * m0
        end

        if j > 1
            neighbor_idx = linear_idx - ni
            A[linear_idx, neighbor_idx] = aym * m0
        end

        if j < nj
            neighbor_idx = linear_idx + ni
            A[linear_idx, neighbor_idx] = ayp * m0
        end

        if k > 1
            neighbor_idx = linear_idx - ni * nj
            A[linear_idx, neighbor_idx] = azm * m0
        end

        if k < nk
            neighbor_idx = linear_idx + ni * nj
            A[linear_idx, neighbor_idx] = azp * m0
        end
    end

    # 対称性チェック
    println("\n対称性を確認中...")

    A_T = transpose(A)
    diff_matrix = A - A_T

    frobenius_norm_diff = norm(diff_matrix, 2)
    frobenius_norm_A = norm(A, 2)
    relative_asymmetry = frobenius_norm_diff / frobenius_norm_A

    max_diff = maximum(abs.(diff_matrix))
    max_A = maximum(abs.(A))

    println()
    @printf("||A - A^T||_F: %.6e\n", frobenius_norm_diff)
    @printf("||A||_F:       %.6e\n", frobenius_norm_A)
    @printf("相対非対称度:   %.6e\n", relative_asymmetry)
    println()
    @printf("max|A - A^T|:  %.6e\n", max_diff)
    @printf("max|A|:        %.6e\n", max_A)
    println()

    # 判定
    if relative_asymmetry < 1e-14
        println("✅ 判定: 係数行列は**対称**です（数値誤差範囲内）")
        println("   → CG法の使用が適切")
    elseif relative_asymmetry < 1e-10
        println("⚠️  判定: 係数行列はほぼ対称ですが、わずかに非対称です")
        println("   → CG法も使用可能だが、BICGSTABの方が理論的に正確")
    else
        println("❌ 判定: 係数行列は**非対称**です")
        println("   → BICGSTABの使用は適切")
    end

    println()
    println("="^80)

    return relative_asymmetry
end

asymmetry = check_matrix_symmetry()

if asymmetry > 1e-10
    println("\n⚠️  注意: Julia版はBICGSTABを使用していますが、行列は非対称です")
    println("   これは非均一格子による正しい動作です")
else
    println("\n✅ Julia版の係数行列も対称です")
    println("   → CG法を使用すれば高速化できる可能性があります")
end
