"""
verify_matrix_symmetry.jl

CalcAX!が生成する行列の対称性を検証するテストスクリプト
"""

using LinearAlgebra
using SparseArrays
using Printf

# プロジェクトのモジュールをロード
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using IHCP_CGM
using IHCP_CGM.Commons
using IHCP_CGM.ThermalProperties

"""
CalcAX!から陽に行列Aを構築する関数

各格子点に対して単位ベクトルを入力し、CalcAX!の出力から行列の列を抽出する
"""
function extract_matrix_from_CalcAX!(
    ni::Int, nj::Int, nk::Int,
    Δh::NTuple{3,Float64},
    Δt::Float64,
    z_faces::Vector{Float64},
    dz_phys::Vector{Float64},
    ρ::Float64
)
    wk = WorkBuffers(ni, nj, nk)

    # generate_guard_cell_gridを使ってガードセル込み配列を生成
    nk_phys = nk - 2
    ZC, ΔZ = IHCP_CGM.generate_guard_cell_grid(nk_phys, dz_phys, z_faces)

    # 一様な熱物性値を設定（対称性確認のため）
    λ_uniform = 50.0  # W/m·K
    cp_uniform = 500.0  # J/kg·K
    fill!(wk.λ, λ_uniform)
    fill!(wk.cp, cp_uniform)

    # マスク設定: 境界セル（k=2, k=nk-1）はmask=0.0
    fill!(wk.mask, 1.0)  # デフォルトは内部セル
    for j in 1:nj, i in 1:ni
        wk.mask[i, j, 1] = 0.0      # 底面ガードセル
        wk.mask[i, j, 2] = 0.0      # 底面物理セル（境界）
        wk.mask[i, j, nk-1] = 0.0   # 表面物理セル（境界）
        wk.mask[i, j, nk] = 0.0     # 表面ガードセル
    end

    # 内部格子点数（境界セルを除く）
    ni_inner = ni - 2
    nj_inner = nj - 2
    # Z方向: ガードセル2個 + 境界セル2個 = 4個除く
    nk_inner = nk - 4
    N = ni_inner * nj_inner * nk_inner

    # 疎行列として構築
    I_idx = Int[]
    J_idx = Int[]
    V_val = Float64[]

    # 格子点インデックス → 線形インデックス変換
    # Z方向は境界セルを除くのでk-3でマッピング
    function grid_to_linear(i, j, k)
        return (i-2) + (j-2)*ni_inner + (k-3)*ni_inner*nj_inner + 1
    end

    # 各格子点で単位ベクトルをテスト（内部セルのみ）
    for k in 3:nk-2  # Z方向の内部セルのみ
        for j in 2:nj-1
            for i in 2:ni-1
                col_idx = grid_to_linear(i, j, k)

                # 単位ベクトルを設定
                fill!(wk.θ, 0.0)
                wk.θ[i, j, k] = 1.0

                # CalcAX!を実行（Ax = bのAx部分を計算、bは無視）
                ax = zeros(Float64, ni, nj, nk)
                IHCP_CGM.CommonSolver.CalcAX!(
                    ax, wk.θ, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, ZC, ΔZ, "sequential"
                )

                # 非ゼロ要素を抽出（内部セルのみ）
                for kk in 3:nk-2  # Z方向の内部セルのみ
                    for jj in 2:nj-1
                        for ii in 2:ni-1
                            val = ax[ii, jj, kk]
                            if abs(val) > 1e-14
                                row_idx = grid_to_linear(ii, jj, kk)
                                push!(I_idx, row_idx)
                                push!(J_idx, col_idx)
                                push!(V_val, val)
                            end
                        end
                    end
                end
            end
        end
    end

    # 疎行列を構築
    A = sparse(I_idx, J_idx, V_val, N, N)
    return A
end


"""
行列の対称性をチェック
"""
function check_symmetry(A::SparseMatrixCSC{Float64, Int})
    N = size(A, 1)

    # A - A^T を計算
    diff = A - A'

    # フロベニウスノルム
    norm_diff = norm(diff)
    norm_A = norm(A)
    rel_error = norm_diff / norm_A

    println("=" ^ 70)
    println("行列対称性検証")
    println("=" ^ 70)
    println("行列サイズ: $(N) × $(N)")
    println("非ゼロ要素数: $(nnz(A))")
    println("||A - A^T||_F: ", norm_diff)
    println("||A||_F: ", norm_A)
    println("相対誤差: ", rel_error)

    if rel_error < 1e-12
        println("✅ 行列は対称です（誤差 < 1e-12）")
        is_symmetric = true
    elseif rel_error < 1e-8
        println("⚠️  行列はほぼ対称です（1e-12 < 誤差 < 1e-8）")
        is_symmetric = true
    else
        println("❌ 行列は非対称です（誤差 > 1e-8）")
        is_symmetric = false
    end

    # 最大非対称要素を表示
    if norm_diff > 1e-14
        println("\n最大非対称要素:")
        for i in 1:min(10, N)
            for j in i+1:N
                diff_val = A[i,j] - A[j,i]
                if abs(diff_val) > 1e-10
                    println("  A[$i,$j] = ", A[i,j], ", A[$j,$i] = ", A[j,i],
                            ", 差分 = ", diff_val)
                end
            end
        end
    end
    println("=" ^ 70)

    return is_symmetric, rel_error
end


"""
メイン実行関数
"""
function main()
    println("\n🔍 CalcAX!行列の対称性検証\n")

    # テストケース1: 小規模格子、均等Z格子
    println("【テストケース1】小規模格子（5×5×5）、均等Z格子")
    ni, nj, nk = 5, 5, 5
    dx, dy = 0.12e-3, 0.12e-3
    dt = 1e-3
    ρ = 7900.0

    # 均等Z格子（物理セル数: nk-2）
    nk_phys = nk - 2
    z_faces_uniform = collect(range(0.0, 0.7e-3, length=nk_phys+1))
    dz_uniform = [z_faces_uniform[k+1] - z_faces_uniform[k] for k in 1:nk_phys]

    A1 = extract_matrix_from_CalcAX!(ni, nj, nk, (dx, dy, 0.0), dt, z_faces_uniform, dz_uniform, ρ)
    is_sym1, err1 = check_symmetry(A1)

    # テストケース2: 小規模格子、非均等Z格子（実際のプロジェクト設定）
    println("\n【テストケース2】小規模格子（5×5×5）、非均等Z格子")

    # Python版と同じZ格子生成（表面側に密集）
    nk_phys = nk - 2
    total_depth = 0.7e-3
    ratio = 1.5

    # 等比数列の和の公式から初項を計算
    dz_first = total_depth * (ratio - 1) / (ratio^nk_phys - 1)
    dz_nonuniform = [dz_first * ratio^(i-1) for i in 1:nk_phys]

    # 境界座標を計算
    z_faces_nonuniform = zeros(nk_phys+1)
    for i in 2:nk_phys+1
        z_faces_nonuniform[i] = z_faces_nonuniform[i-1] + dz_nonuniform[i-1]
    end

    A2 = extract_matrix_from_CalcAX!(ni, nj, nk, (dx, dy, 0.0), dt, z_faces_nonuniform, dz_nonuniform, ρ)
    is_sym2, err2 = check_symmetry(A2)

    # テストケース3: 中規模格子、非均等Z格子
    println("\n【テストケース3】中規模格子（10×10×10）、非均等Z格子")
    ni, nj, nk = 10, 10, 10

    nk_phys = nk - 2
    dz_first = total_depth * (ratio - 1) / (ratio^nk_phys - 1)
    dz_nonuniform = [dz_first * ratio^(i-1) for i in 1:nk_phys]

    # 境界座標を計算
    z_faces_nonuniform = zeros(nk_phys+1)
    for i in 2:nk_phys+1
        z_faces_nonuniform[i] = z_faces_nonuniform[i-1] + dz_nonuniform[i-1]
    end

    A3 = extract_matrix_from_CalcAX!(ni, nj, nk, (dx, dy, 0.0), dt, z_faces_nonuniform, dz_nonuniform, ρ)
    is_sym3, err3 = check_symmetry(A3)

    # 総合結果
    println("\n" * "=" ^ 70)
    println("総合結果")
    println("=" ^ 70)
    @printf("テストケース1（均等Z）:   %s (誤差: %.2e)\n", is_sym1 ? "✅ 対称" : "❌ 非対称", err1)
    @printf("テストケース2（非均等Z）:  %s (誤差: %.2e)\n", is_sym2 ? "✅ 対称" : "❌ 非対称", err2)
    @printf("テストケース3（非均等Z）:  %s (誤差: %.2e)\n", is_sym3 ? "✅ 対称" : "❌ 非対称", err3)
    println("=" ^ 70)

    if is_sym1 && is_sym2 && is_sym3
        println("\n✅ 結論: CalcAX!が生成する行列は対称です")
        println("   → PCGソルバーの使用が理論的に正当化されます")
    else
        println("\n❌ 結論: CalcAX!が生成する行列は非対称です")
        println("   → PBICGSTABソルバーの使用が推奨されます")
    end

    return is_sym1 && is_sym2 && is_sym3
end

# 実行
main()
