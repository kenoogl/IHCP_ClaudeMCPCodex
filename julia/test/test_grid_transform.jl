"""
test_grid_transform.jl

GridTransform.jl のテスト：格子変換とガイドセル初期化
"""

using Test

# テスト用に相対パスからモジュールをインクルード
include("../src/utils/GridTransform.jl")
using .GridTransform
# Phase C実装予定の関数はコメントアウト
# import .GridTransform: BoundaryType, ISOTHERMAL, HEAT_FLUX, CONVECTION
# import .GridTransform: initialize_guard_cells!, compute_z_range, λf
import .GridTransform: convert_to_guard_cell_grid

@testset "GridTransform Tests" begin

  @testset "Test 1: 等間隔格子の変換（3層）" begin
    # IHCP形式の格子定義（等間隔）
    nk = 3
    dz_uniform = 0.1e-3  # 100μm

    dz = fill(dz_uniform, nk)
    dz_b = fill(0.5 * dz_uniform, nk)
    dz_b[1] = Inf  # 底面境界（実際には使用しない）
    dz_t = fill(0.5 * dz_uniform, nk)
    dz_t[nk] = Inf  # 上面境界

    # Heat3ds形式に変換
    Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

    # 配列サイズ確認
    @test length(Z) == nk + 2  # MZ = 5
    @test length(ΔZ) == nk + 1  # MZ - 1 = 4

    # セル中心座標の確認（等間隔の期待値）
    # Z[1]: 下側ガイドセル
    # Z[2]: 最下層セル中心（基準点0.0）
    # Z[3]: 第2層セル中心（0.0 + dz_uniform）
    # Z[4]: 第3層セル中心（0.0 + 2*dz_uniform）
    # Z[5]: 上側ガイドセル

    @test Z[2] ≈ 0.0  # 最下層セル中心（基準点）
    println("Z[1] (下側ガイド): ", Z[1])
    println("Z[2] (最下層): ", Z[2])
    println("Z[3] (第2層): ", Z[3])
    println("Z[4] (第3層): ", Z[4])
    println("Z[5] (上側ガイド): ", Z[5])

    # セル代表幅の確認
    println("\nΔZ:")
    for k in 1:length(ΔZ)
      println("  ΔZ[$k]: ", ΔZ[k])
    end
  end


  @testset "Test 2: ガイドセル配列の初期化" begin
    # Phase C実装予定（initialize_guard_cells!未実装）
    @test_skip "initialize_guard_cells! は Phase C で実装予定"
  end


  @testset "Test 3: λf関数（調和平均 + マスク補正）" begin
    # Phase C実装予定（λf関数未実装）
    @test_skip "λf 関数は Phase C で実装予定"
  end


  @testset "Test 4: compute_z_range（境界条件に応じた計算範囲）" begin
    # Phase C実装予定（compute_z_range, BoundaryType未実装）
    @test_skip "compute_z_range と BoundaryType は Phase C で実装予定"
  end

end
