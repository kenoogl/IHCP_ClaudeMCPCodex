"""
test_grid_transform.jl

GridTransform.jl のテスト：格子変換とガイドセル初期化
"""

using Test

# テスト用に相対パスからモジュールをインクルード
include("../src/utils/GridTransform.jl")
using .GridTransform
import .GridTransform: BoundaryType, ISOTHERMAL, HEAT_FLUX, CONVECTION
import .GridTransform: convert_to_guard_cell_grid, initialize_guard_cells!, compute_z_range, λf

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
    ni, nj, nk = 3, 3, 3

    # 内点データ
    θ_init = fill(300.0, ni, nj, nk)
    λ_init = fill(15.0, ni, nj, nk)

    # ガイドセル配列（サイズ: ni+2, nj+2, nk+2）
    θ = zeros(Float64, ni+2, nj+2, nk+2)
    λ = zeros(Float64, ni+2, nj+2, nk+2)
    mask = zeros(Float64, ni+2, nj+2, nk+2)

    # 初期化
    initialize_guard_cells!(θ, λ, mask, θ_init, λ_init)

    # マスク初期値確認（全て1.0であるべき）
    @test all(mask .== 1.0)

    # 内点データの転送確認
    for k in 1:nk, j in 1:nj, i in 1:ni
      @test θ[i+1, j+1, k+1] ≈ θ_init[i, j, k]
      @test λ[i+1, j+1, k+1] ≈ λ_init[i, j, k]
    end

    println("\n内点データ転送確認: OK")
    println("θ[2,2,2] (内点[1,1,1]): ", θ[2,2,2])
    println("λ[2,2,2] (内点[1,1,1]): ", λ[2,2,2])
  end


  @testset "Test 3: λf関数（調和平均 + マスク補正）" begin
    # 通常の調和平均（両方内点）
    a, b = 15.0, 20.0
    ma, mb = 1.0, 1.0
    result = λf(a, b, ma, mb)
    expected = 2.0 * a * b / (a + b) * (2.0 - div(Int(ma) + Int(mb), 2))
    @test result ≈ expected
    println("\nλf(15, 20, 1, 1) = ", result, " (期待値: ", expected, ")")

    # 一方が境界
    ma, mb = 1.0, 0.0
    result = λf(a, b, ma, mb)
    expected = 2.0 * a * b / (a + b) * (2.0 - div(Int(ma) + Int(mb), 2))
    @test result ≈ expected
    println("λf(15, 20, 1, 0) = ", result, " (期待値: ", expected, ")")

    # 両方境界（λ=0の場合を想定）
    a, b = 0.0, 20.0
    ma, mb = 0.0, 1.0
    result = λf(a, b, ma, mb)
    @test result ≈ 0.0
    println("λf(0, 20, 0, 1) = ", result, " (期待値: 0.0)")
  end


  @testset "Test 4: compute_z_range（境界条件に応じた計算範囲）" begin
    nk = 10  # 計算内点数、配列サイズはnk+2=12

    # パターン1: ISOTHERMAL（下側）, HEAT_FLUX（上側）
    # 典型的なIHCP問題（底面等温、上面熱流束）
    z_range = compute_z_range(nk, ISOTHERMAL, HEAT_FLUX)
    @test z_range == [3, nk+1]
    println("\ncompute_z_range($nk, ISOTHERMAL, HEAT_FLUX) = ", z_range)

    # パターン2: HEAT_FLUX（下側）, ISOTHERMAL（上側）
    z_range = compute_z_range(nk, HEAT_FLUX, ISOTHERMAL)
    @test z_range == [2, nk]
    println("compute_z_range($nk, HEAT_FLUX, ISOTHERMAL) = ", z_range)

    # パターン3: HEAT_FLUX（下側）, CONVECTION（上側）
    # Heat3dsのMode3（底面PCB温度、上面熱伝達、側面断熱）
    z_range = compute_z_range(nk, HEAT_FLUX, CONVECTION)
    @test z_range == [2, nk+1]
    println("compute_z_range($nk, HEAT_FLUX, CONVECTION) = ", z_range)

    # パターン4: ISOTHERMAL（下側）, ISOTHERMAL（上側）
    # 両端固定
    z_range = compute_z_range(nk, ISOTHERMAL, ISOTHERMAL)
    @test z_range == [3, nk]
    println("compute_z_range($nk, ISOTHERMAL, ISOTHERMAL) = ", z_range)

    # パターン5: CONVECTION（下側）, CONVECTION（上側）
    z_range = compute_z_range(nk, CONVECTION, CONVECTION)
    @test z_range == [2, nk+1]
    println("compute_z_range($nk, CONVECTION, CONVECTION) = ", z_range)
  end

end
