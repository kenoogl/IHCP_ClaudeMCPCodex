"""
Phase 1: 熱物性値計算モジュールのテスト

TDD方針に従い、実装前にテストを定義
Python参照データとの数値比較（許容誤差: 1e-12）
"""

using Test
using JSON

# テスト用に相対パスからモジュールをインクルード
include("../src/ThermalProperties.jl")
using .ThermalProperties

# 参照データ読み込み
function load_reference_data()
  json_path = joinpath(@__DIR__, "..", "data", "phase1_reference_data.json")
  return JSON.parsefile(json_path)
end

@testset "Phase 1: ThermalProperties Module" begin

  ref_data = load_reference_data()

  @testset "polyval_numba: 多項式評価" begin
    # ρ計算テスト（498.15K）
    rho_coeffs = Vector{Float64}(ref_data["rho_coeffs"])
    expected_rho = Float64(ref_data["polyval_test_rho"])
    calculated_rho = polyval_numba(rho_coeffs, 498.15)

    @test isapprox(calculated_rho, expected_rho, atol=1e-12)
    println("  ρ(498.15K): $(calculated_rho) ≈ $(expected_rho)")

    # 個別温度での比熱・熱伝導率テスト
    cp_coeffs = Vector{Float64}(ref_data["cp_coeffs"])
    k_coeffs = Vector{Float64}(ref_data["k_coeffs"])

    for (label, test_case) in ref_data["single_temp_tests"]
      temp = test_case["temp"]
      expected_cp = test_case["cp"]
      expected_k = test_case["k"]

      calc_cp = polyval_numba(cp_coeffs, temp)
      calc_k = polyval_numba(k_coeffs, temp)

      @test isapprox(calc_cp, expected_cp, atol=1e-12)
      @test isapprox(calc_k, expected_k, atol=1e-12)

      println("  $(label): cp=$(calc_cp) ≈ $(expected_cp)")
      println("           k=$(calc_k) ≈ $(expected_k)")
    end
  end

  @testset "thermal_properties_calculator: 3D配列計算" begin
    cp_coeffs = Vector{Float64}(ref_data["cp_coeffs"])
    k_coeffs = Vector{Float64}(ref_data["k_coeffs"])

    # Python配列（行優先）をJulia配列（列優先）に変換
    # Python: (ni, nj, nk) = (3, 3, 3)
    # Julia: (nk, nj, ni) = (3, 3, 3) として転置
    test_temp_list = ref_data["test_temp_small"]
    expected_cp_list = ref_data["cp_small"]
    expected_k_list = ref_data["k_small"]

    # 多次元配列の変換（Pythonの[i][j][k]がJuliaで[k,j,i]になる）
    ni, nj, nk = 3, 3, 3
    test_temp = zeros(nk, nj, ni)
    expected_cp = zeros(nk, nj, ni)
    expected_k = zeros(nk, nj, ni)

    for i in 1:ni
      for j in 1:nj
        for k_idx in 1:nk
          test_temp[k_idx, j, i] = test_temp_list[i][j][k_idx]
          expected_cp[k_idx, j, i] = expected_cp_list[i][j][k_idx]
          expected_k[k_idx, j, i] = expected_k_list[i][j][k_idx]
        end
      end
    end

    # 実際の計算実行（これが失敗することを確認する）
    calc_cp, calc_k = thermal_properties_calculator(test_temp, cp_coeffs, k_coeffs)

    # 全要素の比較（許容誤差1e-12）
    @test isapprox(calc_cp, expected_cp, atol=1e-12)
    @test isapprox(calc_k, expected_k, atol=1e-12)

    # 統計情報表示
    println("  配列形状: $(size(test_temp))")
    println("  温度範囲: $(minimum(test_temp)) - $(maximum(test_temp)) K")
    println("  cp範囲: $(minimum(calc_cp)) - $(maximum(calc_cp)) J/(kg·K)")
    println("  k範囲: $(minimum(calc_k)) - $(maximum(calc_k)) W/(m·K)")
    println("  最大誤差(cp): $(maximum(abs.(calc_cp - expected_cp)))")
    println("  最大誤差(k): $(maximum(abs.(calc_k - expected_k)))")
  end

  @testset "境界値テスト" begin
    cp_coeffs = Vector{Float64}(ref_data["cp_coeffs"])
    k_coeffs = Vector{Float64}(ref_data["k_coeffs"])

    # データの最小・最大温度
    sus304_temp = ref_data["sus304_temp"]
    min_temp = Float64(minimum(sus304_temp))
    max_temp = Float64(maximum(sus304_temp))

    # 境界値での計算（外挿警告なし）
    cp_min = polyval_numba(cp_coeffs, min_temp)
    k_min = polyval_numba(k_coeffs, min_temp)
    cp_max = polyval_numba(cp_coeffs, max_temp)
    k_max = polyval_numba(k_coeffs, max_temp)

    @test cp_min > 0
    @test k_min > 0
    @test cp_max > 0
    @test k_max > 0

    println("  境界値($(min_temp)K): cp=$(cp_min), k=$(k_min)")
    println("  境界値($(max_temp)K): cp=$(cp_max), k=$(k_max)")
  end

  @testset "配列形状の一致性" begin
    # 様々なサイズでテスト
    test_sizes = [(2, 2, 2), (4, 3, 5), (10, 1, 1)]
    cp_coeffs = Vector{Float64}(ref_data["cp_coeffs"])
    k_coeffs = Vector{Float64}(ref_data["k_coeffs"])

    for (nk, nj, ni) in test_sizes
      # ランダム温度配列（300-1600K）
      temp_array = 300.0 .+ rand(nk, nj, ni) .* 1300.0

      calc_cp, calc_k = thermal_properties_calculator(temp_array, cp_coeffs, k_coeffs)

      @test size(calc_cp) == size(temp_array)
      @test size(calc_k) == size(temp_array)
      @test all(calc_cp .> 0)
      @test all(calc_k .> 0)

      println("  サイズ$(size(temp_array)): OK")
    end
  end

end

println("\n=== Phase 1 テスト完了 ===")
println("次のステップ: ThermalProperties.jlの実装")