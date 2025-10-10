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
import .ThermalProperties: polyval

# 参照データ読み込み
function load_reference_data()
  json_path = joinpath(@__DIR__, "..", "data", "phase1_reference_data.json")
  return JSON.parsefile(json_path)
end

@testset "Phase 1: ThermalProperties Module" begin

  ref_data = load_reference_data()

  @testset "polyval: 多項式評価" begin
    # ρ計算テスト（498.15K）
    rho_coeffs = Vector{Float64}(ref_data["rho_coeffs"])
    expected_rho = Float64(ref_data["polyval_test_rho"])
    calculated_rho = polyval(rho_coeffs, 498.15)

    @test isapprox(calculated_rho, expected_rho, atol=1e-12)
    println("  ρ(498.15K): $(calculated_rho) ≈ $(expected_rho)")

    # 個別温度での比熱・熱伝導率テスト
    cp_coeffs = Vector{Float64}(ref_data["cp_coeffs"])
    k_coeffs = Vector{Float64}(ref_data["k_coeffs"])

    for (label, test_case) in ref_data["single_temp_tests"]
      temp = test_case["temp"]
      expected_cp = test_case["cp"]
      expected_k = test_case["k"]

      calc_cp = polyval(cp_coeffs, temp)
      calc_k = polyval(k_coeffs, temp)

      @test isapprox(calc_cp, expected_cp, atol=1e-12)
      @test isapprox(calc_k, expected_k, atol=1e-12)

      println("  $(label): cp=$(calc_cp) ≈ $(expected_cp)")
      println("           k=$(calc_k) ≈ $(expected_k)")
    end
  end

  @testset "境界値テスト" begin
    cp_coeffs = Vector{Float64}(ref_data["cp_coeffs"])
    k_coeffs = Vector{Float64}(ref_data["k_coeffs"])

    # データの最小・最大温度
    sus304_temp = ref_data["sus304_temp"]
    min_temp = Float64(minimum(sus304_temp))
    max_temp = Float64(maximum(sus304_temp))

    # 境界値での計算（外挿警告なし）
    cp_min = polyval(cp_coeffs, min_temp)
    k_min = polyval(k_coeffs, min_temp)
    cp_max = polyval(cp_coeffs, max_temp)
    k_max = polyval(k_coeffs, max_temp)

    @test cp_min > 0
    @test k_min > 0
    @test cp_max > 0
    @test k_max > 0

    println("  境界値($(min_temp)K): cp=$(cp_min), k=$(k_min)")
    println("  境界値($(max_temp)K): cp=$(cp_max), k=$(k_max)")
  end


end

println("\n=== Phase 1 テスト完了 ===")
println("次のステップ: ThermalProperties.jlの実装")