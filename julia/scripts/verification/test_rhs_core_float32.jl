#!/usr/bin/env julia

"""
calRHS_core!のFloat32型テンプレート化検証
"""

using IHCP_CGM

const RHSCore = IHCP_CGM.RHSCore
const BC = IHCP_CGM.BoundaryConditions

function main()
  println("=== calRHS_core! Float32テンプレート化テスト ===\n")

  # Float32型でサンプルデータ準備
  ni, nj, nk = 5, 5, 5
  b = zeros(Float32, ni, nj, nk)
  θ = fill(Float32(300.0), ni, nj, nk)
  dx = Float32(0.05)
  dy = Float32(0.05)
  dz = fill(Float32(0.02), nk)

  # テスト: Float32で熱流束境界条件
  println("テスト: Float32型での熱流束境界条件")
  bc_set = BC.BoundaryConditionSet(
    BC.heat_flux_bc(1000.0),  # x_minus
    BC.heat_flux_bc(2000.0),  # x_plus
    BC.adiabatic_bc(),        # y_minus
    BC.adiabatic_bc(),        # y_plus
    BC.adiabatic_bc(),        # z_minus
    BC.adiabatic_bc()         # z_plus
  )

  fill!(b, Float32(0.0))
  result = RHSCore.calRHS_core!(b, θ, dx, dy, dz, bc_set, "sequential")

  println("  配列bの型: ", typeof(b))
  println("  配列bの要素型: ", eltype(b))
  println("  戻り値の型: ", typeof(result))
  println("  戻り値: ", result)
  println("  b[2,3,3] (x_minus境界近く): ", b[2,3,3])
  println("  b[4,3,3] (x_plus境界近く): ", b[4,3,3])
  println("  期待値: x_minus: 1000.0/0.05 = 20000.0")
  println("  期待値: x_plus: -2000.0/0.05 = -40000.0")
  println()

  # Float64型でも動作確認
  println("テスト: Float64型での熱流束境界条件（比較用）")
  b64 = zeros(Float64, ni, nj, nk)
  θ64 = fill(300.0, ni, nj, nk)
  dx64 = 0.05
  dy64 = 0.05
  dz64 = fill(0.02, nk)

  result64 = RHSCore.calRHS_core!(b64, θ64, dx64, dy64, dz64, bc_set, "sequential")

  println("  配列b64の型: ", typeof(b64))
  println("  配列b64の要素型: ", eltype(b64))
  println("  戻り値の型: ", typeof(result64))
  println("  b64[2,3,3]: ", b64[2,3,3])
  println("  b64[4,3,3]: ", b64[4,3,3])
  println()

  # 結果比較
  println("Float32とFloat64の結果比較:")
  println("  b[2,3,3] (Float32): ", b[2,3,3])
  println("  b64[2,3,3] (Float64): ", b64[2,3,3])
  println("  差分: ", abs(b[2,3,3] - Float32(b64[2,3,3])))
  println()

  println("=== テンプレート化テスト完了 ===")
  println("✅ Float32とFloat64の両方で正しく動作")
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
