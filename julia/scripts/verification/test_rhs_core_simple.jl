#!/usr/bin/env julia

"""
calRHS_core!の動作確認テスト
"""

using IHCP_CGM

const RHSCore = IHCP_CGM.RHSCore
const BC = IHCP_CGM.BoundaryConditions

function main()
  println("=== calRHS_core! 動作確認テスト ===\n")

  # サンプルデータ準備
  ni, nj, nk = 5, 5, 5
  b = zeros(Float64, ni, nj, nk)
  θ = fill(300.0, ni, nj, nk)
  dx = 0.05
  dy = 0.05
  dz = fill(0.02, nk)

  # テスト1: 熱流束境界条件
  println("テスト1: 熱流束境界条件")
  bc_set = BC.BoundaryConditionSet(
    BC.heat_flux_bc(1000.0),  # x_minus
    BC.heat_flux_bc(2000.0),  # x_plus
    BC.adiabatic_bc(),        # y_minus
    BC.adiabatic_bc(),        # y_plus
    BC.adiabatic_bc(),        # z_minus
    BC.adiabatic_bc()         # z_plus
  )

  fill!(b, 0.0)
  result = RHSCore.calRHS_core!(b, θ, dx, dy, dz, bc_set, "sequential")
  println("  戻り値の型: ", typeof(result))
  println("  戻り値: ", result)
  println("  b[2,3,3] (x_minus境界近く): ", b[2,3,3])
  println("  b[4,3,3] (x_plus境界近く): ", b[4,3,3])
  println("  期待値: x_minus: 1000.0/0.05 = 20000.0")
  println("  期待値: x_plus: -2000.0/0.05 = -40000.0")
  println()

  # テスト2: 熱伝達境界条件
  println("テスト2: 熱伝達境界条件")
  bc_set2 = BC.BoundaryConditionSet(
    BC.adiabatic_bc(),
    BC.adiabatic_bc(),
    BC.adiabatic_bc(),
    BC.adiabatic_bc(),
    BC.adiabatic_bc(),
    BC.convection_bc(50.0, 350.0)  # z_plus: h=50, T_amb=350
  )

  fill!(b, 0.0)
  result2 = RHSCore.calRHS_core!(b, θ, dx, dy, dz, bc_set2, "sequential")
  println("  戻り値の型: ", typeof(result2))
  println("  戻り値: ", result2)
  println("  b[3,3,4] (z_plus境界近く): ", b[3,3,4])
  println("  期待値: h * (T_amb - T) / dz = 50.0 * (350.0 - 300.0) / 0.02 = 125000.0")
  println()

  # テスト3: 並列実行（スレッド数が1より大きい場合のみ）
  if Threads.nthreads() > 1
    println("テスト3: 並列実行")
    fill!(b, 0.0)
    result3 = RHSCore.calRHS_core!(b, θ, dx, dy, dz, bc_set, "thread")
    println("  戻り値の型: ", typeof(result3))
    println("  戻り値: ", result3)
    println("  b[2,3,3]: ", b[2,3,3])
    println()
  else
    println("スキップ: 並列テスト（スレッド数=1）")
  end

  println("=== テスト完了 ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
