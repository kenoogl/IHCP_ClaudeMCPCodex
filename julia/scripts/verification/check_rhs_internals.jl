#!/usr/bin/env julia

"""
RHSCore関数の内部変数の型安定性チェック
"""

using IHCP_CGM
using InteractiveUtils

const RHSCore = IHCP_CGM.RHSCore
const BC = IHCP_CGM.BoundaryConditions

function main()
  # サンプルデータ準備
  ni, nj, nk = 5, 5, 5
  b = zeros(Float64, ni, nj, nk)
  θ = fill(300.0, ni, nj, nk)
  dx = 0.05
  dy = 0.05
  dz = fill(0.02, nk)

  # 境界条件セット作成
  bc_set = BC.BoundaryConditionSet(
    BC.heat_flux_bc(1000.0),
    BC.heat_flux_bc(2000.0),
    BC.adiabatic_bc(),
    BC.adiabatic_bc(),
    BC.adiabatic_bc(),
    BC.convection_bc(50.0, 300.0)
  )

  println("=== RHS_heat_flux! の型推論 ===")
  @code_warntype RHSCore.RHS_heat_flux!(b, dx, dy, dz, bc_set.x_minus, :x_minus, "sequential")
  println()

  println("=== RHS_convection! の型推論 ===")
  @code_warntype RHSCore.RHS_convection!(b, θ, dx, dy, dz, bc_set.z_plus, :z_plus, "sequential")
  println()

  # Float32でもチェック
  println("=== Float32での型推論チェック ===")
  b32 = zeros(Float32, ni, nj, nk)
  θ32 = fill(Float32(300.0), ni, nj, nk)
  dx32 = Float32(0.05)
  dy32 = Float32(0.05)
  dz32 = fill(Float32(0.02), nk)

  println("\n--- RHS_heat_flux! (Float32) ---")
  @code_warntype RHSCore.RHS_heat_flux!(b32, dx32, dy32, dz32, bc_set.x_minus, :x_minus, "sequential")
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
