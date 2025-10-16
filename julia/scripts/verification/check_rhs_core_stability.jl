#!/usr/bin/env julia

"""
calRHS_core!関数の型安定性チェック
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
    BC.heat_flux_bc(1000.0),  # x_minus
    BC.heat_flux_bc(2000.0),  # x_plus
    BC.adiabatic_bc(),        # y_minus
    BC.adiabatic_bc(),        # y_plus
    BC.adiabatic_bc(),        # z_minus
    BC.convection_bc(50.0, 300.0)   # z_plus
  )

  println("== @code_warntype calRHS_core! ==")
  @code_warntype RHSCore.calRHS_core!(b, θ, dx, dy, dz, bc_set, "sequential")
  println()

  println("== @code_warntype apply_face_bc! ==")
  @code_warntype RHSCore.apply_face_bc!(b, θ, dx, dy, dz, bc_set.x_minus, :x_minus, "sequential")
  println()
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
