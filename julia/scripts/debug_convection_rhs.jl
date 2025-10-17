#!/usr/bin/env julia
"""
対流境界条件のRHS/係数行列デバッグ

Phase 6: RHSと係数行列の値を直接確認
"""

using Printf

# プロジェクトのソースコードをロード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM
using .IHCP_CGM.BoundaryConditions
using .IHCP_CGM.DHCPSolver: set_dhcp_bc_parameters_convection
using .IHCP_CGM.Commons
using .IHCP_CGM.ThermalProperties: set_properties!
using .IHCP_CGM.RHSCore: calRHS_core!
using .IHCP_CGM.CommonSolver: CalcRK!, CalcAX!

function debug_rhs()
  println("="^60)
  println("対流境界条件のRHS/係数行列デバッグ")
  println("="^60)

  # パラメータ設定（小規模）
  ni, nj, nk = 3, 3, 3
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  dt = 1e-3
  rho = 7823.493962874829

  # z方向格子（均等）
  Lz = 0.5e-3
  dz_val = Lz / nk
  dz = fill(dz_val, nk)
  ΔZ = zeros(nk + 2)
  ΔZ[2:nk+1] = dz
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  z_centers = zeros(nk + 2)
  z_centers[1] = 0.0
  z_centers[nk+2] = Lz
  for k in 2:nk+1
    z_centers[k] = (k-1.5) * dz_val
  end

  # 熱物性値係数（SUS304）
  cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

  # 初期温度場（550K均一）
  T_init = fill(550.0, ni, nj, nk)

  # WorkBuffers
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  println("\n【テスト条件】")
  println("格子: $(ni)×$(nj)×$(nk)")
  println("dx=$(dx*1e3) mm, dy=$(dy*1e3) mm, dz=$(dz_val*1e3) mm")
  println("密度: $(rho) kg/m³")
  println("初期温度: 550K（均一）")

  # 境界条件設定（h=1.0, T_∞=300K）
  bc_set = set_dhcp_bc_parameters_convection(1.0, 300.0)
  println("\n【境界条件】")
  println("Z-plus面: 対流 (h=1.0 W/m²·K, T_∞=300K)")
  println("他の面: 断熱")

  # 温度場を work.θ にコピー
  for k in 1:nk, j in 1:nj, i in 1:ni
    work.θ[i+1, j+1, k+1] = T_init[i, j, k]
    work.mask[i+1, j+1, k+1] = 1.0
  end

  # 熱物性値計算
  set_properties!(T_init, work.cp, work.λ, cp_coeffs, k_coeffs)

  # 境界条件適用
  apply_boundary_conditions!(work.θ, work.λ, work.cp, work.mask, bc_set)

  # HC配列生成
  HC = set_BC_coef(bc_set)
  println("\n【HC配列】")
  println("HC = ", HC)

  # RHSコア計算
  calRHS_core!(work.b, work.θ, dx, dy, ΔZ, bc_set, "sequential")

  # 最終RHS計算（DHCPと同じ）
  ddt = inv(dt)
  for k in 2:nk+1, j in 2:nj+1, i in 2:ni+1
    dz_k = ΔZ[k]
    a_p_0 = rho * work.cp[i,j,k] * dx * dy * dz_k * ddt
    work.b[i,j,k] = -( a_p_0 * work.θ[i,j,k] + work.b[i,j,k] )
  end

  println("\n【RHS (work.b)】")
  for k in nk+1:-1:2
    println("k=$(k-1):")
    for j in 2:nj+1
      for i in 2:ni+1
        @printf("  [%d,%d,%d]: %.6e\n", i-1, j-1, k-1, work.b[i,j,k])
      end
    end
  end

  # 残差計算（CalcRK!）
  Δh = (dx, dy, 1.0)
  res0 = CalcRK!(work.pcg_r, work.θ, work.b, work.λ, work.cp, work.mask, rho, Δh, dt, z_centers, ΔZ, HC, "sequential")

  println("\n【初期残差】")
  @printf("res0 = %.6e\n", res0)

  if isnan(res0) || isinf(res0)
    println("❌ 初期残差がNaN/Inf!")
  end

  println("\n【残差ベクトル (work.pcg_r)】")
  for k in nk+1:-1:2
    println("k=$(k-1):")
    for j in 2:nj+1
      for i in 2:ni+1
        @printf("  [%d,%d,%d]: %.6e\n", i-1, j-1, k-1, work.pcg_r[i,j,k])
      end
    end
  end

  # AX計算
  CalcAX!(work.pcg_q, work.θ, Δh, dt, work.λ, work.cp, work.mask, rho, z_centers, ΔZ, HC, "sequential")

  println("\n【AX (work.pcg_q)】")
  for k in nk+1:-1:2
    println("k=$(k-1):")
    for j in 2:nj+1
      for i in 2:ni+1
        @printf("  [%d,%d,%d]: %.6e\n", i-1, j-1, k-1, work.pcg_q[i,j,k])
      end
    end
  end

  println("\n" * "="^60)
end

# メイン実行
if abspath(PROGRAM_FILE) == @__FILE__
  try
    debug_rhs()
  catch e
    println("\n❌ エラーが発生しました: ", e)
    Base.showerror(stdout, e, catch_backtrace())
    println()
    exit(1)
  end
end
