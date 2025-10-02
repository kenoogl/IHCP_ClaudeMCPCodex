"""
Phase 2.2 性能劣化デバッグ
"""

using LinearAlgebra
using Printf
BLAS.set_num_threads(1)

using IHCP_CGM

# 小規模問題で詳細計測
ni, nj, nk = 80, 100, 20
nt = 11
dx, dy, dt = 0.12e-3, 0.12e-3, 1e-3

z_total = 0.7e-3
base = 1.6
dz = zeros(nk)
r_sum = sum(base^(i-1) for i in 1:nk)
dz[1] = z_total / r_sum
for i in 2:nk
  dz[i] = dz[i-1] * base
end
dz_b = copy(dz)
dz_b[1] = Inf
dz_t = copy(dz)
dz_t[end] = Inf

rho = 7900.0
cp_coeffs = [0.0, 0.0, 0.0, 500.0]
k_coeffs = [0.0, 0.0, 0.0, 15.0]
T_initial = fill(300.0, ni, nj, nk)
q_surface = fill(1000.0, ni, nj, nt-1)

println("="^60)
println("Phase 2.2 性能劣化デバッグ")
println("="^60)
println("問題サイズ: $(ni)×$(nj)×$(nk)")
println("時間ステップ: $(nt)")
println()

# コンパイル
println("コンパイル中...")
IHCP_CGM.solve_dhcp!(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
  dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000, verbose=false)

# 詳細計測
println("\n単一ステップの詳細計測:")

# 熱物性値計算
T_test = fill(300.0, ni, nj, nk)
t_thermal = @elapsed begin
  for _ in 1:10
    cp, k = IHCP_CGM.ThermalProperties.thermal_properties_calculator(
      T_test, cp_coeffs, k_coeffs
    )
  end
end
@printf("  熱物性値計算×10: %.1f ms\n", t_thermal * 1000)

# 係数構築
t_build = @elapsed begin
  for _ in 1:10
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = IHCP_CGM.DHCPSolver.build_dhcp_system!(
      T_test, q_surface[:, :, 1], rho, fill(500.0, ni, nj, nk), fill(15.0, ni, nj, nk),
      dx, dy, dz, dz_b, dz_t, dt
    )
  end
end
@printf("  係数構築×10: %.1f ms\n", t_build * 1000)

# 疎行列組み立て
a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = IHCP_CGM.DHCPSolver.build_dhcp_system!(
  T_test, q_surface[:, :, 1], rho, fill(500.0, ni, nj, nk), fill(15.0, ni, nj, nk),
  dx, dy, dz, dz_b, dz_t, dt
)
t_assemble = @elapsed begin
  for _ in 1:10
    A = IHCP_CGM.DHCPSolver.assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)
  end
end
@printf("  疎行列組み立て×10: %.1f ms\n", t_assemble * 1000)

# CG求解（単一）
A = IHCP_CGM.DHCPSolver.assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)
N = ni * nj * nk
x0 = zeros(N)
diag_A = diag(A)
inv_diag = [d != 0.0 ? 1.0 / d : 0.0 for d in diag_A]
Pl = Diagonal(inv_diag)

using IterativeSolvers
t_cg = @elapsed begin
  for _ in 1:10
    x, hist = cg(A, b, Pl=Pl, x=x0, rtol=1e-8, maxiter=1000, log=true)
  end
end
@printf("  CG求解×10: %.1f ms\n", t_cg * 1000)

# フル計算（10ステップ）
println("\nフル計算（10ステップ）:")
t_full = @elapsed begin
  T_all = IHCP_CGM.solve_dhcp!(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000, verbose=false)
end
@printf("  総時間: %.1f ms\n", t_full * 1000)
@printf("  単一ステップあたり: %.1f ms\n", t_full * 1000 / (nt - 1))

println("\n内訳分析:")
@printf("  熱物性値: %.1f ms/step\n", t_thermal * 100 / (nt - 1))
@printf("  係数構築: %.1f ms/step\n", t_build * 100 / (nt - 1))
@printf("  疎行列組立: %.1f ms/step\n", t_assemble * 100 / (nt - 1))
@printf("  CG求解: %.1f ms/step\n", t_cg * 100 / (nt - 1))
total_estimated = (t_thermal + t_build + t_assemble + t_cg) * 100 / (nt - 1)
@printf("  推定合計: %.1f ms/step\n", total_estimated)
@printf("  実測合計: %.1f ms/step\n", t_full * 1000 / (nt - 1))
@printf("  差分（オーバーヘッド）: %.1f ms/step\n",
  t_full * 1000 / (nt - 1) - total_estimated)
