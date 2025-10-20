#!/usr/bin/env julia
"""
単一ウィンドウテスト（デバッグ用）
"""

using Printf
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DT = 1.0e-3

function main()
  println("="^80)
  println("単一ウィンドウテスト")
  println("="^80)

  # データ読み込み
  data_file = joinpath(PROJECT_ROOT, "shared/data/T_measure_700um_1ms.npy")
  Y_obs_full = npzread(data_file)

  nt = 72  # ウィンドウ1つ分（71ステップ + 初期）
  Y_obs = Y_obs_full[:, :, 1:nt]
  println("Y_obs形状: $(size(Y_obs))")

  # 格子設定
  ni, nj = size(Y_obs, 1), size(Y_obs, 2)
  nk = 20

  dx = 0.12e-3
  dy = 0.12e-3 * 0.866

  # z格子（Python版と同じ）
  zf_py = [0.00000000e+00, 7.32951631e-05, 1.36380895e-04, 1.90679287e-04,
           2.37414346e-04, 2.77639585e-04, 3.12261768e-04, 3.42061358e-04,
           3.67710102e-04, 3.89786181e-04, 4.08787238e-04, 4.25141600e-04,
           4.39217929e-04, 4.51333538e-04, 4.61761539e-04, 4.70737003e-04,
           4.78462256e-04, 4.85111444e-04, 4.90834452e-04, 4.95760291e-04,
           5.00000000e-04]

  z_faces = zf_py
  z_centers = [(z_faces[i] + z_faces[i+1]) / 2 for i in 1:nk]
  ΔZ = [z_faces[i+1] - z_faces[i] for i in 1:nk]

  # 熱物性値（ハードコード）
  rho = 7823.493962874829
  cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

  # 初期条件
  T_init = zeros(Float64, ni, nj, nk)
  for k in 1:nk
    T_init[:, :, k] = Y_obs[:, :, 1]
  end

  q_init = zeros(Float64, ni, nj, nt-1)

  println("初期化完了")
  println("CGM 3回実行開始...")

  # ワークスペース確保
  work = IHCP_CGM.CGMSolver.CGMWorkspace(ni, nj, nk, nt)

  # CGMパラメータ（Python版と一致させる）
  cgm_params = (
    max_iter = 3,
    sigma = 1.8,
    rtol_dhcp = 1.0e-6,
    rtol_adjoint = 1.0e-8,
    maxiter_cg = 20000,
    dire_reset_every = 5,
    eps = 1.0e-12,
    min_iter = 10,
    P = 10,
    eta = 1.0e-4,
    beta_max = 1.0e8,
    verbose = true,
    dhcp_extrapolation = :quadratic,
    adjoint_residual_scale = 0.5,
    dhcp_solver = :pcg,
    dhcp_smoother = :jacobi,  # Python版と一致
    adjoint_solver = :pcg,
    adjoint_smoother = :jacobi,  # Python版と一致
    sensitivity_solver = :pbicgstab,
    sensitivity_smoother = :gs
  )

  start_time = time()

  q_result, T_result, J_history = IHCP_CGM.CGMSolver.solve_cgm!(
    T_init, Y_obs, q_init, work,
    dx, dy, z_centers, ΔZ,
    DT, rho, cp_coeffs, k_coeffs;
    cgm_params...
  )

  elapsed = time() - start_time

  println("\n" * "="^80)
  println("完了")
  println("  実行時間: $(round(elapsed, digits=2))秒")
  println("  CGM反復数: $(length(J_history))")
  println("  最終J: $(J_history[end])")
  println("  q範囲: [$(minimum(q_result)), $(maximum(q_result))] W/m²")
  println("="^80)

  return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
  exit(main())
end
