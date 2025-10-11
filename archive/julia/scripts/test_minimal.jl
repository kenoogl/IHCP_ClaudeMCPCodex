#!/usr/bin/env julia
"""
最小限のテストスクリプト（10ステップ）
バグの原因を特定するため
"""

using NPZ
using Printf

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function main()
  println("="^80)
  println("Minimal 10-step test")
  println("="^80)

  # 問題定義
  nt = 10
  dt = 1.0e-3
  ni, nj, nk = 80, 100, 20
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  # CGM parameters
  window_size = 9  # 10ステップに対して9
  overlap = 2
  cgm_iteration = 1
  q_init_value = 0.0

  println("  grid: ni=$ni, nj=$nj, nk=$nk")
  println("  time steps: nt=$nt, dt=$dt s")
  println("  window_size=$window_size, overlap=$overlap")

  # Z方向格子
  z_faces_normalized = range(1.0, stop=0.0, length=nk + 1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  dz = diff(z_faces)
  z_centers = similar(dz)
  z_centers[1] = z_faces[1]
  z_centers[end] = z_faces[end]
  if nk > 2
    z_centers[2:end-1] = (z_faces[2:end-2] .+ z_faces[3:end-1]) ./ 2
  end

  dz_top = zeros(Float64, nk)
  dz_top[end] = Inf
  if nk > 1
    dz_top[1:end-1] = z_centers[2:end] .- z_centers[1:end-1]
  end

  dz_bottom = zeros(Float64, nk)
  dz_bottom[1] = Inf
  if nk > 1
    dz_bottom[2:end] = z_centers[2:end] .- z_centers[1:end-1]
  end

  Z, dZ_guard = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

  # 材料物性値
  rho = 7823.493962874829
  cp_coeffs = Float64[2.0092965857857302e-10, -3.426055709046039e-7, 0.13492793577599277, 469.85285976023886]
  k_coeffs = Float64[4.799122446183406e-12, -8.182993477116119e-9, 0.016176544533897493, 8.117517482517492]

  # ソルバーパラメータ
  rtol_dhcp = 1.0e-6
  maxiter_dhcp = 20000
  rtol_adjoint = 1.0e-8
  maxiter_adjoint = 20000

  # 入力データ読み込み
  println("\nLoading measurement data...")
  data_path = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
  if !isfile(data_path)
    error("Measurement file not found: $data_path")
  end

  T_measure_full = npzread(data_path)
  Y_obs_python = T_measure_full[1:nt, :, :]
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))

  println("  measurement shape: $(size(Y_obs))")

  # 初期温度
  T_measure_init = Y_obs[:, :, 1]
  T_init = repeat(reshape(T_measure_init, ni, nj, 1), outer=(1, 1, nk))
  T_init = convert(Array{Float64,3}, T_init)

  println("  initial temperature range: $(minimum(T_init)) ~ $(maximum(T_init)) K")

  # WorkBuffers
  wk = IHCP_CGM.WorkBuffers(ni + 2, nj + 2, nk + 2)

  # スライディングウィンドウCGM実行
  println("\nRunning sliding window CGM...")
  solve_start = time()

  # 0a3538f時点のAPIに合わせて呼び出し（dz, dz_b, dz_tが必要）
  q_result, windows_info = solve_sliding_window_cgm(
    Y_obs,
    T_init,
    wk,
    dx,
    dy,
    Z,
    dZ_guard,
    dz,
    dz_bottom,
    dz_top,
    dt,
    rho,
    cp_coeffs,
    k_coeffs,
    window_size,
    overlap,
    q_init_value,
    cgm_iteration;
    rtol_dhcp=rtol_dhcp,
    maxiter_dhcp=maxiter_dhcp,
    rtol_adjoint=rtol_adjoint,
    maxiter_adjoint=maxiter_adjoint
  )

  cgm_elapsed = time() - solve_start

  println("\n✅ SUCCESS!")
  println("  elapsed time: $(cgm_elapsed) s")
  println("  q_result shape: $(size(q_result))")
  println("  heat flux range: $(minimum(q_result)) ~ $(maximum(q_result)) W/m^2")

  return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
  try
    exit(main())
  catch err
    println("\n❌ Error: ", err)
    Base.showerror(stdout, err, catch_backtrace())
    println()
    exit(1)
  end
end
