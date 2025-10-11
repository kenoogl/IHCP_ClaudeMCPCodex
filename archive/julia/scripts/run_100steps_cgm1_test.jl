#!/usr/bin/env julia
"""
Julia 100-step CGM validation runner (CGM iteration = 1, with sliding window)

This script reproduces the 100-step inverse heat conduction test.
Uses sliding window approach to match the previous validation results.
"""

using Dates
using Printf
using Statistics
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))

function main()
  println("="^80)
  println("Julia 100-step CGM validation (CGM iteration=1)")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  println("Julia threads: $(Threads.nthreads())")

  total_start = time()

  # Problem definition -----------------------------------------------------
  nt = 100
  dt = 1.0e-3
  ni, nj, nk = 80, 100, 20
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  # CGM parameters (matching previous validation)
  window_size = 71
  overlap = 17
  cgm_iteration = 1
  q_init_value = 0.0

  println("\n[1/6] Grid and material parameters")
  println(@sprintf("  grid: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  spacing: dx=%.3e m, dy=%.3e m", dx, dy))
  println(@sprintf("  time steps: nt=%d, dt=%.3e s", nt, dt))
  println(@sprintf("  CGM: window_size=%d, overlap=%d, iteration=%d", window_size, overlap, cgm_iteration))

  # Build z-direction grid
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

  # Use the same material property values as Python reference
  rho = 7823.493962874829
  cp_coeffs = Float64[2.0092965857857302e-10, -3.426055709046039e-7, 0.13492793577599277, 469.85285976023886]
  k_coeffs = Float64[4.799122446183406e-12, -8.182993477116119e-9, 0.016176544533897493, 8.117517482517492]

  println(@sprintf("  rho (at 498.15K): %.6f kg/m^3", rho))

  # Solver parameters
  rtol_dhcp = 1.0e-6
  maxiter_dhcp = 20000
  rtol_adjoint = 1.0e-8
  maxiter_adjoint = 20000

  # Input preparation ------------------------------------------------------
  println("\n[2/6] Loading measurement data")
  data_path = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
  if !isfile(data_path)
    error("Measurement file not found: $data_path")
  end

  println("  loading measurement data from: $data_path")
  load_start = time()
  T_measure_full = npzread(data_path)
  ni_file = size(T_measure_full, 2)
  nj_file = size(T_measure_full, 3)
  if size(T_measure_full, 1) < nt
    error("Measurement file only has $(size(T_measure_full, 1)) time steps; need $nt")
  end

  # Python format (nt, ni, nj) → Julia format (ni, nj, nt)
  Y_obs_python = T_measure_full[1:nt, :, :]
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))
  Y_obs_python = nothing
  T_measure_full = nothing
  GC.gc()

  println("  measurement subset shape (Julia order): $(size(Y_obs))")
  println(@sprintf("  load time: %.2f s", time() - load_start))

  # Initial temperature from first frame
  T_measure_init = Y_obs[:, :, 1]
  T_init = repeat(reshape(T_measure_init, ni, nj, 1), outer=(1, 1, nk))
  T_init = convert(Array{Float64,3}, T_init)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))

  # WorkBuffers
  wk = IHCP_CGM.WorkBuffers(ni + 2, nj + 2, nk + 2)

  # Sliding window CGM solve -----------------------------------------------
  println("\n[3/6] Running sliding window CGM solver")
  solve_start = time()

  q_result, windows_info = solve_sliding_window_cgm(
    Y_obs,
    T_init,
    wk,
    dx,
    dy,
    Z,
    dZ_guard,
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
  println(@sprintf("  CGM elapsed time: %.2f s", cgm_elapsed))
  println(@sprintf("  iterations completed: %d windows", length(windows_info)))
  println(@sprintf("  Heat flux range: %.4e ~ %.4e W/m^2", minimum(q_result), maximum(q_result)))

  # Forward verification ---------------------------------------------------
  println("\n[4/6] Forward verification (DHCP)")
  verify_start = time()

  wk_verify = IHCP_CGM.WorkBuffers(ni + 2, nj + 2, nk + 2)

  T_verify = IHCP_CGM.solve_dhcp!(
    T_init,
    q_result,
    wk_verify,
    nt,
    rho,
    cp_coeffs,
    k_coeffs,
    dx,
    dy,
    Z,
    dZ_guard,
    dt;
    rtol=rtol_dhcp,
    maxiter=maxiter_dhcp
  )

  elapsed_verify = time() - verify_start

  # Residual analysis
  T_bottom = @view T_verify[:, :, 1, :]  # Bottom surface (k=1)
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))

  println(@sprintf("  Verification time: %.2f s", elapsed_verify))
  println(@sprintf("  RMS residual: %.4e K", rms_error))
  println(@sprintf("  Max residual: %.4e K", max_error))

  # Save results -----------------------------------------------------------
  println("\n[5/6] Saving outputs")
  results_dir = joinpath(PROJECT_ROOT, "shared", "results")
  mkpath(results_dir)

  # Convert to Python memory order for comparison
  q_python = permutedims(q_result, (3, 1, 2))              # (nt-1, ni, nj)
  T_python = permutedims(T_verify, (4, 1, 2, 3))           # (nt, ni, nj, nk)
  Y_obs_python_out = permutedims(Y_obs, (3, 1, 2))        # (nt, ni, nj)

  q_path = joinpath(results_dir, "julia_100steps_cgm1_q.npy")
  T_path = joinpath(results_dir, "julia_100steps_cgm1_T.npy")
  npzwrite(q_path, q_python)
  npzwrite(T_path, T_python)

  # Bundle metadata
  npz_payload = Dict(
    "Y_obs" => Y_obs,
    "T_init" => T_init,
    "dt" => dt,
    "dx" => dx,
    "dy" => dy,
    "dz" => dz,
    "dz_b" => dz_bottom,
    "dz_t" => dz_top,
    "rho" => rho,
    "cp_coeffs" => cp_coeffs,
    "k_coeffs" => k_coeffs,
    "window_size" => window_size,
    "overlap" => overlap,
    "cgm_iteration" => cgm_iteration,
    "q_result" => q_result,
    "T_verify" => T_verify,
    "elapsed_time_cgm" => cgm_elapsed,
    "elapsed_time_verify" => elapsed_verify,
    "rms_error" => rms_error,
    "max_error" => max_error,
    "n_windows" => length(windows_info),
    "timestamp" => Dates.format(Dates.now(), Dates.ISODateTime)
  )

  npz_path = joinpath(results_dir, "julia_100steps_cgm1.npz")
  npzwrite(npz_path, npz_payload)

  println("\nSaved results:")
  println("  heat flux     → $q_path")
  println("  temperature   → $T_path")
  println("  bundle (.npz) → $npz_path")

  # Summary ----------------------------------------------------------------
  println("\n[6/6] Summary")
  total_elapsed = time() - total_start
  println(@sprintf("  Total runtime: %.2f s", total_elapsed))
  println(@sprintf("  CGM share: %.1f%%", 100 * cgm_elapsed / total_elapsed))
  println("="^80)

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
