#!/usr/bin/env julia
"""
Julia 100-step CGM validation runner.

This script reproduces the 100-step inverse heat conduction test documented in
`docs/comparison_test_params.md`. It loads the shared measurement data,
reuses the production CGM solver, and persists results that are directly
comparable with the Python reference implementation.
"""

using Dates
using Printf
using Statistics
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

"""Ensure the solver runs sequentially as required for the comparison test."""
function ensure_single_thread()
  nthreads = Threads.nthreads()
  if nthreads != 1
    @warn "JULIA_NUM_THREADS should be 1 for this validation run" current=nthreads
  end
end

"""
Construct the stretched z-direction grid used by the benchmark.
"""
function build_z_grid(nk::Int, Lz::Float64, stretch_factor::Float64)
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

  return dz, dz_bottom, dz_top, z_faces
end

"""
Load the SUS304 thermal properties and fit third-order polynomials.
"""
function load_material_properties()
  # Use the same polynomial coefficients as the Python reference
  # These were obtained from fitting metal_thermal_properties.csv
  rho_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

  # Use the same density value as Python reference
  rho = 7823.493962874829

  return rho, cp_coeffs, k_coeffs, rho_coeffs
end

"""
Load the first `nt` measurement frames (time-major in the source npy file).
"""
function load_measurement_subset(nt::Int)
  data_path = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
  if !isfile(data_path)
    error("Measurement file not found: $data_path")
  end

  println("  loading measurement data from: $data_path")
  load_start = time()
  T_measure_full = npzread(data_path)
  ni_full = size(T_measure_full, 2)
  nj_full = size(T_measure_full, 3)
  if size(T_measure_full, 1) < nt
    error("Measurement file only has $(size(T_measure_full, 1)) time steps; need $nt")
  end

  subset = Array(T_measure_full[1:nt, :, :])
  if eltype(subset) != Float64
    subset = Float64.(subset)
  end
  T_measure_full = nothing
  GC.gc()

  println("  measurement subset shape (Python order): $(size(subset))")
  println(@sprintf("  load time: %.2f s", time() - load_start))

  return subset, ni_full, nj_full
end

"""
Build the initial temperature volume from the first measurement frame.
"""
function build_initial_temperature(first_frame::AbstractArray{<:Real, 2}, nk::Int)
  base = Float64.(first_frame)
  T_init = repeat(base, 1, 1, nk)
  return Array{Float64, 3}(T_init)
end

"""
Persist results in both Julia-native and Python-friendly layouts.
"""
function save_results(result_dir::String;
    q_julia::Array{Float64,3},
    T_julia::Array{Float64,4},
    Y_obs_julia::Array{Float64,3},
    T_init::Array{Float64,3},
    dz::Vector{Float64},
    dz_bottom::Vector{Float64},
    dz_top::Vector{Float64},
    z_faces::Vector{Float64},
    params::Dict{String,Any})

  mkpath(result_dir)

  # Python expects time-major arrays; convert before writing .npy files.
  q_python = permutedims(q_julia, (3, 1, 2))              # (nt-1, ni, nj)
  T_python = permutedims(T_julia, (4, 1, 2, 3))           # (nt, ni, nj, nk)
  Y_obs_python = permutedims(Y_obs_julia, (3, 1, 2))      # (nt, ni, nj)

  q_path = joinpath(result_dir, "julia_100steps_cgm1_q.npy")
  T_path = joinpath(result_dir, "julia_100steps_cgm1_T.npy")
  npzwrite(q_path, q_python)
  npzwrite(T_path, T_python)

  # Bundle comprehensive metadata for downstream analysis.
  npz_payload = Dict(
    "Y_obs_julia" => Y_obs_julia,
    "Y_obs_python" => Y_obs_python,
    "T_init" => T_init,
    "q_result_julia" => q_julia,
    "q_result_python" => q_python,
    "T_result_julia" => T_julia,
    "T_result_python" => T_python,
    "dz" => dz,
    "dz_bottom" => dz_bottom,
    "dz_top" => dz_top,
    "z_faces" => z_faces,
  )

  merge!(npz_payload, params)

  npz_path = joinpath(result_dir, "julia_100steps_cgm1.npz")
  npzwrite(npz_path, npz_payload)

  println("\nSaved results:")
  println("  heat flux     → $q_path")
  println("  temperature   → $T_path")
  println("  bundle (.npz) → $npz_path")
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("Julia 100-step CGM validation")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  ensure_single_thread()

  total_start = time()

  # Problem definition -----------------------------------------------------
  nt = 100
  dt = 1.0e-3
  ni, nj, nk = 80, 100, 20
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  println("\n[1/6] Grid and material parameters")
  println(@sprintf("  grid: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  spacing: dx=%.3e m, dy=%.3e m", dx, dy))
  println(@sprintf("  time steps: nt=%d, dt=%.3e s", nt, dt))

  dz, dz_bottom, dz_top, z_faces = build_z_grid(nk, Lz, stretch_factor)
  Z, dZ_guard = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

  rho, cp_coeffs, k_coeffs, rho_coeffs = load_material_properties()
  println(@sprintf("  rho (at 498.15K): %.6f kg/m^3", rho))

  # Input preparation ------------------------------------------------------
  println("\n[2/6] Loading measurement data")
  Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)x$(nj_file), expected $(ni)x$(nj)")
  end

  Y_obs = permutedims(Y_obs_python, (2, 3, 1))  # Julia layout (ni, nj, nt)
  Y_obs_python = nothing

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))

  q_init = zeros(Float64, ni, nj, nt - 1)
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  # CGM solve --------------------------------------------------------------
  println("\n[3/6] Running CGM solver (sequential, 1 iteration for comparison)")
  solve_start = time()

  cgm_params = (
    max_iter = 1,
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
  )

  q_result, T_result, J_history = solve_cgm!(
    T_init,
    Y_obs,
    q_init,
    work,
    dx,
    dy,
    Z,
    dZ_guard,
    dt,
    rho,
    cp_coeffs,
    k_coeffs;
    params = cgm_params,
    par = "sequential"
  )

  cgm_elapsed = time() - solve_start
  println(@sprintf("  CGM elapsed time: %.2f s", cgm_elapsed))
  println(@sprintf("  iterations completed: %d", length(J_history)))

  # Residual diagnostics ---------------------------------------------------
  println("\n[4/6] Residual analysis")
  T_bottom = @view T_result[:, :, 1, :]
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))
  println(@sprintf("  RMS residual: %.4e K", rms_error))
  println(@sprintf("  Max residual: %.4e K", max_error))
  println(@sprintf("  Heat flux range: %.4e ~ %.4e W/m^2", minimum(q_result), maximum(q_result)))

  # Save results -----------------------------------------------------------
  println("\n[5/6] Saving outputs")
  results_dir = joinpath(PROJECT_ROOT, "shared", "results")

  metadata = Dict{String, Any}(
    "dt" => dt,
    "dx" => dx,
    "dy" => dy,
    "Z" => Z,
    "dZ_guard" => dZ_guard,
    "rho" => rho,
    "rho_coeffs" => rho_coeffs,
    "cp_coeffs" => cp_coeffs,
    "k_coeffs" => k_coeffs,
    "sigma" => cgm_params.sigma,
    "max_iter" => cgm_params.max_iter,
    "J_history" => collect(J_history),
    "elapsed_cgm" => cgm_elapsed,
    "rms_error" => rms_error,
    "max_error" => max_error,
    "timestamp" => Dates.format(Dates.now(), Dates.ISODateTime)
  )

  save_results(results_dir;
    q_julia = q_result,
    T_julia = T_result,
    Y_obs_julia = Y_obs,
    T_init = T_init,
    dz = dz,
    dz_bottom = dz_bottom,
    dz_top = dz_top,
    z_faces = z_faces,
    params = metadata
  )

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
