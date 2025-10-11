#!/usr/bin/env julia
"""
Julia 10-step medium-size CGM performance test.

中規模グリッド(10×10×10)でフルサイズとの性能比較。
合成データを使用してCGM性能を測定。
"""

using Dates
using Printf
using Statistics

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))

const NT = 10
const DT = 1.0e-3

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

function ensure_single_thread()
  nthreads = Threads.nthreads()
  if nthreads != 1
    @warn "JULIA_NUM_THREADS should be 1 for this validation run" current = nthreads
  end
end

function build_z_grid(nk::Int, Lz::Float64, stretch_factor::Float64)
  z_faces_normalized = range(1.0, stop = 0.0, length = nk + 1)
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

function load_material_properties()
  cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]
  rho = 7823.493962874829  # 定数として扱う
  return rho, cp_coeffs, k_coeffs
end

function generate_synthetic_data(ni::Int, nj::Int, nk::Int, nt::Int)
  """
  合成温度データ生成
  - 初期温度: 300K
  - 時間変化: 線形加熱（0.1K/step）+ ノイズ
  """
  T_init = fill(300.0, ni, nj, nk)
  Y_obs = zeros(Float64, ni, nj, nt)

  for t in 1:nt
    base_temp = 300.0 + 0.1 * (t - 1)
    Y_obs[:, :, t] = base_temp .+ 0.01 .* randn(ni, nj)
  end

  return T_init, Y_obs
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("Julia 10-step medium-size CGM performance test")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  ensure_single_thread()

  total_start = time()

  # Problem definition (中規模: 10×10×10) -----------------------------------
  nt = NT
  dt = DT
  ni, nj, nk = 10, 10, 10  # ←中規模グリッド
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  println("\n[1/5] Grid and material parameters")
  println(@sprintf("  grid: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  spacing: dx=%.3e m, dy=%.3e m", dx, dy))
  println(@sprintf("  time steps: nt=%d, dt=%.3e s", nt, dt))

  dz, dz_bottom, dz_top, z_faces = build_z_grid(nk, Lz, stretch_factor)
  Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

  rho, cp_coeffs, k_coeffs = load_material_properties()
  println(@sprintf("  rho (reference): %.6f kg/m^3", rho))

  # Input preparation (合成データ) -------------------------------------------
  println("\n[2/5] Generating synthetic data")
  T_init, Y_obs = generate_synthetic_data(ni, nj, nk, nt)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))
  println(@sprintf("  observation temperature range: %.2f ~ %.2f K", minimum(Y_obs), maximum(Y_obs)))

  q_init = zeros(Float64, ni, nj, nt - 1)
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  # CGM solve ----------------------------------------------------------------
  println("\n[3/5] Running CGM solver (single window, 1 iteration)")
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
    verbose = true
  )

  q_result, T_result, J_history = solve_cgm!(
    T_init,
    Y_obs,
    q_init,
    work,
    dx,
    dy,
    Z,
    ΔZ,
    dt,
    rho,
    cp_coeffs,
    k_coeffs;
    params = cgm_params
  )

  cgm_elapsed = time() - solve_start
  println(@sprintf("  CGM elapsed time: %.2f s", cgm_elapsed))
  println(@sprintf("  iterations completed: %d", length(J_history)))

  # Residual diagnostics -----------------------------------------------------
  println("\n[4/5] Residual analysis")
  T_bottom = @view T_result[:, :, 1, :]
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))
  println(@sprintf("  RMS residual: %.4e K", rms_error))
  println(@sprintf("  Max residual: %.4e K", max_error))
  println(@sprintf("  Heat flux range: %.4e ~ %.4e W/m^2", minimum(q_result), maximum(q_result)))

  # Summary ------------------------------------------------------------------
  println("\n[5/5] Performance summary")
  total_elapsed = time() - total_start
  println(@sprintf("  Total runtime: %.2f s", total_elapsed))
  println(@sprintf("  CGM share: %.1f%%", 100 * cgm_elapsed / total_elapsed))
  println(@sprintf("  Grid size: %d cells (ni=%d, nj=%d, nk=%d)", ni*nj*nk, ni, nj, nk))
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
