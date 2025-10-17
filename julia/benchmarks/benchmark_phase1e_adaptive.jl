#!/usr/bin/env julia
"""
Phase 1-E Adaptive Tolerance Benchmark (adaptive_tol = true)

Phase 1-Bの最適パラメータ（dhcp_extrapolation=:quadratic, adjoint_residual_scale=0.5）に加えて、
適応的収束判定を有効化（adaptive_tol=true）したベンチマーク。
"""

using Dates
using Printf
using Statistics
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))
const RESULT_PATH = joinpath(PROJECT_ROOT, "shared", "results", "phase1e_adaptive.npz")

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

function build_initial_temperature(first_frame::AbstractArray{<:Real,2}, nk::Int)
  base = Float64.(first_frame)
  T_init = repeat(base, 1, 1, nk)
  return Array{Float64,3}(T_init)
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("Phase 1-E Adaptive Tolerance Benchmark (adaptive_tol = true)")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  ensure_single_thread()

  total_start = time()

  # Problem definition -----------------------------------------------------
  nt = NT
  dt = DT
  ni, nj, nk = 80, 100, 20
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  println("\n[1/5] Grid and material parameters")
  println(@sprintf("  grid: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  spacing: dx=%.3e m, dy=%.3e m", dx, dy))
  println(@sprintf("  time steps: nt=%d, dt=%.3e s", nt, dt))

  dz, dz_bottom, dz_top, z_faces = build_z_grid(nk, Lz, stretch_factor)

  # ガイドセル込み格子を生成
  z_centers, dz_grid = generate_guard_cell_grid(nk, dz, z_faces)

  rho, cp_coeffs, k_coeffs = load_material_properties()
  println(@sprintf("  rho (reference): %.6f kg/m^3", rho))

  # Input preparation ------------------------------------------------------
  println("\n[2/5] Loading measurement data")
  Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)x$(nj_file), expected $(ni)x$(nj)")
  end

  # Phase A: メモリレイアウト最適化 - (nt,ni,nj) → (ni,nj,nt)
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))  # (nt, ni, nj) → (ni, nj, nt)
  Y_obs_python = nothing

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))

  # Phase A: メモリレイアウト最適化 - (nt-1,ni,nj) → (ni,nj,nt-1)
  q_init = zeros(Float64, ni, nj, nt - 1)
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  # CGM solve --------------------------------------------------------------
  println("\n[3/5] Running CGM solver (single window, 1 iteration)")
  println("  Phase 1-B最適パラメータ適用 + adaptive_tol=true（Phase 1-E）")
  solve_start = time()

  # Phase 1-B最適パラメータ + adaptive_tol=true（Phase 1-E）
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
    # Phase 1-B最適化: 初期推定値改善（16.56%改善）
    dhcp_extrapolation = :quadratic,          # DHCP二次外挿（-11.6%）
    adjoint_initial_strategy = :previous,      # Adjoint初期値戦略（デフォルト戦略）
    adjoint_residual_scale = 0.5,             # Adjoint残差スケール（追加-5.3%）
    # Phase 1-E: 適応的収束判定（有効）
    adaptive_tol = true,                       # 適応的収束判定有効化
    tol_min_dhcp = 1.0e-7,                     # DHCP最小tol（推奨）
    tol_max_dhcp = 1.0e-5,                     # DHCP最大tol
    tol_min_adjoint = 1.0e-9,                  # Adjoint最小tol（推奨）
    tol_max_adjoint = 1.0e-5,                  # Adjoint最大tol
    tol_min_sensitivity = 1.0e-7,              # Sensitivity最小tol（推奨）
    tol_max_sensitivity = 1.0e-5,              # Sensitivity最大tol
    kappa_theta = 0.2,                         # 相対変化量スケール係数（推奨）
    warmup_steps = 3                           # 固定tol使用期間（推奨）
  )

  q_result, T_result, J_history = solve_cgm!(
    T_init,
    Y_obs,
    q_init,
    work,  # WorkBuffers追加
    dx,
    dy,
    z_centers,  # ガイドセルグリッド追加
    dz_grid,
    dt,
    rho,
    cp_coeffs,
    k_coeffs;
    params = cgm_params
  )

  cgm_elapsed = time() - solve_start
  println(@sprintf("  CGM elapsed time: %.2f s", cgm_elapsed))
  println(@sprintf("  iterations completed: %d", length(J_history)))

  # Residual diagnostics ---------------------------------------------------
  println("\n[4/5] Residual analysis")
  # Phase A: メモリレイアウト最適化 - T_result[ni,nj,nk,nt] の底面(k=1)
  T_bottom = @view T_result[:, :, 1, :]  # (ni, nj, nt)
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))
  println(@sprintf("  RMS residual: %.4e K", rms_error))
  println(@sprintf("  Max residual: %.4e K", max_error))
  println(@sprintf("  Heat flux range: %.4e ~ %.4e W/m^2", minimum(q_result), maximum(q_result)))

  # Save results -----------------------------------------------------------
  println("\n[5/5] Saving outputs")
  mkpath(dirname(RESULT_PATH))

  # Phase A: メモリレイアウト最適化後の形状
  # q_result: (ni, nj, nt-1), T_result: (ni, nj, nk, nt), Y_obs: (ni, nj, nt)
  metadata = Dict(
    "T" => T_result,
    "q" => q_result,
    "Y_obs" => Y_obs,
    "T_init" => T_init,
    "dt" => dt,
    "dx" => dx,
    "dy" => dy,
    "dz" => dz,
    "dz_bottom" => dz_bottom,
    "dz_top" => dz_top,
    "z_faces" => z_faces,
    "rho" => rho,
    "cp_coeffs" => cp_coeffs,
    "k_coeffs" => k_coeffs,
    "J_history" => collect(J_history),
    "elapsed_cgm" => cgm_elapsed,
    "rms_error" => rms_error,
    "max_error" => max_error,
  )

  npzwrite(RESULT_PATH, metadata)
  println("  saved bundle (.npz) → $RESULT_PATH")

  total_elapsed = time() - total_start
  println("\nSummary")
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
