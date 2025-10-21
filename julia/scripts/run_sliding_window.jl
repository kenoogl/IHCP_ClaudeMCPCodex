#!/usr/bin/env julia
"""
Julia sliding-window CGM driver.

This script mirrors the Python implementation in
`python/validation/run_sliding_window_validation.py`, while reusing
the solver setup proven in `run_10steps_fullsize_test.jl`.

Default configuration:
  • Total time steps (nt): 300
  • Window size: 71 heat-flux steps (needs 72 temperature samples)
  • Overlap: 17 steps
  • dt: 1.0e-3 s
  • Grid: 80 × 100 × 20 (same as run_10steps_fullsize_test.jl)

Outputs:
  shared/results/julia_sliding_window_cgm{N}.npz
    where {N} is the CGM maximum iteration count.
"""

using Dates
using Printf
using Statistics
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))
const RESULT_DIR = joinpath(PROJECT_ROOT, "shared", "results")

const DEFAULT_NT = 300
const DEFAULT_WINDOW_SIZE = 71
const DEFAULT_OVERLAP = 17
const DEFAULT_DT = 1.0e-3
const DEFAULT_SOLVER = :pbicgstab
const DEFAULT_PRECOND = :jacobi
const DEFAULT_CGM_MAX_ITER = 3
const DEFAULT_Q_INIT = 0.0

# ---------------------------------------------------------------------------
# Utility helpers (reused from run_10steps_fullsize_test.jl)
# ---------------------------------------------------------------------------

function parse_command_line_args()
  solver_type = DEFAULT_SOLVER
  precond_type = DEFAULT_PRECOND
  cgm_iter = DEFAULT_CGM_MAX_ITER
  nt = DEFAULT_NT
  window_size = DEFAULT_WINDOW_SIZE
  overlap = DEFAULT_OVERLAP
  dt = DEFAULT_DT
  q_init_value = DEFAULT_Q_INIT

  i = 1
  while i <= length(ARGS)
    arg = ARGS[i]
    if arg == "--help" || arg == "-h"
      println("""
      Usage: julia run_sliding_window.jl [OPTIONS]

      Options:
        --solver SOLVER        DHCP/Adjoint solver (pbicgstab | pcg) [default: pbicgstab]
        --precond PRECOND      Preconditioner (jacobi | gs | none)    [default: jacobi]
        --cgm-iter N           CGM maximum iterations                 [default: 3]
        --nt N                 Number of time steps to use            [default: 300]
        --window N             Sliding-window size (heat-flux steps)  [default: 71]
        --overlap N            Overlap between windows (steps)        [default: 17]
        --dt VALUE             Time step size in seconds              [default: 0.001]
        --q-init VALUE         Initial heat-flux guess (W/m^2)        [default: 0.0]
        -h, --help             Show this help message and exit

      Examples:
        julia run_sliding_window.jl
        julia run_sliding_window.jl --cgm-iter 200 --solver pcg --precond gs
      """)
      exit(0)
    elseif arg == "--solver"
      i += 1
      i > length(ARGS) && error("--solver requires an argument")
      solver_str = lowercase(ARGS[i])
      if solver_str == "pbicgstab"
        solver_type = :pbicgstab
      elseif solver_str == "pcg"
        solver_type = :pcg
      else
        error("Unknown solver type: $(ARGS[i])")
      end
    elseif arg == "--precond"
      i += 1
      i > length(ARGS) && error("--precond requires an argument")
      precond_str = lowercase(ARGS[i])
      if precond_str == "jacobi"
        precond_type = :jacobi
      elseif precond_str == "gs"
        precond_type = :gs
      elseif precond_str == "diagonal"
        precond_type = :diagonal
      elseif precond_str == "none"
        precond_type = :none
      else
        error("Unknown preconditioner: $(ARGS[i])")
      end
    elseif arg == "--cgm-iter"
      i += 1
      i > length(ARGS) && error("--cgm-iter requires an argument")
      cgm_iter = parse(Int, ARGS[i])
    elseif arg == "--nt"
      i += 1
      i > length(ARGS) && error("--nt requires an argument")
      nt = parse(Int, ARGS[i])
    elseif arg == "--window"
      i += 1
      i > length(ARGS) && error("--window requires an argument")
      window_size = parse(Int, ARGS[i])
    elseif arg == "--overlap"
      i += 1
      i > length(ARGS) && error("--overlap requires an argument")
      overlap = parse(Int, ARGS[i])
    elseif arg == "--dt"
      i += 1
      i > length(ARGS) && error("--dt requires an argument")
      dt = parse(Float64, ARGS[i])
    elseif arg == "--q-init"
      i += 1
      i > length(ARGS) && error("--q-init requires an argument")
      q_init_value = parse(Float64, ARGS[i])
    else
      error("Unknown argument: $arg")
    end
    i += 1
  end

  if nt < 2
    error("nt must be >= 2 (got $nt)")
  end
  if window_size < 1
    error("window size must be >= 1 (got $window_size)")
  end
  if overlap < 0
    error("overlap must be >= 0 (got $overlap)")
  end
  if overlap >= window_size
    @warn "Overlap $(overlap) is >= window size $(window_size); results may degenerate"
  end
  if cgm_iter < 1
    error("cgm-iter must be >= 1 (got $cgm_iter)")
  end

  return (
    solver_type = solver_type,
    precond_type = precond_type,
    cgm_iter = cgm_iter,
    nt = nt,
    window_size = window_size,
    overlap = overlap,
    dt = dt,
    q_init_value = q_init_value
  )
end

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

  ΔZ = zeros(nk + 2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  z_centers = zeros(nk + 2)
  z_centers[1] = z_faces[1]
  z_centers[nk+2] = z_faces[nk+1]
  for k in 2:nk+1
    z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
  end

  return z_centers, ΔZ, z_faces, dz
end

function load_material_properties()
  cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]
  rho = 7823.493962874829
  return rho, cp_coeffs, k_coeffs
end

function load_measurement_subset(nt_request::Int)
  data_path = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
  if !isfile(data_path)
    error("Measurement file not found: $data_path")
  end

  println("  loading measurement data from: $data_path")
  flush(stdout)
  load_start = time()
  T_measure_full = npzread(data_path)
  nt_available = size(T_measure_full, 1)
  if nt_request > nt_available
    error("Requested $nt_request time steps but only $nt_available are available")
  end

  subset = Array(T_measure_full[1:nt_request, :, :])
  if eltype(subset) != Float64
    subset = Float64.(subset)
  end
  T_measure_full = nothing
  GC.gc()

  println("  measurement subset shape (Python order): $(size(subset))")
  println(@sprintf("  load time: %.2f s", time() - load_start))
  flush(stdout)

  return subset
end

function build_initial_temperature(first_frame::AbstractArray{<:Real,2}, nk::Int)
  base = Float64.(first_frame)
  T_init = repeat(base, 1, 1, nk)
  return Array{Float64,3}(T_init)
end

# ---------------------------------------------------------------------------
# Sliding-window execution
# ---------------------------------------------------------------------------

function run_sliding_window!(
  Y_obs::Array{Float64,3},
  T_init::Array{Float64,3},
  work::WorkBuffers,
  dx::Float64,
  dy::Float64,
  z_centers::Vector{Float64},
  ΔZ::Vector{Float64},
  dt::Float64,
  rho::Float64,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  window_size::Int,
  overlap::Int,
  q_init_value::Float64,
  cgm_params::NamedTuple
)
  ni, nj, nt = size(Y_obs)
  nk = size(T_init, 3)
  total_flux_steps = nt - 1

  q_global = zeros(Float64, ni, nj, total_flux_steps)
  prev_flux_end = 0
  prev_q_win = nothing
  current_T_start = copy(T_init)
  final_T = similar(T_init)

  window_ids = Int[]
  window_starts = Int[]
  window_ends = Int[]
  window_lengths = Int[]
  window_overlaps = Int[]
  window_iterations = Int[]
  window_J_final = Float64[]
  window_elapsed = Float64[]
  window_q_min = Float64[]
  window_q_max = Float64[]
  window_histories = Vector{Vector{Float64}}()

  total_start = time()
  cgm_time_accum = 0.0
  window_id = 1

  println("\n=== Sliding-window CGM ===")
  println(@sprintf("  total time steps: %d", nt))
  println(@sprintf("  flux steps: %d", total_flux_steps))
  println(@sprintf("  window size: %d", window_size))
  println(@sprintf("  overlap: %d", overlap))
  flush(stdout)

  while prev_flux_end < total_flux_steps
    start_idx = if window_id == 1
      0
    else
      # mirror Python advance: next start = previous start + (length - overlap)
      last_start = window_starts[end]
      last_length = window_lengths[end]
      last_start + max(1, last_length - overlap)
    end

    if start_idx >= total_flux_steps
      break
    end

    max_L = min(window_size, total_flux_steps - start_idx)
    end_idx = start_idx + max_L
    obs_start = start_idx + 1
    obs_end = start_idx + max_L + 1

    @views Y_obs_win = Y_obs[:, :, obs_start:obs_end]

    q_init_win = zeros(Float64, ni, nj, max_L)
    if window_id == 1 || prev_q_win === nothing
      q_init_win .= q_init_value
      println(@sprintf("\n[Window %d] range=[%d,%d] length=%d", window_id, start_idx, end_idx, max_L))
      println(@sprintf("  initial heat flux: constant %.3e W/m^2", q_init_value))
    else
      L_overlap = min(overlap, max_L, size(prev_q_win, 3))
      if L_overlap > 0
        @views q_init_win[:, :, 1:L_overlap] .= prev_q_win[:, :, end-L_overlap+1:end]
      end
      if L_overlap < max_L
        edge = @view prev_q_win[:, :, end]
        for t in L_overlap+1:max_L
          @views q_init_win[:, :, t] .= edge
        end
      end
      println(@sprintf("\n[Window %d] range=[%d,%d] length=%d", window_id, start_idx, end_idx, max_L))
      println(@sprintf("  initial heat flux: inherited (overlap=%d)", L_overlap))
    end

    # Prepare temperature initial condition
    T_start = copy(current_T_start)

    window_time_start = time()
    q_win, T_win, J_hist = solve_cgm!(
      T_start,
      Y_obs_win,
      q_init_win,
      work,
      dx,
      dy,
      z_centers,
      ΔZ,
      dt,
      rho,
      cp_coeffs,
      k_coeffs;
      params = cgm_params
    )
    window_elapsed_time = time() - window_time_start
    cgm_time_accum += window_elapsed_time

    n_iters = length(J_hist)
    J_final = n_iters > 0 ? J_hist[end] : NaN

    # Heat-flux stitching with 0.5/0.5 averaging in overlap
    flux_start = start_idx + 1
    flux_end = start_idx + max_L
    overlap_steps = max(0, prev_flux_end - flux_start + 1)
    overlap_steps = min(overlap_steps, max_L)

    if overlap_steps > 0
      @views begin
        old_slice = q_global[:, :, flux_start:flux_start + overlap_steps - 1]
        new_slice = q_win[:, :, 1:overlap_steps]
        old_slice .= 0.5 .* old_slice .+ 0.5 .* new_slice
      end
    end

    if overlap_steps < max_L
      assign_start = flux_start + overlap_steps
      assign_end = flux_end
      @views q_global[:, :, assign_start:assign_end] .= q_win[:, :, overlap_steps+1:max_L]
    end

    prev_flux_end = max(prev_flux_end, flux_end)
    prev_q_win = copy(q_win)
    @views current_T_start = copy(T_win[:, :, :, end])
    final_T .= current_T_start

    q_min = minimum(q_win)
    q_max = maximum(q_win)

    println(@sprintf("  CGM iterations: %d", n_iters))
    println(@sprintf("  objective: %.6e", J_final))
    println(@sprintf("  elapsed: %.2f s", window_elapsed_time))
    println(@sprintf("  q range: [%.3e, %.3e] W/m^2", q_min, q_max))
    flush(stdout)

    push!(window_ids, window_id)
    push!(window_starts, start_idx)
    push!(window_ends, end_idx)
    push!(window_lengths, max_L)
    push!(window_overlaps, overlap_steps)
    push!(window_iterations, n_iters)
    push!(window_J_final, J_final)
    push!(window_elapsed, window_elapsed_time)
    push!(window_q_min, q_min)
    push!(window_q_max, q_max)
    push!(window_histories, collect(J_hist))

    window_id += 1
  end

  total_elapsed = time() - total_start

  if prev_flux_end < total_flux_steps
    @warn "Sliding-window loop ended before covering all flux steps" covered = prev_flux_end total = total_flux_steps
  end

  summary = (
    window_ids = window_ids,
    window_starts = window_starts,
    window_ends = window_ends,
    window_lengths = window_lengths,
    window_overlaps = window_overlaps,
    window_iterations = window_iterations,
    window_J_final = window_J_final,
    window_elapsed = window_elapsed,
    window_q_min = window_q_min,
    window_q_max = window_q_max,
    window_histories = window_histories,
    cgm_elapsed = cgm_time_accum,
    total_elapsed = total_elapsed
  )

  return q_global, final_T, summary
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("Julia sliding-window CGM runner")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  flush(stdout)

  ensure_single_thread()
  cfg = parse_command_line_args()

  println("\n[Configuration]")
  println("  Solver: $(cfg.solver_type)")
  println("  Preconditioner: $(cfg.precond_type)")
  println("  CGM max iter: $(cfg.cgm_iter)")
  println("  nt: $(cfg.nt)")
  println("  window size: $(cfg.window_size)")
  println("  overlap: $(cfg.overlap)")
  println(@sprintf("  dt: %.3e s", cfg.dt))
  println(@sprintf("  q_init: %.3e W/m^2", cfg.q_init_value))
  flush(stdout)

  total_start = time()

  # Problem definition
  ni, nj, nk = 80, 100, 20
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  println("\n[1/5] Grid and material parameters")
  println(@sprintf("  grid: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  spacing: dx=%.3e m, dy=%.3e m", dx, dy))
  println(@sprintf("  time steps: nt=%d, dt=%.3e s", cfg.nt, cfg.dt))
  flush(stdout)

  z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)
  rho, cp_coeffs, k_coeffs = load_material_properties()
  println(@sprintf("  rho: %.6f kg/m^3", rho))
  flush(stdout)

  println("\n[2/5] Loading measurement data")
  flush(stdout)
  Y_obs_python = load_measurement_subset(cfg.nt)
  ni_file = size(Y_obs_python, 2)
  nj_file = size(Y_obs_python, 3)
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)x$(nj_file), expected $(ni)x$(nj)")
  end

  # (nt, ni, nj) -> (ni, nj, nt)
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))
  Y_obs_python = nothing
  GC.gc()

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))
  flush(stdout)

  println("\n[3/5] Preparing solvers")
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  cgm_params = (
    max_iter = cfg.cgm_iter,
    sigma = 1.8,
    rtol_dhcp = 1.0e-6,
    rtol_adjoint = 1.0e-8,
    maxiter_cg = 20000,
    dire_reset_every = 5,
    eps = 1.0e-12,
    min_iter = min(cfg.cgm_iter, 10),
    P = 10,
    eta = 1.0e-4,
    beta_max = 1.0e8,
    verbose = false,
    dhcp_extrapolation = :quadratic,
    adjoint_residual_scale = 0.5,
    dhcp_solver = cfg.solver_type,
    dhcp_smoother = cfg.precond_type,
    adjoint_solver = cfg.solver_type,
    adjoint_smoother = cfg.precond_type,
    sensitivity_solver = cfg.solver_type,
    sensitivity_smoother = cfg.precond_type
  )

  println("\n[4/5] Running sliding-window CGM")
  flush(stdout)
  q_global, T_final, summary = run_sliding_window!(
    Y_obs,
    T_init,
    work,
    dx,
    dy,
    z_centers,
    ΔZ,
    cfg.dt,
    rho,
    cp_coeffs,
    k_coeffs,
    cfg.window_size,
    cfg.overlap,
    cfg.q_init_value,
    cgm_params
  )

  q_range = extrema(q_global)
  println("\n[5/5] Saving outputs")
  flush(stdout)
  mkpath(RESULT_DIR)
  result_path = joinpath(RESULT_DIR, @sprintf("julia_sliding_window_cgm%d.npz", cfg.cgm_iter))

  metadata = Dict(
    "q_global" => q_global,
    "T_final" => T_final,
    "dt" => cfg.dt,
    "dx" => dx,
    "dy" => dy,
    "dz" => dz,
    "z_faces" => z_faces,
    "z_centers" => z_centers,
    "ΔZ" => ΔZ,
    "rho" => rho,
    "cp_coeffs" => cp_coeffs,
    "k_coeffs" => k_coeffs,
    "nt" => cfg.nt,
    "window_size" => cfg.window_size,
    "overlap" => cfg.overlap,
    "cgm_max_iter" => cfg.cgm_iter,
    "q_init_value" => cfg.q_init_value,
    "window_ids" => collect(summary.window_ids),
    "window_start_idxs" => collect(summary.window_starts),
    "window_end_idxs" => collect(summary.window_ends),
    "window_lengths" => collect(summary.window_lengths),
    "window_overlaps" => collect(summary.window_overlaps),
    "window_iterations" => collect(summary.window_iterations),
    "window_J_final" => collect(summary.window_J_final),
    "window_elapsed" => collect(summary.window_elapsed),
    "window_q_min" => collect(summary.window_q_min),
    "window_q_max" => collect(summary.window_q_max),
    "cgm_elapsed" => summary.cgm_elapsed,
    "total_elapsed" => summary.total_elapsed,
    "q_min" => q_range[1],
    "q_max" => q_range[2]
  )

  npzwrite(result_path, metadata)

  metadata_txt_path = replace(result_path, ".npz" => "_metadata.txt")
  open(metadata_txt_path, "w") do io
    println(io, "solver_type: $(String(cfg.solver_type))")
    println(io, "precond_type: $(String(cfg.precond_type))")
    println(io, "cgm_max_iter: $(cfg.cgm_iter)")
    println(io, "nt: $(cfg.nt)")
    println(io, "window_size: $(cfg.window_size)")
    println(io, "overlap: $(cfg.overlap)")
    println(io, @sprintf("dt: %.6e", cfg.dt))
    println(io, @sprintf("q_init: %.6e", cfg.q_init_value))
    println(io, @sprintf("q_min: %.6e", q_range[1]))
    println(io, @sprintf("q_max: %.6e", q_range[2]))
    println(io, @sprintf("cgm_elapsed: %.2f", summary.cgm_elapsed))
    println(io, @sprintf("total_elapsed: %.2f", summary.total_elapsed))
  end

  println(@sprintf("  saved bundle (.npz) → %s", result_path))
  println(@sprintf("  saved metadata → %s", metadata_txt_path))
  flush(stdout)

  total_elapsed = time() - total_start
  total_windows = length(summary.window_ids)
  total_iters = sum(summary.window_iterations)

  println("\nSummary")
  println(@sprintf("  Total runtime: %.2f s", total_elapsed))
  println(@sprintf("  Sliding-window elapsed: %.2f s", summary.total_elapsed))
  println(@sprintf("  Accumulated CGM time: %.2f s", summary.cgm_elapsed))
  println(@sprintf("  Windows processed: %d", total_windows))
  println(@sprintf("  Total CGM iterations: %d", total_iters))
  println(@sprintf("  Heat-flux range: %.3e ~ %.3e W/m^2", q_range[1], q_range[2]))
  println("="^80)
  flush(stdout)

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
