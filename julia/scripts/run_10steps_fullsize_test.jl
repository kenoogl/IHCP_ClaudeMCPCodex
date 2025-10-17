#!/usr/bin/env julia
"""
Julia 10-step full-size CGM validation runner.

This script mirrors ``python/validation/run_10steps_fullsize_test.py``.
It loads the shared measurement data, runs the CGM solver for a single
full window (nt = 10), and saves heat-flux and temperature results in a
Python-friendly ``.npz`` bundle for downstream comparison.

Usage:
  julia run_10steps_fullsize_test.jl [--solver SOLVER] [--precond PRECOND]

Options:
  --solver SOLVER     Solver type: pbicgstab (default), pcg
  --precond PRECOND   Preconditioner type: jacobi (default), gs, none
"""

using Dates
using Printf
using Statistics
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))
const RESULT_PATH = joinpath(PROJECT_ROOT, "shared", "results", "julia_10steps_fullsize.npz")

const NT = 10
const DT = 1.0e-3

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

"""
    parse_command_line_args() -> (solver_type, precond_type)

コマンドライン引数をパースしてソルバーと前処理の設定を取得

Returns:
  solver_type: ソルバータイプ（Symbol: :pbicgstab, :pcg）
  precond_type: 前処理タイプ（Symbol: :jacobi, :gs, :none）
"""
function parse_command_line_args()
  solver_type = :pbicgstab  # デフォルト
  precond_type = :jacobi    # デフォルト

  i = 1
  while i <= length(ARGS)
    if ARGS[i] == "--help" || ARGS[i] == "-h"
      println("""
      Usage: julia run_10steps_fullsize_test.jl [OPTIONS]

      Options:
        --solver SOLVER       Solver type (default: pbicgstab)
                              Available: pbicgstab, pcg
        --precond PRECOND     Preconditioner type (default: jacobi)
                              Available: jacobi, gs, none
        -h, --help            Show this help message and exit

      Examples:
        julia run_10steps_fullsize_test.jl
        julia run_10steps_fullsize_test.jl --solver pcg --precond gs
      """)
      exit(0)
    elseif ARGS[i] == "--solver"
      if i + 1 > length(ARGS)
        error("--solver requires an argument")
      end
      solver_str = lowercase(ARGS[i + 1])
      if solver_str == "pbicgstab"
        solver_type = :pbicgstab
      elseif solver_str == "pcg"
        solver_type = :pcg
      else
        error("Unknown solver type: $(ARGS[i + 1]). Available: pbicgstab, pcg")
      end
      i += 2
    elseif ARGS[i] == "--precond"
      if i + 1 > length(ARGS)
        error("--precond requires an argument")
      end
      precond_str = lowercase(ARGS[i + 1])
      if precond_str == "jacobi"
        precond_type = :jacobi
      elseif precond_str == "gs"
        precond_type = :gs
      elseif precond_str == "none"
        precond_type = :none
      else
        error("Unknown preconditioner type: $(ARGS[i + 1]). Available: jacobi, gs, none")
      end
      i += 2
    else
      error("Unknown argument: $(ARGS[i])")
    end
  end

  return solver_type, precond_type
end

function ensure_single_thread()
  nthreads = Threads.nthreads()
  if nthreads != 1
    @warn "JULIA_NUM_THREADS should be 1 for this validation run" current = nthreads
  end
end

"""
  build_z_grid(nk, Lz, stretch_factor) -> (z_centers, ΔZ, z_faces, dz)

Python版と同じ格子生成 + ガイドセル込み配列を直接生成

main.jlと同じ方式に統一（convert_to_guard_cell_gridは使用しない）

# Returns
- z_centers: セル中心座標（ガイドセル込み、nk+2個） [m]
- ΔZ: セル幅（ガイドセル込み、nk+2個） [m]
- z_faces: セル境界座標（nk+1個） [m]
- dz: 各層の厚さ（nk個、npz保存用） [m]
"""
function build_z_grid(nk::Int, Lz::Float64, stretch_factor::Float64)
  # Step 1: z_faces生成（Python版と同じ）
  z_faces_normalized = range(1.0, stop = 0.0, length = nk + 1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  # Step 2: dzの計算
  dz = diff(z_faces)

  # Step 3: ΔZ配列（ガイドセル込み、nk+2個）
  ΔZ = zeros(nk + 2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # Step 4: z_centers（ガイドセル込み、nk+2個）
  z_centers = zeros(nk + 2)

  # ガイドセル
  z_centers[1] = z_faces[1]        # 底面ガイドセル
  z_centers[nk+2] = z_faces[nk+1]  # 表面ガイドセル

  # すべての物理セル（k=2からnk+1）を統一した方法で計算
  # ゼロ除算を回避し、test_dhcp_convection.jlと整合性を確保
  for k in 2:nk+1
    z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
  end

  return z_centers, ΔZ, z_faces, dz
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
  flush(stdout)
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
  flush(stdout)

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
  println("Julia 10-step full-size CGM validation")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  flush(stdout)
  ensure_single_thread()

  # コマンドライン引数パース
  solver_type, precond_type = parse_command_line_args()
  println("\n[Configuration]")
  println("  Solver: $solver_type")
  println("  Preconditioner: $precond_type")
  flush(stdout)

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
  flush(stdout)

  # main.jlと同じ直接生成方式
  z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)

  rho, cp_coeffs, k_coeffs = load_material_properties()
  println(@sprintf("  rho (reference): %.6f kg/m^3", rho))
  flush(stdout)

  # Input preparation ------------------------------------------------------
  println("\n[2/5] Loading measurement data")
  flush(stdout)
  Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)x$(nj_file), expected $(ni)x$(nj)")
  end

  # Phase A: メモリレイアウト最適化 - (nt,ni,nj) → (ni,nj,nt)
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))  # (nt, ni, nj) → (ni, nj, nt)
  Y_obs_python = nothing

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))
  flush(stdout)

  # Phase A: メモリレイアウト最適化 - (nt-1,ni,nj) → (ni,nj,nt-1)
  q_init = zeros(Float64, ni, nj, nt - 1)
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  # CGM solve --------------------------------------------------------------
  println("\n[3/5] Running CGM solver (single window, 1 iteration)")
  flush(stdout)
  solve_start = time()

  # Phase 1-B最適パラメータを適用（16.96%改善達成）
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
    adjoint_residual_scale = 0.5,             # Adjoint残差スケール（追加-5.3%）
    # ソルバーと前処理の設定（コマンドライン引数から）
    dhcp_solver = solver_type,
    dhcp_smoother = precond_type,
    adjoint_solver = solver_type,
    adjoint_smoother = precond_type,
    sensitivity_solver = solver_type,
    sensitivity_smoother = precond_type
  )

  q_result, T_result, J_history = solve_cgm!(
    T_init,
    Y_obs,
    q_init,
    work,  # WorkBuffers追加
    dx,
    dy,
    z_centers,  # ガイドセルグリッド追加
    ΔZ,         # セル幅（ガイドセル込み）
    dt,
    rho,
    cp_coeffs,
    k_coeffs;
    params = cgm_params
  )

  cgm_elapsed = time() - solve_start
  println(@sprintf("  CGM elapsed time: %.2f s", cgm_elapsed))
  println(@sprintf("  iterations completed: %d", length(J_history)))
  flush(stdout)

  # Residual diagnostics ---------------------------------------------------
  println("\n[4/5] Residual analysis")
  flush(stdout)
  # Phase A: メモリレイアウト最適化 - T_result[ni,nj,nk,nt] の底面(k=1)
  T_bottom = @view T_result[:, :, 1, :]  # (ni, nj, nt)
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))
  println(@sprintf("  RMS residual: %.4e K", rms_error))
  println(@sprintf("  Max residual: %.4e K", max_error))
  println(@sprintf("  Heat flux range: %.4e ~ %.4e W/m^2", minimum(q_result), maximum(q_result)))
  flush(stdout)

  # Save results -----------------------------------------------------------
  println("\n[5/5] Saving outputs")
  flush(stdout)
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
    "dz" => dz,           # 各層の厚さ（nk個）
    "z_faces" => z_faces, # セル境界座標（nk+1個）
    "z_centers" => z_centers, # セル中心座標（ガイドセル込み、nk+2個）
    "ΔZ" => ΔZ,          # セル幅（ガイドセル込み、nk+2個）
    "rho" => rho,
    "cp_coeffs" => cp_coeffs,
    "k_coeffs" => k_coeffs,
    "J_history" => collect(J_history),
    "elapsed_cgm" => cgm_elapsed,
    "rms_error" => rms_error,
    "max_error" => max_error,
  )

  # ソルバー情報を別途テキストファイルに保存
  metadata_txt_path = replace(RESULT_PATH, ".npz" => "_metadata.txt")
  open(metadata_txt_path, "w") do io
    println(io, "solver_type: $(String(solver_type))")
    println(io, "precond_type: $(String(precond_type))")
  end

  npzwrite(RESULT_PATH, metadata)
  println("  saved bundle (.npz) → $RESULT_PATH")
  flush(stdout)

  total_elapsed = time() - total_start
  println("\nSummary")
  println(@sprintf("  Total runtime: %.2f s", total_elapsed))
  println(@sprintf("  CGM share: %.1f%%", 100 * cgm_elapsed / total_elapsed))
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
