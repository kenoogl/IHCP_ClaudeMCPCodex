#!/usr/bin/env julia
"""
Phase 1-B 追加ベンチマーク: DHCP二次外挿 + Adjoint残差0.1倍

設定:
- DHCP外挿: quadratic
- Adjoint残差スケール: 0.1
- Sensitivity外挿: none
- Adjoint初期値戦略: previous (デフォルト)

目的:
さらなる最適化パラメータの探索
"""

using Dates
using Printf
using Statistics
using NPZ
using JSON

include("../src/IHCP_CGM.jl")
using .IHCP_CGM

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))
const RESULT_DIR = joinpath(PROJECT_ROOT, "shared", "results")

const NT = 10
const DT = 1.0e-3

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

function ensure_single_thread()
  nthreads = Threads.nthreads()
  if nthreads != 1
    @warn "JULIA_NUM_THREADS should be 1 for this benchmark" current = nthreads
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
  rho = 7823.493962874829
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

  println(@sprintf("  load time: %.2f s", time() - load_start))

  return subset, ni_full, nj_full
end

function build_initial_temperature(first_frame::AbstractArray{<:Real,2}, nk::Int)
  base = Float64.(first_frame)
  T_init = repeat(base, 1, 1, nk)
  return Array{Float64,3}(T_init)
end

# ---------------------------------------------------------------------------
# CGM with detailed iteration logging
# ---------------------------------------------------------------------------

"""
CGM実行（詳細ログ付き）

solve_cgm!を呼び出して、各ソルバーの反復数を記録
"""
function run_cgm_with_logging(grid_params, material_params, Y_obs, T_init)
  ni, nj, nk, nt = grid_params.ni, grid_params.nj, grid_params.nk, grid_params.nt
  dx, dy, dt = grid_params.dx, grid_params.dy, grid_params.dt
  Z, ΔZ = grid_params.Z, grid_params.ΔZ
  rho, cp_coeffs, k_coeffs = material_params

  # 初期熱流束
  q_init = zeros(Float64, ni, nj, nt - 1)
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  # CGMパラメータ
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
    verbose = true,  # 詳細ログ出力
    dhcp_extrapolation = :quadratic,
    sensitivity_extrapolation = :none,
    adjoint_initial_strategy = :previous,
    adjoint_residual_scale = 0.1  # ← 0.1倍
  )

  # CGM実行
  solve_start = time()
  q_result, T_result, J_history = solve_cgm!(
    T_init,
    Y_obs,
    q_init,
    work,
    dx, dy,
    Z, ΔZ,
    dt,
    rho,
    cp_coeffs,
    k_coeffs;
    params = cgm_params
  )
  cgm_elapsed = time() - solve_start

  return q_result, T_result, J_history, cgm_elapsed
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("Phase 1-B 追加ベンチマーク: DHCP二次外挿 + Adjoint残差0.1倍")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  ensure_single_thread()

  total_start = time()

  # Problem definition
  nt = NT
  dt = DT
  ni, nj, nk = 80, 100, 20
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  println("\n[1/3] 初期化")
  println(@sprintf("  格子: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  時間ステップ: nt=%d, dt=%.3e s", nt, dt))

  dz, dz_bottom, dz_top, z_faces = build_z_grid(nk, Lz, stretch_factor)
  Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

  rho, cp_coeffs, k_coeffs = load_material_properties()

  # 測定データ読み込み
  println("\n[2/3] 測定データ読み込み")
  Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)x$(nj_file), expected $(ni)x$(nj)")
  end

  Y_obs = permutedims(Y_obs_python, (2, 3, 1))
  Y_obs_python = nothing
  GC.gc()

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)

  # グリッド・物性値パラメータ
  grid_params = (
    ni = ni, nj = nj, nk = nk, nt = nt,
    dx = dx, dy = dy, dt = dt,
    Z = Z, ΔZ = ΔZ
  )
  material_params = (rho, cp_coeffs, k_coeffs)

  # ベンチマーク実行
  println("\n[3/3] ベンチマーク実行")
  println("設定:")
  println("  DHCP外挿: quadratic")
  println("  Adjoint残差スケール: 0.1")
  println("  Sensitivity外挿: none")
  println("  Adjoint初期値戦略: previous")
  println("\n予想実行時間: 約15分")

  # CGM実行（ログ付き）
  q_result, T_result, J_history, cgm_elapsed = run_cgm_with_logging(
    grid_params, material_params, Y_obs, T_init
  )

  # 残差計算
  T_bottom = @view T_result[:, :, 1, :]
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))

  # 結果サマリー
  results = Dict(
    "scenario_name" => "dhcp_quadratic_adjoint_0.1",
    "description" => "DHCP二次外挿 + Adjoint残差0.1倍",
    "dhcp_extrapolation" => "quadratic",
    "sensitivity_extrapolation" => "none",
    "adjoint_initial_strategy" => "previous",
    "adjoint_residual_scale" => 0.1,
    "elapsed_time" => cgm_elapsed,
    "rms_error" => rms_error,
    "max_error" => max_error,
    "J_final" => length(J_history) > 0 ? J_history[end] : NaN,
    "heat_flux_min" => minimum(q_result),
    "heat_flux_max" => maximum(q_result),
    "heat_flux_mean" => mean(q_result)
  )

  println("\n結果サマリー:")
  println(@sprintf("  実行時間: %.2f 秒", cgm_elapsed))
  println(@sprintf("  RMS残差: %.4e K", rms_error))
  println(@sprintf("  最大残差: %.4e K", max_error))
  println(@sprintf("  目的関数: %.4e", results["J_final"]))
  println(@sprintf("  熱流束範囲: %.4e ~ %.4e W/m²", minimum(q_result), maximum(q_result)))

  # 結果保存
  timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
  result_file = joinpath(RESULT_DIR, "adjoint_scale_01_$(timestamp).json")
  mkpath(dirname(result_file))

  # ブランチ情報取得
  branch_name = try
    strip(read(`git rev-parse --abbrev-ref HEAD`, String))
  catch
    "unknown"
  end

  commit_hash = try
    strip(read(`git rev-parse --short HEAD`, String))
  catch
    "unknown"
  end

  output_data = Dict(
    "benchmark_info" => Dict(
      "date" => string(now()),
      "julia_version" => string(VERSION),
      "grid_size" => (ni, nj, nk),
      "time_steps" => nt,
      "dt" => dt,
      "branch" => branch_name,
      "commit" => commit_hash
    ),
    "scenario" => results
  )

  open(result_file, "w") do io
    JSON.print(io, output_data, 2)
  end
  println("\n✅ 結果保存: $result_file")

  # サマリー出力
  total_elapsed = time() - total_start
  println("\n" * "="^80)
  println("ベンチマーク完了")
  println("="^80)
  println(@sprintf("総実行時間: %.2f 秒 (%.2f 分)", total_elapsed, total_elapsed / 60))

  # 比較
  println("\n参考: 既存結果との比較")
  println("-"^80)
  println("ベースライン: 1,063.42秒")
  println("DHCP二次のみ: 940.50秒 (-11.6%)")
  println("DHCP二次 + Adjoint残差0.5倍: 890.47秒 (-16.3%)")
  println("本実験（DHCP二次 + Adjoint残差0.1倍）: $(cgm_elapsed)秒")
  improvement_vs_baseline = (1063.42 - cgm_elapsed) / 1063.42 * 100
  println(@sprintf("  vs ベースライン: %.2f%% 改善", improvement_vs_baseline))
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
