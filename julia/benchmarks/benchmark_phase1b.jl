#!/usr/bin/env julia
"""
Phase 1-B 初期推定値改善ベンチマーク

実験方針:
各改善要素を独立に評価するため、評価対象以外を従来法に固定

ベンチマークシナリオ（全6パターン）:
1. ベースライン: DHCP=none, Adjoint=previous, Sensitivity=none
2. DHCP線形外挿: DHCP=linear, Adjoint=previous, Sensitivity=none
3. DHCP二次外挿: DHCP=quadratic, Adjoint=previous, Sensitivity=none
4. Adjoint改善: DHCP=none, Adjoint=residual, Sensitivity=none
5. Sensitivity線形外挿: DHCP=none, Adjoint=previous, Sensitivity=linear
6. Sensitivity二次外挿: DHCP=none, Adjoint=previous, Sensitivity=quadratic

測定項目:
- 総実行時間
- DHCP平均CG反復数
- Adjoint平均CG反復数
- Sensitivity平均CG反復数
- RMS残差、最大残差
- 目的関数値
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

# ベンチマーク設定（6パターン）
const SCENARIOS = [
  # 1. ベースライン
  (
    name = "baseline",
    description = "ベースライン（Phase 0/1-A相当）",
    dhcp_extrapolation = :none,
    sensitivity_extrapolation = :none,
    adjoint_initial_strategy = :previous
  ),
  # 2. DHCP線形外挿
  (
    name = "dhcp_linear",
    description = "DHCP線形外挿（他は従来法）",
    dhcp_extrapolation = :linear,
    sensitivity_extrapolation = :none,
    adjoint_initial_strategy = :previous
  ),
  # 3. DHCP二次外挿
  (
    name = "dhcp_quadratic",
    description = "DHCP二次外挿（他は従来法）",
    dhcp_extrapolation = :quadratic,
    sensitivity_extrapolation = :none,
    adjoint_initial_strategy = :previous
  ),
  # 4. Adjoint改善
  (
    name = "adjoint_residual",
    description = "Adjoint初期値改善（他は従来法）",
    dhcp_extrapolation = :none,
    sensitivity_extrapolation = :none,
    adjoint_initial_strategy = :residual
  ),
  # 5. Sensitivity線形外挿
  (
    name = "sensitivity_linear",
    description = "Sensitivity線形外挿（他は従来法）",
    dhcp_extrapolation = :none,
    sensitivity_extrapolation = :linear,
    adjoint_initial_strategy = :previous
  ),
  # 6. Sensitivity二次外挿
  (
    name = "sensitivity_quadratic",
    description = "Sensitivity二次外挿（他は従来法）",
    dhcp_extrapolation = :none,
    sensitivity_extrapolation = :quadratic,
    adjoint_initial_strategy = :previous
  ),
]

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
function run_cgm_with_logging(scenario, grid_params, material_params, Y_obs, T_init)
  ni, nj, nk, nt = grid_params.ni, grid_params.nj, grid_params.nk, grid_params.nt
  dx, dy, dt = grid_params.dx, grid_params.dy, grid_params.dt
  z_centers, dz = grid_params.z_centers, grid_params.dz
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
    dhcp_extrapolation = scenario.dhcp_extrapolation,
    sensitivity_extrapolation = scenario.sensitivity_extrapolation,
    adjoint_initial_strategy = scenario.adjoint_initial_strategy
  )

  # CGM実行
  solve_start = time()
  q_result, T_result, J_history = solve_cgm!(
    T_init,
    Y_obs,
    q_init,
    work,
    dx, dy,
    z_centers, dz,
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
# Benchmark runner
# ---------------------------------------------------------------------------

function run_scenario(scenario, grid_params, material_params, Y_obs, T_init)
  println("\n" * "="^80)
  println("シナリオ: $(scenario.name)")
  println("説明: $(scenario.description)")
  println("設定:")
  println("  DHCP外挿: $(scenario.dhcp_extrapolation)")
  println("  Adjoint初期値: $(scenario.adjoint_initial_strategy)")
  println("  Sensitivity外挿: $(scenario.sensitivity_extrapolation)")
  println("="^80)

  ni, nj, nk, nt = grid_params.ni, grid_params.nj, grid_params.nk, grid_params.nt

  # CGM実行（ログ付き）
  q_result, T_result, J_history, cgm_elapsed = run_cgm_with_logging(
    scenario, grid_params, material_params, Y_obs, T_init
  )

  # 残差計算
  T_bottom = @view T_result[:, :, 1, :]
  residual = T_bottom .- Y_obs
  rms_error = sqrt(mean(residual .^ 2))
  max_error = maximum(abs.(residual))

  # 結果サマリー
  results = Dict(
    "scenario_name" => scenario.name,
    "description" => scenario.description,
    "dhcp_extrapolation" => String(scenario.dhcp_extrapolation),
    "sensitivity_extrapolation" => String(scenario.sensitivity_extrapolation),
    "adjoint_initial_strategy" => String(scenario.adjoint_initial_strategy),
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

  return results
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("Phase 1-B 初期推定値改善ベンチマーク")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  println("実験計画: 6シナリオ（ベースライン + 各改善要素を独立評価）")
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

  # ガイドセル込み格子を生成
  z_centers, dz_grid = generate_guard_cell_grid(nk, dz, z_faces)

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
    z_centers = z_centers, dz = dz_grid
  )
  material_params = (rho, cp_coeffs, k_coeffs)

  # ベンチマーク実行
  println("\n[3/3] ベンチマーク実行")
  println("シナリオ数: $(length(SCENARIOS))")
  println("予想実行時間: 約$(length(SCENARIOS) * 18)分")

  all_results = []
  for (idx, scenario) in enumerate(SCENARIOS)
    println("\n" * "█"^80)
    println("進行状況: [$idx/$(length(SCENARIOS))]")
    println("█"^80)

    result = run_scenario(scenario, grid_params, material_params, Y_obs, T_init)
    push!(all_results, result)

    GC.gc()  # メモリ解放

    # 中間進捗
    elapsed_so_far = time() - total_start
    avg_time_per_scenario = elapsed_so_far / idx
    remaining_scenarios = length(SCENARIOS) - idx
    estimated_remaining = avg_time_per_scenario * remaining_scenarios
    println(@sprintf("\n進捗: %.1f%% 完了、残り約%.1f分",
      100.0 * idx / length(SCENARIOS), estimated_remaining / 60))
  end

  # 結果保存
  timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
  result_file = joinpath(RESULT_DIR, "phase1b_benchmark_$(timestamp).json")
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
    "scenarios" => all_results
  )

  open(result_file, "w") do io
    JSON.print(io, output_data, 2)
  end
  println("\n✅ 結果保存: $result_file")

  # サマリー出力
  total_elapsed = time() - total_start
  println("\n" * "="^80)
  println("Phase 1-B ベンチマーク完了")
  println("="^80)
  println(@sprintf("総実行時間: %.2f 秒 (%.2f 分)", total_elapsed, total_elapsed / 60))

  # 比較表
  println("\n" * "="^80)
  println("性能比較サマリー")
  println("="^80)
  println(@sprintf("%-25s | %12s | %12s | %12s | %12s",
    "シナリオ", "実行時間[s]", "改善率[%]", "RMS残差[K]", "目的関数"))
  println("-"^80)

  baseline_time = all_results[1]["elapsed_time"]
  for result in all_results
    improvement = (baseline_time - result["elapsed_time"]) / baseline_time * 100
    println(@sprintf("%-25s | %12.2f | %+11.2f%% | %12.4e | %12.4e",
      result["scenario_name"],
      result["elapsed_time"],
      improvement,
      result["rms_error"],
      result["J_final"]
    ))
  end
  println("="^80)

  # 改善効果の分析
  println("\n改善効果分析:")
  println("-"^80)
  dhcp_linear = all_results[2]["elapsed_time"]
  dhcp_quadratic = all_results[3]["elapsed_time"]
  adjoint_residual = all_results[4]["elapsed_time"]
  sens_linear = all_results[5]["elapsed_time"]
  sens_quadratic = all_results[6]["elapsed_time"]

  println(@sprintf("DHCP線形外挿: %.2f%% 改善", (baseline_time - dhcp_linear) / baseline_time * 100))
  println(@sprintf("DHCP二次外挿: %.2f%% 改善", (baseline_time - dhcp_quadratic) / baseline_time * 100))
  println(@sprintf("Adjoint初期値改善: %.2f%% 改善", (baseline_time - adjoint_residual) / baseline_time * 100))
  println(@sprintf("Sensitivity線形外挿: %.2f%% 改善", (baseline_time - sens_linear) / baseline_time * 100))
  println(@sprintf("Sensitivity二次外挿: %.2f%% 改善", (baseline_time - sens_quadratic) / baseline_time * 100))
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
