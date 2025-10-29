#!/usr/bin/env julia
"""
DHCP（直接熱伝導問題）ソルバー単体テスト

run_10steps_fullsize_test.jlから抽出したDHCPソルバーのみの基本テスト。
コアソルバーの性能と精度を確認するための軽量テストスクリプト。

Usage:
  julia test_dhcp_solver.jl [OPTIONS]

Options:
  --solver SOLVER     Solver type: pbicgstab (default), pcg
  --precond PRECOND   Preconditioner type: diagonal (default), gs, none
  --nt STEPS          Number of time steps (default: 10)
  --basesize SIZE     FLoops basesize parameter (default: 600)
  --ni NI             X-direction grid count (default: 80)
  --nj NJ             Y-direction grid count (default: 100)
  --nk NK             Z-direction grid count (default: 20)
  --synthetic         Enable synthetic test mode (no measurement data required)

Examples:
  # Default (measurement data mode)
  julia test_dhcp_solver.jl

  # Synthetic mode with custom grid size
  julia test_dhcp_solver.jl --synthetic --ni 40 --nj 50 --nk 10

  # Large-scale synthetic test
  JULIA_NUM_THREADS=4 julia test_dhcp_solver.jl --synthetic --ni 160 --nj 200 --nk 40 --basesize 600
"""

using Dates
using Printf
using Statistics
using NPZ

include("../src/IHCP_CGM.jl")
using .IHCP_CGM
using .IHCP_CGM.Commons: set_backend_config

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

"""
    parse_command_line_args() -> (solver_type, precond_type, nt, basesize, ni, nj, nk, synthetic)

コマンドライン引数をパースしてソルバー、前処理、ステップ数、basesize、格子数、モードを取得

Returns:
  solver_type: ソルバータイプ（Symbol: :pbicgstab, :pcg）
  precond_type: 前処理タイプ（Symbol: :diagonal, :gs, :none）
  nt: タイムステップ数（Int）
  basesize: FLoops basesizeパラメータ（Int）
  ni, nj, nk: 格子数（Int）
  synthetic: 合成テストモードフラグ（Bool）
"""
function parse_command_line_args()
  solver_type = :pbicgstab  # デフォルト
  precond_type = :diagonal  # デフォルト（Python版互換）
  nt = 10                   # デフォルト
  basesize = 600            # 最適値（スレッド数×basesize最適化実験の結果）

  # 新規追加: 格子数とモード
  ni = 80
  nj = 100
  nk = 20
  synthetic = false

  i = 1
  while i <= length(ARGS)
    if ARGS[i] == "--help" || ARGS[i] == "-h"
      println("""
      Usage: julia test_dhcp_solver.jl [OPTIONS]

      Options:
        --solver SOLVER       Solver type (default: pbicgstab)
                              Available: pbicgstab, pcg
        --precond PRECOND     Preconditioner type (default: diagonal)
                              Available: diagonal, gs, none
        --nt STEPS            Number of time steps (default: 10)
        --basesize SIZE       FLoops basesize parameter (default: 600)
        --ni NI               X-direction grid count (default: 80)
        --nj NJ               Y-direction grid count (default: 100)
        --nk NK               Z-direction grid count (default: 20)
        --synthetic           Enable synthetic test mode (no measurement data)
        -h, --help            Show this help message and exit

      Examples:
        julia test_dhcp_solver.jl
        julia test_dhcp_solver.jl --solver pcg --precond gs --nt 20
        julia test_dhcp_solver.jl --synthetic --ni 40 --nj 50 --nk 10
        julia test_dhcp_solver.jl --solver pbicgstab --precond gs --basesize 600 --synthetic
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
      if precond_str == "diagonal"
        precond_type = :diagonal
      elseif precond_str == "gs"
        precond_type = :gs
      elseif precond_str == "none"
        precond_type = :none
      else
        error("Unknown preconditioner type: $(ARGS[i + 1]). Available: diagonal, gs, none")
      end
      i += 2
    elseif ARGS[i] == "--nt"
      if i + 1 > length(ARGS)
        error("--nt requires an argument")
      end
      nt = parse(Int, ARGS[i + 1])
      if nt < 1
        error("--nt must be positive: $(ARGS[i + 1])")
      end
      i += 2
    elseif ARGS[i] == "--basesize"
      if i + 1 > length(ARGS)
        error("--basesize requires an argument")
      end
      basesize = parse(Int, ARGS[i + 1])
      if basesize < 1
        error("--basesize must be positive: $(ARGS[i + 1])")
      end
      i += 2
    elseif ARGS[i] == "--ni"
      if i + 1 > length(ARGS)
        error("--ni requires an argument")
      end
      ni = parse(Int, ARGS[i + 1])
      if ni < 1
        error("--ni must be positive: $(ARGS[i + 1])")
      end
      i += 2
    elseif ARGS[i] == "--nj"
      if i + 1 > length(ARGS)
        error("--nj requires an argument")
      end
      nj = parse(Int, ARGS[i + 1])
      if nj < 1
        error("--nj must be positive: $(ARGS[i + 1])")
      end
      i += 2
    elseif ARGS[i] == "--nk"
      if i + 1 > length(ARGS)
        error("--nk requires an argument")
      end
      nk = parse(Int, ARGS[i + 1])
      if nk < 1
        error("--nk must be positive: $(ARGS[i + 1])")
      end
      i += 2
    elseif ARGS[i] == "--synthetic"
      synthetic = true
      i += 1
    else
      error("Unknown argument: $(ARGS[i])")
    end
  end

  return solver_type, precond_type, nt, basesize, ni, nj, nk, synthetic
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

"""
    prepare_synthetic_test(ni, nj, nk, nt) -> (T_init, Y_obs)

合成テストモード用の初期化

# 引数
- ni, nj, nk: 格子数
- nt: タイムステップ数

# 戻り値
- T_init: 初期温度場 (ni, nj, nk) [K]
- Y_obs: 観測データ（合成） (ni, nj, nt) [K]

# 温度設定
一様温度（T_base = 573.15K）
"""
function prepare_synthetic_test(ni::Int, nj::Int, nk::Int, nt::Int)
  println("\n[2/5] Preparing synthetic test data")
  flush(stdout)

  # 一様温度
  T_base = 573.15  # 300°C in Kelvin
  T_init = fill(T_base, ni, nj, nk)
  Y_obs = fill(T_base, ni, nj, nt)

  println("  Mode: Synthetic test")
  println(@sprintf("  Initial temperature: %.2f K (uniform)", T_base))
  println(@sprintf("  Grid size: %d × %d × %d (N=%d)", ni, nj, nk, ni*nj*nk))
  flush(stdout)

  return T_init, Y_obs
end

"""
    prepare_measurement_test(ni, nj, nk, nt) -> (T_init, Y_obs)

測定データモード用の初期化（既存の処理を関数化）

# サイズ検証
測定データのサイズ（80×100）と一致する必要あり
"""
function prepare_measurement_test(ni::Int, nj::Int, nk::Int, nt::Int)
  println("\n[2/5] Loading measurement data")
  flush(stdout)

  Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)

  # サイズ検証
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)×$(nj_file), expected $(ni)×$(nj)")
  end

  # メモリレイアウト最適化 - (nt,ni,nj) → (ni,nj,nt)
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))
  Y_obs_python = nothing
  GC.gc()

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))
  flush(stdout)

  return T_init, Y_obs
end

"""
    analyze_residuals(T_result, Y_obs, synthetic)

残差分析（モード依存）

# 合成モード
- 温度統計のみ表示
- 残差計算なし（Y_obsが意味を持たない）

# 測定データモード
- RMS残差、最大残差を計算
- 測定データとの詳細比較
"""
function analyze_residuals(T_result::AbstractArray{<:Real,4},
                          Y_obs::AbstractArray{<:Real,3},
                          synthetic::Bool)
  T_bottom = @view T_result[:, :, 1, :]  # (ni, nj, nt)

  if synthetic
    # 合成モード: 簡易分析
    println("  Mode: Synthetic test")
    println(@sprintf("  Mean temperature: %.2f K", mean(T_result)))
    println(@sprintf("  Temperature range: %.2f ~ %.2f K", minimum(T_result), maximum(T_result)))
    println(@sprintf("  Temperature change: %.4e K (max-min)", maximum(T_result) - minimum(T_result)))
    println("  Note: No residual analysis (synthetic mode)")
  else
    # 測定データモード: 詳細分析
    residual = T_bottom .- Y_obs
    rms_error = sqrt(mean(residual .^ 2))
    max_error = maximum(abs.(residual))
    mean_temp = mean(T_result)

    println("  Mode: Measurement data")
    println(@sprintf("  RMS residual: %.4e K", rms_error))
    println(@sprintf("  Max residual: %.4e K", max_error))
    println(@sprintf("  Mean temperature: %.2f K", mean_temp))
    println(@sprintf("  Temperature range: %.2f ~ %.2f K", minimum(T_result), maximum(T_result)))
  end
  flush(stdout)
end

# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

function main()
  println("="^80)
  println("DHCP（直接熱伝導問題）ソルバー単体テスト")
  println("="^80)
  println("Project root: $PROJECT_ROOT")
  flush(stdout)

  # コマンドライン引数パース（拡張版）
  solver_type, precond_type, nt, basesize, ni, nj, nk, synthetic = parse_command_line_args()

  # Backend設定（並列化パラメータ）
  set_backend_config(basesize=basesize)

  println("\n[Configuration]")
  println("  Solver: $solver_type")
  println("  Preconditioner: $precond_type")
  println("  Time steps: $nt")
  println("  Grid size: $ni × $nj × $nk (N=$(ni*nj*nk))")
  println("  FLoops basesize: $basesize")
  println("  Julia threads: $(Threads.nthreads())")
  println("  Test mode: $(synthetic ? "Synthetic" : "Measurement data")")
  flush(stdout)

  total_start = time()

  # Problem definition -----------------------------------------------------
  dt = 1.0e-3
  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  Lz = 0.5e-3
  stretch_factor = 3.0

  println("\n[1/5] Grid and material parameters")
  println(@sprintf("  grid: ni=%d, nj=%d, nk=%d", ni, nj, nk))
  println(@sprintf("  spacing: dx=%.3e m, dy=%.3e m", dx, dy))
  println(@sprintf("  time steps: nt=%d, dt=%.3e s", nt, dt))
  flush(stdout)

  z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)

  rho, cp_coeffs, k_coeffs = load_material_properties()
  println(@sprintf("  rho (reference): %.6f kg/m^3", rho))
  flush(stdout)

  # Input preparation ------------------------------------------------------
  # モード分岐（重要な変更点）
  if synthetic
    T_init, Y_obs = prepare_synthetic_test(ni, nj, nk, nt)
  else
    T_init, Y_obs = prepare_measurement_test(ni, nj, nk, nt)
  end

  # DHCP初期化
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  # 既知の熱流束（簡易テスト用にゼロを仮定）
  q_test = zeros(Float64, ni, nj, nt - 1)
  println("\n[3/5] 設定された熱流束パターン")
  println(@sprintf("  q_test shape: %s", size(q_test)))
  println(@sprintf("  q_test range: %.4e ~ %.4e W/m^2", minimum(q_test), maximum(q_test)))
  flush(stdout)

  # DHCP solve --------------------------------------------------------------
  println("\n[4/5] Running DHCP solver")
  flush(stdout)
  solve_start = time()

  # solve_dhcp!の引数:
  # T_initial, q_surface, wk, nt, ρ, cp_coeffs, k_coeffs, dx, dy, ZC, dz, dt
  # 注意: dzはガイドセル込み配列（ΔZ）を渡す必要がある
  # キーワード引数: rtol, maxiter, verbose, solver, smoother
  T_result, iter_counts, residuals = solve_dhcp!(
    T_init,           # T_initial (ni, nj, nk)
    q_test,           # q_surface (ni, nj, nt-1)
    work,             # WorkBuffers
    nt,               # タイムステップ数
    rho,              # 密度
    cp_coeffs,        # 比熱係数
    k_coeffs,         # 熱伝導率係数
    dx,               # x方向格子幅
    dy,               # y方向格子幅
    z_centers,        # セル中心座標（ガイドセル込み、nk+2個）
    ΔZ,               # セル幅（ガイドセル込み、nk+2個）
    dt;               # 時間ステップ幅
    rtol = 1.0e-6,
    maxiter = 20000,
    verbose = true,
    solver = solver_type,
    smoother = precond_type
  )

  dhcp_elapsed = time() - solve_start
  println(@sprintf("  DHCP elapsed time: %.2f s", dhcp_elapsed))
  flush(stdout)

  # Residual diagnostics ---------------------------------------------------
  println("\n[5/5] Residual analysis")
  flush(stdout)
  analyze_residuals(T_result, Y_obs, synthetic)

  total_elapsed = time() - total_start
  println("\nSummary")
  println(@sprintf("  Total runtime: %.2f s", total_elapsed))
  println(@sprintf("  DHCP share: %.1f%%", 100 * dhcp_elapsed / total_elapsed))
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
