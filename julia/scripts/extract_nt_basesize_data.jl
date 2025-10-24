#!/usr/bin/env julia
"""
nt-basesize測定結果からデータを抽出してCSV出力

入力:
  - shared/results/nt_basesize/nt*_bs*.log (新規測定、12ファイル)
  - shared/results/threads_basesize/step1_t4_bs*.log (既存nt=10データ、4ファイル)

出力:
  - shared/results/nt_basesize/measurement_results.csv
    列: nt, basesize, threads, dhcp_time, total_time, solver, precond

使い方:
  julia julia/scripts/extract_nt_basesize_data.jl
"""

using Printf
using DataFrames
using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
    extract_dhcp_time(log_path::String) -> Union{Float64, Nothing}

ログファイルからDHCP実行時間を抽出

パターン: "  DHCP elapsed time: XX.XX s"
"""
function extract_dhcp_time(log_path::String)
  if !isfile(log_path)
    @warn "Log file not found: $log_path"
    return nothing
  end

  for line in eachline(log_path)
    m = match(r"DHCP elapsed time:\s+([\d\.]+)\s+s", line)
    if m !== nothing
      return parse(Float64, m.captures[1])
    end
  end

  @warn "DHCP elapsed time not found in: $log_path"
  return nothing
end

"""
    extract_total_time(log_path::String) -> Union{Float64, Nothing}

ログファイルから総実行時間を抽出

パターン: "  Total runtime: XX.XX s"
"""
function extract_total_time(log_path::String)
  if !isfile(log_path)
    return nothing
  end

  for line in eachline(log_path)
    m = match(r"Total runtime:\s+([\d\.]+)\s+s", line)
    if m !== nothing
      return parse(Float64, m.captures[1])
    end
  end

  return nothing
end

"""
    extract_new_measurements() -> DataFrame

新規測定データ（nt=30,50,100）を抽出
"""
function extract_new_measurements()
  nt_values = [30, 50, 100]
  basesize_values = [400, 500, 600, 700]
  threads = 4
  solver = "pbicgstab"
  precond = "gs"

  results = DataFrame(
    nt = Int[],
    basesize = Int[],
    threads = Int[],
    dhcp_time = Float64[],
    total_time = Float64[],
    solver = String[],
    precond = String[]
  )

  log_dir = joinpath(PROJECT_ROOT, "shared", "results", "nt_basesize")

  println("Extracting new measurements (nt=30,50,100)...")
  for nt in nt_values
    for basesize in basesize_values
      log_path = joinpath(log_dir, "nt$(nt)_bs$(basesize).log")

      dhcp_time = extract_dhcp_time(log_path)
      total_time = extract_total_time(log_path)

      if dhcp_time !== nothing && total_time !== nothing
        push!(results, (nt, basesize, threads, dhcp_time, total_time, solver, precond))
        println(@sprintf("  nt=%3d, bs=%4d: DHCP=%.2fs, Total=%.2fs", nt, basesize, dhcp_time, total_time))
      else
        @warn "  nt=$nt, basesize=$basesize: Data extraction failed"
      end
    end
  end

  return results
end

"""
    extract_existing_measurements() -> DataFrame

既存測定データ（nt=10）を抽出
"""
function extract_existing_measurements()
  nt = 10
  basesize_values = [400, 500, 600, 700]
  threads = 4
  solver = "pbicgstab"
  precond = "gs"

  results = DataFrame(
    nt = Int[],
    basesize = Int[],
    threads = Int[],
    dhcp_time = Float64[],
    total_time = Float64[],
    solver = String[],
    precond = String[]
  )

  log_dir = joinpath(PROJECT_ROOT, "shared", "results", "threads_basesize")

  println("\nExtracting existing measurements (nt=10)...")
  for basesize in basesize_values
    log_path = joinpath(log_dir, "step1_t4_bs$(basesize).log")

    dhcp_time = extract_dhcp_time(log_path)
    total_time = extract_total_time(log_path)

    if dhcp_time !== nothing && total_time !== nothing
      push!(results, (nt, basesize, threads, dhcp_time, total_time, solver, precond))
      println(@sprintf("  nt=%3d, bs=%4d: DHCP=%.2fs, Total=%.2fs", nt, basesize, dhcp_time, total_time))
    else
      @warn "  nt=$nt, basesize=$basesize: Data extraction failed"
    end
  end

  return results
end

"""
    main()

メイン処理: データ抽出とCSV出力
"""
function main()
  println("="^80)
  println("nt-basesize measurement data extraction")
  println("="^80)

  # 新規測定データを抽出
  new_data = extract_new_measurements()

  # 既存測定データを抽出
  existing_data = extract_existing_measurements()

  # 統合
  all_data = vcat(existing_data, new_data)

  # nt, basesizeでソート
  sort!(all_data, [:nt, :basesize])

  # 出力
  output_path = joinpath(PROJECT_ROOT, "shared", "results", "nt_basesize", "measurement_results.csv")
  CSV.write(output_path, all_data)

  println("\n" * "="^80)
  println("Data extraction completed!")
  println("="^80)
  println("Total measurements: $(nrow(all_data))")
  println("Output file: $output_path")
  println("")

  # サマリー表示
  println("Summary by nt:")
  for nt_val in sort(unique(all_data.nt))
    subset = filter(row -> row.nt == nt_val, all_data)
    println(@sprintf("  nt=%3d: %d measurements", nt_val, nrow(subset)))
  end

  println("\nSummary by basesize:")
  for bs_val in sort(unique(all_data.basesize))
    subset = filter(row -> row.basesize == bs_val, all_data)
    println(@sprintf("  basesize=%4d: %d measurements", bs_val, nrow(subset)))
  end

  println("\n" * "="^80)

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
