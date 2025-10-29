#!/usr/bin/env julia
"""
Step2測定結果（threads × basesize）からデータを抽出してCSV出力

入力:
  - shared/results/threads_basesize/step2_t*_bs*.log (18ファイル)

出力:
  - shared/results/threads_basesize/step2_results.csv
    列: threads, basesize, cgm_time, total_time, solver, precond

使い方:
  julia julia/scripts/extract_step2_data.jl
"""

using Printf
using DataFrames
using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
    extract_cgm_time(log_path::String) -> Union{Float64, Nothing}

ログファイルからCGM実行時間を抽出

パターン: "  CGM elapsed time: XX.XX s"
"""
function extract_cgm_time(log_path::String)
  if !isfile(log_path)
    @warn "Log file not found: $log_path"
    return nothing
  end

  for line in eachline(log_path)
    m = match(r"CGM elapsed time:\s+([\d\.]+)\s+s", line)
    if m !== nothing
      return parse(Float64, m.captures[1])
    end
  end

  @warn "CGM elapsed time not found in: $log_path"
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
    extract_step2_measurements() -> DataFrame

Step2測定データ（threads × basesize）を抽出
"""
function extract_step2_measurements()
  threads_values = [2, 4, 8]
  basesize_values = [400, 500, 600, 1000, 2000, 5000]
  solver = "pbicgstab"
  precond = "gs"

  results = DataFrame(
    threads = Int[],
    basesize = Int[],
    cgm_time = Float64[],
    total_time = Float64[],
    solver = String[],
    precond = String[]
  )

  log_dir = joinpath(PROJECT_ROOT, "shared", "results", "threads_basesize")

  println("Extracting Step2 measurements...")
  for threads in threads_values
    for basesize in basesize_values
      log_path = joinpath(log_dir, "step2_t$(threads)_bs$(basesize).log")

      cgm_time = extract_cgm_time(log_path)
      total_time = extract_total_time(log_path)

      if cgm_time !== nothing && total_time !== nothing
        push!(results, (threads, basesize, cgm_time, total_time, solver, precond))
        println(@sprintf("  threads=%d, bs=%4d: CGM=%.2fs, Total=%.2fs", threads, basesize, cgm_time, total_time))
      else
        @warn "  threads=$threads, basesize=$basesize: Data extraction failed"
      end
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
  println("Step2 (threads × basesize) measurement data extraction")
  println("="^80)

  # 測定データを抽出
  data = extract_step2_measurements()

  # threads, basesizeでソート
  sort!(data, [:threads, :basesize])

  # 出力
  output_path = joinpath(PROJECT_ROOT, "shared", "results", "threads_basesize", "step2_results.csv")
  CSV.write(output_path, data)

  println("\n" * "="^80)
  println("Data extraction completed!")
  println("="^80)
  println("Total measurements: $(nrow(data))")
  println("Output file: $output_path")
  println("")

  # サマリー表示
  println("Summary by threads:")
  for t_val in sort(unique(data.threads))
    subset = filter(row -> row.threads == t_val, data)
    println(@sprintf("  threads=%d: %d measurements", t_val, nrow(subset)))
  end

  println("\nSummary by basesize:")
  for bs_val in sort(unique(data.basesize))
    subset = filter(row -> row.basesize == bs_val, data)
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
