#!/usr/bin/env julia
"""
show_sliding_window.jl

スライディングウィンドウアルゴリズムの分割状況を表示

Usage:
  julia julia/scripts/show_sliding_window.jl --nt <nt> --window <window> --overlap <overlap>

Example:
  julia julia/scripts/show_sliding_window.jl --nt 10 --window 5 --overlap 2
  julia julia/scripts/show_sliding_window.jl --nt 100 --window 35 --overlap 10
  julia julia/scripts/show_sliding_window.jl --nt 100 --window 60 --overlap 15
"""

using Printf

# test_sliding_window_algorithm.jlのアルゴリズムを利用
function calculate_windows(nt::Int, window_size::Int, overlap::Int)
  windows = []
  start_idx = 0
  window_id = 1  # Julia式: 1始まり
  safety_counter = 0
  safety_limit = nt * 5

  total_flux_steps = nt - 1

  while start_idx < total_flux_steps
    safety_counter += 1
    if safety_counter > safety_limit
      @warn "Safety break"
      break
    end

    max_L = min(window_size, total_flux_steps - start_idx)
    end_idx = start_idx + max_L

    push!(windows, (
      id = window_id,
      start_idx = start_idx + 1,  # Julia式: 1始まり
      end_idx = end_idx + 1,       # Julia式: 1始まり
      max_L = max_L
    ))

    step = max(1, max_L - overlap)
    start_idx += step
    window_id += 1
  end

  return windows
end

function display_windows(nt::Int, window_size::Int, overlap::Int)
  println("="^80)
  println("スライディングウィンドウ分割状況")
  println("="^80)
  println()
  println("パラメータ:")
  println("  nt (総ステップ数): $nt")
  println("  window (ウィンドウサイズ): $window_size")
  println("  overlap (オーバーラップ): $overlap")
  println("  熱流束ステップ数 (nt-1): $(nt-1)")
  println()

  windows = calculate_windows(nt, window_size, overlap)

  println("-"^80)
  println("ウィンドウ詳細:")
  println("-"^80)
  @printf("%-4s %-8s %-8s %-8s %-8s %s\n", "ID", "開始", "終了", "長さ", "ステップ", "種別")
  println("-"^80)

  main_count = 0
  short_count = 0

  for (i, w) in enumerate(windows)
    step_val = i < length(windows) ? windows[i+1].start_idx - w.start_idx : 0
    
    # 主要ウィンドウか短いウィンドウかを判定
    kind = if w.max_L >= overlap
      main_count += 1
      "主要 ●"
    else
      short_count += 1
      "短い ○"
    end
    
    @printf("%-4d %-8d %-8d %-8d %-8s %s\n",
      w.id, w.start_idx, w.end_idx, w.max_L,
      step_val > 0 ? string(step_val) : "N/A", kind)
  end

  println()
  println("-"^80)
  println("統計:")
  println("-"^80)
  println("  総ウィンドウ数: $(length(windows))個")
  println("  主要ウィンドウ (length >= overlap): $(main_count)個")
  println("  短いウィンドウ (length < overlap): $(short_count)個")
  println("  主要ウィンドウ割合: $(@sprintf("%.1f%%", 100.0 * main_count / length(windows)))")
  println()

  # カバー率確認（1始まりインデックスで確認）
  covered = falses(nt-1)
  for w in windows
    for idx in w.start_idx:w.end_idx-1
      if 1 <= idx <= nt-1
        covered[idx] = true
      end
    end
  end

  coverage = count(covered) / (nt-1) * 100
  println("  カバー率: $(@sprintf("%.1f%%", coverage))")

  if coverage < 100.0
    uncovered = findall(.!covered)
    println("  ⚠️  未カバー領域あり: $uncovered")
  else
    println("  ✅ 全領域カバー済み")
  end

  println()
  println("="^80)
end

# コマンドライン引数パース
function parse_args()
  nt = nothing
  window = nothing
  overlap = nothing

  i = 1
  while i <= length(ARGS)
    if ARGS[i] == "--nt" && i < length(ARGS)
      nt = parse(Int, ARGS[i+1])
      i += 2
    elseif ARGS[i] == "--window" && i < length(ARGS)
      window = parse(Int, ARGS[i+1])
      i += 2
    elseif ARGS[i] == "--overlap" && i < length(ARGS)
      overlap = parse(Int, ARGS[i+1])
      i += 2
    elseif ARGS[i] == "--help" || ARGS[i] == "-h"
      println("Usage: julia julia/scripts/show_sliding_window.jl --nt NT --window WINDOW --overlap OVERLAP")
      println()
      println("Options:")
      println("  --nt       総時間ステップ数")
      println("  --window   ウィンドウサイズ")
      println("  --overlap  オーバーラップ")
      println()
      println("Example:")
      println("  julia julia/scripts/show_sliding_window.jl --nt 10 --window 5 --overlap 2")
      println("  julia julia/scripts/show_sliding_window.jl --nt 100 --window 35 --overlap 10")
      exit(0)
    else
      i += 1
    end
  end

  if isnothing(nt) || isnothing(window) || isnothing(overlap)
    println("エラー: 必須引数が不足しています")
    println("使用方法: julia julia/scripts/show_sliding_window.jl --nt NT --window WINDOW --overlap OVERLAP")
    println("ヘルプ: julia julia/scripts/show_sliding_window.jl --help")
    exit(1)
  end

  return nt, window, overlap
end

# メイン
function main()
  nt, window_size, overlap = parse_args()
  
  # 検証
  if nt < 2
    println("エラー: nt は 2 以上である必要があります")
    return 1
  end
  
  if window_size < 1
    println("エラー: window は 1 以上である必要があります")
    return 1
  end
  
  if overlap < 0
    println("エラー: overlap は 0 以上である必要があります")
    return 1
  end

  display_windows(nt, window_size, overlap)
  
  return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
  exit(main())
end
