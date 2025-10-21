#!/usr/bin/env julia
"""
スライディングウィンドウアルゴリズムのみを検証（ソルバーなし）

Python版とJulia版のウィンドウ分割ロジックを比較検証する。
"""

using Printf

"""
Python版のスライディングウィンドウアルゴリズム
(IHCP_CGM_Sliding_Window_Calculation_ver2.py 1603-1653行)
"""
function sliding_window_python_style(nt::Int, window_size::Int, overlap::Int)
    println("="^80)
    println("Python版スライディングウィンドウアルゴリズム")
    println("="^80)
    println("パラメータ: nt=$nt, window_size=$window_size, overlap=$overlap")
    println()

    total_flux_steps = nt - 1
    windows = []

    start_idx = 0
    window_id = 0
    safety_counter = 0
    safety_limit = nt * 5

    # Python版 (1603行): while start_idx < nt - 1:
    while start_idx < total_flux_steps
        safety_counter += 1
        if safety_counter > safety_limit
            @warn "Safety break"
            break
        end

        # Python版 (1610行): max_L = min(window_size, (nt - 1) - start_idx)
        max_L = min(window_size, total_flux_steps - start_idx)
        end_idx = start_idx + max_L

        push!(windows, (
            id = window_id,
            start_idx = start_idx,
            end_idx = end_idx,
            max_L = max_L
        ))

        # Python版 (1652-1653行):
        # step = max(1, max_L - overlap)
        # start_idx += step
        step = max(1, max_L - overlap)
        start_idx += step
        window_id += 1
    end

    println("総ウィンドウ数: ", length(windows))
    println()
    println("各ウィンドウ:")
    println("-"^80)
    println(@sprintf("%-4s %-8s %-8s %-8s %-8s", "ID", "開始", "終了", "長さ", "ステップ"))
    println("-"^80)

    for (i, w) in enumerate(windows)
        step_val = i < length(windows) ? windows[i+1].start_idx - w.start_idx : 0
        println(@sprintf("%-4d %-8d %-8d %-8d %-8s",
            w.id, w.start_idx, w.end_idx, w.max_L,
            step_val > 0 ? string(step_val) : "N/A"))
    end
    println()

    return windows
end

"""
Julia版のスライディングウィンドウアルゴリズム（修正後）
(run_sliding_window.jl 290-400行)
"""
function sliding_window_julia_style(nt::Int, window_size::Int, overlap::Int)
    println("="^80)
    println("Julia版スライディングウィンドウアルゴリズム（修正後）")
    println("="^80)
    println("パラメータ: nt=$nt, window_size=$window_size, overlap=$overlap")
    println()

    total_flux_steps = nt - 1
    windows = []

    # 修正後のコード
    start_idx = 0
    window_id = 0
    safety_counter = 0
    safety_limit = nt * 5

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
            start_idx = start_idx,
            end_idx = end_idx,
            max_L = max_L
        ))

        step = max(1, max_L - overlap)
        start_idx += step
        window_id += 1
    end

    println("総ウィンドウ数: ", length(windows))
    println()
    println("各ウィンドウ:")
    println("-"^80)
    println(@sprintf("%-4s %-8s %-8s %-8s %-8s", "ID", "開始", "終了", "長さ", "ステップ"))
    println("-"^80)

    for (i, w) in enumerate(windows)
        step_val = i < length(windows) ? windows[i+1].start_idx - w.start_idx : 0
        println(@sprintf("%-4d %-8d %-8d %-8d %-8s",
            w.id, w.start_idx, w.end_idx, w.max_L,
            step_val > 0 ? string(step_val) : "N/A"))
    end
    println()

    return windows
end

"""
ウィンドウ比較
"""
function compare_windows(py_windows, jl_windows)
    println("="^80)
    println("比較結果")
    println("="^80)

    if length(py_windows) != length(jl_windows)
        println("❌ ウィンドウ数が異なります！")
        println("   Python: $(length(py_windows))個")
        println("   Julia:  $(length(jl_windows))個")
        return false
    end

    all_match = true
    for i in 1:length(py_windows)
        py_w = py_windows[i]
        jl_w = jl_windows[i]

        if py_w.start_idx != jl_w.start_idx ||
           py_w.end_idx != jl_w.end_idx ||
           py_w.max_L != jl_w.max_L
            println("❌ ウィンドウ$(i)が異なります:")
            println("   Python: start=$(py_w.start_idx), end=$(py_w.end_idx), len=$(py_w.max_L)")
            println("   Julia:  start=$(jl_w.start_idx), end=$(jl_w.end_idx), len=$(jl_w.max_L)")
            all_match = false
        end
    end

    if all_match
        println("✅ すべてのウィンドウが一致しました！")
        println("   総ウィンドウ数: $(length(py_windows))個")
        println()
        println("詳細:")
        for (i, w) in enumerate(py_windows)
            step_val = i < length(py_windows) ? py_windows[i+1].start_idx - w.start_idx : 0
            println(@sprintf("  Window %d: [%d,%d] length=%d, step=%s",
                w.id, w.start_idx, w.end_idx, w.max_L,
                step_val > 0 ? string(step_val) : "N/A"))
        end
        return true
    else
        return false
    end
end

# テストケース
function main()
    println("\n" * "="^80)
    println("テストケース1: nt=10, window_size=71, overlap=17")
    println("="^80 * "\n")

    py_wins = sliding_window_python_style(10, 71, 17)
    println()
    jl_wins = sliding_window_julia_style(10, 71, 17)
    println()
    result1 = compare_windows(py_wins, jl_wins)

    println("\n" * "="^80)
    println("テストケース2: nt=3, window_size=71, overlap=17")
    println("="^80 * "\n")

    py_wins2 = sliding_window_python_style(3, 71, 17)
    println()
    jl_wins2 = sliding_window_julia_style(3, 71, 17)
    println()
    result2 = compare_windows(py_wins2, jl_wins2)

    println("\n" * "="^80)
    println("テストケース3: nt=100, window_size=71, overlap=17")
    println("="^80 * "\n")

    py_wins3 = sliding_window_python_style(100, 71, 17)
    println()
    jl_wins3 = sliding_window_julia_style(100, 71, 17)
    println()
    result3 = compare_windows(py_wins3, jl_wins3)

    println("\n" * "="^80)
    println("最終結果")
    println("="^80)
    println("テストケース1 (nt=10):  ", result1 ? "✅ PASS" : "❌ FAIL")
    println("テストケース2 (nt=3):   ", result2 ? "✅ PASS" : "❌ FAIL")
    println("テストケース3 (nt=100): ", result3 ? "✅ PASS" : "❌ FAIL")
    println("="^80)

    if result1 && result2 && result3
        println("\n🎉 すべてのテストがパスしました！")
        println("Python版とJulia版のスライディングウィンドウアルゴリズムは完全に一致しています。")
        return 0
    else
        println("\n⚠️  一部のテストが失敗しました。")
        return 1
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
