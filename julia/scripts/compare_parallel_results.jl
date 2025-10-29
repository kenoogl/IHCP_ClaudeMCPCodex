#!/usr/bin/env julia
"""
compare_parallel_results.jl

並列化実装の数値精度検証スクリプト

1スレッド（基準）と8スレッド（並列化）の計算結果を比較し、
数値精度の劣化がないことを確認する。

使用方法:
  julia --project=julia julia/scripts/compare_parallel_results.jl
"""

using NPZ

println("="^80)
println("並列化実装の数値精度検証")
println("="^80)
println()

# データ読み込み
println("[1/4] Loading results...")
data1 = npzread("shared/results/julia_sliding_window_cgm2_1thread.npz")
data8 = npzread("shared/results/julia_sliding_window_cgm2_8threads.npz")
println("✅ Files loaded successfully")
println()

# データ構造確認
println("[2/4] Checking data structure...")
println("  1-thread data keys: ", keys(data1))
println("  8-thread data keys: ", keys(data8))

q1 = data1["q_global"]
q8 = data8["q_global"]

println("  q1 size: ", size(q1))
println("  q8 size: ", size(q8))
@assert size(q1) == size(q8) "Size mismatch between 1-thread and 8-thread results!"
println("✅ Data structure validated")
println()

# 誤差計算
println("[3/4] Computing errors...")
abs_diff = abs.(q1 .- q8)
abs_error_max = maximum(abs_diff)
abs_error_mean = sum(abs_diff) / length(abs_diff)

rel_diff = abs.(q1 .- q8) ./ (abs.(q1) .+ 1e-12)
rel_error_max = maximum(rel_diff)
rel_error_mean = sum(rel_diff) / length(rel_diff)

println("  Absolute error:")
println("    Maximum: ", abs_error_max, " W/m²")
println("    Mean   : ", abs_error_mean, " W/m²")
println()
println("  Relative error:")
println("    Maximum: ", rel_error_max)
println("    Mean   : ", rel_error_mean)
println()

# 統計情報
println("  Heat flux statistics (1-thread):")
println("    Min: ", minimum(q1), " W/m²")
println("    Max: ", maximum(q1), " W/m²")
println("    Mean: ", sum(q1) / length(q1), " W/m²")
println()
println("  Heat flux statistics (8-threads):")
println("    Min: ", minimum(q8), " W/m²")
println("    Max: ", maximum(q8), " W/m²")
println("    Mean: ", sum(q8) / length(q8), " W/m²")
println()

# 判定
println("[4/4] Validation...")
threshold = 1e-10
if rel_error_max < threshold
    println("✅ PASS: Numerical precision preserved!")
    println("   Maximum relative error (", rel_error_max, ") < threshold (", threshold, ")")
    exit(0)
else
    println("❌ FAIL: Numerical precision degradation detected!")
    println("   Maximum relative error (", rel_error_max, ") >= threshold (", threshold, ")")

    # 詳細分析
    println()
    println("Detailed analysis:")
    max_idx = argmax(rel_diff)
    println("  Worst point index: ", max_idx)
    println("  1-thread value: ", q1[max_idx], " W/m²")
    println("  8-thread value: ", q8[max_idx], " W/m²")
    println("  Absolute difference: ", abs_diff[max_idx], " W/m²")
    println("  Relative difference: ", rel_diff[max_idx])

    exit(1)
end
