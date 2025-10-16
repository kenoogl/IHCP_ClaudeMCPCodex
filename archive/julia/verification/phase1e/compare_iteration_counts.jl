#!/usr/bin/env julia
"""
Phase 1-E: 標準版 vs 調整版の反復回数比較

標準版（κ=0.2, warmup=3）と調整版（κ=0.1, warmup=4）の
各ソルバーのステップごとの反復回数を比較表示する。

標準版のデータは既知（レポートから）、調整版は実測値を記録する。
"""

using Printf

println("="^80)
println("Phase 1-E: 反復回数比較 - 標準版 vs 調整版")
println("="^80)

# 標準版（κ=0.2, warmup=3）のデータ（既知）
println("\n## 標準版（κ=0.2, warmup=3）の反復回数")
println("出典: shared/results/performance_phase1e_adaptive_tolerance.md")

# DHCP Solver
println("\n### DHCP Solver")
println("| 時間ステップ | ベースライン | 標準版 | 削減率 | 備考 |")
println("|-------------|------------|--------|--------|------|")
dhcp_baseline = [397, 316, 326, 417, 480, 429, 396, 685, 551]
dhcp_standard = [397, 316, 293, 540, 474, 733, 315, 279, 539]
for (i, (b, s)) in enumerate(zip(dhcp_baseline, dhcp_standard))
    t = i + 1
    change = 100 * (s - b) / b
    marker = abs(change) > 20 ? (change > 0 ? "⚠️ 悪化" : "🏆 改善") : ""
    println(@sprintf("| t=%d | %d | %d | %+.1f%% | %s |", t, b, s, change, marker))
end
total_b = sum(dhcp_baseline)
total_s = sum(dhcp_standard)
println(@sprintf("| **合計** | **%d** | **%d** | **%+.2f%%** | |", total_b, total_s, 100*(total_s-total_b)/total_b))

# Adjoint Solver
println("\n### Adjoint Solver")
println("| 時間ステップ | ベースライン | 標準版 | 削減率 | 備考 |")
println("|-------------|------------|--------|--------|------|")
adj_baseline = [638, 685, 653, 681, 601, 618, 564, 699, 923]  # t=1～9
adj_standard = [577, 718, 539, 595, 470, 460, 612, 658, 664]
for (i, (b, s)) in enumerate(zip(adj_baseline, adj_standard))
    t = i
    change = 100 * (s - b) / b
    marker = abs(change) > 20 ? (change > 0 ? "⚠️ 悪化" : "🏆 改善") : ""
    println(@sprintf("| t=%d | %d | %d | %+.1f%% | %s |", t, b, s, change, marker))
end
total_b_adj = sum(adj_baseline)
total_s_adj = sum(adj_standard)
println(@sprintf("| **合計** | **%d** | **%d** | **%+.2f%%** | |", total_b_adj, total_s_adj, 100*(total_s_adj-total_b_adj)/total_b_adj))

# Sensitivity Solver
println("\n### Sensitivity Solver")
println("| 時間ステップ | ベースライン | 標準版 | 削減率 | 備考 |")
println("|-------------|------------|--------|--------|------|")
sens_baseline = [397, 540, 455, 463, 392, 449, 450, 455, 491]
sens_standard = [430, 638, 360, 578, 339, 396, 401, 458, 483]
for (i, (b, s)) in enumerate(zip(sens_baseline, sens_standard))
    t = i + 1
    change = 100 * (s - b) / b
    marker = abs(change) > 20 ? (change > 0 ? "⚠️ 悪化" : "🏆 改善") : ""
    println(@sprintf("| t=%d | %d | %d | %+.1f%% | %s |", t, b, s, change, marker))
end
total_b_sens = sum(sens_baseline)
total_s_sens = sum(sens_standard)
println(@sprintf("| **合計** | **%d** | **%d** | **%+.2f%%** | |", total_b_sens, total_s_sens, 100*(total_s_sens-total_b_sens)/total_b_sens))

println("\n" * "="^80)
println("標準版の特徴:")
println("  - DHCP: t=5,7で悪化（+29.5%, +70.9%）⚠️")
println("  - Adjoint: t=2で悪化（+4.8%）、t=7でわずかに悪化（+8.5%）")
println("  - Sensitivity: t=3,5で悪化（+18.1%, +24.8%）⚠️")
println("="^80)

println("\n## 調整版（κ=0.1, warmup=4）の期待効果")
println("\n調整の狙い:")
println("  - より保守的な緩和（κ=0.2 → 0.1）")
println("  - より長いwarmup期間（3 → 4ステップ）")
println("  → 標準版で見られた悪化を抑制")
println("\n期待される改善:")
println("  1. DHCP t=5,7の悪化を20-50%抑制")
println("  2. Sensitivity t=3,5の悪化を30-60%抑制")
println("  3. 全体的により安定した収束")
println("\n実測データを取得するには:")
println("  julia --project=julia julia/benchmarks/benchmark_phase1e_tuned.jl")
println("="^80)
