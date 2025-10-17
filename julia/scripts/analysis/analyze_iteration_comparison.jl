using Printf

println("=== Phase 1-E vs 体積積分形式 反復回数比較 ===\n")

# Phase 1-E反復回数（ドキュメントから）
phase1e_dhcp = [397, 316, 293, 540, 474, 733, 315, 279, 539]     # t=2-10
phase1e_adjoint = [664, 658, 612, 460, 470, 595, 539, 718, 577]  # t=9-1（後退）
phase1e_sensitivity = [430, 638, 360, 578, 339, 396, 401, 458, 483]  # t=2-10

# 体積積分形式反復回数（ログから抽出）
volume_dhcp = [11, 12, 11, 12, 11, 12, 12, 12, 12]  # t=2-10
volume_adjoint = [9, 11, 11, 11, 11, 12, 13, 13, 12]  # t=9-1（後退）
volume_sensitivity = [12, 11, 10, 11, 10, 11, 10, 11, 10]  # t=2-10

# DHCP比較
println("### DHCP Solver（前進時間積分）")
println("| 時間ステップ | Phase 1-E | 体積積分形式 | 削減率 |")
println("|-------------|----------|------------|--------|")
total_phase1e_dhcp = 0
total_volume_dhcp = 0
for (i, t) in enumerate(2:10)
  global total_phase1e_dhcp, total_volume_dhcp
  reduction = (1 - volume_dhcp[i] / phase1e_dhcp[i]) * 100
  @printf("| t=%d | %d | %d | **%.1f%%** |\n", t, phase1e_dhcp[i], volume_dhcp[i], reduction)
  total_phase1e_dhcp += phase1e_dhcp[i]
  total_volume_dhcp += volume_dhcp[i]
end
total_reduction_dhcp = (1 - total_volume_dhcp / total_phase1e_dhcp) * 100
@printf("| **合計** | **%d** | **%d** | **%.1f%%** 🏆 |\n", total_phase1e_dhcp, total_volume_dhcp, total_reduction_dhcp)
@printf("\n平均反復回数: Phase 1-E=%.1f, 体積積分=%.1f\n\n", total_phase1e_dhcp/9, total_volume_dhcp/9)

# Adjoint比較
println("### Adjoint Solver（後退時間積分）")
println("| 時間ステップ | Phase 1-E | 体積積分形式 | 削減率 |")
println("|-------------|----------|------------|--------|")
total_phase1e_adj = 0
total_volume_adj = 0
for (i, t) in enumerate(9:-1:1)
  global total_phase1e_adj, total_volume_adj
  reduction = (1 - volume_adjoint[i] / phase1e_adjoint[i]) * 100
  @printf("| t=%d | %d | %d | **%.1f%%** |\n", t, phase1e_adjoint[i], volume_adjoint[i], reduction)
  total_phase1e_adj += phase1e_adjoint[i]
  total_volume_adj += volume_adjoint[i]
end
total_reduction_adj = (1 - total_volume_adj / total_phase1e_adj) * 100
@printf("| **合計** | **%d** | **%d** | **%.1f%%** 🏆 |\n", total_phase1e_adj, total_volume_adj, total_reduction_adj)
@printf("\n平均反復回数: Phase 1-E=%.1f, 体積積分=%.1f\n\n", total_phase1e_adj/9, total_volume_adj/9)

# Sensitivity比較
println("### Sensitivity Solver（前進時間積分）")
println("| 時間ステップ | Phase 1-E | 体積積分形式 | 削減率 |")
println("|-------------|----------|------------|--------|")
total_phase1e_sens = 0
total_volume_sens = 0
for (i, t) in enumerate(2:10)
  global total_phase1e_sens, total_volume_sens
  reduction = (1 - volume_sensitivity[i] / phase1e_sensitivity[i]) * 100
  @printf("| t=%d | %d | %d | **%.1f%%** |\n", t, phase1e_sensitivity[i], volume_sensitivity[i], reduction)
  total_phase1e_sens += phase1e_sensitivity[i]
  total_volume_sens += volume_sensitivity[i]
end
total_reduction_sens = (1 - total_volume_sens / total_phase1e_sens) * 100
@printf("| **合計** | **%d** | **%d** | **%.1f%%** 🏆 |\n", total_phase1e_sens, total_volume_sens, total_reduction_sens)
@printf("\n平均反復回数: Phase 1-E=%.1f, 体積積分=%.1f\n\n", total_phase1e_sens/9, total_volume_sens/9)

# 総合
println("### 総合比較")
total_phase1e = total_phase1e_dhcp + total_phase1e_adj + total_phase1e_sens
total_volume = total_volume_dhcp + total_volume_adj + total_volume_sens
total_reduction = (1 - total_volume / total_phase1e) * 100
@printf("- Phase 1-E総反復回数: %d回\n", total_phase1e)
@printf("- 体積積分形式総反復回数: %d回\n", total_volume)
@printf("- **削減率: %.1f%%** 🏆\n", total_reduction)
@printf("- **削減倍率: %.1f倍** 🚀\n", total_phase1e / total_volume)
