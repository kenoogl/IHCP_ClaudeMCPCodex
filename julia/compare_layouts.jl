"""
compare_layouts.jl

旧レイアウト版と新レイアウト版の数値比較
"""

using NPZ
using Statistics
using Printf

println("="^70)
println("メモリレイアウト変更前後の計算結果比較")
println("="^70)

# 旧レイアウト版の結果読み込み
old_data = npzread("layout_old_result.npz")
q_old = old_data["q_global"]
println("\n旧レイアウト版（main）:")
println("  形状: $(size(q_old))  # 期待: (nt-1, ni, nj) = (49, 20, 20)")
println("  範囲: [$(minimum(q_old)), $(maximum(q_old))] W/m²")

# 新レイアウト版の結果読み込み
new_data = npzread("layout_new_result.npz")
q_new = new_data["q_global"]
println("\n新レイアウト版（tuning2）:")
println("  形状: $(size(q_new))  # 期待: (ni, nj, nt-1) = (20, 20, 49)")
println("  範囲: [$(minimum(q_new)), $(maximum(q_new))] W/m²")

# 配列形状を揃える（旧 → 新）
println("\n配列形状変換中...")
q_old_reshaped = permutedims(q_old, (2, 3, 1))  # (nt-1,ni,nj) → (ni,nj,nt-1)
println("  旧レイアウト変換後: $(size(q_old_reshaped))")

# 数値比較
println("\n" * "="^70)
println("数値精度比較")
println("="^70)

diff = q_old_reshaped .- q_new
abs_diff = abs.(diff)

max_abs_diff = maximum(abs_diff)
mean_abs_diff = mean(abs_diff)
rms_diff = sqrt(mean(diff.^2))

# 相対誤差
max_val = max(maximum(abs.(q_old_reshaped)), maximum(abs.(q_new)))
max_rel_diff = max_abs_diff / max_val
mean_rel_diff = mean_abs_diff / max_val

println("\n絶対誤差:")
println("  最大値: $(max_abs_diff) W/m²")
println("  平均値: $(mean_abs_diff) W/m²")
println("  RMS値:  $(rms_diff) W/m²")

println("\n相対誤差:")
println("  最大値: $(max_rel_diff) ($(max_rel_diff * 100)%)")
println("  平均値: $(mean_rel_diff) ($(mean_rel_diff * 100)%)")

# 判定基準
threshold_abs = 1e-10  # 絶対誤差閾値
threshold_rel = 1e-12  # 相対誤差閾値

println("\n" * "="^70)
println("判定結果")
println("="^70)

if max_abs_diff < threshold_abs && max_rel_diff < threshold_rel
  println("\n✅ 完全一致: メモリレイアウト変更は数値計算に影響なし")
  println("  最大絶対誤差 $(max_abs_diff) < $(threshold_abs)")
  println("  最大相対誤差 $(max_rel_diff) < $(threshold_rel)")
elseif max_rel_diff < 1e-6
  println("\n✅ 実用上一致: 相対誤差 < 0.0001%")
  println("  最大相対誤差: $(max_rel_diff * 100)%")
else
  println("\n⚠️ 差異検出: さらに調査が必要")
  println("  最大絶対誤差: $(max_abs_diff)")
  println("  最大相対誤差: $(max_rel_diff * 100)%")
end

# 詳細統計（10点サンプル）
println("\n" * "="^70)
println("詳細比較（サンプル10点）")
println("="^70)

sample_indices = [
  (1, 1, 1),
  (10, 10, 10),
  (20, 20, 25),
  (5, 15, 40),
  (15, 5, 30),
  (1, 20, 49),
  (20, 1, 49),
  (10, 10, 1),
  (10, 10, 49),
  (20, 20, 49)
]

println("\n  (i, j, t)  |   旧値    |   新値    |   差分    |  相対誤差")
println("-"^70)

for (i, j, t) in sample_indices
  old_val = q_old_reshaped[i, j, t]
  new_val = q_new[i, j, t]
  diff_val = old_val - new_val
  rel_err = abs(diff_val) / max(abs(old_val), abs(new_val), 1e-10)

  @printf("  (%2d,%2d,%2d) | %9.3e | %9.3e | %9.3e | %9.3e\n",
          i, j, t, old_val, new_val, diff_val, rel_err)
end

println("\n" * "="^70)
println("検証完了")
println("="^70)
