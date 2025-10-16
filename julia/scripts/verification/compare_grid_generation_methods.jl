"""
compare_grid_generation_methods.jl

2つのz方向格子生成方法を比較：
1. Python版と同じ方法（nk個）→ convert_to_guard_cell_grid() でガイドセル追加
2. main.jlの直接方法（nk+2個を直接生成）

目的: 両方が同じ結果を生成することを確認
"""

using Printf

include("../../src/utils/GridTransform.jl")
using .GridTransform

# Python版の格子生成（run_10steps_fullsize_test.jlと同じ）
function method1_indirect(nk::Int, Lz::Float64, stretch_factor::Float64)
  """方法1: Python版 + convert_to_guard_cell_grid()"""

  # Step 1: Python版と同じ格子生成（nk個のみ）
  z_faces_normalized = range(1.0, stop=0.0, length=nk+1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  dz = diff(z_faces)
  z_centers = similar(dz)  # nk個
  z_centers[1] = z_faces[1]
  z_centers[end] = z_faces[end]
  if nk > 2
    z_centers[2:end-1] = (z_faces[2:end-2] .+ z_faces[3:end-1]) ./ 2
  end

  # Step 2: dz_top, dz_bottomの計算
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

  # Step 3: convert_to_guard_cell_grid()でガイドセル追加
  z_centers_guard, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

  return z_faces, z_centers, dz, z_centers_guard, dz_grid
end

# main.jlの直接方法（修正後）
function method2_direct(nk::Int, Lz::Float64, stretch_factor::Float64)
  """方法2: main.jlの直接方法（修正後）"""

  # Step 1: z_faces生成（Python版と同じ）
  z_faces_normalized = range(1.0, stop=0.0, length=nk+1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  # Step 2: dzの計算
  dz = diff(z_faces)

  # Step 3: ΔZ配列（ガイドセル込み、nk+2個）
  ΔZ = zeros(nk+2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # Step 4: z_centers（ガイドセル込み、nk+2個）
  z_centers = zeros(nk+2)

  # ガイドセル
  z_centers[1] = z_faces[1]      # 底面ガイドセル
  z_centers[nk+2] = z_faces[nk+1]  # 表面ガイドセル

  # 底面物理セル（境界に重なる）
  z_centers[2] = z_faces[1]

  # 表面物理セル（境界に重なる）
  z_centers[nk+1] = z_faces[nk+1]

  # 中間物理セルのみ（k=3からnkまで）
  if nk > 2
    for k in 3:nk
      z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
    end
  end

  return z_faces, dz, z_centers, ΔZ
end

# テスト実行
println("="^80)
println("2つのz方向格子生成方法の比較")
println("="^80)

nk = 20
Lz = 0.5e-3
stretch_factor = 3.0

println("\n【共通設定】")
println(@sprintf("  物理セル数（nk）: %d", nk))
println(@sprintf("  z方向長さ（Lz）: %.3e m", Lz))
println(@sprintf("  伸張係数: %.1f", stretch_factor))

# 方法1: Python版 + convert_to_guard_cell_grid()
println("\n" * "-"^80)
println("方法1: Python版（nk個）+ convert_to_guard_cell_grid()")
println("-"^80)
z_faces1, z_centers_py, dz1, z_centers_guard1, dz_grid1 = method1_indirect(nk, Lz, stretch_factor)

println(@sprintf("  z_faces数: %d", length(z_faces1)))
println(@sprintf("  z_centers（Python版）数: %d", length(z_centers_py)))
println(@sprintf("  z_centers_guard数: %d", length(z_centers_guard1)))
println(@sprintf("  dz_grid数: %d", length(dz_grid1)))
println(@sprintf("  底面ガイドセル: %.6e", z_centers_guard1[1]))
println(@sprintf("  底面物理セル: %.6e", z_centers_guard1[2]))
println(@sprintf("  第2層物理セル: %.6e", z_centers_guard1[3]))
println(@sprintf("  表面物理セル: %.6e", z_centers_guard1[nk+1]))
println(@sprintf("  表面ガイドセル: %.6e", z_centers_guard1[nk+2]))

# 方法2: main.jlの直接方法
println("\n" * "-"^80)
println("方法2: main.jlの直接方法（修正後）")
println("-"^80)
z_faces2, dz2, z_centers2, ΔZ2 = method2_direct(nk, Lz, stretch_factor)

println(@sprintf("  z_faces数: %d", length(z_faces2)))
println(@sprintf("  z_centers数: %d", length(z_centers2)))
println(@sprintf("  ΔZ数: %d", length(ΔZ2)))
println(@sprintf("  底面ガイドセル: %.6e", z_centers2[1]))
println(@sprintf("  底面物理セル: %.6e", z_centers2[2]))
println(@sprintf("  第2層物理セル: %.6e", z_centers2[3]))
println(@sprintf("  表面物理セル: %.6e", z_centers2[nk+1]))
println(@sprintf("  表面ガイドセル: %.6e", z_centers2[nk+2]))

# 比較
println("\n" * "="^80)
println("比較結果")
println("="^80)

# z_facesの一致確認
println("\n【z_faces（境界面座標）】")
max_diff_faces = maximum(abs.(z_faces1 .- z_faces2))
println(@sprintf("  最大絶対誤差: %.6e", max_diff_faces))
if max_diff_faces < 1e-14
  println("  ✓ 完全一致")
else
  println("  ⚠️ 不一致")
end

# z_centersの一致確認（ガイドセル込み）
println("\n【z_centers（ガイドセル込みセル中心座標）】")
max_diff_centers = maximum(abs.(z_centers_guard1 .- z_centers2))
println(@sprintf("  最大絶対誤差: %.6e", max_diff_centers))
if max_diff_centers < 1e-14
  println("  ✓ 完全一致")
else
  println("  ⚠️ 不一致")
  println("\n  詳細な誤差:")
  for k in 1:(nk+2)
    diff = abs(z_centers_guard1[k] - z_centers2[k])
    if diff > 1e-14
      println(@sprintf("    z_centers[%2d]: 方法1=%.6e, 方法2=%.6e, 誤差=%.6e",
        k, z_centers_guard1[k], z_centers2[k], diff))
    end
  end
end

# dzの一致確認
println("\n【dz（セル代表幅）】")
# dz_grid1はnk+1個、ΔZ2はnk+2個なので比較範囲を調整
# dz_grid1[k] vs ΔZ2[k] for k=2:nk
if length(dz_grid1) == nk + 1 && length(ΔZ2) == nk + 2
  # dz_grid1は[2:nk+1]が物理セルの代表幅
  # ΔZ2は[2:nk+1]が物理セルの幅
  # 比較対象: dz_grid1 vs ΔZ2[2:nk+2]? または dz_grid1[2:nk] vs ΔZ2[2:nk+1]?

  # まずは直接比較できる部分を確認
  println(@sprintf("  dz_grid1の長さ: %d", length(dz_grid1)))
  println(@sprintf("  ΔZ2の長さ: %d", length(ΔZ2)))

  # GridTransformの定義: dz[k] = 0.5*(z_centers[k+1] - z_centers[k-1]) (k=2～nk)
  # main.jlの定義: ΔZ[k] = セル幅（物理セル）

  # これらは異なる定義なので、直接比較は難しい
  println("  注意: dz_grid（GridTransform）とΔZ（main.jl）は異なる定義")
  println("  → dz_gridは中心間距離ベース、ΔZはセル幅ベース")
  println("  → 直接比較は困難、物理的な意味の確認が必要")
else
  println("  ⚠️ 配列サイズが異なります")
end

# Python版z_centersとの一致確認
println("\n【Python版z_centers（nk個）との比較】")
let
  local max_err = 0.0
  for k in 1:nk
    py_val = z_centers_py[k]
    m1_val = z_centers_guard1[k+1]  # ガイドセル分シフト
    m2_val = z_centers2[k+1]        # ガイドセル分シフト

    diff_m1 = abs(m1_val - py_val)
    diff_m2 = abs(m2_val - py_val)

    max_err = max(max_err, diff_m1, diff_m2)

    if diff_m1 > 1e-14 || diff_m2 > 1e-14
      println(@sprintf("  物理セル %2d: Python=%.6e, 方法1=%.6e, 方法2=%.6e",
        k, py_val, m1_val, m2_val))
      println(@sprintf("              誤差: 方法1=%.6e, 方法2=%.6e", diff_m1, diff_m2))
    end
  end
  println(@sprintf("  最大絶対誤差: %.6e", max_err))
  global max_diff_py = max_err
end
if max_diff_py < 1e-14
  println("  ✓ 両方法ともPython版と完全一致")
else
  println("  ⚠️ Python版と不一致")
end

# 結論
println("\n" * "="^80)
println("結論")
println("="^80)

if max_diff_centers < 1e-14 && max_diff_py < 1e-14
  println("✅ 両方法は完全に一致しており、Python版とも一致します")
  println("   → どちらの方法を使っても正しい結果が得られます")
else
  println("⚠️ 不一致が検出されました")
  if max_diff_py < 1e-14
    println("   → Python版とは一致しているので、両方法とも正しい可能性があります")
    println("   → ただし、内部表現が異なる可能性があります")
  else
    println("   → Python版との不一致があります")
    println("   → 実装を見直す必要があります")
  end
end

println("\n" * "="^80)
