"""
check_z_grid_logic.jl

Python版とJulia版のz方向格子生成ロジックの比較検証

目的:
1. 現在のJulia版実装の問題点を明確化
2. Python版との対応関係を確認
3. 正しい実装案を検証
"""

using Printf

function python_z_grid_generation(nz::Int, Lz::Float64, stretch_factor::Float64)
  """Python版の格子生成ロジック（nz個の物理セル）"""

  # Grid faces plotting for Z direction
  z_faces = range(1, 0, length=nz+1)  # Generate nz+1 grid surfaces

  # Irregular grid generation from regular grid faces
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1)) .* (exp.(stretch_factor .* z_faces) .- 1)

  z_centers = zeros(nz)

  # The Bottom-most grid center coordinate
  z_centers[1] = z_faces[1]

  # The top-most grid center coordinate
  z_centers[end] = z_faces[end]

  # Intermediate grids
  if nz > 2
    z_centers[2:end-1] = (z_faces[2:end-2] .+ z_faces[3:end-1]) ./ 2
  end

  # Boundary component Calculation
  dz = diff(z_faces)

  return z_faces, z_centers, dz
end

function julia_current_implementation(nk::Int, dz::Vector{Float64})
  """現在のJulia版実装（ガイドセル込み、nk+2個）"""

  # セル境界座標（z_faces）を計算
  z_faces = zeros(nk + 1)
  z_faces[1] = 0.0  # 底面
  for k in 1:nk
    z_faces[k+1] = z_faces[k] + dz[k]
  end

  # ΔZ配列（ガイドセル込み）
  ΔZ = zeros(nk+2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # セル中心座標（z_centers）を計算
  z_centers = zeros(nk+2)
  z_centers[1] = z_faces[1]
  z_centers[2] = z_faces[1]      # 底面セル中心=底面境界
  z_centers[end-1] = z_faces[end]  # 表面セル中心=表面境界
  z_centers[end]   = z_faces[end]

  # ⚠️ 問題のループ: z_centers[2]とz_centers[nk+1]を上書きしてしまう
  for k in 2:(nk+1)
    z_centers[k] = (z_faces[k] + z_faces[k-1]) * 0.5
  end

  return z_faces, z_centers, ΔZ
end

function julia_corrected_implementation(nk::Int, dz::Vector{Float64})
  """修正版Julia実装（ガイドセル込み、nk+2個）"""

  # セル境界座標（z_faces）を計算
  z_faces = zeros(nk + 1)
  z_faces[1] = 0.0  # 底面
  for k in 1:nk
    z_faces[k+1] = z_faces[k] + dz[k]
  end

  # ΔZ配列（ガイドセル込み）
  ΔZ = zeros(nk+2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # セル中心座標（z_centers）を計算
  z_centers = zeros(nk+2)

  # ガイドセル
  z_centers[1] = z_faces[1]      # 底面ガイドセル
  z_centers[nk+2] = z_faces[nk+1]  # 表面ガイドセル

  # 底面物理セル（境界に重なる）
  z_centers[2] = z_faces[1]

  # 表面物理セル（境界に重なる）
  z_centers[nk+1] = z_faces[nk+1]

  # ✅ 修正: 中間物理セルのみ（k=3からnkまで）
  if nk > 2
    for k in 3:nk
      z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
    end
  end

  return z_faces, z_centers, ΔZ
end

# テストケース: nz=20（Python版）、nk=20（Julia版）
println("="^80)
println("Z方向格子生成ロジックの検証")
println("="^80)

nz = 20
Lz = 0.5e-3
stretch_factor = 3.0

println("\n【Python版】物理セル数: $nz")
py_z_faces, py_z_centers, py_dz = python_z_grid_generation(nz, Lz, stretch_factor)

println(@sprintf("  z_faces数: %d（境界面）", length(py_z_faces)))
println(@sprintf("  z_centers数: %d（物理セル中心）", length(py_z_centers)))
println(@sprintf("  dz数: %d（各層の厚さ）", length(py_dz)))
println(@sprintf("  底面境界: %.6e", py_z_faces[1]))
println(@sprintf("  表面境界: %.6e", py_z_faces[end]))
println(@sprintf("  底面セル中心: %.6e", py_z_centers[1]))
println(@sprintf("  表面セル中心: %.6e", py_z_centers[end]))
println(@sprintf("  第2層セル中心: %.6e", py_z_centers[2]))

nk_with_guide = nz + 2
println(string("\n【Julia版（現在の実装）】物理セル数: ", nz, "、ガイドセル込み: ", nk_with_guide))
jl_z_faces, jl_z_centers, jl_ΔZ = julia_current_implementation(nz, py_dz)

println(@sprintf("  z_faces数: %d（境界面）", length(jl_z_faces)))
println(@sprintf("  z_centers数: %d（ガイドセル込み）", length(jl_z_centers)))
println(@sprintf("  ΔZ数: %d（ガイドセル込み）", length(jl_ΔZ)))
println(@sprintf("  底面境界: %.6e", jl_z_faces[1]))
println(@sprintf("  表面境界: %.6e", jl_z_faces[end]))
println(@sprintf("  底面ガイドセル: %.6e", jl_z_centers[1]))
println(@sprintf("  底面物理セル: %.6e", jl_z_centers[2]))
println(@sprintf("  第2層物理セル: %.6e", jl_z_centers[3]))
println(@sprintf("  表面物理セル: %.6e", jl_z_centers[end-1]))
println(@sprintf("  表面ガイドセル: %.6e", jl_z_centers[end]))

println(string("\n【Julia版（修正版）】物理セル数: ", nz, "、ガイドセル込み: ", nk_with_guide))
jl_corr_z_faces, jl_corr_z_centers, jl_corr_ΔZ = julia_corrected_implementation(nz, py_dz)

println(@sprintf("  z_faces数: %d（境界面）", length(jl_corr_z_faces)))
println(@sprintf("  z_centers数: %d（ガイドセル込み）", length(jl_corr_z_centers)))
println(@sprintf("  ΔZ数: %d（ガイドセル込み）", length(jl_corr_ΔZ)))
println(@sprintf("  底面境界: %.6e", jl_corr_z_faces[1]))
println(@sprintf("  表面境界: %.6e", jl_corr_z_faces[end]))
println(@sprintf("  底面ガイドセル: %.6e", jl_corr_z_centers[1]))
println(@sprintf("  底面物理セル: %.6e", jl_corr_z_centers[2]))
println(@sprintf("  第2層物理セル: %.6e", jl_corr_z_centers[3]))
println(@sprintf("  表面物理セル: %.6e", jl_corr_z_centers[end-1]))
println(@sprintf("  表面ガイドセル: %.6e", jl_corr_z_centers[end]))

# 比較検証
println("\n" * "="^80)
println("Python版とJulia版の対応関係検証")
println("="^80)

println("\n【物理セルの対応】")
println(@sprintf("%-30s | %-20s | %-20s | %-20s", "位置", "Python [idx]", "Julia現在 [idx]", "Julia修正 [idx]"))
println("-"^80)

# 底面セル
py_idx = 1
jl_idx = 2
println(@sprintf("%-30s | %.6e [%2d] | %.6e [%2d] | %.6e [%2d]",
  "底面物理セル",
  py_z_centers[py_idx], py_idx,
  jl_z_centers[jl_idx], jl_idx,
  jl_corr_z_centers[jl_idx], jl_idx))

# 第2層セル
py_idx = 2
jl_idx = 3
println(@sprintf("%-30s | %.6e [%2d] | %.6e [%2d] | %.6e [%2d]",
  "第2層物理セル",
  py_z_centers[py_idx], py_idx,
  jl_z_centers[jl_idx], jl_idx,
  jl_corr_z_centers[jl_idx], jl_idx))

# 表面セル
py_idx = nz
jl_idx = nz + 1
println(@sprintf("%-30s | %.6e [%2d] | %.6e [%2d] | %.6e [%2d]",
  "表面物理セル",
  py_z_centers[py_idx], py_idx,
  jl_z_centers[jl_idx], jl_idx,
  jl_corr_z_centers[jl_idx], jl_idx))

# 誤差計算
println("\n" * "="^80)
println("誤差解析")
println("="^80)

# 現在の実装との誤差
println("\n【現在のJulia実装 vs Python版】物理セル中心座標の誤差")
let
  local max_err = 0.0
  for i in 1:nz
    py_val = py_z_centers[i]
    jl_val = jl_z_centers[i+1]  # ガイドセル分シフト
    err = abs(jl_val - py_val)
    max_err = max(max_err, err)
    if err > 1e-10
      println(@sprintf("  物理セル %2d: Python=%.6e, Julia現在=%.6e, 誤差=%.6e", i, py_val, jl_val, err))
    end
  end
  println(@sprintf("  最大絶対誤差: %.6e", max_err))
  global max_err_current = max_err
end

# 修正版との誤差
println("\n【修正版Julia実装 vs Python版】物理セル中心座標の誤差")
let
  local max_err = 0.0
  for i in 1:nz
    py_val = py_z_centers[i]
    jl_val = jl_corr_z_centers[i+1]  # ガイドセル分シフト
    err = abs(jl_val - py_val)
    max_err = max(max_err, err)
    if err > 1e-10
      println(@sprintf("  物理セル %2d: Python=%.6e, Julia修正=%.6e, 誤差=%.6e", i, py_val, jl_val, err))
    end
  end
  println(@sprintf("  最大絶対誤差: %.6e", max_err))
  global max_err_corrected = max_err
end

# 結論
println("\n" * "="^80)
println("結論")
println("="^80)

if max_err_current > 1e-10
  println("⚠️ 現在のJulia実装には問題があります！")
  println("   ループが底面・表面セルの中心座標を上書きしています。")
else
  println("✓ 現在のJulia実装は正しいです。")
end

if max_err_corrected < 1e-10
  println("✓ 修正版Julia実装はPython版と一致します。")
else
  println("⚠️ 修正版Julia実装にも誤差があります。")
end

println("\n【推奨事項】")
if max_err_current > 1e-10
  println("1. main.jlの格子生成ループを修正する必要があります")
  println("2. 修正箇所: main.jl:179-181")
  println("3. 修正内容: ループを `for k in 3:nk` に変更")
end

println("\n" * "="^80)
