"""
verify_run_10steps_grid.jl

run_10steps_fullsize_test.jlの格子生成が正しいか検証
"""

using Printf

# run_10steps_fullsize_test.jlからbuild_z_grid関数を抽出
function build_z_grid(nk::Int, Lz::Float64, stretch_factor::Float64)
  # Step 1: z_faces生成（Python版と同じ）
  z_faces_normalized = range(1.0, stop = 0.0, length = nk + 1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  # Step 2: dzの計算
  dz = diff(z_faces)

  # Step 3: ΔZ配列（ガイドセル込み、nk+2個）
  ΔZ = zeros(nk + 2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # Step 4: z_centers（ガイドセル込み、nk+2個）
  z_centers = zeros(nk + 2)

  # ガイドセル
  z_centers[1] = z_faces[1]        # 底面ガイドセル
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

  return z_centers, ΔZ, z_faces, dz
end

# Python版の格子生成（検証用）
function python_z_grid_generation(nk::Int, Lz::Float64, stretch_factor::Float64)
  z_faces_normalized = range(1.0, stop=0.0, length=nk+1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1)

  dz = diff(z_faces)
  z_centers_py = zeros(nk)
  z_centers_py[1] = z_faces[1]
  z_centers_py[end] = z_faces[end]
  if nk > 2
    z_centers_py[2:end-1] = (z_faces[2:end-2] .+ z_faces[3:end-1]) ./ 2
  end

  return z_faces, z_centers_py, dz
end

# テスト実行
println("="^80)
println("run_10steps_fullsize_test.jl 格子生成検証")
println("="^80)

nk = 20
Lz = 0.5e-3
stretch_factor = 3.0

println("\n【設定】")
println(@sprintf("  物理セル数（nk）: %d", nk))
println(@sprintf("  z方向長さ（Lz）: %.3e m", Lz))
println(@sprintf("  伸張係数: %.1f", stretch_factor))

# Python版
println("\n" * "-"^80)
println("Python版格子生成")
println("-"^80)
py_z_faces, py_z_centers, py_dz = python_z_grid_generation(nk, Lz, stretch_factor)
println(@sprintf("  z_faces数: %d", length(py_z_faces)))
println(@sprintf("  z_centers数: %d", length(py_z_centers)))
println(@sprintf("  底面物理セル: %.6e", py_z_centers[1]))
println(@sprintf("  第2層物理セル: %.6e", py_z_centers[2]))
println(@sprintf("  表面物理セル: %.6e", py_z_centers[end]))

# run_10steps版
println("\n" * "-"^80)
println("run_10steps_fullsize_test.jl版格子生成")
println("-"^80)
z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)
println(@sprintf("  z_faces数: %d", length(z_faces)))
println(@sprintf("  z_centers数: %d（ガイドセル込み）", length(z_centers)))
println(@sprintf("  ΔZ数: %d（ガイドセル込み）", length(ΔZ)))
println(@sprintf("  底面ガイドセル: %.6e", z_centers[1]))
println(@sprintf("  底面物理セル: %.6e", z_centers[2]))
println(@sprintf("  第2層物理セル: %.6e", z_centers[3]))
println(@sprintf("  表面物理セル: %.6e", z_centers[nk+1]))
println(@sprintf("  表面ガイドセル: %.6e", z_centers[nk+2]))

# 比較
println("\n" * "="^80)
println("Python版との比較")
println("="^80)

println("\n【物理セル中心座標】")
let
  local max_err = 0.0
  for k in 1:nk
    py_val = py_z_centers[k]
    jl_val = z_centers[k+1]  # ガイドセル分シフト

    err = abs(jl_val - py_val)
    max_err = max(max_err, err)

    if err > 1e-14
      println(@sprintf("  物理セル %2d: Python=%.6e, Julia=%.6e, 誤差=%.6e",
        k, py_val, jl_val, err))
    end
  end
  println(@sprintf("  最大絶対誤差: %.6e", max_err))

  if max_err < 1e-14
    println("  ✅ Python版と完全一致")
  else
    println("  ❌ Python版と不一致")
  end
  global result = (max_err < 1e-14)
end

println("\n" * "="^80)
println("結論")
println("="^80)

if result
  println("✅ run_10steps_fullsize_test.jlの格子生成はPython版と完全に一致しています")
  println("   → 修正成功")
else
  println("❌ Python版との不一致があります")
  println("   → 修正が必要")
end

println("\n" * "="^80)
