"""
grid_helpers.jl

格子生成ヘルパー関数
"""

"""
    generate_guard_cell_grid(nk::Int, dz::Vector{Float64}, z_faces::Vector{Float64}) -> (z_centers, dz_grid)

Python版と同じ方法でZ方向格子を生成し、ガイドセル込み配列を直接生成（main.jlと同じ方式）

# Arguments
- `nk::Int`: 物理セル数（Z方向）
- `dz::Vector{Float64}`: セル幅配列 (nk,) [m]
- `z_faces::Vector{Float64}`: セル境界座標 (nk+1,) [m]

# Returns
- `z_centers::Vector{Float64}`: セル中心z座標（ガイドセル含む） (nk+2,) [m]
- `dz_grid::Vector{Float64}`: セル代表幅（ガイドセル含む） (nk+2,) [m]

# Notes
- z_centers[1]: 底面ガイドセル
- z_centers[2]～z_centers[nk+1]: 計算内点（物理セル）
- z_centers[nk+2]: 表面ガイドセル
- dz_grid[1]およびdz_grid[nk+2]はガイドセル用
"""
function generate_guard_cell_grid(nk::Int, dz::Vector{Float64}, z_faces::Vector{Float64})
  # 入力検証
  @assert length(dz) == nk "dz length must equal nk"
  @assert length(z_faces) == nk + 1 "z_faces length must equal nk + 1"

  # dz_grid（ガイドセル込み、nk+2個）
  dz_grid = zeros(Float64, nk+2)
  dz_grid[2:nk+1] = dz[1:nk]
  dz_grid[1] = dz_grid[2]
  dz_grid[nk+2] = dz_grid[nk+1]

  # z_centers（ガイドセル込み、nk+2個）
  z_centers = zeros(Float64, nk+2)

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

  return z_centers, dz_grid
end
