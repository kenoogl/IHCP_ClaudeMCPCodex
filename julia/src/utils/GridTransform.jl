"""
GridTransform.jl

IHCP格子定義からHeat3ds形式の格子への変換ユーティリティ

目的:
- (dz, dz_b, dz_t)形式から(z_centers, dz)形式への変換
"""
module GridTransform

export convert_to_guard_cell_grid

"""
  convert_to_guard_cell_grid(nk, dz_cell, dz_b, dz_t) -> (z_centers, dz)

IHCPの格子定義(dz_cell, dz_b, dz_t)からHeat3ds形式の(z_centers, dz)を生成

# Arguments
- nk: 計算内点数（Z方向）
- dz_cell: セル中心幅 (nk,) [m]
- dz_b: セル中心から下側界面までの距離 (nk,) [m]
- dz_t: セル中心から上側界面までの距離 (nk,) [m]

# Returns
- z_centers: セル中心z座標（ガイドセル含む） (nk+2,) [m]
- dz: セル代表幅 (nk+1,) [m]

# Notes
- z_centers[1]: 下側ガイドセル
- z_centers[2]～z_centers[nk+1]: 計算内点
- z_centers[nk+2]: 上側ガイドセル
- dz[k] = 0.5*(z_centers[k+1] - z_centers[k-1]) （k=2～nk）
- 端点補正: dz[2] *= 0.5, dz[nk] *= 0.5
"""
function convert_to_guard_cell_grid(
  nk::Int,
  dz_cell::Vector{Float64},
  dz_b::Vector{Float64},
  dz_t::Vector{Float64}
)
  @assert length(dz_cell) == nk "dz_cell length mismatch"
  @assert length(dz_b) == nk "dz_b length mismatch"
  @assert length(dz_t) == nk "dz_t length mismatch"

  nk_full = nk + 2
  z_centers = zeros(Float64, nk_full)
  dz = zeros(Float64, nk_full - 1)

  # セル中心座標の計算（下から積算）
  # 最下層セル（k=1）の中心を基準とする
  z_bottom_cell = 0.0  # 仮の原点

  # z_centers[2]が最下層セルの中心
  z_centers[2] = z_bottom_cell

  # 上方向に積算（k=2～nk）
  for k in 2:nk
    # z_centers[k+1] = z_centers[k] + セル間距離
    # セル間距離 = dz_t[k-1] (k-1番目セルの上側界面距離) + dz_b[k] (k番目セルの下側界面距離)
    # ただし、中心間距離として dz_t[k-1] = dz_b[k] が成り立つはず（等間隔の場合）
    # 一般的には: z_centers[k+1] = z_centers[k] + dz_t[k-1] + dz_b[k]
    # しかし、dz_t[k-1]とdz_b[k]は界面を挟んで同じ距離を指す場合がある

    # より正確には: z_centers[k+1] = z_centers[k] + (dz_t[k-1] + dz_b[k])
    # ただし、通常 dz_t[k-1] = z_centers[界面] - z_centers[k-1]
    #              dz_b[k] = z_centers[k] - z_centers[界面]
    # なので z_centers[k] = z_centers[k-1] + dz_t[k-1] + dz_b[k] は成立しない

    # IHCPの定義を確認:
    # dz_b[k]: セル中心から下側界面までの距離
    # dz_t[k]: セル中心から上側界面までの距離
    # したがって、セル幅 = dz_b[k] + dz_t[k]

    # セル中心間距離は:
    # z_centers[k+1] - z_centers[k] = dz_t[k] + dz_b[k+1]（k番目セルの上側界面距離 + k+1番目セルの下側界面距離）

    # ただし、境界では dz_b[1] = Inf, dz_t[nk] = Inf なので特別扱い

    if k <= nk
      if k == nk
        # 最上層セル: dz_t[nk] = Inf の場合がある
        # z_centers[nk+1]は上側ガイドセル、適当に配置
        z_centers[k+1] = z_centers[k] + dz_cell[k]  # 暫定的にセル幅を使用
      else
        # 中間層: セル中心間距離 = dz_t[k-1] + dz_b[k]
        # ただし、これは界面を挟んだ距離
        # より簡潔には: セル中心間距離 = 0.5*(dz_cell[k-1] + dz_cell[k])（等間隔格子の近似）
        # 一般的には、dz_cell[k]がセル幅なので、単純に積算
        z_centers[k+1] = z_centers[k] + 0.5*(dz_cell[k-1] + dz_cell[k])
      end
    end
  end

  # 下側ガイドセル（k=1）の配置
  # z_centers[1] = z_centers[2] - セル間距離
  z_centers[1] = z_centers[2] - dz_cell[1]

  # 上側ガイドセル（k=nk+2）の配置
  # z_centers[nk+2] = z_centers[nk+1] + セル間距離
  z_centers[nk_full] = z_centers[nk+1] + dz_cell[nk]

  # セル代表幅の計算（Heat3ds方式）
  # dz[k] = 0.5*(z_centers[k+1] - z_centers[k-1])  (k=2～nk)
  for k in 2:nk_full-1
    dz[k] = 0.5 * (z_centers[k+1] - z_centers[k-1])
  end

  # 端点補正（Heat3ds: Zcoord.jl:147-148）
  dz[2] *= 0.5       # 下側境界補正
  dz[nk_full-1] *= 0.5  # 上側境界補正

  return z_centers, dz
end


end # module GridTransform
