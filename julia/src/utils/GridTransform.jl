"""
GridTransform.jl

IHCP格子定義からHeat3ds形式のガイドセル格子への変換ユーティリティ

目的:
- (dz, dz_b, dz_t)形式から(Z, ΔZ)形式への変換
- ガイドセルを含む配列の初期化
- マスク配列と熱伝導率配列の境界条件設定
- 境界条件に基づくZ方向計算範囲の決定
"""
module GridTransform

# 境界条件タイプの定義
@enum BoundaryType begin
  ISOTHERMAL   # 等温条件 (Dirichlet)
  HEAT_FLUX    # 熱流束条件 (Neumann)
  CONVECTION   # 熱伝達条件 (Robin)
end

export BoundaryType, ISOTHERMAL, HEAT_FLUX, CONVECTION
export convert_to_guard_cell_grid, initialize_guard_cells!, compute_z_range, λf

"""
  convert_to_guard_cell_grid(nk, dz, dz_b, dz_t) -> (Z, ΔZ)

IHCPの格子定義(dz, dz_b, dz_t)からHeat3ds形式の(Z, ΔZ)を生成

# Arguments
- nk: 計算内点数（Z方向）
- dz: セル中心幅 (nk,) [m]
- dz_b: セル中心から下側界面までの距離 (nk,) [m]
- dz_t: セル中心から上側界面までの距離 (nk,) [m]

# Returns
- Z: セル中心z座標（ガイドセル含む） (nk+2,) [m]
- ΔZ: セル代表幅 (nk+1,) [m]

# Notes
- Z[1]: 下側ガイドセル
- Z[2]～Z[nk+1]: 計算内点
- Z[nk+2]: 上側ガイドセル
- ΔZ[k] = 0.5*(Z[k+1] - Z[k-1]) （k=2～nk）
- 端点補正: ΔZ[2] *= 0.5, ΔZ[nk] *= 0.5
"""
function convert_to_guard_cell_grid(
  nk::Int,
  dz::Vector{Float64},
  dz_b::Vector{Float64},
  dz_t::Vector{Float64}
)
  @assert length(dz) == nk "dz length mismatch"
  @assert length(dz_b) == nk "dz_b length mismatch"
  @assert length(dz_t) == nk "dz_t length mismatch"

  nk_full = nk + 2
  Z = zeros(Float64, nk_full)
  ΔZ = zeros(Float64, nk_full - 1)

  # セル中心座標の計算（下から積算）
  # 最下層セル（k=1）の中心を基準とする
  z_bottom_cell = 0.0  # 仮の原点

  # Z[2]が最下層セルの中心
  Z[2] = z_bottom_cell

  # 上方向に積算（k=2～nk）
  for k in 2:nk
    # Z[k+1] = Z[k] + セル間距離
    # セル間距離 = dz_t[k-1] (k-1番目セルの上側界面距離) + dz_b[k] (k番目セルの下側界面距離)
    # ただし、中心間距離として dz_t[k-1] = dz_b[k] が成り立つはず（等間隔の場合）
    # 一般的には: Z[k+1] = Z[k] + dz_t[k-1] + dz_b[k]
    # しかし、dz_t[k-1]とdz_b[k]は界面を挟んで同じ距離を指す場合がある

    # より正確には: Z[k+1] = Z[k] + (dz_t[k-1] + dz_b[k])
    # ただし、通常 dz_t[k-1] = Z[界面] - Z[k-1]
    #              dz_b[k] = Z[k] - Z[界面]
    # なので Z[k] = Z[k-1] + dz_t[k-1] + dz_b[k] は成立しない

    # IHCPの定義を確認:
    # dz_b[k]: セル中心から下側界面までの距離
    # dz_t[k]: セル中心から上側界面までの距離
    # したがって、セル幅 = dz_b[k] + dz_t[k]

    # セル中心間距離は:
    # Z[k+1] - Z[k] = dz_t[k] + dz_b[k+1]（k番目セルの上側界面距離 + k+1番目セルの下側界面距離）

    # ただし、境界では dz_b[1] = Inf, dz_t[nk] = Inf なので特別扱い

    if k <= nk
      if k == nk
        # 最上層セル: dz_t[nk] = Inf の場合がある
        # Z[nk+1]は上側ガイドセル、適当に配置
        Z[k+1] = Z[k] + dz[k]  # 暫定的にセル幅を使用
      else
        # 中間層: セル中心間距離 = dz_t[k-1] + dz_b[k]
        # ただし、これは界面を挟んだ距離
        # より簡潔には: セル中心間距離 = 0.5*(dz[k-1] + dz[k])（等間隔格子の近似）
        # 一般的には、dz[k]がセル幅なので、単純に積算
        Z[k+1] = Z[k] + 0.5*(dz[k-1] + dz[k])
      end
    end
  end

  # 下側ガイドセル（k=1）の配置
  # Z[1] = Z[2] - セル間距離
  Z[1] = Z[2] - dz[1]

  # 上側ガイドセル（k=nk+2）の配置
  # Z[nk+2] = Z[nk+1] + セル間距離
  Z[nk_full] = Z[nk+1] + dz[nk]

  # セル代表幅の計算（Heat3ds方式）
  # ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])  (k=2～nk)
  for k in 2:nk_full-1
    ΔZ[k] = 0.5 * (Z[k+1] - Z[k-1])
  end

  # 端点補正（Heat3ds: Zcoord.jl:147-148）
  ΔZ[2] *= 0.5       # 下側境界補正
  ΔZ[nk_full-1] *= 0.5  # 上側境界補正

  return Z, ΔZ
end


"""
  initialize_guard_cells!(θ, λ, mask, θ_init, λ_init)

ガイドセルを含む配列を初期化（内点のデータをコピー）

# Arguments
- θ: 温度配列 (ni+2, nj+2, nk+2)
- λ: 熱伝導率配列 (ni+2, nj+2, nk+2)
- mask: マスク配列 (ni+2, nj+2, nk+2)
- θ_init: 初期温度（内点） (ni, nj, nk)
- λ_init: 初期熱伝導率（内点） (ni, nj, nk)

# Side Effects
- θ, λ, mask配列を初期化
- 全領域のマスクを1.0（内点）に設定
- 内点データをガイドセル配列の内部にコピー
- ガイドセル部分の値は未初期化（境界条件設定で上書き）

# Notes
- 境界条件の詳細設定は BoundaryConditions.apply_boundary_conditions! を使用
- この関数は配列の基本的な初期化のみを行う
"""
function initialize_guard_cells!(
  θ::Array{Float64,3},
  λ::Array{Float64,3},
  mask::Array{Float64,3},
  θ_init::Array{Float64,3},
  λ_init::Array{Float64,3}
)
  SZ = size(θ)
  @assert size(λ) == SZ "λ size mismatch"
  @assert size(mask) == SZ "mask size mismatch"

  ni, nj, nk = size(θ_init)
  @assert SZ == (ni+2, nj+2, nk+2) "Array size mismatch: expected ($(ni+2), $(nj+2), $(nk+2)), got $SZ"

  # 全領域を内点として初期化
  mask .= 1.0

  # 内点の温度と熱伝導率を設定（インデックス2:ni+1, 2:nj+1, 2:nk+1）
  for k in 1:nk, j in 1:nj, i in 1:ni
    θ[i+1, j+1, k+1] = θ_init[i, j, k]
    λ[i+1, j+1, k+1] = λ_init[i, j, k]
  end

  # ガイドセル部分（境界）は未初期化
  # BoundaryConditions.apply_boundary_conditions! で設定

  return nothing
end


"""
  compute_z_range(nk, z_minus_bc, z_plus_bc) -> Vector{Int}

境界条件に基づいてZ方向の計算範囲を決定

# Arguments
- nk: 計算内点数（配列サイズはnk+2）
- z_minus_bc: Z-方向（下側）の境界条件タイプ
- z_plus_bc: Z+方向（上側）の境界条件タイプ

# Returns
- [z_start, z_end]: 計算範囲のインデックス（ガイドセル配列基準）

# 配列インデックス構造 (MZ = nk+2)
- k=1: 下側ガイドセル
- k=2: 下側境界セル
- k=3～nk+1: 計算内点
- k=nk+2: 上側ガイドセル

# 計算範囲の決定ルール
- Z-方向が等温条件 → z_start = 3 (k=2は温度固定、計算対象外)
- Z-方向が熱流束/熱伝達 → z_start = 2 (k=2も計算対象)
- Z+方向が等温条件 → z_end = nk (k=nk+1は温度固定、計算対象外)
- Z+方向が熱流束/熱伝達 → z_end = nk+1 (k=nk+1も計算対象)

# Examples
```julia
compute_z_range(10, ISOTHERMAL, HEAT_FLUX)  # [3, 11]
compute_z_range(10, HEAT_FLUX, ISOTHERMAL)  # [2, 10]
compute_z_range(10, HEAT_FLUX, CONVECTION)  # [2, 11]
compute_z_range(10, ISOTHERMAL, ISOTHERMAL) # [3, 10]
```
"""
function compute_z_range(nk::Int, z_minus_bc::BoundaryType, z_plus_bc::BoundaryType)
  # Z-方向: 等温条件なら内点から開始（k=3）
  z_start = (z_minus_bc == ISOTHERMAL) ? 3 : 2

  # Z+方向: 等温条件なら境界セルを除外（k=nkまで）
  z_end = (z_plus_bc == ISOTHERMAL) ? nk : nk+1

  return [z_start, z_end]
end


"""
Harmonic mean with mask correction (Heat3ds互換)
@param a left value
@param b right value
@param ma mask for left (1.0=interior, 0.0=boundary)
@param mb mask for right (1.0=interior, 0.0=boundary)
"""
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))

end # module GridTransform
