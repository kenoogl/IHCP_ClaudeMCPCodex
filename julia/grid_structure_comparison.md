# Z方向格子構造の比較：Heat3ds vs IHCP

**作成日**: 2025年10月3日
**目的**: ガイドセル方式（Heat3ds）への移行のため、両実装の格子定義を詳細比較

---

## 1. Heat3dsの格子構造（ガイドセル方式）

### 1.1 データ構造

```julia
# 配列サイズ: [ni+2, nj+2, nk+2]（ガイドセル含む）
θ[i, j, k]  # 温度配列
λ[i, j, k]  # 熱伝導率配列
mask[i, j, k]  # マスク配列（1.0=内点、0.0=境界）

# 計算範囲（内点のみ）
for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
    # ガイドセルを除く計算
end
```

### 1.2 Z方向格子配列

```julia
Z[k]::Vector{Float64}     # サイズ: (nk+2,)
ΔZ[k]::Vector{Float64}    # サイズ: (nk+1,)

# Z[k]: k番目のセル中心のz座標
# - Z[1]: 下側ガイドセル
# - Z[2]～Z[nk+1]: 計算内点
# - Z[nk+2]: 上側ガイドセル

# ΔZ[k]: k番目のセルの代表幅（体積分散用）
# 定義: ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])  (k=2～nk)
# 端点補正: ΔZ[2] *= 0.5, ΔZ[nk] *= 0.5
```

**生成関数（Zcoord.jl:137-151）**：
```julia
function genZ!(Z::Vector{Float64}, ΔZ::Vector{Float64}, SZ, ox, dz)
    mz = SZ[3]
    Zcase2!(Z, SZ)  # 不等間隔格子定義

    for k in 2:mz-1
        ΔZ[k] = 0.5*(Z[k+1] - Z[k-1])
    end

    ΔZ[2] = 0.5*ΔZ[2]      # 下側境界補正
    ΔZ[mz-1] = 0.5*ΔZ[mz-1]  # 上側境界補正
end
```

### 1.3 係数計算（調和平均 + マスク補正）

**λf関数（common.jl:22）**：
```julia
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))

# ma, mb: マスク値（1.0=内点、0.0=境界）
# 補正係数:
#   ma=1, mb=1 (両方内点): 2.0 - 1 = 1.0 → 通常の調和平均
#   ma=1, mb=0 (一方境界): 2.0 - 0 = 2.0 → 調和平均を2倍
#   ma=0, mb=0 (両方境界): 2.0 - 0 = 2.0 → （稀）
```

**Z方向係数計算（heat3d_NonUniform.jl:549-557）**：
```julia
# 下側界面
λ_b = λ[i, j, k-1]
m_b = mask[i, j, k-1]
azm = λf(λ0, λ_b, m0, m_b) / (dz_b * ΔZ[k])

# 上側界面
λ_t = λ[i, j, k+1]
m_t = mask[i, j, k+1]
azp = λf(λ0, λ_t, m0, m_t) / (dz_t * ΔZ[k])
```

ここで：
- `dz_b = Z[k] - Z[k-1]` （セル中心間距離）
- `dz_t = Z[k+1] - Z[k]`
- `ΔZ[k]`: セル体積の代表幅

---

## 2. IHCPの格子構造（現行実装）

### 2.1 データ構造

```julia
# 配列サイズ: [ni, nj, nk]（ガイドセルなし）
T[i, j, k]  # 温度配列
k_array[i, j, k]  # 熱伝導率配列

# 計算範囲（全領域）
for k_idx in 1:nk, j in 1:nj, i in 1:ni
    # if文で境界判定
end
```

### 2.2 Z方向格子配列

```julia
dz[k]::Vector{Float64}     # サイズ: (nk,)
dz_b[k]::Vector{Float64}   # サイズ: (nk,)
dz_t[k]::Vector{Float64}   # サイズ: (nk,)

# dz[k]: k番目のセル中心幅
# dz_b[k]: k番目のセル中心から下側界面までの距離
# dz_t[k]: k番目のセル中心から上側界面までの距離

# 境界条件の表現:
# - dz_b[1] = Inf  （底面熱流束境界）
# - dz_t[nk] = Inf （上面熱流束境界）
```

**等間隔格子の例（test_dhcp_solver.jl:48-52）**：
```julia
dz = fill(0.1e-3, nk)
dz_b = fill(0.1e-3, nk)
dz_b[1] = Inf
dz_t = fill(0.1e-3, nk)
dz_t[end] = Inf
```

### 2.3 係数計算（if文で境界処理）

**Z方向係数計算（DHCPSolver.jl:172-189）**：
```julia
# 下側界面
if k_idx == 1
    a_b[p] = 0.0  # 底面熱流束境界（後でRHSに加算）
else
    k_b = k[i, j, k_idx-1]
    k_p = k[i, j, k_idx]
    k_harmonic = 2.0 * k_p * k_b / (k_p + k_b)
    dz_b_k = dz_b[k_idx]
    a_b[p] = k_harmonic * dx * dy / dz_b_k
end

# 上側界面
if k_idx == nk
    a_t[p] = 0.0  # 上面断熱境界
else
    k_t = k[i, j, k_idx+1]
    k_p = k[i, j, k_idx]
    k_harmonic = 2.0 * k_p * k_t / (k_p + k_t)
    dz_t_k = dz_t[k_idx]
    a_t[p] = k_harmonic * dx * dy / dz_t_k
end
```

---

## 3. 格子構造の対応関係

### 3.1 配列サイズの対応

| 項目 | Heat3ds | IHCP |
|------|---------|------|
| 温度配列 | `[ni+2, nj+2, nk+2]` | `[ni, nj, nk]` |
| 計算内点数 | `ni × nj × nk` | `ni × nj × nk` |
| ガイドセル | あり（1層ずつ） | なし |
| 境界処理 | マスク配列+λf関数 | if文 |

### 3.2 Z方向格子配列の対応

| Heat3ds | IHCP | 説明 |
|---------|------|------|
| `Z[k]` | - | セル中心z座標（ガイドセル含む） |
| `ΔZ[k]` | `dz[k]` | セル代表幅 |
| `Z[k] - Z[k-1]` | `dz_b[k]` | 下側セル中心間距離 |
| `Z[k+1] - Z[k]` | `dz_t[k]` | 上側セル中心間距離 |

**重要な違い**：
- Heat3ds: `Z[k]`から界面距離を計算
- IHCP: `dz_b[k]`, `dz_t[k]`を直接保持

### 3.3 係数計算の対応

| 項目 | Heat3ds | IHCP |
|------|---------|------|
| 調和平均 | `λf(a, b, ma, mb)` | `2.0*k_p*k_b/(k_p+k_b)` |
| マスク補正 | `2.0-div(ma+mb,2)` | なし |
| 境界処理 | `mask=0.0`, `λ=0.0` | `if k_idx == 1` |
| 体積補正 | `/ ΔZ[k]` | 係数に含まれる |

---

## 4. ガイドセル方式への移行設計

### 4.1 新しい格子配列の定義

```julia
# 新しいデータ構造（Method 1: ガイドセル方式）
struct GridInfo
    # 内点数（計算領域）
    ni::Int
    nj::Int
    nk::Int

    # 配列サイズ（ガイドセル含む）
    ni_full::Int  # = ni + 2
    nj_full::Int  # = nj + 2
    nk_full::Int  # = nk + 2

    # XY方向格子（等間隔）
    dx::Float64
    dy::Float64

    # Z方向格子（不等間隔）
    Z::Vector{Float64}    # サイズ: (nk_full,) = セル中心z座標（ガイドセル含む）
    ΔZ::Vector{Float64}   # サイズ: (nk_full-1,) = セル代表幅
end
```

### 4.2 Z方向格子の生成関数

```julia
"""
IHCPの(dz, dz_b, dz_t)からHeat3ds形式の(Z, ΔZ)を生成
"""
function convert_to_guard_cell_grid(
    nk::Int,
    dz::Vector{Float64},
    dz_b::Vector{Float64},
    dz_t::Vector{Float64}
)
    nk_full = nk + 2
    Z = zeros(Float64, nk_full)
    ΔZ = zeros(Float64, nk_full - 1)

    # セル中心座標の計算（下から積算）
    Z[1] = 0.0  # 下側ガイドセル（仮）
    Z[2] = dz_b[1]  # 最下層セル中心

    for k in 2:nk
        Z[k+1] = Z[k] + dz_t[k-1]  # または Z[k] + dz_b[k]
    end

    Z[nk_full] = Z[nk_full-1] + dz_t[nk]  # 上側ガイドセル

    # セル代表幅の計算
    for k in 2:nk_full-1
        ΔZ[k] = 0.5 * (Z[k+1] - Z[k-1])
    end

    # 端点補正
    ΔZ[2] *= 0.5
    ΔZ[nk_full-1] *= 0.5

    return Z, ΔZ
end
```

### 4.3 マスク配列とλ配列の初期化

```julia
"""
境界条件に基づいてマスク配列とλ配列を初期化
"""
function initialize_guard_cells!(
    θ::Array{Float64,3},
    λ::Array{Float64,3},
    mask::Array{Float64,3},
    bc_params  # 境界条件パラメータ
)
    SZ = size(θ)

    # 全領域を内点として初期化
    mask .= 1.0

    # ガイドセルをマスク=0に設定
    # X方向
    mask[1, :, :] .= 0.0
    mask[SZ[1], :, :] .= 0.0

    # Y方向
    mask[:, 1, :] .= 0.0
    mask[:, SZ[2], :] .= 0.0

    # Z方向
    mask[:, :, 1] .= 0.0      # 下側（熱流束境界）
    mask[:, :, SZ[3]] .= 0.0  # 上側（断熱境界）

    # 熱流束境界条件: λ=0.0に設定
    λ[:, :, 1] .= 0.0   # 下側ガイドセル
    λ[:, :, SZ[3]] .= 0.0  # 上側ガイドセル（断熱）

    # 内点のλは熱物性値で設定（別途）
end
```

### 4.4 係数計算（マトリックスフリー版）

```julia
"""
Heat3ds風のマトリックスフリーA*x演算
"""
function apply_dhcp_operator!(
    ap::Array{Float64,3},
    p::Array{Float64,3},
    λ::Array{Float64,3},
    mask::Array{Float64,3},
    dx::Float64, dy::Float64,
    Z::Vector{Float64},
    ΔZ::Vector{Float64},
    rho::Float64,
    cp::Array{Float64,3},
    dt::Float64
)
    SZ = size(p)

    # 計算範囲: 内点のみ（ガイドセル除く）
    for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i, j, k]
        m0 = mask[i, j, k]

        # Z方向の界面距離
        dz_b = Z[k] - Z[k-1]
        dz_t = Z[k+1] - Z[k]

        # X方向係数（調和平均 + マスク補正）
        λ_w = λ[i-1, j, k]
        m_w = mask[i-1, j, k]
        axm = λf(λ0, λ_w, m0, m_w) * dy * ΔZ[k] / dx

        λ_e = λ[i+1, j, k]
        m_e = mask[i+1, j, k]
        axp = λf(λ0, λ_e, m0, m_e) * dy * ΔZ[k] / dx

        # Y方向係数
        λ_s = λ[i, j-1, k]
        m_s = mask[i, j-1, k]
        aym = λf(λ0, λ_s, m0, m_s) * dx * ΔZ[k] / dy

        λ_n = λ[i, j+1, k]
        m_n = mask[i, j+1, k]
        ayp = λf(λ0, λ_n, m0, m_n) * dx * ΔZ[k] / dy

        # Z方向係数
        λ_b = λ[i, j, k-1]
        m_b = mask[i, j, k-1]
        azm = λf(λ0, λ_b, m0, m_b) * dx * dy / (dz_b * ΔZ[k])

        λ_t = λ[i, j, k+1]
        m_t = mask[i, j, k+1]
        azp = λf(λ0, λ_t, m0, m_t) * dx * dy / (dz_t * ΔZ[k])

        # 対角項（時間項 + 空間項）
        vol = dx * dy * ΔZ[k]
        dd = rho * cp[i,j,k] * vol / dt + (axp + axm + ayp + aym + azp + azm)

        # 7点ステンシル演算
        ss = ( axp * p[i+1, j  , k  ] + axm * p[i-1, j  , k  ]
             + ayp * p[i  , j+1, k  ] + aym * p[i  , j-1, k  ]
             + azp * p[i  , j  , k+1] + azm * p[i  , j  , k-1] )

        # マスクを適用（境界点は計算しない）
        ap[i, j, k] = (dd * p[i, j, k] - ss) * m0
    end
end

# λf関数（Heat3ds互換）
function λf(a::Float64, b::Float64, ma::Float64, mb::Float64)
    if a + b == 0.0
        return 0.0
    end
    correction = 2.0 - div(ma + mb, 2)
    return 2.0 * a * b / (a + b) * correction
end
```

---

## 5. 移行に伴う主要変更点

### 5.1 データ構造の変更

| 項目 | 変更前 | 変更後 |
|------|--------|--------|
| 温度配列 | `T[ni, nj, nk]` | `T[ni+2, nj+2, nk+2]` |
| 熱伝導率 | `k[ni, nj, nk]` | `λ[ni+2, nj+2, nk+2]` |
| マスク配列 | なし | `mask[ni+2, nj+2, nk+2]` |
| Z格子 | `(dz, dz_b, dz_t)` | `(Z, ΔZ)` |

### 5.2 インデックス範囲の変更

```julia
# 変更前
for k_idx in 1:nk, j in 1:nj, i in 1:ni
    # 境界をif文で判定
end

# 変更後
for k in 2:nk+1, j in 2:nj+1, i in 2:ni+1
    # ガイドセルを除く範囲
end
```

### 5.3 境界条件の設定方法

**変更前（if文）**：
```julia
if k_idx == 1
    a_b[p] = 0.0  # 底面熱流束
end
```

**変更後（マスク+λ=0）**：
```julia
# 初期化時に設定
mask[:, :, 1] .= 0.0  # 下側ガイドセル
λ[:, :, 1] .= 0.0

# 係数計算では自動的にゼロになる
azm = λf(λ0, λ_b, m0, m_b) / (dz_b * ΔZ[k])  # λ_b=0.0 → azm=0.0
```

---

## 6. テスト戦略

### 6.1 段階的検証

1. **Phase A: 格子変換関数のテスト**
   - `convert_to_guard_cell_grid`の動作確認
   - 等間隔格子で`Z`, `ΔZ`が正しく生成されるか

2. **Phase B: 単一時間ステップのDHCP**
   - 等間隔格子、小規模問題（5×5×3）
   - 既存実装との数値一致確認

3. **Phase C: 不等間隔格子のDHCP**
   - Z方向不等間隔格子での動作確認
   - Python参照データとの比較

4. **Phase D: 全ソルバーの統合**
   - Adjoint, CGM, Sliding Windowの順に移行
   - 全505テストの再実行

### 6.2 数値精度の確認

- **許容誤差**: 相対誤差 < 1e-12（機械精度レベル）
- **境界条件**: 熱流束が正しく適用されるか
- **マスク効果**: 境界点が計算から除外されるか

---

## 7. 実装スケジュール

| Phase | 作業内容 | 期間 | リスク |
|-------|---------|------|--------|
| **Phase A** | 格子変換関数実装・テスト | 1週 | 低 |
| **Phase B** | DHCP単一ステップ（等間隔） | 2週 | 中 |
| **Phase C** | DHCP不等間隔格子 | 1週 | 中 |
| **Phase D** | Adjoint移行 | 2週 | 高 |
| **Phase E** | CGM移行 | 1週 | 中 |
| **Phase F** | Sliding Window移行 | 1週 | 低 |
| **Phase G** | 全テスト統合 | 1週 | 低 |
| **合計** | | **9週** | |

---

## 8. まとめ

### 8.1 主要な技術的課題

1. **Z方向格子の変換**
   - `(dz, dz_b, dz_t)` → `(Z, ΔZ)` の正確な変換
   - 不等間隔格子での界面距離の計算

2. **マスク配列の境界条件適用**
   - 熱流束境界: `mask=0.0`, `λ=0.0`
   - 断熱境界: `mask=0.0`, `λ=0.0`

3. **既存テストとの互換性**
   - 505個全テストを再実行
   - 数値精度が維持されるか確認

### 8.2 期待される効果

1. **性能向上**: 15-50%の高速化（疎行列組み立て削減）
2. **コードの簡潔化**: if文削減、統一的な境界処理
3. **保守性向上**: Heat3dsとの共通化

---

**次のステップ**: Phase Aの実装開始（格子変換関数とテスト）
