# ガイドセル方式実装サマリー

**作成日**: 2025年10月3日
**ブランチ**: tuning2
**目的**: Heat3ds形式のガイドセル方式をIHCPに導入し、マトリックスフリーBiCGstab実装の基盤を構築

---

## 1. 実装完了項目

### 1.1 GridTransform.jl モジュール

**ファイルパス**: `julia/src/utils/GridTransform.jl`

**実装機能**:

#### (1) BoundaryType列挙型定義
```julia
@enum BoundaryType begin
  ISOTHERMAL   # 等温条件 (Dirichlet)
  HEAT_FLUX    # 熱流束条件 (Neumann)
  CONVECTION   # 熱伝達条件 (Robin)
end
```

#### (2) convert_to_guard_cell_grid関数
```julia
function convert_to_guard_cell_grid(nk, dz, dz_b, dz_t) -> (Z, ΔZ)
```
- **目的**: IHCP格子定義`(dz, dz_b, dz_t)`からHeat3ds形式`(Z, ΔZ)`への変換
- **入力**:
  - `nk`: 計算内点数
  - `dz`: セル中心幅 (nk,)
  - `dz_b`: セル中心から下側界面までの距離 (nk,)
  - `dz_t`: セル中心から上側界面までの距離 (nk,)
- **出力**:
  - `Z`: セル中心z座標（ガイドセル含む） (nk+2,)
  - `ΔZ`: セル代表幅 (nk+1,)
- **状態**: ⚠️ 格子座標計算ロジックの検証が必要

#### (3) initialize_guard_cells!関数
```julia
function initialize_guard_cells!(θ, λ, mask, θ_init, λ_init)
```
- **目的**: ガイドセルを含む配列の初期化
- **機能**:
  - マスク配列を全て1.0（内点）に初期化
  - 内点データ`(ni,nj,nk)`をガイドセル配列`(ni+2,nj+2,nk+2)`の内部にコピー
  - ガイドセル部分は未初期化（境界条件設定で上書き）
- **状態**: ✅ 実装完了

#### (4) compute_z_range関数
```julia
function compute_z_range(nk, z_minus_bc, z_plus_bc) -> Vector{Int}
```
- **目的**: 境界条件に基づいてZ方向の計算範囲を決定
- **入力**:
  - `nk`: 計算内点数（配列サイズはnk+2）
  - `z_minus_bc`: Z-方向（下側）の境界条件タイプ
  - `z_plus_bc`: Z+方向（上側）の境界条件タイプ
- **出力**: `[z_start, z_end]` 計算範囲のインデックス
- **決定ルール**:
  - Z-方向が等温条件 → `z_start = 3` （k=2は温度固定、計算対象外）
  - Z-方向が熱流束/熱伝達 → `z_start = 2` （k=2も計算対象）
  - Z+方向が等温条件 → `z_end = nk` （k=nk+1は温度固定、計算対象外）
  - Z+方向が熱流束/熱伝達 → `z_end = nk+1` （k=nk+1も計算対象）
- **例**:
  ```julia
  compute_z_range(10, ISOTHERMAL, HEAT_FLUX)  # [3, 11]
  compute_z_range(10, HEAT_FLUX, ISOTHERMAL)  # [2, 10]
  compute_z_range(10, HEAT_FLUX, CONVECTION)  # [2, 11]
  ```
- **状態**: ✅ 実装完了

#### (5) λf関数（調和平均 + マスク補正）
```julia
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))
```
- **目的**: Heat3ds互換の調和平均関数（マスク補正付き）
- **マスク補正係数**:
  - `ma=1, mb=1` (両方内点): `2.0 - 1 = 1.0` → 通常の調和平均
  - `ma=1, mb=0` (一方境界): `2.0 - 0 = 2.0` → 調和平均を2倍
  - `ma=0, mb=0` (両方境界): `2.0 - 0 = 2.0` → 2倍
- **状態**: ✅ 実装完了（Heat3ds common.jl:22と同一形式）

---

### 1.2 test_grid_transform.jl テストスイート

**ファイルパス**: `julia/test/test_grid_transform.jl`

**実装テストケース**:

#### Test 1: 等間隔格子の変換（3層）
- **目的**: `convert_to_guard_cell_grid`の基本動作確認
- **検証項目**:
  - 配列サイズ: `Z` (nk+2=5), `ΔZ` (nk+1=4)
  - セル中心座標の妥当性
  - セル代表幅の計算
- **状態**: ⚠️ 実行待ち（格子座標計算ロジックの検証が必要）

#### Test 2: ガイドセル配列の初期化
- **目的**: `initialize_guard_cells!`の動作確認
- **検証項目**:
  - マスク配列が全て1.0に初期化されるか
  - 内点データ`(3,3,3)`がガイドセル配列`(5,5,5)`の内部に正しくコピーされるか
- **状態**: ✅ テスト作成完了

#### Test 3: λf関数（調和平均 + マスク補正）
- **目的**: λf関数の動作確認
- **検証項目**:
  - 通常の調和平均（両方内点）
  - 一方が境界の場合の補正
  - λ=0の場合の挙動
- **状態**: ✅ テスト作成完了

#### Test 4: compute_z_range（境界条件に応じた計算範囲）
- **目的**: `compute_z_range`の全パターン検証
- **検証項目**:
  - パターン1: `(ISOTHERMAL, HEAT_FLUX)` → `[3, 11]`
  - パターン2: `(HEAT_FLUX, ISOTHERMAL)` → `[2, 10]`
  - パターン3: `(HEAT_FLUX, CONVECTION)` → `[2, 11]`
  - パターン4: `(ISOTHERMAL, ISOTHERMAL)` → `[3, 10]`
  - パターン5: `(CONVECTION, CONVECTION)` → `[2, 11]`
- **状態**: ✅ テスト作成完了

---

## 2. 配列構造の整理

### 2.1 Heat3ds形式（ガイドセル方式）

```julia
# 配列サイズ
MX = NX + 2  # X方向（ガイドセル含む）
MY = NY + 2  # Y方向（ガイドセル含む）
MZ = NZ + 2  # Z方向（ガイドセル含む）

# 配列定義
θ[MX, MY, MZ]     # 温度配列
λ[MX, MY, MZ]     # 熱伝導率配列
mask[MX, MY, MZ]  # マスク配列（1.0=内点、0.0=境界）

# Z方向格子配列
Z[MZ]      # セル中心z座標（ガイドセル含む）
ΔZ[MZ-1]   # セル代表幅
```

### 2.2 計算範囲（インデックス）

**X, Y方向**:
```julia
for j in 2:MY-1, i in 2:MX-1
  # 計算内点（ガイドセル除く）
end
```

**Z方向**:
```julia
# z_rangeで動的に決定
for k in z_range[1]:z_range[2]
  # 境界条件に応じた計算範囲
end
```

**配列インデックス構造**（Z方向、MZ = nk+2）:
```
k=1:           下側ガイドセル
k=2:           下側境界セル（等温なら計算対象外、熱流束なら計算対象）
k=3～nk+1:     計算内点
k=nk+2:        上側ガイドセル
```

---

## 3. 境界条件とz_rangeの関係

### 3.1 Heat3dsの例（heat3d_nu.jl:211-212）

```julia
# Mode3: Z-方向がPCB温度（等温相当）、Z+方向が熱伝達
z_range[1] = 3           # 等温条件のため内点から開始
z_range[2] = SZ[3]-1     # 熱伝達条件のため境界セルまで計算
```

### 3.2 全パターン一覧（nk=10の場合）

| Z-方向 | Z+方向 | z_start | z_end | 計算範囲 | 用途例 |
|--------|--------|---------|-------|----------|--------|
| ISOTHERMAL | HEAT_FLUX | 3 | 11 | [3:11] | IHCP典型問題（底面等温、上面熱流束） |
| HEAT_FLUX | ISOTHERMAL | 2 | 10 | [2:10] | 逆問題（底面熱流束、上面等温） |
| HEAT_FLUX | CONVECTION | 2 | 11 | [2:11] | Heat3ds Mode3（PCB、上面熱伝達） |
| ISOTHERMAL | ISOTHERMAL | 3 | 10 | [3:10] | 両端固定問題 |
| CONVECTION | CONVECTION | 2 | 11 | [2:11] | 両端熱伝達 |

---

## 4. マスク配列と境界条件

### 4.1 マスク配列の役割

```julia
mask[i, j, k] = 1.0  # 内点（計算対象）
mask[i, j, k] = 0.0  # 境界点（計算対象外、または特別扱い）
```

### 4.2 境界条件の実装（Heat3ds boundary_conditions.jl）

**等温条件**:
```julia
mask[boundary] = 0.0      # 境界点を計算対象外に
λ[boundary] = λ[interior]  # 内点側の値と同じ
θ[boundary] = T_specified  # 指定温度に固定
```

**熱流束条件**:
```julia
mask[boundary] = 0.0  # 境界点を計算対象外に
λ[boundary] = 0.0     # 熱伝導率をゼロに（自動的に伝導項がゼロ）
```

**熱伝達条件**:
```julia
mask[boundary] = 0.0       # 境界点を計算対象外に
λ[boundary] = 0.0          # 熱伝導率をゼロに
θ[boundary] = T_ambient    # 周囲温度に設定
```

### 4.3 λf関数の境界処理

```julia
# 係数計算（Z方向上側界面の例）
λ_t = λ[i, j, k+1]
m_t = mask[i, j, k+1]
azp = λf(λ0, λ_t, m0, m_t) * dx * dy / (dz_t * ΔZ[k])

# λ[i, j, k+1] = 0.0（熱流束境界）の場合:
# λf(λ0, 0.0, 1.0, 0.0) = 0.0 → azp = 0.0
# つまり、境界での伝導項が自動的にゼロになる
```

---

## 5. 次のステップ

### 5.1 短期目標（Phase A: 格子変換機能の検証）

- [ ] **convert_to_guard_cell_gridの格子座標計算ロジックを検証**
  - 等間隔格子でZ, ΔZが正しく生成されるか
  - 不等間隔格子での動作確認
  - セル中心間距離の計算方法の見直し（現在の実装は暫定的）

- [ ] **test_grid_transform.jlの実行**
  - 全4テストケースの動作確認
  - 数値精度の検証

### 5.2 中期目標（Phase B: マトリックスフリーDHCPソルバー）

- [ ] **apply_dhcp_operator!関数の実装**
  - Heat3ds CalcAX!を参考にしたマトリックスフリーA*x演算
  - ガイドセル配列`(ni+2, nj+2, nk+2)`対応
  - z_rangeを使った計算範囲制御

- [ ] **solve_dhcp_matrix_free!関数の実装**
  - LinearOperatorsパッケージの活用
  - BiCGstab法の統合（IterativeSolvers.jl）
  - 前処理器の検討（対角前処理）

- [ ] **小規模問題での数値検証**
  - 既存solve_dhcp!との結果比較
  - 数値精度: 相対誤差 < 1e-12

### 5.3 長期目標（Phase C～G: 全ソルバーの移行）

- [ ] **Phase C**: Adjoint随伴ソルバーのマトリックスフリー化
- [ ] **Phase D**: CGM共役勾配法の統合
- [ ] **Phase E**: Sliding Windowソルバーの統合
- [ ] **Phase F**: 全505テストの再実行と検証
- [ ] **Phase G**: 性能評価とチューニング

---

## 6. 期待される効果

### 6.1 性能向上

**現状のボトルネック**（tuning2プロファイル結果より）:
- 疎行列CG法が全体の**90%以上**を占める
- 疎行列組み立て（`assemble_dhcp_matrix`）のコスト

**期待される改善**:
- 疎行列組み立て削減 → **15-50%の高速化**
- SIMD最適化（if文削減、連続メモリアクセス）→ 追加10-20%
- 合計: **25-70%の性能向上**

### 6.2 コードの簡潔化

**現行方式**（if文による境界処理）:
```julia
if k_idx == 1
  a_b[p] = 0.0  # 底面熱流束
else
  k_b = k[i, j, k_idx-1]
  k_p = k[i, j, k_idx]
  a_b[p] = 2.0 * k_p * k_b / (k_p + k_b) * dx * dy / dz_b[k_idx]
end
```

**ガイドセル方式**（マスク配列による統一処理）:
```julia
# 境界条件設定時に λ[boundary] = 0.0 に設定済み
λ_b = λ[i, j, k-1]
m_b = mask[i, j, k-1]
azm = λf(λ0, λ_b, m0, m_b) * dx * dy / (dz_b * ΔZ[k])
# λ_b = 0.0 なら自動的に azm = 0.0
```

**利点**:
- if文が不要（SIMD最適化に有利）
- 境界条件の種類に依存しない統一的な処理
- Heat3dsとの共通化（保守性向上）

---

## 7. 実装状況まとめ

| 項目 | 状態 | 備考 |
|------|------|------|
| **GridTransform.jl** | ✅ 実装完了 | λfをcommon.jl形式に修正済み |
| **BoundaryType定義** | ✅ 実装完了 | ISOTHERMAL, HEAT_FLUX, CONVECTION |
| **convert_to_guard_cell_grid** | ⚠️ 検証待ち | 格子座標計算ロジックの見直しが必要 |
| **initialize_guard_cells!** | ✅ 実装完了 | |
| **compute_z_range** | ✅ 実装完了 | 全5パターン対応 |
| **λf関数** | ✅ 実装完了 | Heat3ds互換 |
| **test_grid_transform.jl** | ✅ 作成完了 | 4テストケース、実行待ち |
| **マトリックスフリーDHCP** | ⏳ 未着手 | 次フェーズ |

---

## 8. 技術的課題

### 8.1 格子座標計算ロジック（convert_to_guard_cell_grid）

**現状の実装**:
```julia
# セル中心間距離の計算（暫定的）
Z[k+1] = Z[k] + 0.5*(dz[k-1] + dz[k])
```

**問題点**:
- IHCP格子定義`(dz, dz_b, dz_t)`との対応が不明確
- 不等間隔格子での動作が未検証

**対策**:
- 等間隔格子での動作確認（テスト実行）
- Python参照データとの比較
- Heat3ds Zcoord.jl（genZ!関数）との対応確認

### 8.2 境界条件モジュールの統合

**現状**:
- GridTransform.jl内で独自にBoundaryType定義
- Heat3dsのBoundaryConditions.jlは別プロジェクト

**今後の方針**:
- IHCP専用のBoundaryConditionsモジュールを作成するか
- GridTransform.jlに境界条件設定機能を統合するか
- 検討が必要

---

**作成者**: Claude Code
**最終更新**: 2025年10月3日
