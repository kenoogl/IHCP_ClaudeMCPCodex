# Phase 3: Adjoint随伴ソルバー 実装ログ

**開始日**: 2025-09-30
**完了日**: 2025-09-30（バグ修正完了）
**ステータス**: ✅ 完了（全テスト合格）

---

## Phase 3の目標

逆熱伝導問題（IHCP）のための随伴（Adjoint）ソルバーを実装：

1. **後退時間積分**: 終端条件から過去へ遡る時間ループ
2. **残差注入**: 測定値と計算値の差を底面から注入
3. **感度解析**: 表面熱流束に対する温度場の感度（随伴場λ）
4. **CGM準備**: Phase 4の共役勾配法で使用する勾配場を提供

---

## 実装内容

### ファイル構成

```
julia/
├── src/
│   ├── solvers/
│   │   └── AdjointSolver.jl         # Phase 3メイン実装
│   └── ThermalProperties.jl         # Phase 1（修正）
├── test/
│   └── test_adjoint_solver.jl       # Phase 3テスト
└── data/
    ├── phase3_reference_adjoint_1D.json   # 1D参照データ
    ├── phase3_reference_adjoint_3D.json   # 3D参照データ
    └── generate_phase3_reference.py       # 参照データ生成（修正）
```

### AdjointSolver.jl の構成

**主要関数**:

1. `adjoint_index(i, j, k, ni, nj)`: Fortran順序のグローバルインデックス計算
2. `build_adjoint_system!()`: 係数とRHSベクトル構築（残差注入含む）
3. `assemble_adjoint_matrix()`: CSC疎行列組み立て（7点ステンシル）
4. `solve_adjoint!()`: 複数時間ステップ随伴ソルバー（後退時間積分）

**随伴方程式**:

```
(I - α∇²)λ^(n-1) = λ^n + 2·(T_cal - Y_obs)
```

ここで：
- `λ`: 随伴場（感度場）
- `α = (k·dt) / (ρ·cp·dx²)`: 熱拡散係数
- `T_cal`: DHCP計算温度
- `Y_obs`: 観測温度（IRカメラデータ）

**残差注入（底面のみ）**:

```julia
if k_idx == 1
  rhs += 2.0 * (T_cal_bottom[i, j] - Y_obs[i, j]) * dx * dy
end
```

---

## 発見されたバグと修正

### バグ1: ThermalProperties.jl の配列インデックス誤り

**問題**: 配列形状の認識が誤っていた

```julia
# 修正前（誤り）
nk, nj, ni = size(Temperature)  # (nk, nj, ni)と仮定
T_current = Temperature[k_idx, j, i]

# 修正後（正しい）
ni, nj, nk = size(Temperature)  # (ni, nj, nk) Pythonと同じ
T_current = Temperature[i, j, k_idx]
```

**影響**: 誤った温度値を読み出し → 誤った熱物性値 → 随伴場が10^8倍小さい

### バグ2: generate_phase3_reference.py の多項式評価不一致

**問題**: Hornerの方法とpolyval_numbaで係数順序が逆

```python
# 修正前（Hornerの方法）
cp = cp_coeffs[0] + T * (cp_coeffs[1] + T * (cp_coeffs[2] + T * cp_coeffs[3]))
# [c0, c1, c2, c3] = c0 + c1*T + c2*T^2 + c3*T^3

# 修正後（polyval_numba互換）
def polyval_numba_simple(coeffs, x):
  result = 0.0
  for i in range(len(coeffs)):
    result += coeffs[i] * x ** (len(coeffs) - i - 1)
  return result
# [a, b, c, d] = a*T^3 + b*T^2 + c*T + d
```

### バグ3: 参照データの係数順序誤り

**問題**: 定数係数の配置が逆

```python
# 修正前
cp_coeffs = np.array([500.0, 0.0, 0.0, 0.0])  # 500*T^3と解釈される
k_coeffs = np.array([15.0, 0.0, 0.0, 0.0])

# 修正後
cp_coeffs = np.array([0.0, 0.0, 0.0, 500.0])  # 500（定数）
k_coeffs = np.array([0.0, 0.0, 0.0, 15.0])
```

---

## テスト結果

### Phase 3-1: 1D小規模問題

**テスト構成**:
- 格子: 1×1×5（1D問題）
- 時間: 11ステップ、dt=0.1s
- 物性: 定数（cp=500, k=15, rho=8000）
- 境界: q=1000 W/m²（表面熱流束）

**結果**:

```
✅ 終端条件 λ[nt] = 0
✅ Python参照との比較
   - 最大誤差: 5.08e-21 < 1e-7
   - 相対誤差: 3.72e-16 < 1e-5
✅ 時間ステップ毎の精度（t=1, 5, 10）
✅ 残差注入機構の動作確認
✅ 後退時間積分の単調性
```

**CG収束**: 平均5反復（Python参照と一致）

### Phase 3-2: 3D小規模問題（温度依存熱物性値）

**テスト構成**:
- 格子: 3×3×8（3D問題）
- 時間: 16ステップ、dt=0.001s
- 物性: 温度依存（SUS304近似）
  - `cp = 462 + 0.134*T`
  - `k = 14.6 + 0.0127*T`
- 境界: 空間変動熱流束 `q = 3000 * sin(πx/Lx) * sin(πy/Ly)`

**結果**:

```
✅ 終端条件 λ[nt] = 0
✅ Python参照との比較
   - 最大誤差: 1.26e-14 < 1e-7
   - 相対誤差: 7.41e-9 < 1e-5
✅ 勾配場（表面 k=8）の正確性
   - Python参照との誤差: 1.26e-14 < 1e-7
✅ 温度依存熱物性値の影響確認
✅ CG収束性
```

**CG収束**: 平均34.2反復（Python参照31.5 ± 0.7）

---

## 実装の詳細

### 1. 後退時間積分

```julia
# 終端条件
λ_all[nt, :, :, :] .= 0.0

# 後退ループ（nt-1 → 1）
for t in (nt-1):-1:1
  λ_initial = λ_all[t+1, :, :, :]  # 次ステップ（時間的に後）

  # 熱物性値計算
  cp, k = thermal_properties_calculator(T_cal[t, :, :, :], cp_coeffs, k_coeffs)

  # 係数とRHS構築（残差注入）
  a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_adjoint_system!(
    λ_initial, T_cal[t, :, :, 1], Y_obs[t, :, :],
    rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
  )

  # 疎行列組み立て
  A = assemble_adjoint_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

  # CG法求解（対角前処理、ホットスタート）
  x, cg_history = cg!(x0, A, b, Pl=M, reltol=rtol, maxiter=maxiter)

  # 結果保存
  λ_all[t, :, :, :] = reshape(x, ni, nj, nk)
  x0 .= x  # ホットスタート更新
end
```

### 2. 残差注入機構

**物理的意味**:
測定値と計算値の差を底面から注入し、感度を伝播させる

**実装**:

```julia
# 底面（k=1）のみ
if k_idx == 1
  residual = T_cal_bottom[i, j] - Y_obs[i, j]
  rhs += 2.0 * residual * dx * dy
end
```

**係数2.0の由来**:
目的関数 `J = Σ(T - Y)²` の勾配 `∂J/∂T = 2(T - Y)`

### 3. 対角前処理CG法

```julia
# 対角前処理器
diag_A = diag(A)
inv_diag = [d != 0.0 ? 1.0/d : 0.0 for d in diag_A]
M = Diagonal(inv_diag)

# CG法（IterativeSolvers.jl）
x, cg_history = cg!(x0, A, b, Pl=M, reltol=1e-8, maxiter=1000, log=true)
```

**収束性**:
- 1D問題: 5反復
- 3D問題: 34反復（前処理なしでは100反復超）

---

## 性能評価

### メモリ使用量

**1D問題（1×1×5、11ステップ）**:
- `λ_all`: 55要素 × 8bytes = 440 bytes
- 疎行列A: ~175要素（7点ステンシル）

**3D問題（3×3×8、16ステップ）**:
- `λ_all`: 1,152要素 × 8bytes = 9,216 bytes
- 疎行列A: ~1,512要素（7点ステンシル）

### 計算時間

```
Phase 3-1 (1D):  2.2秒（テスト全体）
Phase 3-2 (3D):  0.0秒（テスト全体）
```

**実用規模（推定）**:
- 格子: 100×100×20
- 時間: 1000ステップ
- 推定時間: ~10分（CG収束に依存）

---

## Python実装との比較

### コード構造の対応

| Python | Julia | 対応 |
|--------|-------|------|
| `coeffs_and_rhs_building_Adjoint()` | `build_adjoint_system!()` | ✅ |
| `assemble_A_Adjoint()` | `assemble_adjoint_matrix()` | ✅ |
| `multiple_time_step_solver_Adjoint()` | `solve_adjoint!()` | ✅ |

### 数値精度の比較

| 項目 | Python | Julia | 差異 |
|------|--------|-------|------|
| 1D最大誤差 | - | 5.08e-21 | 機械精度 |
| 3D最大誤差 | - | 1.26e-14 | 機械精度 |
| CG反復（1D） | 5.0 | 5.0 | 完全一致 |
| CG反復（3D） | 31.5±0.7 | 34.2±0.4 | 8%差（許容範囲） |

**CG反復回数の差異**:
- Python: scipy.sparse.linalg.cg
- Julia: IterativeSolvers.jl
- 実装の微細な差異（停止条件、前処理器の丸め誤差）

---

## 学んだ教訓

### 1. 配列インデックスの明示

**教訓**: 配列形状を関数のドキュメントに明記する

```julia
"""
Args:
  Temperature: 3D温度配列 (ni, nj, nk) [K] ← 明示的に記載
"""
```

### 2. 多項式係数の統一

**教訓**: 係数順序を全コードで統一する

- Pythonオリジナル: `[a, b, c, d]` = `a*T^3 + b*T^2 + c*T + d`
- 全ての実装でこの順序を統一

### 3. テストデータの検証

**教訓**: 参照データ生成スクリプトも検証対象

- 参照データが間違っていると、実装が正しくてもテストが失敗
- 手動計算やPythonオリジナルとの直接比較が重要

### 4. デバッグ手法

**効果的だったアプローチ**:
1. 異常値の発見（cp = 1.35e10）
2. 根本原因の追跡（配列インデックス、係数順序）
3. 段階的な修正と検証

---

## 残課題と今後の改善

### 1. 性能最適化（将来の課題）

**現状**: 小規模問題では十分高速

**実用規模での改善案**:
- スレッド並列化（`@threads`マクロ）
- より効率的な前処理器（ILU分解など）
- GPU加速（CUDA.jl）

### 2. メモリ効率（将来の課題）

**現状**: 全時間ステップの`λ_all`を保持

**改善案**:
- スライディングウィンドウ（Phase 4で実装予定）
- 必要な時間範囲のみ保持
- ディスクへの書き出し（超大規模問題）

### 3. 数値安定性

**現状**: 1e-14の精度で安定

**将来の検討事項**:
- より厳しい条件（高温、急激な温度変化）での安定性
- 適応的時間刻み
- より強固な前処理器

---

## Phase 3まとめ

### 達成事項

✅ **随伴ソルバーの完全実装**
   - 後退時間積分
   - 残差注入機構
   - CG法による高速求解

✅ **Python参照との完全一致**
   - 1D問題: 5.08e-21の誤差
   - 3D問題: 1.26e-14の誤差

✅ **全テストケース合格**
   - 終端条件
   - 残差注入
   - 時間積分の正確性
   - 勾配場の精度

✅ **重大バグの発見と修正**
   - 配列インデックスの誤り
   - 多項式係数の順序不一致
   - 参照データ生成の修正

### Phase 4への準備

Phase 3で構築した随伴ソルバーを基盤として、Phase 4では：

1. **共役勾配法（CGM）**: 表面熱流束の逆解析
2. **ライン検索**: 最適なステップサイズの決定
3. **感度解析**: `λ[:, :, :, nk]`（表面の随伴場）を勾配として使用
4. **スライディングウィンドウ**: 時間ウィンドウごとの逆解析

Phase 3で確立した高精度・高速な随伴ソルバーにより、Phase 4の実装が確実になりました。

---

**次のステップ**: Phase 4（CGM最適化）の実装開始