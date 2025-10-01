# Phase 3 Adjointソルバー バグ修正レポート

**日付**: 2025-09-30
**対象**: Julia Adjoint随伴ソルバー実装
**問題**: 随伴場が10^8倍小さい値を返す致命的バグ
**結果**: 全テスト合格（1e-7以内の精度達成）

---

## 問題の概要

### 初期状態（バグ発生時）

Phase 3テスト実行結果：

```
Julia Adjoint結果:
  λ_all範囲: [-5.05e-13, 1.90e-14]  ← 異常に小さい

Python参照データ:
  λ_ref範囲: [-1.36e-5, 5.14e-7]    ← 期待値

最大誤差: 1.36e-5（基準 < 1e-7を超過）
相対誤差: 0.999...（ほぼ完全不一致）
```

**症状**: Julia実装の随伴場が約10^8倍小さい → 実装に根本的な問題

---

## バグ調査プロセス

### 1. AdjointSolver.jl の詳細確認

- 残差注入の係数（2.0）：正しい
- 底面インデックス（k=1）へのアクセス：正しい
- RHS構築の符号：正しい
- 時間ステップのインデックス：正しい

→ **Adjoint実装自体には問題なし**

### 2. thermal_properties_calculator の検証

手動デバッグで異常値を発見：

```julia
# 期待値
cp_manual = 500.0 J/(kg·K)
k_manual = 15.0 W/(m·K)

# Julia実装の実際の出力
cp[1,1,1] = 1.35e10 J/(kg·K)  ← 異常に大きい！
k[1,1,1] = 4.05e8 W/(m·K)     ← 異常に大きい！
```

→ **熱物性値計算に重大な問題**

### 3. 配列インデックスの調査

**ThermalProperties.jl（修正前）**:

```julia
nk, nj, ni = size(Temperature)  # ← 誤った順序！

for i in 1:ni
  for j in 1:nj
    for k_idx in 1:nk
      T_current = Temperature[k_idx, j, i]  # ← 誤ったアクセス
```

しかし、`Temperature`は`(ni, nj, nk)`形状で渡されているため、
`Temperature[k_idx, j, i]`は間違った要素を読み出していた。

**Python実装（正しい）**:

```python
ni, nj, nk = Temperature.shape  # (ni, nj, nk)
for i in prange(ni):
    for j in range(nj):
        for k_ijk in range(nk):
            T_current = Temperature[i, j, k_ijk]
```

### 4. polyval_numba の係数順序問題

**Phase 3参照データ生成スクリプト（修正前）**:

```python
# generate_phase3_reference.py (49行)
cp[i, j, k_idx] = cp_coeffs[0] + T_val * (cp_coeffs[1] + ...)  # Hornerの方法
```

これは`cp_coeffs = [c0, c1, c2, c3]`で`c0 + c1*T + c2*T^2 + c3*T^3`を計算。

しかし、オリジナルの`polyval_numba`（60行）では：

```python
result += coeffs[i] * x ** (len(coeffs) - i - 1)
```

これは`coeffs = [a, b, c, d]`で`a*x^3 + b*x^2 + c*x + d`を計算。

**係数の順序が逆！**

参照データの係数`[500, 0, 0, 0]`は、Hornerの方法では`500`（定数）だが、
`polyval_numba`では`500*T^3`として解釈される。

---

## バグ修正内容

### 修正1: ThermalProperties.jl の配列インデックス修正

**修正前**:

```julia
function thermal_properties_calculator(
  Temperature::Array{Float64, 3},
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64}
)
  nk, nj, ni = size(Temperature)  # ← 誤り！

  for i in 1:ni
    for j in 1:nj
      for k_idx in 1:nk
        T_current = Temperature[k_idx, j, i]  # ← 誤ったアクセス
```

**修正後**:

```julia
function thermal_properties_calculator(
  Temperature::Array{Float64, 3},
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64}
)
  ni, nj, nk = size(Temperature)  # ← 修正: Pythonと同じ (ni, nj, nk)

  for k_idx in 1:nk
    for j in 1:nj
      for i in 1:ni
        T_current = Temperature[i, j, k_idx]  # ← 正しいアクセス
```

**変更点**:
- 配列形状の認識を`(ni, nj, nk)`に修正（Pythonと同じ）
- ループ順序を変更（メモリアクセスの最適化）
- インデックスアクセスを`Temperature[i, j, k_idx]`に修正

### 修正2: generate_phase3_reference.py の多項式評価修正

**修正前**:

```python
def thermal_properties_calculator_simple(T, cp_coeffs, k_coeffs):
  for i in range(ni):
    for j in range(nj):
      for k_idx in range(nk):
        T_val = T[i, j, k_idx]
        # Hornerの方法（係数順序が逆）
        cp[i, j, k_idx] = cp_coeffs[0] + T_val * (cp_coeffs[1] + ...)
```

**修正後**:

```python
def polyval_numba_simple(coeffs, x):
  """多項式評価（Pythonオリジナル57-61行と同じ）"""
  result = 0.0
  for i in range(len(coeffs)):
    result += coeffs[i] * x ** (len(coeffs) - i - 1)
  return result

def thermal_properties_calculator_simple(T, cp_coeffs, k_coeffs):
  for i in range(ni):
    for j in range(nj):
      for k_idx in range(nk):
        T_val = T[i, j, k_idx]
        # polyval_numbaと同じ多項式評価
        cp[i, j, k_idx] = polyval_numba_simple(cp_coeffs, T_val)
```

### 修正3: 係数順序の修正

**修正前**:

```python
cp_coeffs = np.array([500.0, 0.0, 0.0, 0.0])  # [c0, c1, c2, c3]
k_coeffs = np.array([15.0, 0.0, 0.0, 0.0])
```

**修正後**:

```python
# 係数順序: [a, b, c, d] for a*T^3 + b*T^2 + c*T + d
# 定数の場合は [0, 0, 0, cp_const]
cp_coeffs = np.array([0.0, 0.0, 0.0, 500.0])
k_coeffs = np.array([0.0, 0.0, 0.0, 15.0])
```

**1D問題（定数物性）**:
- 修正前: `[500, 0, 0, 0]` → `500*T^3` = 1.35e10 (T=300K)
- 修正後: `[0, 0, 0, 500]` → `500` = 500.0

**3D問題（温度依存）**:
- 修正前: `[462, 0.134, 0, 0]` → `462*T^3 + 0.134*T^2`
- 修正後: `[0, 0, 0.134, 462]` → `0.134*T + 462`

---

## テスト結果

### Phase 3-1: 1D小規模問題

**修正前**:

```
λ_all範囲: [-5.05e-13, 1.90e-14]
最大誤差: 1.36e-5 > 1e-7  ← 失敗
相対誤差: 0.999             ← 失敗
```

**修正後**:

```
============================================================
Phase 3テスト: Adjoint 1D小規模問題
============================================================
格子: 1×1×5
時間: nt=11, dt=0.1 s

Julia Adjoint求解完了
  平均CG反復回数: 5.0
  λ_all範囲: [-1.3644e-5, 5.1353e-7]  ← Python参照と一致！

誤差統計（全時間ステップ）:
  最大絶対誤差: 5.08e-21  ← 1e-7以内 ✓
  平均絶対誤差: 1.09e-21
  相対誤差: 3.72e-16      ← 1e-5以内 ✓

Test Summary: 全8テスト合格
```

### Phase 3-2: 3D小規模問題（温度依存熱物性値）

```
============================================================
Phase 3テスト: Adjoint 3D小規模問題（温度依存熱物性値）
============================================================
格子: 3×3×8
時間: nt=16, dt=0.001 s

Julia Adjoint求解完了
  平均CG反復回数: 34.2 ± 0.41
  λ_all範囲: [-1.7031e-6, 1.3668e-6]

誤差統計（全時間ステップ）:
  最大絶対誤差: 1.26e-14  ← 1e-7以内 ✓
  平均絶対誤差: 2.57e-15
  相対誤差: 7.41e-9       ← 1e-5以内 ✓

勾配場（表面 k=8）:
  Python参照との誤差: 1.26e-14  ← 1e-7以内 ✓

Test Summary: 全5テスト合格
```

---

## 影響範囲の確認

### 1. Phase 1（熱物性値計算）への影響

Phase 1のテストは**定数係数を使用していない**ため、バグの影響は受けていない。

```julia
# Phase 1では実際の多項式係数を使用
cp_coeffs = [2.01e-10, -3.43e-7, 0.135, 469.85]  # 全ての項が非ゼロ
k_coeffs = [4.80e-12, -8.18e-9, 0.0162, 8.12]
```

このため、Phase 1のテストは既に合格していた。

### 2. Phase 2（DHCPソルバー）への影響

Phase 2でも実際の多項式係数を使用しているため、バグの影響なし。

### 3. 今後の注意事項

**配列形状の一貫性**:
- Pythonとの互換性のため、全ての3D配列は`(ni, nj, nk)`形状を使用
- コメントで配列形状を明記

**多項式係数の順序**:
- 全てのコードで`[a, b, c, d]` = `a*T^3 + b*T^2 + c*T + d`形式を統一
- テストデータ生成時は`polyval_numba`と同じ実装を使用

---

## まとめ

### バグの根本原因

1. **ThermalProperties.jl**: 配列インデックスの順序ミス
   - `(nk, nj, ni)`を仮定していたが、実際は`(ni, nj, nk)`
   - 誤った温度値を読み出し → 誤った熱物性値

2. **generate_phase3_reference.py**: 多項式評価の不一致
   - Hornerの方法（定数項が先）とpolyval_numba（高次項が先）の混在
   - テストデータの係数順序が誤り

### 修正の効果

- **精度**: 最大誤差が1.36e-5 → 5.08e-21（約15桁改善）
- **整合性**: Python参照実装と完全一致（機械精度内）
- **信頼性**: 全テストケース合格（1D、3D両方）

### 学んだ教訓

1. **配列インデックスの明示**: 配列形状をコメントで明記する重要性
2. **多項式評価の統一**: 係数順序を全コードで統一する必要性
3. **テストデータの検証**: 参照データ生成スクリプトも検証対象
4. **デバッグの手法**: 異常値から根本原因を追跡する重要性

---

## 次のステップ

Phase 3が完了したので、次は**Phase 4: CGM最適化**に進みます。

**Phase 4の内容**:
- 共役勾配法（CGM）による表面熱流束の逆解析
- ライン検索アルゴリズム
- 感度解析と勾配計算
- スライディングウィンドウ実装

Phase 3で構築した随伴ソルバーを基盤として、完全な逆解析システムを実装します。