# セッションサマリー: スライディングウィンドウ発散問題調査

**日時**: 2025-10-21
**ブランチ**: sliding-window-validation
**状態**: 🔴 Julia版発散問題未解決

---

## 📋 実施内容

### 1. オフバイワンエラー修正
**ファイル**: `julia/scripts/run_sliding_window_validation.jl:283`

**修正内容**:
```julia
# 修正前（❌ 72ステップ取得）
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]

# 修正後（✅ 71ステップ）
Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L]
```

**検証結果**:
- ✅ ウィンドウ長は71ステップに修正完了
- ❌ しかし発散問題は解決せず

---

## 🔍 発見した問題

### Python版 vs Julia版の比較結果

| 項目 | Python版 | Julia版 | 差異 |
|------|----------|---------|------|
| **完了状況** | ✅ 23ウィンドウ完了 | ❌ ウィンドウ1で停止 | - |
| **総実行時間** | 251秒 | タイムアウト | - |
| **t=2反復数** | 32回 | 96回 | **3.0倍** |
| **t=10反復数** | 17回 | 96回 | **5.6倍** |
| **t=2残差** | 0.169 | 0.365 | **2.2倍** |
| **t=10残差** | 0.013 | 0.162 | **12.9倍** |
| **t=62反復数** | 16回 | 1260回 | **78倍** |
| **t=63反復数** | 16回 | 17256回 | **1078倍** |

### 詳細な初期残差履歴（最初の10ステップ）

#### Python版（正常）
```
[t=1/71] Iteration= 36 : Res_0= 0.3628
[t=2/71] Iteration= 32 : Res_0= 0.1690
[t=3/71] Iteration= 28 : Res_0= 0.0855
[t=4/71] Iteration= 27 : Res_0= 0.0484
[t=5/71] Iteration= 27 : Res_0= 0.0311
[t=6/71] Iteration= 25 : Res_0= 0.0226
[t=7/71] Iteration= 22 : Res_0= 0.0181
[t=8/71] Iteration= 18 : Res_0= 0.0155
[t=9/71] Iteration= 17 : Res_0= 0.0138
[t=10/71] Iteration= 17 : Res_0= 0.0126
```
→ **残差が順調に減少**（0.363 → 0.013、28倍減少）

#### Julia版（発散）
```
[t=2/71] Iteration= 96 : Res_0= 0.3654
[t=3/71] Iteration= 99 : Res_0= 0.0688  ← 異常に小さい
[t=4/71] Iteration= 96 : Res_0= 0.1661
[t=5/71] Iteration= 97 : Res_0= 0.0936
[t=6/71] Iteration= 96 : Res_0= 0.1853
[t=7/71] Iteration= 96 : Res_0= 0.1884
[t=8/71] Iteration= 97 : Res_0= 0.1329
[t=9/71] Iteration= 97 : Res_0= 0.1722
[t=10/71] Iteration= 96 : Res_0= 0.1624
```
→ **残差が高止まり**（0.365 → 0.162、2.2倍減少のみ）

---

## 🎯 根本原因の仮説

### ❌ 配列スライスの問題ではない
- ウィンドウ長71ステップは正しく設定済み
- 問題は別の箇所にある

### ✅ 最も可能性が高い原因

**Julia版は最初から係数行列Aまたは右辺ベクトルbが間違っている**

#### 根拠:
1. **最初のステップから反復数が3-5倍多い**
   - 係数行列の条件数が悪い可能性

2. **残差が高止まり**
   - 右辺ベクトルbまたは初期値が不適切

3. **t=3の異常パターン**
   - Julia版のt=3だけ残差0.0688（他は0.09-0.19）
   - 計算エラーまたはインデックスミスの可能性

---

## 🔧 次セッションで確認すべき項目

### 優先度1: 時間ステップのインデックス確認

**ファイル**:
- `julia/src/solvers/DHCPSolver.jl`
- `python/validation/run_sliding_window_validation.py`

**確認内容**:
1. **時間ループの開始インデックス**
   - Python: `for t in range(1, nt)` → 1から開始
   - Julia: `for t in 2:nt` → 2から開始
   - 境界条件適用のタイミング

2. **初期温度場の扱い**
   - T_old, T_newの初期化
   - t=1のスキップ有無

### 優先度2: 熱流束境界条件の適用方法

**確認内容**:
1. **Z-plus境界での熱流束設定**
   - Python版とJulia版の実装差異
   - 配列インデックスの1-based/0-based変換

2. **CalcAX!関数内の境界条件**
   - マトリクスフリー実装での境界適用タイミング

### 優先度3: 初期温度場とq_initの設定

**確認内容**:
1. **T_init形状と値の検証**
   - Python: `(ni, nj, nk)`
   - Julia: `(ni, nj, nk)` (同じ)
   - 値の一致確認

2. **q_init形状と値の検証**
   - Python: `(nt-1, ni, nj)` = `(70, 80, 100)`
   - Julia: `(ni, nj, nt-1)` = `(80, 100, 70)` ← **次元順序が逆？**

---

## 📂 重要ファイル

### ログファイル
- `shared/results/python_phase1_fixed.log` (2801行) - Python版完全成功
- `shared/results/julia_phase1_FINAL.log` (112行) - Julia版発散で停止

### 修正済みファイル
- `julia/scripts/run_sliding_window_validation.jl:283` - Y_obsスライス修正済み

### 使用中の反復法
- **DHCP & Adjoint**: PCG（Preconditioned Conjugate Gradient、前処理なし）
- **Sensitivity**: PBiCGSTAB（前処理なし）

---

## 🚀 次セッションのクイックスタート

### 1. 状態確認
```bash
git status
git log --oneline -5
```

### 2. 時間ステップインデックス比較
```bash
# Python版の時間ループ確認
grep -A 10 "for t in range" python/validation/run_sliding_window_validation.py

# Julia版の時間ループ確認
grep -A 10 "for t in" julia/src/solvers/DHCPSolver.jl
```

### 3. q_init配列形状確認
```bash
# Python版
grep -B 5 -A 5 "q_init.*shape\|q_init.*zeros" python/validation/run_sliding_window_validation.py

# Julia版
grep -B 5 -A 5 "q_init.*zeros" julia/scripts/run_sliding_window_validation.jl
```

---

## 📊 補足情報

### Python版ウィンドウ別サマリー（参考）
```
ウィンドウ1: [0,71]    DHCP 3717回  Adjoint 14227回  35.03秒
ウィンドウ2: [54,125]  DHCP 3410回  Adjoint 14055回  32.96秒
...
全23ウィンドウ完了、総251秒
```

### Julia版（発散状況）
```
ウィンドウ1のみ実行:
t=1-61:  正常（反復96-179回）
t=62:    発散開始（1260回）
t=63:    完全発散（17256回、137秒）
→ タイムアウト
```

---

**最終更新**: 2025-10-21 08:00
**次セッション担当者へ**: このファイルを最初に読んでから作業開始してください
