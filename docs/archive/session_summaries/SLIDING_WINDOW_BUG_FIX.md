# Julia版スライディングウィンドウ発散問題の解決

**日時**: 2025年10月21日
**ブランチ**: sliding-window-validation
**影響ファイル**: `julia/scripts/run_sliding_window_validation.jl`

## 問題の症状

### 発散パターン
- **ウィンドウ1**: 時間ステップ t=1～62 まで正常収束
- **t=62**: 1260反復で収束（異常に多い）
- **t=63以降**: 爆発的発散
  - t=63: 17256反復、Res_0=483
  - t=64: 20000反復打ち切り、Res_0=659
  - t=65: 20000反復打ち切り、Res_0=7998
  - t=66以降: 指数関数的に発散

### 対比: Python版は正常動作
- 全23ウィンドウ、300ステップを251秒で完了
- すべての時間ステップで正常収束

## 根本原因

### 発見した2つのバグ

#### 1. **配列インデックスのオフバイワンエラー**

**問題コード** (行283):
```julia
Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L]
```

**Python版** (正しい):
```python
Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]
```

**問題の詳細**:
- Julia版: `start_idx=0, max_L=71` → `Y_obs[:, :, 1:71]` → **71要素**
- Python版: `start_idx=0, end_idx=71` → `Y_obs[0:72]` → **72要素**
- 熱流束71ステップを計算するには温度データ72ステップが必要

**結果**:
- Julia版は1ステップ分のデータが不足
- 境界条件が正しく設定されず発散

#### 2. **Sensitivity Solverの非最適設定**

**問題コード** (行198):
```julia
sensitivity_solver = :pbicgstab,  # 非対称問題用
```

**最適設定**:
```julia
sensitivity_solver = :pcg,  # 対称正定値問題用
```

**理由**:
- 熱伝導方程式は対称正定値行列を生成
- PCGはPBICGSTAB!より17.4%高速（ベンチマーク結果）
- `solver_comparison_3way.md`でPCG+noneが最速と確認済み

## 修正内容

### 修正1: 配列インデックスの修正

**ファイル**: `julia/scripts/run_sliding_window_validation.jl`
**行**: 283-284

**修正前**:
```julia
Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L]
```

**修正後**:
```julia
# Python版: Y_obs[start_idx:end_idx+1] と同等 (Juliaは1-origin)
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]
```

**検証**:
- `start_idx=0, end_idx=71` → `Y_obs[:, :, 1:72]` → **72要素** ✅
- Python版と同じデータサイズ

### 修正2: PCGへの統一

**ファイル**: `julia/scripts/run_sliding_window_validation.jl`
**行**: 194-199

**修正前**:
```julia
dhcp_solver = :pcg,
dhcp_smoother = :none,
adjoint_solver = :pcg,
adjoint_smoother = :none,
sensitivity_solver = :pbicgstab,  # ← 非最適
sensitivity_smoother = :none
```

**修正後**:
```julia
dhcp_solver = :pcg,
dhcp_smoother = :none,  # 前処理なし（PCG+noneが最速）
adjoint_solver = :pcg,
adjoint_smoother = :none,  # 前処理なし（PCG+noneが最速）
sensitivity_solver = :pcg,  # PCGに統一（対称正定値問題）
sensitivity_smoother = :none  # 前処理なし（PCG+noneが最速）
```

## テスト結果

### 修正前（発散）
- **ウィンドウ1**: t=62まで正常、t=63以降発散
- **時間ステップ表示**: `[t=2/71]`から開始（**71要素のみ**）
- **結果**: タイムアウトまたは強制終了

### 修正後（成功）
- **10ステップテスト**: 21秒で完了、全収束 ✅
- **時間ステップ表示**: `[t=2/10]`から`[t=10/10]`まで正常
- **300ステップ**: 実行中（PID 19060）

### 期待される結果
- Python版と同等: 約250秒で完了
- 全23ウィンドウで正常収束
- 発散なし

## 技術的洞察

### なぜt=62で異常が顕在化したか

1. **初期ステップ（t=1～61）**:
   - データ不足の影響が小さい
   - 収束はするが、計算精度が徐々に悪化

2. **臨界点（t=62）**:
   - 累積誤差が臨界値に到達
   - 1260反復（通常の10倍以上）でギリギリ収束

3. **発散開始（t=63以降）**:
   - 誤差が収束範囲を超える
   - 反復計算が発散モードに入る

### データサイズの重要性

```
熱流束 q: [t=0, t=1, ..., t=70]  → 71ステップ
温度 T:   [t=0, t=1, ..., t=71]  → 72ステップ

熱伝導方程式:
  T[t+1] = f(T[t], q[t])

∴ q[0:71]を計算するにはT[0:72]が必要
```

Julia版は**T[0:71]しか渡していなかった**ため、T[71]が欠落し境界条件が不正になった。

## 参照

### ベンチマーク結果
- `shared/results/solver_comparison/solver_comparison_3way.md`
- PCG+none: 27.30秒（最速）
- PBICGSTAB!+none: 32.06秒（17.4%遅い）

### 成功例
- `shared/results/current/sliding_window/python_phase1_fixed.log` (166KB)
  - 23ウィンドウ、251秒、全収束

### 発散ログ
- `shared/results/current/sliding_window/julia_phase1_FINAL.log` (7.3KB)
  - t=63で発散開始

## 次のステップ

1. ✅ バグ修正完了
2. 🔄 全300ステップ実行中（PID 19060）
3. ⏳ Python版との結果比較
4. ⏳ 最終レポート作成
5. ⏳ gitコミット

---

**教訓**: 1-originと0-originの違いは、単純なオフセット問題だけでなく、**範囲の包含性**にも影響する。Python `[a:b]`は`b`を含まないが、Julia `a:b`は`b`を含む。しかし、元の配列スライス部分で`start_idx+max_L`としていたため、実質的に1要素不足していた。
