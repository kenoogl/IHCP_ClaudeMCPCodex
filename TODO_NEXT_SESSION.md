# 次セッション作業ガイド

**作成日時**: 2025年10月21日 14:30
**ブランチ**: sliding-window-validation
**最新コミット**: `d292900` - "fix: Julia版スライディングウィンドウの配列インデックスバグ修正"

## ✅ 完了した作業（このセッション）

### 🐛 Julia版スライディングウィンドウ発散問題の解決

**問題症状**:
- t=1～62: 正常収束
- t=62: 異常（1260反復）
- t=63以降: 爆発的発散（17256反復→20000打ち切り）

**根本原因**:
1. **配列インデックスのオフバイワンエラー**
   - `Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L]`
   - 71要素のみ（本来72必要）
   - 熱流束71ステップ計算には温度データ72ステップ必要

2. **Sensitivity Solver非最適設定**
   - PBICGSTAB!使用（非対称問題用）
   - PCGが17.4%高速（対称正定値問題）

**修正内容**:
```julia
# 修正1: 配列スライス (行283-284)
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]  # 72要素

# 修正2: PCGに統一 (行198)
sensitivity_solver = :pcg  # PCGに統一
```

**テスト結果**:
- ✅ 10ステップテスト: 21秒で完了、全収束
- ⚠️ 300ステップテスト: プロセス異常（9分以上データ読み込み中で停止）

**詳細レポート**: `SLIDING_WINDOW_BUG_FIX.md`

### 📊 Python版成功データ

**実行完了**: `shared/results/current/sliding_window/python_phase1_fixed.log` (166KB)
- 総ウィンドウ数: 23
- 総実行時間: 251.38秒
- 全ウィンドウで正常収束 ✅
- q範囲: [-5.611e+04, 2.843e+02] W/m²

## 🚧 次セッションで実施すべきタスク

### 1. Julia版実行異常の調査（優先度: 最高）

**現象**:
- PID 19060が9分以上データ読み込み中で停止
- ログ: "データ読み込み中..." 以降進まず

**調査項目**:
```bash
# プロセス確認
ps aux | grep julia | grep sliding

# ログ確認
tail -100 shared/results/julia_phase1_FIXED.log

# 小規模テスト再実行
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3 --nt 72 --window 71
```

**可能性**:
1. NPZ読み込みの問題（大容量ファイル1.1GB）
2. メモリ不足
3. Julia環境の問題
4. 別のバグ

### 2. 修正版の動作確認

**ステップバイステップ**:

#### Step 1: 最小テスト（10ステップ、overlap無し）
```bash
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --nt 10 --window 10 --overlap 0
```
- 期待: 20秒程度で完了
- 確認: 全収束、発散なし

#### Step 2: ウィンドウ1テスト（72ステップ）
```bash
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 3 --nt 72 --window 71 --overlap 0
```
- 期待: 60秒程度で完了
- 確認: `[t=1/71]`から`[t=71/71]`まで正常

#### Step 3: フルテスト（300ステップ）
```bash
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
```
- 期待: 250秒程度で完了（Python版と同等）
- 確認: 23ウィンドウ、全収束

### 3. Python版との比較分析

**比較項目**:
- 実行時間
- ウィンドウ数
- 各ウィンドウの収束状況
- 熱流束範囲
- CGM反復回数

**Python版データ**: `shared/results/current/sliding_window/python_phase1_fixed.log`

### 4. NPZ保存エラーの修正（低優先度）

**エラーメッセージ**:
```
❌ Error: ArgumentError("cannot reinterpret `Any` as `UInt8`, type `Any` is not a bits type")
```

**原因**: `windows_info`のDict型がNPZ保存と互換性なし

**修正方針**:
- Dictを保存可能な形式に変換
- または別ファイル（JSON）に保存

## 📁 重要なファイル

### 新規作成（このセッション）
- `SLIDING_WINDOW_BUG_FIX.md` - バグ修正詳細レポート

### 修正済み
- `julia/scripts/run_sliding_window_validation.jl` - 配列スライス修正、PCG統一

### 参照ログ
- `shared/results/current/sliding_window/python_phase1_fixed.log` (166KB) - Python成功例
- `shared/results/current/sliding_window/julia_phase1_FINAL.log` (7.3KB) - Julia発散例
- `shared/results/julia_phase1_FIXED.log` - 修正版（異常終了）

### ベンチマーク
- `shared/results/solver_comparison/solver_comparison_3way.md` - PCG+noneが最速

## 🔬 技術メモ

### 配列インデックスの正しい理解

**Python** (0-origin):
```python
# nt=300, window_size=71
start_idx = 0
end_idx = start_idx + max_L = 71
Y_obs_win = Y_obs[0:72, :, :]  # [start:end+1] → 72要素
```

**Julia** (1-origin):
```julia
# nt=300, window_size=71
start_idx = 0
end_idx = start_idx + max_L = 71
Y_obs_win = Y_obs[:, :, 1:72]  # [start+1:end+1] → 72要素
```

**なぜ72要素必要か**:
```
熱流束 q: [t=0, t=1, ..., t=70]  → 71ステップ
温度 T:   [t=0, t=1, ..., t=71]  → 72ステップ

熱伝導方程式: T[t+1] = f(T[t], q[t])
∴ q[0:70]を計算するにはT[0:71]が必要
```

### PCG vs PBICGSTAB!

**ベンチマーク結果** (10ステップ、80×100×20格子):
- PCG+none: 27.30秒 ⭐ **最速**
- PBICGSTAB!+none: 32.06秒（17.4%遅い）
- PBICGSTAB!+GS: 33.62秒（23.2%遅い）

**理由**:
- 熱伝導方程式 → 対称正定値行列
- PCG → 対称正定値問題専用（最適）
- PBICGSTAB! → 非対称問題用（過剰）

## 🐛 既知の問題

1. **Julia版300ステップ実行が異常に遅い**
   - データ読み込みで9分以上停止
   - 要調査・修正

2. **NPZ保存エラー**
   - windows_info (Dict型)が保存不可
   - 低優先度（計算自体は成功）

## 🔄 推奨作業順序

1. ✅ Julia版実行異常の原因調査
2. ✅ 最小テスト（10ステップ）実行
3. ✅ ウィンドウ1テスト（72ステップ）実行
4. ✅ 成功したらフルテスト（300ステップ）実行
5. ✅ Python版との詳細比較
6. ✅ 比較レポート作成
7. ✅ NPZ保存エラー修正（時間があれば）
8. ✅ gitコミット・プッシュ

---

**次セッション開始時のコマンド**:

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5
cat TODO_NEXT_SESSION.md
cat SLIDING_WINDOW_BUG_FIX.md

# 実行異常の調査
ps aux | grep julia | grep sliding
tail -100 shared/results/julia_phase1_FIXED.log

# 最小テスト実行
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --nt 10 --window 10 --overlap 0 2>&1 | head -200
```

**期待される次のマイルストーン**:
- [ ] Julia版300ステップ正常完了
- [ ] Python版と同等の結果（±5%以内）
- [ ] 発散なし、全ウィンドウ収束
- [ ] 実行時間250秒前後
