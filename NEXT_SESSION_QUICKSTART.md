# 次セッション クイックスタート

## 📌 現在地

**ブランチ**: `sliding-window-validation`
**最新コミット**: `f0e7fd4` - 次セッション用TODO作成
**日時**: 2025年10月21日

---

## ⚡ 即座に実行すべきコマンド

### 1. 環境確認（5秒）
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -3
```

### 2. 問題の核心を確認（10秒）
```bash
# Python版の最初の10ステップ（正常）
echo "=== Python版 ==="
grep "^\[t=" shared/results/python_phase1_fixed.log | head -10

# Julia版は発散のためログなし
echo ""
echo "=== Julia版: 発散のため再実行が必要 ==="
```

### 3. 小規模テスト実行（30秒）
```bash
# 小さなウィンドウで問題を再現
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --window 10 --nt 20 \
  > shared/results/julia_debug_10steps.log 2>&1 &

# ログ監視
tail -f shared/results/julia_debug_10steps.log
# Ctrl+Cで停止
```

---

## 🎯 最優先タスク（30分以内）

### タスク1: ウィンドウ1の最初の比較（10分）

**目的**: PythonとJuliaで最初のステップから違いがあるかを確認

**実行**:
```bash
# Python版: ウィンドウ1の最初の10ステップ
grep "^\[t=" shared/results/python_phase1_fixed.log | \
  awk '{print $1, "iters="$4, "Res_0="$7}' | head -10

# Julia版: 小規模テスト実行後に同じ形式で出力
grep "^\[t=" shared/results/julia_debug_10steps.log | \
  awk '{print $1, "Res_0="$7}' | head -10
```

**期待**: 最初のステップで既に残差が違うなら、初期化の問題。同じなら途中からの問題。

### タスク2: インデックス確認（10分）

**疑いのある箇所**:
```julia
# julia/scripts/run_sliding_window_validation.jl:283
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]
```

**確認**:
1. `start_idx`の初期値は0（Julia: 0-indexed風に使用）
2. しかし配列スライスは1-indexedなので `start_idx+1` が必要
3. Python版との対応を確認

### タスク3: 配列次元の確認（10分）

**Python版**:
- Y_obs形状: `(nt, ni, nj)`
- スライス: `Y_obs[start:end+1, :, :]`

**Julia版**:
- Y_obs形状: `(ni, nj, nt)`
- スライス: `Y_obs[:, :, start+1:end+1]`

**確認**: 次元の順序が正しく変換されているか

---

## 🔍 問題の詳細

### 観察された事実

1. **Python版は完全に正常**
   - 反復数: 16-17回/ステップ
   - 残差: 3.6e-1 → 5.3e-3
   - 実行時間: 251秒

2. **Julia版は発散**
   - t=61: Res_0=51.93（100倍大きい！）
   - t=64: Res_0=659.49（発散開始）
   - t=71: Res_0=1.88e21（完全発散）

3. **シングルウィンドウテストは成功**
   - Julia版の基本実装は正しい
   - 問題はスライディングウィンドウ固有

### 推測される原因（優先度順）

1. **配列インデックスのズレ（90%）**
   - PythonとJuliaで配列次元が逆
   - スライディングウィンドウのインデックス計算が間違っている

2. **ウィンドウ間継承の問題（8%）**
   - 前のウィンドウの温度場が正しく継承されていない
   - しかし、ウィンドウ1でも問題があるなら関係ない

3. **数値精度の問題（2%）**
   - 可能性は低い（シングルウィンドウは成功）

---

## 📂 重要なファイル

### 読むべきファイル
1. `TODO_NEXT_SESSION.md` - 詳細な状況説明
2. `julia/scripts/run_sliding_window_validation.jl:273-320` - スライディングウィンドウループ
3. `python/validation/run_sliding_window_validation.py:102-211` - Python版の同じ部分

### 結果ファイル
- `shared/results/python_sliding_window_cgm3.npz` - Python成功結果（19MB）
- `shared/results/python_phase1_fixed.log` - Pythonログ（2801行）

---

## 🚨 デバッグのヒント

### 最初に確認すべきこと

1. **Y_obsのスライス**
   ```julia
   # Julia版: これが正しいか？
   Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]

   # start_idx=0, end_idx=70の場合
   # → Y_obs[:, :, 1:71] (71ステップ)
   ```

2. **ウィンドウの範囲**
   ```julia
   # ウィンドウ1: start_idx=0, max_L=71, end_idx=71
   # → Y_obsは1:72が必要？（nt+1個）
   ```

3. **温度場の時間次元**
   - T_cal: `(ni, nj, nk, nt)` - 71ステップなら72個の時刻（0を含む）
   - Y_obs: `(ni, nj, nt)` - 71ステップなら72個の時刻

### 簡単なテスト

```julia
# julia REPLで確認
using NPZ
Y_obs = npzread("shared/data/T_measure_700um_1ms.npy")
println(size(Y_obs))  # (18143, 80, 100)

# ウィンドウ1のスライス
start_idx = 0
end_idx = 70
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]
println(size(Y_obs_win))  # (80, 100, 71) のはず
```

---

## 📝 次回への引き継ぎ

**もし今回のセッションで解決できなかった場合**:

1. 小規模テストのログを保存
2. Python版とJulia版の最初の10ステップを比較した結果をメモ
3. 疑わしいコード箇所を特定
4. 次回は特定箇所のデバッグから開始

**もし解決した場合**:

1. 修正内容をコミット
2. Julia版Phase 1を完全実行
3. Python版との結果比較
4. Phase 2（CGM 20000回）へ進む

---

**作成日時**: 2025年10月21日
**所要時間**: このセッションで約2時間作業
**次回予想時間**: 30分-1時間（問題特定と修正）
