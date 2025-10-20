# 次セッションTODO

**日時**: 2025年10月21日 03:44
**ブランチ**: sliding-window-validation
**状態**: Julia版スライディングウィンドウ発散問題調査中

## 🎯 現在の状況

### 問題の特定
Julia版の`run_sliding_window_validation.jl`でスライディングウィンドウ計算を実行すると、**ウィンドウ1のt=63から発散**が始まる。

### 重要な発見 ✅

1. **ウィンドウ1単体では正常動作**
   - `test_window1_comparison.py`: Python版完了、全収束 ✅
   - `test_window1_comparison.jl`: Julia版CGM反復1まで完了、全収束 ✅
   - **結論**: 単一ウィンドウの計算ロジックは正しい

2. **スライディングウィンドウスクリプトでの発散**
   - `run_sliding_window_validation.jl`実行時のみ発散
   - 発散箇所: ウィンドウ1、t=63/72から（残差483→20000回反復でも収束せず）
   - Python版は同じ条件で正常動作

3. **Y_obsスライスの形状問題を発見**
   - テストスクリプト作成時、Julia版が72ステップ取得していた（正しくは71）
   - 修正済み: `Y_obs_win = Y_obs_full_julia[:, :, start_idx+1:end_idx+1]`

## 🔥 最も疑わしいバグ: Y_obsスライスのオフバイワンエラー

### 問題箇所

**ファイル**: `julia/scripts/run_sliding_window_validation.jl:281-283`

```julia
max_L = min(window_size, (nt - 1) - start_idx)
end_idx = start_idx + max_L
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]
```

### バグの詳細

**ウィンドウ1の場合**:
- `nt = 300` (Y_obsの全時間ステップ数)
- `window_size = 71`
- `start_idx = 0`
- `max_L = min(71, 299) = 71`
- `end_idx = 0 + 71 = 71`
- `Y_obs_win = Y_obs[:, :, 1:72]` → **72ステップ** ❌

**正しくは**:
- `Y_obs_win = Y_obs[:, :, 1:71]` → **71ステップ** ✅
- または `end_idx = 70` にすべき

### 影響

- Y_obsが72ステップ、q_initが70ステップで**ステップ数不一致**
- CGMソルバー内でインデックスエラーや配列アクセス違反
- 結果として数値的に不安定になり発散

### 修正案

**オプション1**: スライスを明示的に指定
```julia
Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L+1]
```

**オプション2**: end_idxの計算を修正
```julia
end_idx = start_idx + max_L - 1  # 71ではなく70
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]  # 1:71
```

## 🔍 次セッションの最優先タスク

### タスク1: Y_obsスライスの修正と検証 🚨最優先（5分）

1. **修正実施**
   ```julia
   # julia/scripts/run_sliding_window_validation.jl:283を修正
   # 修正前
   Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]

   # 修正後（オプション1推奨）
   Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L+1]
   ```

2. **デバッグ出力追加**
   ```julia
   if verbose
       println("\n=== ウィンドウ $window_id デバッグ ===")
       println("start_idx=$start_idx, end_idx=$end_idx, max_L=$max_L")
       println("Y_obs_win形状: $(size(Y_obs_win))")
       println("期待: (ni=$ni, nj=$nj, nt=$(max_L+1))")
       println("q_init_win形状: $(size(q_init_win))")
       println("期待: (ni=$ni, nj=$nj, nt=$max_L)")
   end
   ```

3. **再実行**
   ```bash
   julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
   ```

### タスク2: 修正確認後、完全比較（10分）

修正が成功したら：

1. **Julia版ウィンドウ1テストを完走させる**
   ```bash
   timeout 900 julia --project=julia test_window1_comparison.jl
   ```

2. **Python版と比較**
   ```bash
   python compare_window1_results.py
   ```

3. **スライディングウィンドウ全体を実行**
   ```bash
   julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
   ```

### タスク3: 結果検証（5分）

- 発散が解消されたか確認
- Python版との結果比較
- 最終レポート作成

## 📁 作成済みファイル

### テストスクリプト
1. `test_window1_comparison.py` - Python版ウィンドウ1単体テスト ✅
2. `test_window1_comparison.jl` - Julia版ウィンドウ1単体テスト ✅
3. `compare_window1_results.py` - 結果比較スクリプト ✅

### 結果ファイル
- `shared/results/window1_python_test.npz` - Python版ウィンドウ1結果 ✅
- `shared/results/julia_window1_test.log` - Julia版実行ログ（途中まで） ✅
- `shared/results/python_sliding_window_cgm3.npz` - Python版全体結果 ✅

## 🔬 技術的詳細

### 発散パターン（修正前のJulia版）
```
[t=62/72] converged: Iteration= 1260 : Res_0= 63.54
[t=63/72] converged: Iteration= 17256 : Res_0= 483.41
[t=64/72] not converged: Iteration=20000 : Res_0= 659.49
```

### 正常パターン（期待される動作）
```
[t=62/71] converged: Iteration= 16 : Res_0= 0.00560
[t=63/71] converged: Iteration= 16 : Res_0= 0.00557
[t=64/71] converged: Iteration= 16 : Res_0= 0.00554
```

## 📊 比較データ

### Python版（正常）
- ウィンドウ1: J = [190008, 178540, 168294]
- q範囲: [-2.12e+04, 2.94e+02] W/m²
- Y_obs_win形状: (71, 80, 100) → 71ステップ ✅
- q_init形状: (70, 80, 100) → 70ステップ ✅

### Julia版（修正前・異常）
- ウィンドウ1 CGM反復0: J = 1.88059e+05
- **Y_obs_win形状: (80, 100, 72)** → 72ステップ ❌
- q_init形状: (80, 100, 70) → 70ステップ ✅
- **ステップ数不一致！**

## 🚀 実行コマンド

### デバッグ実行（修正後）
```bash
# 修正を反映して実行
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3 2>&1 | tee shared/results/julia_fixed.log
```

### ウィンドウ1単体テスト（修正確認用）
```bash
# Julia版（タイムアウト延長）
timeout 900 julia --project=julia test_window1_comparison.jl

# 比較実行
python compare_window1_results.py
```

## 🎓 学び

1. **配列インデックスの落とし穴**
   - `end_idx = start_idx + max_L` は `max_L+1`個の要素を含む
   - Julia: `arr[1:71]` = 71要素
   - 常に `length(start:end)` = `end - start + 1` を確認すること

2. **ステップ数の一貫性**
   - 温度観測: `nt`ステップ (t=0, t=1, ..., t=nt-1)
   - 熱流束: `nt-1`ステップ (時間ステップ間の量)
   - Y_obsとq_initのステップ数は必ず整合性を確認

3. **テスト駆動デバッグの有効性**
   - 単体テストで問題を切り分け
   - 小さなテストケースで動作確認してから全体を実行

## 🔄 次回セッション開始手順

1. **git状態確認**
   ```bash
   cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
   git status
   git log --oneline -5
   ```

2. **このファイルを読む**
   ```bash
   cat TODO_NEXT_SESSION.md
   ```

3. **最優先修正実施**
   - `julia/scripts/run_sliding_window_validation.jl:283`のバグ修正
   - デバッグ出力追加

4. **テスト実行**
   ```bash
   julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
   ```

5. **問題解消確認**
   - 発散が起きないか確認
   - Python版との結果比較

## 📝 Gitコミット準備

修正後、以下をコミット：

```bash
# 修正ファイル
git add julia/scripts/run_sliding_window_validation.jl

# テストスクリプト
git add test_window1_comparison.py
git add test_window1_comparison.jl
git add compare_window1_results.py

# ドキュメント
git add TODO_NEXT_SESSION.md

# コミット
git commit -m "fix: Julia版スライディングウィンドウのY_obsスライス修正

**問題**:
Y_obs_winが72ステップ取得していた（正しくは71ステップ）
→ q_init（70ステップ）とステップ数不一致
→ t=63以降で発散

**修正内容**:
- run_sliding_window_validation.jl:283を修正
- Y_obs_winスライスを明示的に指定
- デバッグ出力を追加

**検証**:
- 単体テスト（test_window1_comparison.jl）で正常動作確認
- Python版との比較で整合性確認

fixes #Julia版発散問題
"
```

---
**作成日時**: 2025-10-21 03:44
**作成者**: Claude Code
**セッション**: sliding-window-validation
**推定修正時間**: 5-10分
