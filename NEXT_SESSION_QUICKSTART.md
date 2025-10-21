# 次セッション クイックスタートガイド

## ⚡ 即座に実行するコマンド

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# 1. 実態確認（最優先）
ps aux | grep -E "(python|julia)" | grep -v grep

# 2. git状態確認
git status
git log --oneline -5

# 3. 最新セッションのサマリー確認
ls -t SESSION_SUMMARY*.md | head -1 | xargs cat
```

## 🎯 最新セッションの成果（2025-10-21）

✅ **バグのあるスクリプト削除完了**:
   - `julia/scripts/test_single_window.jl` (UndefVarError)
   - `test_window1_comparison.jl` (API不一致)
   - `julia/scripts/run_sliding_window_validation.jl` (重複)

✅ **関連ログファイル6個削除**

✅ **現在の動作確認済みスクリプト**:
   - `julia/scripts/run_sliding_window.jl` ⭐ メイン実行スクリプト
   - `julia/scripts/parse_sliding_window_log.py` ⭐ 統計抽出ツール

## 📋 次セッションのタスク

### 🔍 現状把握（必須）

```bash
# 削除されたファイルの確認
git status --short

# 残っているスクリプト確認
ls -1 julia/scripts/*.jl

# 未追跡ファイルの確認
git status | grep "^??"
```

### Phase 1: クリーンアップとコミット

```bash
# 削除をステージング
git add -u

# 未追跡の新規ファイルを確認してから追加
git add julia/scripts/run_sliding_window.jl
git add julia/scripts/parse_sliding_window_log.py
git add python/scripts/

# コミット
git commit -m "chore: バグのあるスクリプト削除、動作確認済みスクリプトのみ保持

削除:
- julia/scripts/test_single_window.jl (UndefVarError)
- test_window1_comparison.jl (API不一致)
- julia/scripts/run_sliding_window_validation.jl (run_sliding_window.jlで代替)
- 関連ログファイル6個

追加:
- julia/scripts/run_sliding_window.jl (動作確認済み)
- julia/scripts/parse_sliding_window_log.py (統計抽出ツール)
- python/scripts/ (Python版補助スクリプト)

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

### Phase 2: 小規模テスト（所要時間: 5分）

```bash
# PCG + none で10ステップ実行
julia --project=julia julia/scripts/run_sliding_window.jl \
  --solver pcg --precond none --cgm-iter 3 --nt 10 \
  2>&1 | tee shared/results/julia_sw_10steps_pcg_none_CLEAN.log

# 統計抽出
python julia/scripts/parse_sliding_window_log.py \
  shared/results/julia_sw_10steps_pcg_none_CLEAN.log
```

### Phase 3: フルスケールテスト（所要時間: 5分）

```bash
# PCG + none で300ステップ実行（必要に応じて）
julia --project=julia julia/scripts/run_sliding_window.jl \
  --solver pcg --precond none --cgm-iter 3 --nt 300 \
  2>&1 | tee shared/results/julia_sw_300steps_pcg_none_CLEAN.log

# 統計抽出
python julia/scripts/parse_sliding_window_log.py \
  shared/results/julia_sw_300steps_pcg_none_CLEAN.log
```

## 📊 期待される結果

| 項目 | Python版 | Julia版（PCG+none） |
|------|----------|---------------------|
| 実行時間 | ~250秒 | ~200秒（予想） |
| ウィンドウ数 | 4-5 | 4-5 |
| DHCP反復数 | 16-36回/step | 同等 |

## 🚨 重要な注意事項

1. **system-reminderは無視する**
   - 古い履歴情報（実際には動いていない）
   - 必ず `ps aux` で実態確認

2. **実行前にプロセス確認**
   ```bash
   ps aux | grep -E "(python|julia)" | grep -v grep
   ```

3. **削除されたスクリプトを参照しない**
   - ❌ `julia/scripts/run_sliding_window_validation.jl`
   - ❌ `julia/scripts/test_single_window.jl`
   - ❌ `test_window1_comparison.jl`
   - ✅ `julia/scripts/run_sliding_window.jl` のみ使用

4. **ログファイルは必ず tee で保存**
   - 統計抽出に必要
   - 後から分析可能

## 📁 重要ファイル

**動作確認済みスクリプト**:
- `julia/scripts/run_sliding_window.jl` - メイン実行スクリプト ⭐
- `julia/scripts/parse_sliding_window_log.py` - 統計抽出ツール ⭐

**セッションサマリー**:
- `SESSION_SUMMARY_2025-10-21_BUGFIX.md` - 最新（バグ修正）
- `SESSION_SUMMARY_2025-10-21_CODEX_SUCCESS.md` - スクリプト作成成功
- その他のサマリー複数

**分析レポート**:
- `PYTHON_JULIA_10STEPS_CGM3_COMPARISON.md` - Python-Julia比較
- `SOLVER_PRECOND_COMPARISON_REPORT_10steps.md` - ソルバー・前処理比較

## ✅ 実行チェックリスト

- [ ] 実態確認（ps aux）
- [ ] git status確認（削除ファイル3個、新規ファイル多数）
- [ ] Phase 1: クリーンアップとコミット
- [ ] Phase 2: 10ステップテスト実行
- [ ] Phase 2: 統計抽出
- [ ] Phase 3: 300ステップテスト実行（任意）
- [ ] Phase 3: 統計抽出
- [ ] 結果分析・ドキュメント更新

## 🗂️ 現在のgit状態

```
M  NEXT_SESSION_QUICKSTART.md
D  julia/scripts/run_sliding_window_validation.jl
D  julia/scripts/test_single_window.jl
D  test_window1_comparison.jl
?? PYTHON_JULIA_10STEPS_CGM3_COMPARISON.md
?? SESSION_SUMMARY_2025-10-21_*.md (複数)
?? SOLVER_PRECOND_COMPARISON_REPORT_10steps.md
?? julia/scripts/parse_sliding_window_log.py
?? julia/scripts/run_sliding_window.jl
?? python/scripts/
```

---

**準備完了！** 次セッションではgitクリーンアップ後、クリーンな環境でテスト実行を行います。
