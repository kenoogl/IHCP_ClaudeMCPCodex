# 次セッション クイックスタートガイド

**作成日**: 2025年10月17日 21:42
**ブランチ**: `tuning7-volume-integral-complete`
**最終コミット**: 4d2ceaf

---

## 🚀 30秒でセッションを開始

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -3
```

**期待される出力**:
```
On branch tuning7-volume-integral-complete
Your branch is up to date with 'origin/tuning7-volume-integral-complete'.

4d2ceaf docs: セッション完了、体積積分形式で28倍高速化達成
c2cff55 fix: 大規模格子での無限ループ問題 - Z方向格子生成の整合性を修正
6010f7f feat: run_10steps_fullsize_test.jlにソルバー・前処理指定機能を追加
```

---

## 📋 前セッションの成果（一言サマリー）

**無限ループ問題を解決し、Phase 1-Eから28倍高速化を達成！**
- 実行時間: 841秒 → **29.6秒**（96.5%改善）
- RMS誤差: 7.744 K → **0.591 K**（92.4%改善）

---

## 🎯 次セッションの最優先タスク

### Task 1: Phase 1-Eとの詳細比較（5分）

```bash
julia --project=julia -e '
  using NPZ
  phase1e = npzread("shared/results/phase1e_adaptive.npz")
  current = npzread("shared/results/julia_10steps_fullsize.npz")

  println("=== Phase 1-E vs 体積積分形式 ===")
  println("実行時間: ", phase1e["elapsed_cgm"], " vs ", current["elapsed_cgm"])
  println("RMS誤差: ", phase1e["rms_error"], " vs ", current["rms_error"])
  println("最大誤差: ", phase1e["max_error"], " vs ", current["max_error"])
'
```

### Task 2: 動作確認テスト（2分）

```bash
# 小規模格子で動作確認
julia --project=julia julia/scripts/test_dhcp_convection.jl 10 5
```

---

## 📝 Claudeへの指示例

```
TODO_NEXT_SESSION.mdを確認してください。

前セッションで体積積分形式の無限ループ問題を解決し、
28倍の高速化を達成しました（841秒→29.6秒）。

Phase 1-Eとの詳細比較を行い、なぜこれほど高速化したのかを
調査してください。まず結果ファイルを読み込んで比較しましょう。
```

---

## 🔗 重要なファイル

- **詳細ガイド**: `TODO_NEXT_SESSION.md` ← 必読！
- **調査レポート**: `INVESTIGATION_RESULTS.md`
- **実行結果**: `shared/results/julia_10steps_fullsize.npz`
- **実行ログ**: `run_10steps_gs_fixed.log`

---

**準備完了！次セッションをお楽しみください** 🎉
