# 次回セッション開始手順

**作業**: Z方向格子変数名変更（第5弾・最終）

## 📊 現在の進捗

**完了**: 21/29ファイル（72.4%）
- ✅ 第1弾（3ab7b12）: ソースコード4ファイル
- ✅ 第2弾（841e4be）: ソースコード6ファイル
- ✅ 第3弾（2c96f9f）: テストコード3ファイル
- ✅ 第4弾（6303b27）: スクリプト8ファイル

**残り**: 8ファイル（27.6%）- ベンチマークのみ

## 🎯 次回の作業内容

### 第5弾（最終）: ベンチマーク8ファイル

以下のファイルで `Z` → `z_centers`、`ΔZ` → `dz` に変更：

1. `julia/benchmarks/tuning4_allocation_reduction.jl`
2. `julia/benchmarks/benchmark_phase1b.jl`
3. `julia/benchmarks/benchmark_adjoint_scale_01.jl`
4. `julia/benchmarks/benchmark_adjoint_scale_02.jl`
5. `julia/benchmarks/benchmark_combined_improvement.jl`
6. `julia/benchmarks/benchmark_phase1e_baseline.jl`
7. `julia/benchmarks/benchmark_phase1e_adaptive.jl`
8. `julia/benchmarks/benchmark_phase1e_tuned.jl`

## 📝 開始コマンド

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5
cat VARIABLE_RENAME_PROGRESS.md
```

## 🔄 作業の依頼方法

Claudeに以下のように依頼：

> 「NEXT_SESSION.mdを読んで、第5弾（ベンチマーク8ファイル）の変更作業を開始してください」

## 📌 重要事項

- **変更パターン**: `Z, ΔZ` → `z_centers, dz_grid`
- **変更方法**: 手動で慎重に（置換ミスを防ぐため）
- **確認**: 変更後にgrep検索で漏れチェック

---

**作成日時**: 2025-10-16 15:40
**ブランチ**: tuning7
**最新コミット**: abc846f
