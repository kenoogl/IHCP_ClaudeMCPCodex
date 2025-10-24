# 次セッションで実施するタスク

**作成日**: 2025年10月24日
**最終更新**: 2025年10月24日 19:00 (中規模ウィンドウ計画策定完了)
**ブランチ**: parallelization
**現在の状態**: 中規模ウィンドウ比較計画策定完了、実行準備完了
**最新コミット**: 未コミット（2つの計画書を新規作成）

---

## 📊 このセッションで完了したこと

### ✅ 実施内容

1. ✅ **TODO_NEXT_SESSION.mdの確認**
   - 前セッションの状況確認完了
   - Python-Julia CGM比較最終レポート確認済み

2. ✅ **未実施の性能測定項目の整理**
   - Phase 5.2 Step 3 Phase 2（大ウィンドウ）が未完了と確認
   - 性能改善提案書の未実施項目を整理

3. ✅ **大ウィンドウPython-Julia比較計画の策定**
   - 計画書作成: `docs/plans/large_window_python_julia_comparison_plan.md`
   - **重大な問題発見**: nt=10ではwindow=71に対して不足
   - 修正版: nt=100に変更、実行時間30-45時間に増加

4. ✅ **中規模ウィンドウ比較計画の策定（推奨版）**
   - 計画書作成: `docs/plans/medium_window_comparison_plan.md` ⭐
   - **最適パラメータ**: nt=100, window=60, overlap=15
   - **ウィンドウ数**: 2個（100%カバー）
   - **実行時間**: 12-14時間（現実的！）

---

## 🎯 次セッションで実施すべきタスク

### 優先度A: 中規模ウィンドウ比較の実行（最推奨）

#### 背景

中規模ウィンドウ計画（`docs/plans/medium_window_comparison_plan.md`）は以下の理由で最推奨：

1. **現実的な実行時間**: 12-14時間（1日で完了可能）
2. **適切なウィンドウ構成**: 2個のウィンドウで100%カバー
3. **段階的評価**: 小（5）→中（60）→大（71）のスケーラビリティ評価
4. **Phase 5.2の実質完了**: 小+中ウィンドウで十分な知見を取得

#### 実行計画

**Phase 1: Julia版basesize最適化**（5-6時間）

```bash
# Test 1-1: basesize=500
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 500 \
  2>&1 | tee shared/results/medium_window_basesize500_cgm3.log

# Test 1-2: basesize=1000
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 1000 \
  2>&1 | tee shared/results/medium_window_basesize1000_cgm3.log

# Test 1-3: basesize=2000
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 2000 \
  2>&1 | tee shared/results/medium_window_basesize2000_cgm3.log
```

**Phase 2: Python-Julia比較**（6-8時間）

```bash
# Python版
export NUMBA_NUM_THREADS=4
python3 python/scripts/run_sliding_window.py \
  --nt 100 --cgm-iter 3 --window 60 --overlap 15 \
  --output python_medium_window_cgm3 \
  2>&1 | tee shared/results/python_medium_window_cgm3.log

# Julia版（最適basesize使用）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize <最適値> \
  2>&1 | tee shared/results/julia_medium_window_cgm3.log
```

**期待される成果**:
- 中規模ウィンドウでのbasesize最適値特定
- Python-Julia完全比較（性能・数値安定性）
- 小→中スケーラビリティ評価
- **Phase 5.2の実質完了**

---

### 優先度B: プロジェクト完了宣言（代替案）

Phase 5.2を現状（85%完了）で終了し、プロジェクト完了とする選択肢：

**理由**:
- 小ウィンドウ測定完了済み
- Python-Julia比較最終レポート完成
- 主要な知見（Julia版推奨、Python版発散）は取得済み

**次のアクション**:
- 論文・報告書作成
- 研究成果の公表

---

### 優先度C: 性能改善継続（Phase 1実装）

中規模ウィンドウ比較を省略し、性能改善に進む選択肢：

**提案1: 前処理改善**（最高優先度）
- ILU(0)、マルチグリッド前処理の実装
- 期待効果: CG反復30-50%削減
- 実装期間: 2-3週間

**提案4: 初期推定値改善**
- 前ステップ解利用、線形外挿
- 期待効果: CG反復15-25%削減
- 実装期間: 1週間

---

## 📁 重要なファイルの場所

### 新規作成された計画書（未コミット）

```
docs/plans/large_window_python_julia_comparison_plan.md   # 大ウィンドウ計画（修正版）
docs/plans/medium_window_comparison_plan.md                # 中規模ウィンドウ計画（推奨版）⭐
```

### 既存の重要なレポート

```
docs/reports/python_julia_cgm_comparison_final.md         # Python-Julia比較最終レポート
docs/reports/PHASE5_2_COMPLETION_SUMMARY.md                # Phase 5.2完了報告書（85%）
docs/plans/performance_improvement_proposals.md             # 性能改善提案書
```

### 実行結果ログ（既存）

```
shared/results/julia_cgm_iter3_basesize500.log   # Julia CGM 3回
shared/results/python_cgm_iter3.log               # Python CGM 3回
shared/results/julia_cgm_iter10_basesize500.log  # Julia CGM 10回
shared/results/python_cgm_iter10_rerun.log       # Python CGM 10回
```

---

## 🔧 注意事項

### Git状態

**現在のブランチ**: parallelization
**未コミット変更**: 2ファイル（新規作成）
- `docs/plans/large_window_python_julia_comparison_plan.md`
- `docs/plans/medium_window_comparison_plan.md`

**コミット前の作業**:
1. 新規計画書をgit add
2. 適切なコミットメッセージで commit
3. リモートに push

**推奨コミットメッセージ**:
```
docs: 中規模・大ウィンドウPython-Julia比較計画策定

**新規作成**:
- `large_window_python_julia_comparison_plan.md`: 大ウィンドウ計画（修正版、nt=100）
- `medium_window_comparison_plan.md`: 中規模ウィンドウ計画（推奨版、nt=100, window=60）

**主要な発見**:
- 当初の大ウィンドウ計画（nt=10, window=71）は実行不可能と判明
- nt=100に修正すると実行時間30-45時間に増加
- 代替案として中規模ウィンドウ（window=60）を提案（12-14時間）

**推奨事項**:
- 中規模ウィンドウ計画の採用を最推奨（現実的な実行時間）
- Phase 5.2の実質完了が可能

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>
```

---

## 🚀 次セッション開始時の手順

### 1. 実態確認（最優先）

```bash
# 現在のディレクトリ確認
pwd

# git状態確認
git status
git log --oneline -5

# ブランチ確認
git branch

# 未コミットファイル確認
ls -lht docs/plans/*.md | head -5
```

### 2. 新規計画書の確認

```bash
# 中規模ウィンドウ計画確認（推奨版）
cat docs/plans/medium_window_comparison_plan.md

# 大ウィンドウ計画確認（参考）
cat docs/plans/large_window_python_julia_comparison_plan.md
```

### 3. 方針決定

以下から選択：

**オプション1（最推奨）**: 中規模ウィンドウ比較実行
- Phase 1から開始（Julia basesize最適化、5-6時間）
- Phase 2に進行（Python-Julia比較、6-8時間）
- 合計12-14時間で Phase 5.2実質完了

**オプション2**: プロジェクト完了宣言
- Phase 5.2を85%完了として終了
- 論文・報告書作成に進む

**オプション3**: 性能改善継続
- 前処理改善（提案1）から開始
- 実装期間2-3週間

---

## 📝 重要な技術的知見（再確認）

### Python-Julia比較結果

**実行速度**:
- Python版が3.3~3.6倍高速

**数値安定性**:
- **Julia版が圧倒的に安定**（熱流束が約2400倍小さい）
- Python版は発散傾向（熱流束10^8オーダー）

**最終結論**:
- **Julia版を実用推奨** ✅
- Python版は使用非推奨 ❌

### basesize最適化結果

| 測定対象 | 最適basesize | 実行時間 |
|---------|-------------|---------|
| DHCP単体 (1ステップ) | 1000 | 6.7秒 |
| CGM全体 (10ステップ) | 1000 | 19.5秒 |
| 小ウィンドウ (window=5) | 500 | 34.34秒 |

**重要**: 最適basesizeは問題構造に依存

---

## 📊 プロジェクト全体の完成度

### 完了項目

- ✅ **Phase 1-6**: 505テスト全合格
- ✅ **Phase 5.2 Step 1-3 Phase 1**: basesize最適化完了（85%）
- ✅ **Python-Julia比較**: 小ウィンドウ（CGM 3回/10回）完了
- ✅ **最終レポート**: Python-Julia比較最終レポート完成
- ✅ **Python版CGMバグ修正**: 完了

### 未完了項目

- ⏸️ **Phase 5.2 Step 3 Phase 2**: 大ウィンドウ測定（85%→100%）
- 🎯 **代替案**: 中規模ウィンドウ測定（推奨）
- 📅 **将来**: 性能改善（Phase 1-3実装）

---

## 💡 次セッションでの推奨アクション

### 最優先推奨（★★★★★）

**中規模ウィンドウ比較計画の実行**

1. 新規計画書をgitコミット・プッシュ
2. Phase 1を開始（Julia basesize最適化、5-6時間）
3. Phase 2に進行（Python-Julia比較、6-8時間）
4. Phase 5.2完了レポート作成
5. プロジェクト完了宣言

**期待される成果**:
- Phase 5.2の実質完了（小+中ウィンドウ測定完了）
- 小→中スケーラビリティの完全評価
- basesize最適化の完全検証
- 次フェーズへの円滑な移行

**実行可能性**: 1日（12-14時間）で完了可能 ✅

---

**次セッション開始時は、このファイルを必ず最初に読んでください。**

**推奨事項**: 中規模ウィンドウ比較計画の採用を強く推奨します。
