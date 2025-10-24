# 次セッションで実施するタスク

**作成日**: 2025年10月24日
**最終更新**: 2025年10月24日 20:30 (中規模ウィンドウ計画v2.0修正完了)
**ブランチ**: parallelization
**現在の状態**: 中規模ウィンドウ比較計画v2.0完成、パラメータ検証済み、実行準備完了
**最新コミット**: 9d74586 (中規模ウィンドウ計画書v2.0修正、コミット・プッシュ済み)

---

## 📊 このセッションで完了したこと

### ✅ 実施内容

1. ✅ **TODO_NEXT_SESSION.mdの確認**
   - 前セッションの状況確認完了
   - 中規模ウィンドウ計画（v1.0）の確認

2. ✅ **スライディングウィンドウアルゴリズムの検証**
   - **重大な問題発見**: window=60, overlap=15では**17個のウィンドウ**が生成される
   - 正しいパラメータの特定: window=35, overlap=10で**14個のウィンドウ**（主要5個）
   - カバー率100%確認済み

3. ✅ **中規模ウィンドウ比較計画書の修正（v1.0 → v2.0）**
   - ファイル: `docs/plans/medium_window_comparison_plan.md`
   - **パラメータ修正**: window=60→35, overlap=15→10
   - **ウィンドウ構成修正**: 2個→14個（主要5個 + 短いもの9個）
   - 全コマンド例を更新

4. ✅ **Git操作完了**
   - コミット: 9d74586 「docs: 中規模ウィンドウ計画書修正、パラメータ検証によりv2.0へ更新」
   - プッシュ: parallelizationブランチに正常プッシュ

### 🔍 発見した重要な知見

**スライディングウィンドウアルゴリズムの仕様**:
- 最後のウィンドウが`max_L < overlap`になると、`step = max(1, max_L - overlap) = 1`となる
- この場合、1ステップずつ進む短いウィンドウが大量に生成される
- これはアルゴリズムの仕様であり、バグではない

**最適パラメータ（nt=100）**:
- window=35, overlap=10: 主要5個（適切） ✅
- window=60, overlap=15: 主要4個 + 短い13個（不適切） ❌

---

## 🎯 次セッションで実施すべきタスク

### 優先度A: 中規模ウィンドウ比較の実行（最推奨★★★★★）

#### 背景

中規模ウィンドウ計画v2.0（`docs/plans/medium_window_comparison_plan.md`）は以下の理由で最推奨：

1. **現実的な実行時間**: 12-14時間（1日で完了可能）
2. **適切なウィンドウ構成**: 主要5個のウィンドウで100%カバー（短いウィンドウ9個は自動生成）
3. **段階的評価**: 小（5）→中（35）のスケーラビリティ評価
4. **Phase 5.2の実質完了**: 小+中ウィンドウで十分な知見を取得
5. **パラメータ検証済み**: アルゴリズム検証により正確性を確認

#### 実行計画

**Phase 1: Julia版basesize最適化**（5-6時間）

```bash
# Test 1-1: basesize=500（小ウィンドウ最適値）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 35 --overlap 10 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 500 \
  2>&1 | tee shared/results/medium_window_basesize500_cgm3.log

# Test 1-2: basesize=1000（Step 2最適値）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 35 --overlap 10 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 1000 \
  2>&1 | tee shared/results/medium_window_basesize1000_cgm3.log

# Test 1-3: basesize=2000（中規模用推定値）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 35 --overlap 10 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 2000 \
  2>&1 | tee shared/results/medium_window_basesize2000_cgm3.log
```

**Phase 2: Python-Julia比較**（6-8時間）

```bash
# Python版 CGM 3回
export NUMBA_NUM_THREADS=4
python3 python/scripts/run_sliding_window.py \
  --nt 100 --cgm-iter 3 --window 35 --overlap 10 \
  --output python_medium_window_cgm3 \
  2>&1 | tee shared/results/python_medium_window_cgm3.log

# Julia版 CGM 3回（最適basesize使用）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 35 --overlap 10 --cgm-iter 3 \
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
- 中規模ウィンドウは12-14時間の実行時間が必要

**次のアクション**:
- 論文・報告書作成
- 研究成果の公表
- 実用実装への移行

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

### 修正済み計画書（コミット済み）

```
docs/plans/medium_window_comparison_plan.md (v2.0) ⭐ 最新版
  - window=35, overlap=10
  - 14個のウィンドウ（主要5個 + 短いもの9個）
  - アルゴリズム検証済み、100%カバー確認済み
```

### 既存の計画書（参考）

```
docs/plans/large_window_python_julia_comparison_plan.md
  - window=71, overlap=17（大ウィンドウ）
  - 実行時間: 30-45時間（非推奨）
```

### 既存の重要なレポート

```
docs/reports/python_julia_cgm_comparison_final.md         # Python-Julia比較最終レポート
docs/reports/PHASE5_2_COMPLETION_SUMMARY.md                # Phase 5.2完了報告書（85%）
docs/plans/performance_improvement_proposals.md             # 性能改善提案書
```

### 実行結果ログ（既存）

```
shared/results/julia_cgm_iter3_basesize500.log   # Julia CGM 3回（小ウィンドウ）
shared/results/python_cgm_iter3.log               # Python CGM 3回（小ウィンドウ）
shared/results/julia_cgm_iter10_basesize500.log  # Julia CGM 10回
shared/results/python_cgm_iter10_rerun.log       # Python CGM 10回
```

---

## 🔧 注意事項

### Git状態

**現在のブランチ**: parallelization
**最新コミット**: 9d74586
**リモート状態**: 最新（プッシュ済み）
**作業ツリー**: クリーン（未コミット変更なし）

**最新コミット内容**:
```
9d74586 docs: 中規模ウィンドウ計画書修正、パラメータ検証によりv2.0へ更新
```

### 中規模ウィンドウパラメータ（v2.0）

**正しいパラメータ**:
```
nt = 100
window = 35
overlap = 10
```

**ウィンドウ構成**:
- 総ウィンドウ数: 14個
- 主要ウィンドウ: 5個（長さ10～35ステップ）
  - Window 1: [0, 35] (35ステップ)
  - Window 2: [25, 60] (35ステップ、10ステップオーバーラップ)
  - Window 3: [50, 85] (35ステップ、10ステップオーバーラップ)
  - Window 4: [75, 99] (24ステップ、10ステップオーバーラップ)
  - Window 5: [89, 99] (10ステップ、10ステップオーバーラップ)
- 短いウィンドウ: 9個（長さ1～9ステップ、アルゴリズム仕様による自動生成）
- カバー率: 100%確認済み

**誤ったパラメータ（v1.0、修正済み）**:
```
nt = 100
window = 60  ← 17個のウィンドウを生成（過剰）
overlap = 15
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

# 最新の計画書確認
ls -lht docs/plans/*.md | head -5
```

### 2. 中規模ウィンドウ計画書v2.0の確認

```bash
# 中規模ウィンドウ計画確認（v2.0、推奨版）
cat docs/plans/medium_window_comparison_plan.md

# パラメータ確認（重要！）
grep -A5 "window_size = " docs/plans/medium_window_comparison_plan.md
```

### 3. 方針決定

以下から選択：

**オプション1（最推奨★★★★★）**: 中規模ウィンドウ比較実行
- Phase 1から開始（Julia basesize最適化、5-6時間）
- Phase 2に進行（Python-Julia比較、6-8時間）
- 合計12-14時間で Phase 5.2実質完了
- **パラメータ**: nt=100, window=35, overlap=10（検証済み）

**オプション2**: プロジェクト完了宣言
- Phase 5.2を85%完了として終了
- 論文・報告書作成に進む
- 実用実装への移行

**オプション3**: 性能改善継続
- 前処理改善（提案1）から開始
- 実装期間2-3週間
- ILU(0)、マルチグリッド前処理

---

## 📝 重要な技術的知見（再確認）

### Python-Julia比較結果（小ウィンドウ）

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
| 中ウィンドウ (window=35) | **測定予定** | **測定予定** |

**重要**: 最適basesizeは問題構造に依存

### スライディングウィンドウアルゴリズムの仕様

**重要な発見**:
- 最後のウィンドウが短くなると、1ステップずつ進む短いウィンドウが自動生成される
- これはアルゴリズムの仕様であり、エラーではない
- 主要ウィンドウ数と短いウィンドウ数を区別する必要がある

**パラメータ設計の指針**:
- `(nt - 1) % (window - overlap)` が小さくなるようにパラメータを選択
- 小ウィンドウ: 5個中4個が主要（80%）
- 中ウィンドウ: 14個中5個が主要（36%）← これでも適切

---

## 📊 プロジェクト全体の完成度

### 完了項目

- ✅ **Phase 1-6**: 505テスト全合格
- ✅ **Phase 5.2 Step 1-3 Phase 1**: basesize最適化完了（85%）
- ✅ **Python-Julia比較**: 小ウィンドウ（CGM 3回/10回）完了
- ✅ **最終レポート**: Python-Julia比較最終レポート完成
- ✅ **Python版CGMバグ修正**: 完了
- ✅ **中規模ウィンドウ計画書v2.0**: 作成・検証・コミット完了

### 未完了項目

- ⏸️ **Phase 5.2 Step 3 Phase 2**: 中規模ウィンドウ測定（85%→100%）
  - 実行時間: 12-14時間（1日で完了可能）
  - パラメータ: nt=100, window=35, overlap=10（検証済み）
- 📅 **将来**: 性能改善（Phase 1-3実装）
  - 前処理改善、初期推定値改善など

---

## 💡 次セッションでの推奨アクション

### 最優先推奨（★★★★★）

**中規模ウィンドウ比較計画v2.0の実行**

1. 計画書v2.0の確認（パラメータ: nt=100, window=35, overlap=10）
2. Phase 1を開始（Julia basesize最適化、5-6時間）
   - basesize=500, 1000, 2000の3パターン測定
3. Phase 2に進行（Python-Julia比較、6-8時間）
   - Python版とJulia版（最適basesize）の完全比較
4. 中規模ウィンドウ比較レポート作成
5. Phase 5.2完了宣言

**期待される成果**:
- Phase 5.2の実質完了（小+中ウィンドウ測定完了）
- 小→中スケーラビリティの完全評価（主要4個→5個）
- basesize最適化の完全検証
- 次フェーズへの円滑な移行

**実行可能性**: 1日（12-14時間）で完了可能 ✅

---

**次セッション開始時は、このファイルを必ず最初に読んでください。**

**推奨事項**: 中規模ウィンドウ比較計画v2.0の採用を強く推奨します。パラメータは検証済みです。

**重要**: コマンド実行時は必ずパラメータを確認してください（window=35, overlap=10）。
