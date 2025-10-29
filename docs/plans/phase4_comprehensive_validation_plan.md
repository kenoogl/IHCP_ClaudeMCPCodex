# Phase 4: 包括的な数値精度検証計画

**作成日**: 2025年10月22日
**ステータス**: 実施中

---

## 🎯 目的

BiCGSTABソルバーの4つの前処理器（none, diagonal, jacobi, gs）すべてに対して、並列化実装の数値精度検証を実施する。

---

## 📋 検証計画

### 前処理器一覧

1. **none**: 前処理なし（恒等変換）
2. **diagonal**: 対角スケーリング
3. **jacobi**: 加重Jacobi法（5回反復）
4. **gs**: Red-Black Gauss-Seidel法（5回反復）

### 検証条件

| パラメータ | 値 |
|-----------|---|
| nt | 10 |
| window | 5 |
| overlap | 2 |
| cgm-iter | 2 |
| スレッド数（基準） | 1 |
| スレッド数（並列化） | 8 |

### 測定項目

各前処理器について以下を測定：

1. **実行時間**
   - 1スレッド実行時間
   - 8スレッド実行時間
   - スピードアップ比
   - 並列化効率

2. **数値精度**
   - 最大相対誤差
   - 平均相対誤差
   - 最大絶対誤差

3. **収束性能**
   - 反復回数（DHCP, Adjoint, Sensitivity）
   - CGM目的関数値

---

## 🔬 実施手順

### Step 1: データ生成（各前処理器）

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# 前処理器リスト
PRECONDITIONERS=("none" "diagonal" "jacobi" "gs")

for precond in "${PRECONDITIONERS[@]}"; do
  echo "=== Processing preconditioner: $precond ==="

  # 1スレッド実行
  JULIA_NUM_THREADS=1 julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 10 --window 5 --overlap 2 --cgm-iter 2 \
    --precond $precond

  # 結果保存
  mv shared/results/julia_sliding_window_cgm2.npz \
     shared/results/julia_sliding_window_cgm2_1thread_${precond}.npz

  # 8スレッド実行
  JULIA_NUM_THREADS=8 julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 10 --window 5 --overlap 2 --cgm-iter 2 \
    --precond $precond

  # 結果保存
  mv shared/results/julia_sliding_window_cgm2.npz \
     shared/results/julia_sliding_window_cgm2_8threads_${precond}.npz
done
```

### Step 2: 精度検証（各前処理器）

```bash
# 精度検証スクリプト実行
for precond in "${PRECONDITIONERS[@]}"; do
  echo "=== Validating preconditioner: $precond ==="
  julia --project=julia julia/scripts/compare_parallel_results.jl \
    shared/results/julia_sliding_window_cgm2_1thread_${precond}.npz \
    shared/results/julia_sliding_window_cgm2_8threads_${precond}.npz \
    > shared/results/validation_${precond}.txt
done
```

### Step 3: 結果集計とレポート作成

実行時間、精度、収束性能を表形式でまとめる。

---

## 📊 期待される結果

### 仮説

1. **none前処理**
   - 反復回数: 最多
   - 実行時間: 中程度（反復回数多いが前処理オーバーヘッドなし）
   - 精度: 高い（並列化の影響を受けにくい）
   - スピードアップ: 良好（並列化オーバーヘッド最小）

2. **diagonal前処理**
   - 反復回数: 中程度
   - 実行時間: 短い（反復回数減少、前処理軽量）
   - 精度: 高い（対角スケーリングのみ、並列化安全）
   - スピードアップ: 最良

3. **jacobi前処理**
   - 反復回数: 少ない
   - 実行時間: 中程度（反復回数最少だが前処理オーバーヘッドあり）
   - 精度: 高い（修正済み、2バッファ方式）
   - スピードアップ: 良好

4. **gs前処理**
   - 反復回数: 最少
   - 実行時間: 場合による（前処理オーバーヘッド最大）
   - 精度: 要確認（Red-Black順序付けによる並列化）
   - スピードアップ: 要確認

### 成功基準

すべての前処理器で：
- **相対誤差 < 1e-6**（実用上十分な精度）
- **スピードアップ > 2.0倍**（8スレッドで）
- **反復回数が1スレッドと同等**（±10%以内）

---

## 📝 実施状況

### ✅ 完了

- [x] jacobi前処理での検証（Phase 4-1 ~ 4-4）
  - 1スレッド: 173.27秒
  - 8スレッド: 52.73秒（3.29倍高速化）
  - 最大相対誤差: 1.19e-7 ✅

### 🔄 進行中

- [ ] none前処理での検証
- [ ] diagonal前処理での検証
- [ ] gs前処理での検証
- [ ] 全前処理器の結果比較レポート

---

## 🚨 既知の課題

### Jacobi前処理器

✅ **修正済み**（b1e4c66コミット）
- 問題: 並列実行時のレースコンディション
- 修正: 2バッファ方式の正しい実装

### GS（Gauss-Seidel）前処理器

⚠️ **要確認**
- Red-Black順序付けが正しく実装されているか
- 並列実行時の収束性能
- 数値精度の保持

---

## 📂 生成されるファイル

```
shared/results/
├── julia_sliding_window_cgm2_1thread_none.npz
├── julia_sliding_window_cgm2_8threads_none.npz
├── julia_sliding_window_cgm2_1thread_diagonal.npz
├── julia_sliding_window_cgm2_8threads_diagonal.npz
├── julia_sliding_window_cgm2_1thread_jacobi.npz  (既存)
├── julia_sliding_window_cgm2_8threads_jacobi.npz (既存)
├── julia_sliding_window_cgm2_1thread_gs.npz
├── julia_sliding_window_cgm2_8threads_gs.npz
├── validation_none.txt
├── validation_diagonal.txt
├── validation_jacobi.txt
└── validation_gs.txt
```

---

## 🔗 関連ドキュメント

- 並列化実装計画: `docs/plans/parallelization_implementation_plan.md`
- 次セッションガイド: `TODO_NEXT_SESSION.md`
- 比較スクリプト: `julia/scripts/compare_parallel_results.jl`

---

**次のステップ**: 残り3つの前処理器（none, diagonal, gs）での検証を実施
