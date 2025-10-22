# 次のセッション開始ガイド

**作成日時**: 2025年10月23日
**セッション状態**: Phase 1-4完了（全前処理器検証完了）
**現在のブランチ**: parallelization
**最新コミット**: 34a115e

---

## 🎯 現在の状態

### ✅ 完了した作業（Phase 1-4）

#### Phase 1: 全スクリプトの並列化対応 ✅

1. **Phase 1-1**: スレッド数表示の追加
2. **Phase 1-2**: --parオプションの追加
3. **Phase 1-3**: CGMSolver.jlのparパラメータ伝播

#### Phase 2: デフォルト設定の変更 ✅

全ソルバーと主要スクリプトのデフォルトを`par="thread"`に変更

#### Phase 3: テスト実施と問題修正 ✅

1. **Phase 3-1**: 動作確認テスト
   - 1スレッド実行: ✅ 正常動作（65.16秒）
   - 8スレッド実行: ⚠️ 反復回数異常増加を検出（87回/step）

2. **Phase 3-2**: 問題調査と修正
   - **問題**: Jacobi前処理器の並列実行時レースコンディション
   - **修正**: 2バッファ方式を正しく実装
   - **結果**: 反復回数が正常化（87回→21.8回/step）

3. **Phase 3-3**: 修正後の動作確認
   - 8スレッド + jacobi前処理: ✅ 正常動作（30.61秒、2.13倍高速化）

#### Phase 4: 全前処理器の並列化検証 ✅ **完了**

**検証実施**（4前処理器 × 2スレッド設定 = 8実行）:

1. **none前処理器**
   - 1スレッド: 94.80秒、8スレッド: 37.12秒
   - スピードアップ: 2.63倍、並列化効率: 32.8%
   - 精度: ❌ 最大相対誤差3.68%（絶対誤差0.05 W/m²は小）

2. **diagonal前処理器**
   - 1スレッド: 147.62秒、8スレッド: 50.76秒
   - スピードアップ: 2.98倍、並列化効率: 37.2%
   - 精度: ❌ 最大相対誤差6.68%（絶対誤差0.11 W/m²は小）

3. **jacobi前処理器** ⭐ **推奨**
   - 1スレッド: 170.75秒、8スレッド: 50.30秒
   - スピードアップ: **3.39倍**、並列化効率: **42.4%** 🏆
   - 精度: ✅ 最大相対誤差1.19e-7（機械精度レベル）

4. **gs前処理器** ⭐ **推奨**
   - 1スレッド: 97.03秒、8スレッド: **34.00秒** 🏆 最速
   - スピードアップ: 2.94倍、並列化効率: 36.8%
   - 精度: ✅ 最大相対誤差1.12e-7（機械精度レベル）

**検証結果の詳細レポート**: `docs/reports/phase4_preconditioner_validation_report.md`

---

## 📊 コミット履歴

```
34a115e docs: Phase 4全前処理器（none, diagonal, jacobi, gs）並列化検証完了
d410a08 docs: 次セッション用ガイド更新（Phase 4包括的検証準備完了）
87ebe49 docs: Phase 4包括的な精度検証計画とスクリプト作成
54a9414 docs: 次セッション用ガイド更新（Phase 1-3完了、Jacobi修正完了）
b1e4c66 fix: Jacobi前処理器の並列実行時レースコンディションを修正
```

---

## 📈 性能測定結果まとめ

### 全前処理器比較（CGM反復2回、nt=10）

| 前処理器 | 1スレッド | 8スレッド | スピードアップ | 並列化効率 | 精度判定 | 推奨度 |
|---------|----------|----------|--------------|-----------|---------|-------|
| **jacobi** | 170.75s | 50.30s | **3.39倍** 🏆 | **42.4%** | ✅ 1.19e-7 | ⭐⭐⭐ |
| **gs** | 97.03s | **34.00s** 🏆 | 2.94倍 | 36.8% | ✅ 1.12e-7 | ⭐⭐⭐ |
| **none** | 94.80s | 37.12s | 2.63倍 | 32.8% | ❌ 3.68% | ⚠️ |
| **diagonal** | 147.62s | 50.76s | 2.98倍 | 37.2% | ❌ 6.68% | ⚠️ |

**結論**:
- **最高スピードアップ**: jacobi（3.39倍）
- **最速実行時間**: gs（34.00秒、8スレッド）
- **精度合格**: jacobi、gs（機械精度レベル）
- **実用推奨**: jacobi前処理器（精度と性能の最良バランス）

---

## 🚀 次に実施すべき作業

### Phase 5: 性能ベンチマーク（推奨、オプション）

**目的**: スケーラビリティの定量評価

#### 測定内容
1. スレッド数スケーリング（1, 2, 4, 8スレッド）
2. 問題サイズ依存性（小規模・中規模）
3. 並列化効率の定量評価

#### 実行例（jacobi前処理器）

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# 小規模問題（nt=10, window=5, overlap=2, CGM反復2）
for threads in 1 2 4 8; do
  echo "=== $threads threads ==="
  JULIA_NUM_THREADS=$threads julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 10 --window 5 --overlap 2 --cgm-iter 2 \
    --precond jacobi \
    2>&1 | tee "shared/results/bench_small_${threads}threads.log"
done

# 中規模問題（nt=100, window=30, overlap=10, CGM反復2）
for threads in 1 2 4 8; do
  echo "=== $threads threads ==="
  JULIA_NUM_THREADS=$threads julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 100 --window 30 --overlap 10 --cgm-iter 2 \
    --precond jacobi \
    2>&1 | tee "shared/results/bench_medium_${threads}threads.log"
done
```

#### 期待される結果
- 理想的なスケーリング: 8スレッドで8倍（効率100%）
- 実測見込み: 8スレッドで3-4倍（効率37-50%）
- 問題サイズ大で効率向上（過去の傾向）

---

### 代替案: 本番適用準備

Phase 5をスキップして本番適用する場合：

1. **デフォルト設定の確認**
   - 前処理器: jacobi（推奨）
   - スレッド数: 8（環境に応じて調整）
   - 並列化モード: thread

2. **ドキュメント整備**
   - ユーザーガイド作成
   - 性能ベンチマーク結果の公開

3. **本番データでの試験運用**
   - 大規模データセットでの動作確認
   - メモリ使用量のモニタリング

---

## 📂 重要なファイル

### 実装済みファイル
- `julia/src/solvers/CommonSolver.jl` ⭐ 全前処理器実装済み
- `julia/scripts/run_sliding_window.jl` ⭐ 並列化対応完了
- `julia/scripts/compare_parallel_results_with_args.jl` ⭐ 精度検証スクリプト

### Phase 4生成ファイル
- `docs/reports/phase4_preconditioner_validation_report.md` ⭐ **包括的レポート**
- `shared/results/validation_{none,diagonal,jacobi,gs}.txt` - 各前処理器検証結果
- `shared/results/julia_sliding_window_cgm2_metadata.txt` - CGM2メタデータ

### データファイル（未コミット、ローカル保存）
```
shared/results/
├── julia_sliding_window_cgm2_1thread_none.npz
├── julia_sliding_window_cgm2_8threads_none.npz
├── julia_sliding_window_cgm2_1thread_diagonal.npz
├── julia_sliding_window_cgm2_8threads_diagonal.npz
├── julia_sliding_window_cgm2_1thread_jacobi.npz
├── julia_sliding_window_cgm2_8threads_jacobi.npz
├── julia_sliding_window_cgm2_1thread_gs.npz
└── julia_sliding_window_cgm2_8threads_gs.npz

（各ファイル約1.8MB、計14.4MB）
```

### ログファイル（未コミット、ローカル保存）
```
shared/results/
├── log_1thread_none.txt
├── log_8threads_none.txt
├── log_1thread_diagonal.txt
├── log_8threads_diagonal.txt
├── log_1thread_jacobi.txt
├── log_8threads_jacobi.txt
├── log_1thread_gs.txt
└── log_8threads_gs.txt
```

---

## ⚠️ 重要な注意事項

### スレッド数は明示的に指定（必須）

```bash
# 正しい（8スレッドで並列実行）
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_sliding_window.jl

# 間違い（1スレッドになる）
julia --project=julia julia/scripts/run_sliding_window.jl
```

### 前処理器の選択

#### 推奨: jacobi前処理器
- **精度**: ✅ 機械精度レベル（1.19e-7）
- **性能**: ✅ 最高スピードアップ（3.39倍）
- **安定性**: ✅ レースコンディション修正済み
- **用途**: 全般的な用途に推奨

#### 推奨: gs前処理器
- **精度**: ✅ 機械精度レベル（1.12e-7）
- **性能**: ✅ 最速実行時間（34秒）
- **収束性**: ✅ 最少反復回数
- **用途**: 最速実行が必要な場合

#### 検討対象: none/diagonal前処理器
- 相対誤差が閾値超過だが、絶対誤差は小さい
- 実用上の問題は要検討
- 特殊用途以外は非推奨

### Jacobi前処理器の修正（重要）

✅ **修正済み**（b1e4c66コミット）
- 並列実行時のレースコンディションを2バッファ方式で解消
- 反復回数が正常化（87回→21.8回/step）
- Phase 4検証で高精度と高性能を確認

---

## 🔗 次のセッション開始手順

### 1. 環境確認

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5
```

**期待される出力**:
```
On branch parallelization
Your branch is up to date with 'origin/parallelization'.
...

34a115e docs: Phase 4全前処理器（none, diagonal, jacobi, gs）並列化検証完了
d410a08 docs: 次セッション用ガイド更新（Phase 4包括的検証準備完了）
...
```

### 2. このファイル確認

```bash
cat TODO_NEXT_SESSION.md
```

### 3. Phase 4完了状態の確認

```bash
# レポート確認
cat docs/reports/phase4_preconditioner_validation_report.md

# データファイル確認
ls -lh shared/results/julia_sliding_window_cgm2_*thread*.npz

# 検証結果確認
ls -lh shared/results/validation_*.txt
```

### 4. 次の作業を選択

#### オプションA: Phase 5実施（推奨）
スケーラビリティベンチマークの実施

#### オプションB: 本番適用準備
ドキュメント整備と本番データでの試験運用

#### オプションC: 精度調査（オプション）
none/diagonal前処理器の精度問題の詳細調査

---

## 📝 主要な成果物

### Phase 4で作成されたドキュメント

1. **包括的検証レポート** ⭐
   - ファイル: `docs/reports/phase4_preconditioner_validation_report.md`
   - 内容: 全4前処理器の詳細な性能・精度比較
   - サイズ: 11KB

2. **検証結果ファイル**
   - `shared/results/validation_none.txt`
   - `shared/results/validation_diagonal.txt`
   - `shared/results/validation_gs.txt`
   - ※ `validation_jacobi.txt`は前セッションで生成済み

3. **メタデータファイル**
   - `shared/results/julia_sliding_window_cgm2_metadata.txt`
   - 最新CGM実行の詳細パラメータ

---

## 🎓 学んだこと（Phase 4から）

### 並列化の知見

1. **前処理器による性能差**
   - jacobi: 最高のスピードアップ（3.39倍）
   - gs: 最速の絶対実行時間（34秒）
   - none/diagonal: 並列化効率が相対的に低い（32-37%）

2. **精度と性能のトレードオフ**
   - 機械精度を維持できる前処理器: jacobi, gs
   - 相対誤差が大きい前処理器: none, diagonal
   - ただし絶対誤差はすべて実用上小さい

3. **並列化効率の向上傾向**
   - CGM反復1回: 26.6%（jacobi）
   - CGM反復2回: 42.4%（jacobi）
   - 問題サイズが大きいほど効率向上

### 実装上の教訓

1. **レースコンディションの重要性**
   - Jacobi前処理器の2バッファ方式が成功の鍵
   - 並列実行時の共有データアクセスに注意

2. **検証の重要性**
   - 性能だけでなく精度検証も必須
   - 相対誤差と絶対誤差の両方を評価

3. **スケーラビリティ**
   - 理想的な線形スケーリングは困難（効率37-42%）
   - 実用上は3倍程度のスピードアップで十分有効

---

**次のセッション開始時**: 必ずこのファイルを読んでから作業を開始してください。

**最重要タスク**: Phase 4完了、次はPhase 5（性能ベンチマーク）または本番適用準備
