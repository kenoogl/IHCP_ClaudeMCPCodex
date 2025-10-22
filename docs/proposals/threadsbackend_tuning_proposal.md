# ThreadsBackendチューニング提案

**作成日時**: 2025年10月23日
**対象**: FLoops.jl の ThreadsBackend パラメータ最適化
**目的**: 並列化効率の改善（現在36% → 目標50-60%）

---

## 🎯 現状の問題

### 現在の実装（commons.jl:39-41）
```julia
function get_backend(par::String)
    return (par == "thread") ? ThreadedEx() : SequentialEx()
end
```

**問題点**:
- `ThreadedEx()` = デフォルトパラメータ使用
- `basesize` 未指定 → FLoopsの自動判定
- `stride` 未指定 → デフォルト動作
- **結果**: 8スレッドで効率36%（Phase 5測定）

---

## 📚 ThreadsBackend パラメータ解説

### ThreadedEx() の内部動作
```julia
# デフォルト（現在の実装）
ThreadedEx()

# 完全な形
ThreadedEx(basesize = 1, stride = 1, nthreads = nothing)
```

### パラメータの意味

#### 1. `basesize`（チャンクサイズ）
**定義**: 各スレッドが一度に処理する最小反復回数

```julia
# 例: basesize=1000
ThreadedEx(basesize=1000)

# 効果:
- スレッドあたり最低1,000回の反復を保証
- オーバーヘッド削減
- キャッシュ局所性向上
```

**現在の推定値**:
```
総反復: 137,592回
8スレッド: 17,199回/スレッド（自動分割）
```

#### 2. `stride`（メモリアクセスパターン）
**定義**: 連続するスレッドの反復間隔

```julia
# stride=1（デフォルト、連続ブロック）
Thread 1: 反復 1, 2, 3, ..., N/8
Thread 2: 反復 N/8+1, N/8+2, ...

# stride=8（インターリーブ）
Thread 1: 反復 1, 9, 17, 25, ...
Thread 2: 反復 2, 10, 18, 26, ...
```

**効果**:
- stride=1: キャッシュ局所性◎、False Sharing△
- stride=大: キャッシュ局所性△、False Sharing◎

---

## 🔬 最適化戦略

### 戦略1: basesize チューニング（推奨）

#### 理論的な最適値
```
目標: オーバーヘッド < 5%

現在の問題サイズ:
  総反復: 137,592回
  8スレッド: 17,199回/スレッド

最適なbasesize推定:
  basesize = 5,000 ~ 20,000

理由:
- 小さすぎる（<1,000）: オーバーヘッド増大
- 適切（5,000-20,000）: バランス良好
- 大きすぎる（>50,000）: 負荷不均衡
```

#### 提案実装
```julia
function get_backend(par::String; basesize::Int=10000)
    if par == "thread"
        return ThreadedEx(basesize=basesize)
    else
        return SequentialEx()
    end
end
```

#### 期待される効果
```
basesize=10,000の場合:
  チャンク数: 137,592 ÷ 10,000 ≈ 14チャンク
  8スレッドでの分配: 各スレッドに1-2チャンク

予測:
- オーバーヘッド削減: 20-30%
- 効率向上: 36% → 45-50%
```

### 戦略2: stride チューニング（False Sharing対策）

#### 問題の診断
```julia
# False Sharingが疑われる場合
# 例: 連続する(i,j,k)への同時書き込み

@floop backend for k in 2:nk-1, j in 2:nj-1, i in 2:ni-1
    array[i,j,k] = f(i,j,k)
end

# i方向が連続 → キャッシュライン競合の可能性
```

#### 提案（要実験的検証）
```julia
# False Sharing軽減のためのstride
function get_backend_with_stride(par::String)
    if par == "thread"
        # キャッシュラインサイズ（64byte）考慮
        # Float64 = 8byte → 8要素/キャッシュライン
        return ThreadedEx(basesize=10000, stride=8)
    else
        return SequentialEx()
    end
end
```

**注意**: strideは慎重に調整が必要
- キャッシュ局所性とのトレードオフ
- 実測による検証必須

---

## 📊 測定計画

### Phase 5.1: basesize スイープ測定

```bash
# 小規模問題でbasesize最適化
for basesize in 1 1000 5000 10000 20000 50000; do
  echo "=== basesize=$basesize ==="
  JULIA_NUM_THREADS=8 julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 10 --window 5 --overlap 2 --cgm-iter 2 --precond gs \
    --basesize $basesize \
    2>&1 | tee "shared/results/bench_gs_8t_basesize${basesize}.log"
done
```

**測定項目**:
- 実行時間
- スレッドあたりのチャンク数
- 並列化効率

**期待される結果**:
```
basesize=1 (デフォルト): 33.99秒、効率36%
basesize=5,000:         ~30秒、効率42%（予測）
basesize=10,000:        ~28秒、効率48%（予測）
basesize=20,000:        ~27秒、効率50%（予測）
basesize=50,000:        ~29秒、効率46%（予測、負荷不均衡）
```

### Phase 5.2: stride 検証（オプション）

```bash
# stride効果の検証
for stride in 1 4 8 16; do
  echo "=== stride=$stride ==="
  JULIA_NUM_THREADS=8 julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 10 --window 5 --overlap 2 --cgm-iter 2 --precond gs \
    --basesize 10000 --stride $stride \
    2>&1 | tee "shared/results/bench_gs_8t_stride${stride}.log"
done
```

---

## 💻 実装提案

### 修正1: commons.jl

```julia
"""
並列動作バックエンドを返す

# キーワード引数
- `basesize::Int`: チャンクサイズ（デフォルト: 10000）
- `stride::Int`: メモリアクセスストライド（デフォルト: 1）
"""
function get_backend(par::String; basesize::Int=10000, stride::Int=1)
    if par == "thread"
        return ThreadedEx(basesize=basesize, stride=stride)
    else
        return SequentialEx()
    end
end
```

### 修正2: 各ソルバーの呼び出し側

#### オプションA: デフォルト値使用（推奨）
```julia
# 変更不要（commons.jlのデフォルトが適用される）
backend = get_backend(par)
```

#### オプションB: 明示的指定
```julia
# 細かい調整が必要な場合
backend = get_backend(par, basesize=15000, stride=1)
```

### 修正3: コマンドライン引数追加（run_sliding_window.jl）

```julia
# パーサー追加
@add_arg_table! s begin
    "--basesize"
        help = "ThreadsBackend basesize (chunk size)"
        arg_type = Int
        default = 10000
    "--stride"
        help = "ThreadsBackend stride (memory access pattern)"
        arg_type = Int
        default = 1
end

# バックエンド生成
backend = get_backend(args["par"],
                      basesize=args["basesize"],
                      stride=args["stride"])
```

---

## 📈 期待される効果

### 保守的な予測
```
現在（basesize=デフォルト）:
  8スレッド効率: 36%
  実行時間: 33.99秒

最適化後（basesize=10,000）:
  8スレッド効率: 45-50%（+9-14ポイント）
  実行時間: 27-30秒（約15-20%高速化）
```

### 楽観的な予測
```
最適化後 + stride調整:
  8スレッド効率: 50-55%
  実行時間: 25-28秒（約20-25%高速化）
```

---

## ⚠️ 注意事項

### 1. 問題サイズ依存性
```
小規模問題（137,592反復）:
  basesize=10,000 → 14チャンク → 効果大

大規模問題（1,100,736反復、空間格子拡大）:
  basesize=10,000 → 110チャンク → 効果さらに大
```

### 2. 前処理器による違い
```
GS前処理器:
  現在: 効率36%
  予測: 効率45-50%（basesize最適化）

Jacobi前処理器:
  現在: 効率42%（Phase 4データ）
  予測: 効率50-55%（basesize最適化）
```

### 3. ハードウェア依存性
```
CPU: Apple M2 Pro（8コア）
L1キャッシュ: 128KB/コア
L2キャッシュ: 不明
メモリバンド幅: 200GB/s

→ basesizeの最適値はハードウェアに依存
→ 実測による調整が必須
```

---

## 🎯 推奨アクション

### Step 1: basesize スイープ測定（必須）
1. commons.jlにbasesize引数を追加
2. Phase 5.1測定スクリプト実行
3. 最適basesizeを特定

**優先度**: 🔴 高（即座に効果が期待できる）
**工数**: 1-2時間（実装30分、測定1-1.5時間）

### Step 2: 最適basesize適用（必須）
1. 測定結果から最適値を決定
2. デフォルト値として設定
3. 全ソルバーに適用

**優先度**: 🔴 高
**工数**: 30分

### Step 3: stride 検証（オプション）
1. False Sharingが疑われる場合のみ実施
2. Phase 5.2測定スクリプト実行
3. 効果があれば適用

**優先度**: 🟡 中（basesize最適化後に判断）
**工数**: 1-2時間

---

## 📝 参考資料

### FLoops.jl ドキュメント
- https://juliafolds.github.io/FLoops.jl/stable/
- ThreadsBackend API: https://juliafolds.github.io/FLoops.jl/stable/reference/api/#FLoops.ThreadedEx

### 関連する並列化の知見
- キャッシュラインサイズ: 64バイト（一般的）
- Float64配列: 8バイト/要素 → 8要素/キャッシュライン
- False Sharing回避: 隣接要素への同時アクセスを避ける

---

## 🎓 まとめ

### 現状の課題
- デフォルトパラメータ使用で効率36%
- basesize未指定によるオーバーヘッド

### 解決策
- basesize=10,000程度に設定
- オーバーヘッド削減とキャッシュ効率向上

### 期待効果
- 効率: 36% → 45-50%（+9-14ポイント）
- 実行時間: 15-20%高速化

### 次のステップ
1. commons.jl修正（basesize引数追加）
2. Phase 5.1測定（basesize スイープ）
3. 最適値適用

**これは空間格子サイズに依存しない改善策であり、
現行格子（80×100×20）でも即座に効果が期待できる！**
