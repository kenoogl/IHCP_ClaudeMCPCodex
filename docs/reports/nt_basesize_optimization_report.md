# nt値とbasesize最適化レポート

**作成日**: 2025年10月25日
**測定対象**: DHCP単体ソルバー
**格子サイズ**: 80×100×20 (N=160,000セル)

---

## エグゼクティブサマリー

nt値（タイムステップ数: 10/30/50/100）とFLoops basesizeパラメータ（400/500/600/700）の組み合わせ16パターンを体系的に測定し、nt値に応じた最適basesize設定を特定しました。

### 主要な発見

> **注意**: このセクションは測定完了後に実データで更新してください

1. **全体最適値**: [測定後に記入]
2. **nt依存性**: [測定後に記入]
3. **basesize依存性**: [測定後に記入]
4. **推奨設定**: [測定後に記入]

---

## 測定結果（完全版）

### 測定データ一覧

> **データソース**: `shared/results/nt_basesize/measurement_results.csv`

| nt | basesize=400 | basesize=500 | basesize=600 | basesize=700 |
|----|-------------|-------------|-------------|-------------|
| **10** | [TBD] | [TBD] | [TBD] | [TBD] |
| **30** | [TBD] | [TBD] | [TBD] | [TBD] |
| **50** | [TBD] | [TBD] | [TBD] | [TBD] |
| **100** | [TBD] | [TBD] | [TBD] | [TBD] |

**単位**: 秒（DHCP実行時間）

### nt値別の最適basesize

| nt | 最適basesize | DHCP時間 | 備考 |
|----|-------------|----------|------|
| **10** | [TBD] | [TBD] | [測定後に記入] |
| **30** | [TBD] | [TBD] | [測定後に記入] |
| **50** | [TBD] | [TBD] | [測定後に記入] |
| **100** | [TBD] | [TBD] | [測定後に記入] |

### basesize別の最適nt

| basesize | 最速nt | DHCP時間 | 備考 |
|----------|--------|----------|------|
| **400** | [TBD] | [TBD] | [測定後に記入] |
| **500** | [TBD] | [TBD] | [測定後に記入] |
| **600** | [TBD] | [TBD] | [測定後に記入] |
| **700** | [TBD] | [TBD] | [測定後に記入] |

---

## 詳細分析

### 1. nt値の影響

#### タイムステップ数とスケーリング

> **分析**: 測定完了後、nt値に対する実行時間のスケーリング特性を分析

```
実行時間（秒）
 |
 |  [測定後にグラフ追加]
 |
 |_____________________________ → nt値
    10    30    50    100
```

**期待される傾向**:
- nt値に対してほぼ線形にスケール（反復ソルバーのため）
- 各nt値で最適basesizeは異なる可能性

### 2. basesizeの影響

#### basesizeと実行時間の関係（nt別）

> **分析**: 測定完了後、各nt値でのbasesize依存性を分析

```
実行時間（秒）
 |
 |  [測定後にグラフ追加]
 |
 |_____________________________ → basesize
    400   500   600   700
```

**期待される傾向**:
- U字カーブ（小さすぎても大きすぎても遅い）
- 最適basesizeはnt値に応じて変化する可能性

### 3. nt × basesizeの相互作用

#### 相互作用マトリックス（高速化率）

> **基準**: nt=10, basesize=400の実行時間を1.00×とした相対値

| nt | basesize=400 | basesize=500 | basesize=600 | basesize=700 |
|----|-------------|-------------|-------------|-------------|
| **10** | 1.00× | [TBD] | [TBD] | [TBD] |
| **30** | [TBD] | [TBD] | [TBD] | [TBD] |
| **50** | [TBD] | [TBD] | [TBD] | [TBD] |
| **100** | [TBD] | [TBD] | [TBD] | [TBD] |

---

## 推奨設定

### 実用推奨値

> **注意**: 測定完了後に実データに基づいて更新

#### nt値に応じた推奨basesize

```julia
# nt値に応じた最適basesize設定
function get_optimal_basesize(nt::Int)
  if nt <= 10
    return 600  # [測定後に更新]
  elseif nt <= 30
    return 600  # [測定後に更新]
  elseif nt <= 50
    return 600  # [測定後に更新]
  else
    return 600  # [測定後に更新]
  end
end
```

#### 汎用設定（4スレッド固定）

```julia
# 最速構成
JULIA_NUM_THREADS=4
const OPTIMAL_BASESIZE = 600  # [測定後に更新]
```

### 避けるべき設定

> **注意**: 測定完了後に実データに基づいて更新

```julia
# ❌ 避けるべき設定
const BAD_BASESIZE_1 = [TBD]  # [理由を記入]
const BAD_BASESIZE_2 = [TBD]  # [理由を記入]
```

---

## 理論的考察

### basesizeの最適値とnt値の関係

**理論式**:
```
最適basesize ≈ f(総要素数, スレッド数, nt)
```

> **分析**: 測定完了後、最適basesizeのnt依存性を理論的に考察

**仮説**:
1. nt値が大きいほど、キャッシュ効率が重要になる → basesize小さめ?
2. nt値が小さいほど、並列化オーバーヘッド削減が重要 → basesize大きめ?

**実測値との比較**: [測定後に記入]

---

## 測定環境

- **CPU**: Apple Silicon（4コア想定）
- **Julia**: 1.x
- **FLoops**: ThreadedEx backend
- **スレッド数**: 4固定
- **問題サイズ**: 80×100×20 = 160,000セル
- **ソルバー**: PBICGSTAB
- **前処理**: Gauss-Seidel

---

## 再現手順

### 1. 測定実行

```bash
# 新規測定（nt=30,50,100）
chmod +x run_nt_basesize_measurements.sh
./run_nt_basesize_measurements.sh

# 中断から再開する場合
./run_nt_basesize_measurements.sh --resume
```

### 2. データ抽出

```bash
# CSVファイル生成（nt=10の既存データも含む）
julia julia/scripts/extract_nt_basesize_data.jl
```

出力: `shared/results/nt_basesize/measurement_results.csv`

### 3. レポート更新

1. CSVデータを確認
2. このMarkdownファイルの [TBD] 箇所を実データで更新
3. グラフ・考察セクションを完成

---

## データファイル

### 測定ログ（16ファイル）

**新規測定** (nt=30,50,100):
```
shared/results/nt_basesize/nt30_bs400.log
shared/results/nt_basesize/nt30_bs500.log
shared/results/nt_basesize/nt30_bs600.log
shared/results/nt_basesize/nt30_bs700.log
shared/results/nt_basesize/nt50_bs400.log
shared/results/nt_basesize/nt50_bs500.log
shared/results/nt_basesize/nt50_bs600.log
shared/results/nt_basesize/nt50_bs700.log
shared/results/nt_basesize/nt100_bs400.log
shared/results/nt_basesize/nt100_bs500.log
shared/results/nt_basesize/nt100_bs600.log
shared/results/nt_basesize/nt100_bs700.log
```

**既存データ** (nt=10):
```
shared/results/threads_basesize/step1_t4_bs400.log
shared/results/threads_basesize/step1_t4_bs500.log
shared/results/threads_basesize/step1_t4_bs600.log
shared/results/threads_basesize/step1_t4_bs700.log
```

### 集約データ

```
shared/results/nt_basesize/measurement_results.csv  # 全16測定の集約
shared/results/nt_basesize/progress.txt             # 進捗記録
shared/results/nt_basesize/completed.txt            # 完了リスト
```

---

## 次のステップ

### 追加検証の推奨

1. **より大きな問題サイズ**: 200×200×40格子での検証
2. **スライディングウィンドウ**: 実際のCGM計算での最適basesize検証
3. **他のソルバー**: PCGでの検証
4. **メモリ使用量**: 各設定でのメモリプロファイル

---

## 結論

> **注意**: 測定完了後に記入

**最適設定**:
- [測定後に記入]

**重要な知見**:
1. [測定後に記入]
2. [測定後に記入]
3. [測定後に記入]

---

## 変更履歴

| 日付 | バージョン | 変更内容 |
|------|-----------|---------|
| 2025-10-25 | 0.1 | テンプレート作成 |
| [TBD] | 1.0 | 測定完了、実データ反映 |

---

**レポート作成者**: Claude Code
**測定実施**: [測定完了後に日付記入]
