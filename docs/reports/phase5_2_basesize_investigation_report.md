# Phase 5.2 Basesize効果検証レポート

**日付**: 2025年10月23日
**目的**: FLoopsのbasesizeパラメータが性能に与える影響を検証
**ブランチ**: parallelization

---

## 実行概要

Julia版IHCP-CGMソルバーにおいて、FLoopsの`basesize`パラメータ（並列処理の最小チャンクサイズ）が性能に与える影響を検証しました。

---

## Step 2: 10ステップフルサイズテスト結果

### テスト条件
- **格子**: 80×100×20 (N=160,000)
- **時間ステップ**: 10ステップ
- **CGM反復**: 1回
- **ソルバー**: PBICGSTAB
- **前処理**: Gauss-Seidel

### 結果サマリー

| 設定 | スレッド数 | 並列モード | basesize | 実行時間 | 順位 |
|------|-----------|-----------|----------|---------|------|
| **最速** | 4 | thread | **1000** | **19.45秒** | 🥇 1位 |
| 良好 | 4 | sequential | 1000 | 32.04秒 | 2位 |
| 良好 | 1 | thread | 1000 | 44.57秒 | 3位 |
| 遅い | 4 | thread | 1 | 295.26秒 | 4位 |
| **最悪** | 1 | thread | **1** | **8054.46秒** | ❌ 5位 |
| - | 4 | thread | 10000 | 28.66秒 | 参考 |

### 重要な発見

#### 1. basesize=1の致命的問題 ⚠️

**Thread=1 + basesize=1の組み合わせで壊滅的な性能劣化**:
- 実行時間: **8054秒 (134分)**
- basesize=1000比で **180倍遅い** (44.57秒 → 8054秒)
- 実用不可能なレベル

**原因推定**:
- チャンクサイズが小さすぎることによるオーバーヘッド
- スレッド数1の場合、タスクスケジューリングのコストが支配的
- 160,000要素を1要素ずつ処理するため、関数呼び出しが膨大

#### 2. 最適basesize値

**basesize=1000が最適**:
- Thread=4環境で **19.45秒** を達成
- basesize=1比で **15.2倍高速化** (295.26秒 → 19.45秒)
- basesize=10000 (28.66秒) よりも高速

**理由**:
- チャンクサイズ1000 ≈ 格子サイズ(80×100=8000)の12.5%
- タスク数とスケジューリングコストのバランスが最適
- キャッシュ効率とスレッド並列性の両立

#### 3. 並列モードの影響

**sequentialモードでも高速**:
- Thread=4 + sequential + basesize=1000: **32.04秒**
- threadモード(19.45秒)の1.6倍程度
- 並列化なしでも実用的な速度

---

## Step 3: スライディングウィンドウテスト結果

### テスト条件
- **格子**: 80×100×20
- **時間ステップ**: 10ステップ
- **ウィンドウサイズ**: 5ステップ
- **オーバーラップ**: 2ステップ
- **CGM反復**: 3回

### 結果サマリー

| basesize | 並列モード | 実行時間 | 順位 |
|----------|-----------|---------|------|
| **500** | thread | **34.34秒** | 🥇 1位 |
| 100 | thread | 43.85秒 | 2位 |
| 1000 | thread | 50.30秒 | 3位 |
| 10000 | thread | 84.64秒 | 4位 |
| 1000 | sequential | 95.53秒 | 5位 |
| 10 | thread | 136.11秒 | 6位 |

### 発見事項

#### 1. 小ウィンドウではbasesize=500が最適

**basesize=500が最速**:
- 実行時間: **34.34秒**
- basesize=1000 (50.30秒) より **1.46倍高速**
- basesize=10000 (84.64秒) より **2.46倍高速**

**理由推定**:
- 小ウィンドウ(5ステップ)では問題サイズが小さい
- より細かい並列粒度が有効
- スライディングウィンドウの反復処理に適合

#### 2. basesize=10の性能劣化

**basesize=10で大幅に遅化**:
- 実行時間: 136.11秒
- basesize=500比で **3.96倍遅い**
- チャンクサイズが小さすぎる影響

---

## 総合結論

### 1. basesize=1は絶対に避けるべき ⚠️

- Thread=1環境で **180倍の性能劣化**
- Thread=4環境でも **15倍遅い**
- デフォルト値として危険

### 2. 推奨basesize値

| 問題サイズ | 推奨basesize | 理由 |
|-----------|-------------|------|
| **フルサイズ (N=160,000)** | **1000** | 最速 (19.45秒) |
| **小ウィンドウ** | **500** | 最速 (34.34秒) |
| **一般的な値** | **500-1000** | 安全かつ高速 |

### 3. 実装への反映

**現在の実装**:
```julia
# julia/src/solvers/CommonSolver.jl
basesize = get(kwargs, :basesize, 1000)  # デフォルト1000
```

**推奨事項**:
- ✅ デフォルト値1000は適切
- ✅ ユーザーが変更可能な設計も維持
- ⚠️ basesize=1は警告を出すべき

### 4. 性能改善の実績

**Step 2 (10ステップフルサイズ)**:
- 最悪ケース: 8054秒 (basesize=1, thread=1)
- 最良ケース: 19.45秒 (basesize=1000, thread=4)
- **改善倍率: 414倍** 🚀

**Step 3 (スライディングウィンドウ)**:
- 最悪ケース: 136.11秒 (basesize=10)
- 最良ケース: 34.34秒 (basesize=500)
- **改善倍率: 3.96倍** 🚀

---

## 次のステップ

### Phase 5.2完了 ✅

**達成事項**:
- [x] Step 2: フルサイズテストでbasesize効果を検証
- [x] Step 3: 小ウィンドウテストでbasesize効果を検証
- [x] 最適basesize値の特定 (500-1000)
- [x] basesize=1の致命的問題を発見

### Phase 5.2 Step 4: 大ウィンドウ測定（保留中）

**保留中のテスト**:
- ウィンドウサイズ71、オーバーラップ17での測定
- basesize=500での性能確認

**状態**: 実行途中でキル（未完了）

### Phase 5.3への移行検討

**次の課題**:
1. Python版との性能比較
2. 大ウィンドウでの最適化検証
3. 総合的な性能レポート作成

---

## 付録: 詳細データ

### Step 2詳細ログ

すべての実行ログは以下に保存:
- `shared/results/step2_fullsize_basesize1.log`
- `shared/results/step2_fullsize_basesize1000.log`
- `shared/results/step2_fullsize_basesize10000.log`
- `shared/results/step2_fullsize_sequential_basesize1000.log`
- `shared/results/step2_fullsize_thread1_basesize1.log`
- `shared/results/step2_fullsize_thread1_basesize1000.log`

### Step 3詳細ログ

すべての実行ログは以下に保存:
- `shared/results/step3_sliding_small_basesize10.log`
- `shared/results/step3_sliding_small_basesize100.log`
- `shared/results/step3_sliding_small_basesize500.log`
- `shared/results/step3_sliding_small_basesize1000.log`
- `shared/results/step3_sliding_small_basesize10000.log`
- `shared/results/step3_sliding_small_sequential.log`

---

**作成者**: Claude Code
**検証環境**: macOS 14.6.0, Julia 1.10, 4スレッド
**最終更新**: 2025年10月23日 21:50
