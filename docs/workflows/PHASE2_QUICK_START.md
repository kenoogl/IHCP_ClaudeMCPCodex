# Phase 2（並列化）クイックスタートガイド

**作成日時**: 2025年10月14日 02:00
**Phase 1完了状態**: 累積改善 -27.34%（1,072秒 → 779秒）

---

## 🎯 Phase 2の目標

### 目標性能

- **現在の実行時間**: 779.48秒（13.0分）
- **Phase 2目標**: 30-50%高速化
- **Phase 2完了後**: **390-545秒（6.5-9.1分）**
- **総合累積改善**: **49-64%**（ベースライン1,072秒比）
- **総合目標達成**: **50-60%改善に到達確実** 🎯

---

## 📋 Phase 1の成果（前提条件）

### 累積改善: -27.34%

| Phase | 改善効果 | 主な内容 |
|-------|---------|---------|
| Phase 0 | -0.34% | 配列アロケーション削減 |
| Phase 1-A | -0.15% | 型安定化 |
| **Phase 1-B** | **-16.56%** | 初期推定値改善（最大貢献）🏆 |
| **Phase 1-E調整版** | **-12.95%** | 適応的収束判定（第2の貢献）🎉 |

### 最適パラメータ（Phase 1完了時）

```julia
cgm_params = (
    # Phase 1-B最適化
    dhcp_extrapolation = :quadratic,          # DHCP二次外挿
    adjoint_residual_scale = 0.5,             # Adjoint残差スケール

    # Phase 1-E調整版（採用）
    adaptive_tol = true,                       # 適応的収束判定有効化
    kappa_theta = 0.1,                         # より保守的な緩和 ✨
    warmup_steps = 4,                          # より長いwarmup期間 ✨
    tol_min_dhcp = 1.0e-7,
    tol_min_adjoint = 1.0e-9,
    tol_min_sensitivity = 1.0e-7,
    tol_max_dhcp = 1.0e-5,
    tol_max_adjoint = 1.0e-5,
    tol_max_sensitivity = 1.0e-5,
)
```

---

## 🚀 Phase 2実装計画

### Week 1: 並列化基礎実装

#### 1. 並列化対象の選定（優先順位）

**高優先度（必須）**:
1. **PBICGSTAB!内部ループ**
   - CalcAX!, CalcRK!の並列化
   - 内積計算（Fdot1, Fdot2）の並列化
   - ベクトル演算（BiCG1!, Triad!, BICG2!）の並列化

2. **Preconditioner（RB-SOR）**
   - Red-Black着色による並列化
   - 既にマルチカラー構造あり（活用可能）

**中優先度**:
3. **CGMループ外の並列化**
   - 時間ステップ間の並列化（スライディングウィンドウ）
   - データ読み込みの並列化

#### 2. 並列化戦略

**既存の並列バックエンド活用**:
- FLoops.jl + ThreadsX.jl が既に実装済み
- `par="thread"` パラメータで並列化切り替え可能
- 既存の並列構造を活用・最適化

**スレッド分割戦略**:
- z方向（k軸）での分割が最適（データ局所性が高い）
- スレッド数: 4-8スレッド（Apple Silicon M-series）
- 負荷分散: 格子点数に基づく動的分割

#### 3. データ競合の回避

**リダクション操作**:
- `@reduce`を活用（FLoops.jl標準機能）
- 内積計算の並列リダクション
- 残差ノルム計算の並列リダクション

**ワーク配列管理**:
- スレッドローカルバッファの確保
- false sharingの回避（パディング）

### Week 2: 性能最適化とベンチマーク

#### 1. スレッド数チューニング

```bash
# 1スレッド（ベースライン）
JULIA_NUM_THREADS=1 julia --project=julia scripts/run_10steps_fullsize_test.jl

# 4スレッド
JULIA_NUM_THREADS=4 julia --project=julia scripts/run_10steps_fullsize_test.jl

# 8スレッド
JULIA_NUM_THREADS=8 julia --project=julia scripts/run_10steps_fullsize_test.jl
```

#### 2. キャッシュ効率化

- データアクセスパターンの最適化
- メモリアライメントの確認
- NUMA対応（必要に応じて）

#### 3. ベンチマーク計画

| シナリオ | スレッド数 | 期待実行時間 | 期待改善 |
|---------|-----------|-------------|---------|
| ベースライン | 1 | 779秒 | - |
| Phase 2-A | 4 | 545-585秒 | 25-30% |
| Phase 2-B | 8 | 390-470秒 | 40-50% |

---

## 📊 成功基準

### 必須条件（Phase 2完了）

1. **性能改善**: 30%以上の高速化 ✅
2. **精度維持**: RMS/Max残差の変化なし ✅
3. **テスト通過**: 505項目全通過 ✅
4. **スケーラビリティ**: 4-8スレッドで効率的に並列化 ✅

### 目標達成条件（プロジェクト完了）

1. **総合改善**: 50-60%達成 🎯
2. **実行時間**: 430-536秒（7-9分） 🎯
3. **論文品質**: 性能評価レポート完成 🎯

---

## 🔍 既存の並列化実装の確認

### 現在の並列化状況

```julia
# julia/src/solvers/CommonSolver.jl
function CalcAX!(...)
    backend = get_backend(par)  # "sequential" or "thread"
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        # 計算処理
    end
end
```

**確認項目**:
1. `par="thread"` での実行確認
2. スレッド数による性能スケーリング測定
3. ボトルネック箇所の特定（プロファイリング）

---

## 📁 Phase 2で作成するファイル

### 実装ファイル（予定）

- `julia/src/solvers/ParallelConfig.jl` - 並列化設定管理
- `julia/scripts/benchmark_parallel.jl` - 並列化ベンチマーク

### レポートファイル（予定）

- `shared/results/performance_phase2_parallel.md` - Phase 2総合レポート
- `docs/parallel_implementation_guide.md` - 並列化実装ガイド
- `docs/parallel_performance_analysis.md` - 性能分析

---

## 🎯 Phase 2開始時のチェックリスト

### 事前確認

- [ ] Phase 1完了状態の確認（779.48秒、-27.34%）
- [ ] テスト全通過の確認（505項目）
- [ ] git状態の確認（working tree clean）
- [ ] ブランチ作成（tuning8 または phase2）

### 環境確認

- [ ] Julia 1.12.0インストール確認
- [ ] FLoops.jl, ThreadsX.jlインストール確認
- [ ] スレッド数設定確認（JULIA_NUM_THREADS）
- [ ] CPU情報確認（Apple Silicon M-series）

### ベースライン測定

- [ ] シングルスレッド（par="sequential"）ベンチマーク
- [ ] マルチスレッド（par="thread"）現状ベンチマーク
- [ ] プロファイリング実行（ボトルネック特定）

---

## 📖 参考ドキュメント

### Phase 1関連

- `docs/workflows/SESSION_STATE.md` - 最新のセッション状態
- `shared/results/performance_phase1e_tuned.md` - Phase 1-E調整版レポート
- `docs/performance_improvement_proposals.md` - 全Phase計画

### 並列化関連

- FLoops.jl公式ドキュメント: https://juliafolds.github.io/FLoops.jl/
- ThreadsX.jl公式ドキュメント: https://github.com/tkf/ThreadsX.jl
- Julia並列計算ガイド: https://docs.julialang.org/en/v1/manual/multi-threading/

---

## 🚀 Phase 2開始コマンド

### ブランチ作成

```bash
# tuning7から新ブランチ作成
git checkout -b tuning8
# または
git checkout -b phase2
```

### ベンチマーク実行（ベースライン）

```bash
# シングルスレッド（現状確認）
JULIA_NUM_THREADS=1 julia --project=julia scripts/run_10steps_fullsize_test.jl
```

---

**Phase 2を開始する準備が整いました！** 🚀

次のセッションで「Phase 2開始」と伝えてください。

---

**作成者**: Claude Code
**作成日時**: 2025年10月14日 02:00
