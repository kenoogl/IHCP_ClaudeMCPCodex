# Julia版並列化実装計画（修正版）

**作成日**: 2025年10月22日
**システム環境**: 8コアCPU、現在1スレッドで動作
**目的**: 既存FLoopsインフラを活用してマルチスレッド並列化を有効化し、性能を向上させる

---

## 1. 現状分析

### 1.1 既存インフラ（✅ 導入済み）

- **FLoops.jl** v0.2.2: データ並列ループフレームワーク
- **ThreadsX.jl** v0.1.12: スレッド並列リダクション
- **@floop マクロ**: 48箇所で既に使用済み
- **par切り替え機構**: `get_backend(par)` で "sequential" / "thread" を切り替え可能

### 1.2 並列化済みの箇所

#### 各ソルバー（既に@floopで実装済み）
- **CommonSolver.jl**: `PBiCGSTAB!()`, `PCG!()`, `CalcRK!()`, `CalcAX!()`, 前処理器
- **DHCPSolver.jl**: RHS構築（`calRHS!`）- 3箇所
- **AdjointSolver.jl**: RHS構築（`calRHS!`）- 3箇所
- **SensitivitySolver.jl**: RHS構築（`calRHS!`）- 3箇所
- **RHSCore.jl**: 境界条件適用（X±, Y±, Z±の6面）

### 1.3 問題点（🔴 未活用）

- **Julia起動時**: 1スレッドで実行中
- **デフォルト設定**: 全ソルバーで `par="sequential"`
  - `solve_cgm!()`: CGMSolver.jl line 260
  - `solve_dhcp!()`: DHCPSolver.jl line 205
  - `solve_adjoint_mf!()`: AdjointSolver.jl line 227

---

## 2. 並列化実装計画（2段階アプローチ）

### Phase 1: 全スクリプトの並列化対応修正

**目標**: 全てのJuliaスクリプトでスレッド数を明示的に表示し、並列化を実行時に選択可能にする

#### 実装内容

**タスク1-1**: スレッド数表示の追加（全スクリプト）

対象ファイル:
```
julia/scripts/
├── run_sliding_window.jl
├── run_10steps_fullsize_test.jl
├── check_type_stability.jl
├── test_dhcp_convection.jl
└── その他全てのスクリプト
```

追加コード（各スクリプトの冒頭、main処理開始前）:
```julia
# 並列化情報の表示
println("=" ^ 80)
println("Julia parallel execution info:")
println("  Available threads: $(Threads.nthreads())")
println("  Parallelization mode: $(par)")
println("=" ^ 80)
println()
```

**タスク1-2**: コマンドライン引数の追加（主要スクリプト）

対象: `run_sliding_window.jl`, `run_10steps_fullsize_test.jl`

追加する引数:
```julia
--par [sequential|thread]    並列化モード選択 [default: sequential]
```

引数パース例:
```julia
elseif arg == "--par"
  i += 1
  i > length(ARGS) && error("--par requires an argument")
  par_str = lowercase(ARGS[i])
  if par_str == "sequential"
    par = "sequential"
  elseif par_str == "thread"
    par = "thread"
  else
    error("Unknown par mode: $(ARGS[i]). Use 'sequential' or 'thread'")
  end
```

ヘルプメッセージ更新:
```julia
--par MODE             Parallelization mode (sequential | thread) [default: sequential]
```

**タスク1-3**: SlidingWindowSolver.jl の修正

`sliding_window_CGM_q_saving_mf!()` 関数:
- `par` パラメータを `solve_cgm!()` に伝播
- 現在 `par="sequential"` がハードコードされている可能性を確認・修正

---

### Phase 2: デフォルト設定の変更

**目標**: 明示的に指定しない場合、デフォルトで並列化を有効にする

#### 実装内容

各ソルバーのデフォルト引数を `par="thread"` に変更:

**変更対象ファイルと行番号**:

1. **julia/src/solvers/CGMSolver.jl**: line 260
   ```julia
   # 修正前
   function solve_cgm!(...; par::String="sequential", ...)

   # 修正後
   function solve_cgm!(...; par::String="thread", ...)
   ```

2. **julia/src/solvers/DHCPSolver.jl**: line 205
   ```julia
   # 修正前
   function solve_dhcp!(...; par::String="sequential", ...)

   # 修正後
   function solve_dhcp!(...; par::String="thread", ...)
   ```

3. **julia/src/solvers/AdjointSolver.jl**: line 227
   ```julia
   # 修正前
   function solve_adjoint_mf!(...; par::String="sequential", ...)

   # 修正後
   function solve_adjoint_mf!(...; par::String="thread", ...)
   ```

4. **julia/src/solvers/SensitivitySolver.jl**: 該当する関数
   ```julia
   # 同様に修正
   ```

5. **julia/src/solvers/SlidingWindowSolver.jl**: 該当する関数
   ```julia
   # 同様に修正
   ```

#### 下位互換性の確保

既存テストコードでは明示的に `par="sequential"` を指定することで、シングルスレッド動作を維持:

```julia
# テストコード例
solve_dhcp!(...; par="sequential")  # 明示的にシングルスレッド指定
```

---

## 3. テスト計画（Phase 1 & 2 完了後）

### 3.1 実行方法

#### スレッド数の明示的指定

```bash
# 8スレッドで実行
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par thread

# 4スレッドで実行
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par thread

# 1スレッド（比較用）
JULIA_NUM_THREADS=1 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par sequential
```

### 3.2 テスト項目

#### テスト1: 動作確認（必須）

**目的**: 並列化実装が正しく動作することを確認

**手順**:
1. Phase 1完了後: `--par thread` オプションが正しく動作するか確認
2. Phase 2完了後: デフォルトで並列化が有効になることを確認
3. スレッド数表示が正しく出力されることを確認

**期待結果**:
```
================================================================================
Julia parallel execution info:
  Available threads: 8
  Parallelization mode: thread
================================================================================
```

#### テスト2: 数値精度検証（必須）

**目的**: 並列化による数値精度の劣化がないことを確認

**手順**:
1. 小規模問題（nt=10, window=5, overlap=2, CGM=3）を実行
2. 1スレッド vs 8スレッドで結果を比較
3. 熱流束・温度場の相対誤差を計算

**成功基準**:
- 相対誤差 < 1e-10（機械精度レベル）

#### テスト3: 性能ベンチマーク（推奨）

**目的**: 並列化の効果を定量的に測定

**測定設定**:

| 問題サイズ | nt | window | overlap | CGM反復 |
|-----------|-----|--------|---------|---------|
| 小規模     | 10  | 5      | 2       | 3       |
| 中規模     | 100 | 30     | 10      | 3       |
| フルサイズ | 300 | 71     | 17      | 3       |

**測定項目**:
1. 全体実行時間
2. ソルバー別累積時間（DHCP / Adjoint / Sensitivity）
3. スケーラビリティ（1 vs 2 vs 4 vs 8スレッド）

**成功基準**:

| スレッド数 | 目標スピードアップ | 並列化効率 |
|-----------|------------------|-----------|
| 1         | 1.00x (基準)      | 100%      |
| 2         | 1.80x 以上       | 90% 以上  |
| 4         | 3.40x 以上       | 85% 以上  |
| 8         | 5.60x 以上       | 70% 以上  |

---

## 4. 実装チェックリスト

### Phase 1: 全スクリプトの並列化対応

- [ ] **タスク1-1**: スレッド数表示の追加
  - [ ] `run_sliding_window.jl`
  - [ ] `run_10steps_fullsize_test.jl`
  - [ ] `check_type_stability.jl`
  - [ ] `test_dhcp_convection.jl`
  - [ ] その他全スクリプト

- [ ] **タスク1-2**: コマンドライン引数の追加
  - [ ] `run_sliding_window.jl` に `--par` オプション追加
  - [ ] `run_10steps_fullsize_test.jl` に `--par` オプション追加
  - [ ] ヘルプメッセージ更新

- [ ] **タスク1-3**: SlidingWindowSolver.jl の修正
  - [ ] `par` パラメータの伝播確認・修正

### Phase 2: デフォルト設定の変更

- [ ] `CGMSolver.jl` line 260: `par="thread"`
- [ ] `DHCPSolver.jl` line 205: `par="thread"`
- [ ] `AdjointSolver.jl` line 227: `par="thread"`
- [ ] `SensitivitySolver.jl`: 該当箇所修正
- [ ] `SlidingWindowSolver.jl`: 該当箇所修正

### テスト（Phase 1 & 2 完了後）

- [ ] **テスト1**: 動作確認
  - [ ] `--par thread` オプション動作確認
  - [ ] スレッド数表示確認
  - [ ] デフォルト動作確認（Phase 2後）

- [ ] **テスト2**: 数値精度検証
  - [ ] 小規模問題実行（1スレッド）
  - [ ] 小規模問題実行（8スレッド）
  - [ ] 結果比較・相対誤差計算
  - [ ] 精度基準クリア確認（< 1e-10）

- [ ] **テスト3**: 性能ベンチマーク
  - [ ] 小規模問題（1, 2, 4, 8スレッド）
  - [ ] スピードアップ計算
  - [ ] 並列化効率計算
  - [ ] 結果レポート作成

---

## 5. リスク管理

### リスク1: 数値精度の劣化

**対策**: テスト2で厳密な精度検証を実施、相対誤差 < 1e-10 を維持

### リスク2: 期待効果が得られない

**対策**: テスト3で効果を測定、不十分な場合は原因分析・追加最適化を検討

### リスク3: メモリ競合・デッドロック

**対策**: FLoopsは安全（データ並列のみ）、各スレッドは独立領域にアクセス

---

## 6. 期待される成果

### Phase 1完了時

✅ **機能**: 全スクリプトで並列化モードを選択可能
✅ **可視性**: スレッド数・並列化モードが明示的に表示
✅ **柔軟性**: 実行時に並列化ON/OFFを切り替え可能

### Phase 2完了時

✅ **利便性**: デフォルトで並列化が有効（明示的指定不要）
✅ **下位互換**: 既存テストは `par="sequential"` で動作維持

### テスト完了時

✅ **性能**: 8コアで5-7倍高速化達成（目標）
✅ **精度**: 数値精度を維持（相対誤差 < 1e-10）
✅ **信頼性**: 全テスト合格

### 数値例（推定）

現状（1スレッド）でフルサイズ計算が1000秒の場合:
- **並列化後（8スレッド）**: 約150-200秒（**800秒短縮**）

---

## 7. 実装の注意事項

### スレッド数の明示的指定（必須）

```bash
# 必ずJULIA_NUM_THREADSを明示的に指定
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_sliding_window.jl --par thread

# 環境変数なしで実行すると1スレッドになる
julia --project=julia julia/scripts/run_sliding_window.jl --par thread  # ← これは1スレッド！
```

### テスト時の推奨手順

1. **Phase 1実装 → 動作確認**
2. **Phase 2実装 → 動作確認**
3. **Phase 1 & 2完了後 → 全テスト実行**

この順序により、問題が発生した場合の切り分けが容易になります。

---

**計画書バージョン**: 2.0（修正版）
**最終更新**: 2025年10月22日
