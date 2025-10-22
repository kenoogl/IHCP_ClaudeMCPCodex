# 次のセッション開始ガイド

**作成日時**: 2025年10月22日
**セッション状態**: 並列化作業ブランチ作成完了、実装準備完了
**現在のブランチ**: parallelization
**ベースブランチ**: main (b8bb0e6)

---

## 🎯 現在の状態

### ✅ 完了した作業

1. **Python版随伴問題のバグ調査完了**
   - 熱流束が約10^11倍小さい異常を特定
   - 原因候補: RHS符号の不一致（Julia版は符号反転あり、Python版はなし）
   - 詳細: `docs/reports/python_adjoint_bug_investigation.md`
   - **決定**: 問題は棚上げし、並列化を優先

2. **Julia版並列化実装計画策定完了**
   - FLoops/ThreadsXが既に導入済み（48箇所で使用）
   - 現状: 1スレッドで動作、デフォルトpar="sequential"
   - システム: 8コア利用可能
   - 計画書: `docs/plans/parallelization_implementation_plan.md`

3. **ブランチ操作完了**
   - ✅ sliding-window-validation を main にマージ完了
   - ✅ main をリモートにプッシュ完了
   - ✅ 新しい parallelization ブランチ作成・プッシュ完了

---

## 🚀 次に実施すべき作業

### Phase 1: 全スクリプトの並列化対応

#### タスク1-1: スレッド数表示の追加

**対象ファイル**: `julia/scripts/` 内の全スクリプト

主要スクリプト:
- `run_sliding_window.jl`
- `run_10steps_fullsize_test.jl`
- `check_type_stability.jl`
- `test_dhcp_convection.jl`
- その他全て

**追加コード** (各スクリプトのmain処理開始前):
```julia
# 並列化情報の表示
println("=" ^ 80)
println("Julia parallel execution info:")
println("  Available threads: $(Threads.nthreads())")
println("  Parallelization mode: $(par)")
println("=" ^ 80)
println()
```

#### タスク1-2: --par オプションの追加

**対象**: `run_sliding_window.jl`, `run_10steps_fullsize_test.jl`

**引数パース追加**:
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

**ヘルプメッセージ更新**:
```julia
--par MODE             Parallelization mode (sequential | thread) [default: sequential]
```

#### タスク1-3: SlidingWindowSolver.jl の修正

`julia/src/solvers/SlidingWindowSolver.jl` の `sliding_window_CGM_q_saving_mf!()` 関数:
- `par` パラメータを `solve_cgm!()` に正しく伝播させる
- 現在のデフォルト値やハードコードを確認・修正

---

### Phase 2: デフォルト設定の変更

各ソルバーのデフォルトを `par="thread"` に変更:

1. **julia/src/solvers/CGMSolver.jl** line 260
   ```julia
   function solve_cgm!(...; par::String="thread", ...)
   ```

2. **julia/src/solvers/DHCPSolver.jl** line 205
   ```julia
   function solve_dhcp!(...; par::String="thread", ...)
   ```

3. **julia/src/solvers/AdjointSolver.jl** line 227
   ```julia
   function solve_adjoint_mf!(...; par::String="thread", ...)
   ```

4. **julia/src/solvers/SensitivitySolver.jl** 該当箇所
   ```julia
   function solve_sensitivity_mf!(...; par::String="thread", ...)
   ```

5. **julia/src/solvers/SlidingWindowSolver.jl** 該当箇所
   ```julia
   function sliding_window_CGM_q_saving_mf!(...; par::String="thread", ...)
   ```

---

### Phase 1&2完了後: テスト実施

#### テスト1: 動作確認

```bash
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par thread
```

**確認項目**:
- スレッド数表示が正しく出力される
- `--par thread` オプションが正しく動作する
- エラーなく実行完了する

#### テスト2: 数値精度検証

**目的**: 並列化による数値精度の劣化がないことを確認

**実行**:
```bash
# 1スレッド（基準）
JULIA_NUM_THREADS=1 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par sequential

# 8スレッド（並列化）
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par thread
```

**成功基準**:
- 熱流束・温度場の相対誤差 < 1e-10（機械精度レベル）

#### テスト3: 性能ベンチマーク

**測定設定**:

| 問題サイズ | nt  | window | overlap | CGM反復 |
|-----------|-----|--------|---------|---------|
| 小規模     | 10  | 5      | 2       | 3       |
| 中規模     | 100 | 30     | 10      | 3       |
| フルサイズ | 300 | 71     | 17      | 3       |

**実行例**:
```bash
# 各スレッド数で実行
for threads in 1 2 4 8; do
  echo "=== Testing with $threads threads ==="
  JULIA_NUM_THREADS=$threads julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt 10 --window 5 --overlap 2 --cgm-iter 3 --par thread
done
```

**成功基準**:

| スレッド数 | 目標スピードアップ | 並列化効率 |
|-----------|------------------|-----------|
| 1         | 1.00x (基準)      | 100%      |
| 2         | 1.80x 以上       | 90% 以上  |
| 4         | 3.40x 以上       | 85% 以上  |
| 8         | 5.60x 以上       | 70% 以上  |

---

## 📂 重要なファイル

### 計画書・レポート
- `docs/plans/parallelization_implementation_plan.md` - 並列化実装計画（詳細版）
- `docs/reports/python_adjoint_bug_investigation.md` - Python版バグ調査報告

### 修正対象ソースコード

#### スクリプト
- `julia/scripts/run_sliding_window.jl`
- `julia/scripts/run_10steps_fullsize_test.jl`
- `julia/scripts/check_type_stability.jl`
- `julia/scripts/test_dhcp_convection.jl`

#### ソルバー
- `julia/src/solvers/CGMSolver.jl`
- `julia/src/solvers/DHCPSolver.jl`
- `julia/src/solvers/AdjointSolver.jl`
- `julia/src/solvers/SensitivitySolver.jl`
- `julia/src/solvers/SlidingWindowSolver.jl`

---

## 🔧 環境情報

- **システム**: 8コアCPU
- **現在のスレッド数**: 1（デフォルト）
- **Julia**: FLoops v0.2.2, ThreadsX v0.1.12 導入済み
- **ブランチ**: parallelization
- **ベース**: main (b8bb0e6)

---

## 📝 実装チェックリスト

### Phase 1: 全スクリプトの並列化対応

- [ ] **タスク1-1**: スレッド数表示の追加
  - [ ] run_sliding_window.jl
  - [ ] run_10steps_fullsize_test.jl
  - [ ] check_type_stability.jl
  - [ ] test_dhcp_convection.jl
  - [ ] その他全スクリプト

- [ ] **タスク1-2**: --par オプションの追加
  - [ ] run_sliding_window.jl
  - [ ] run_10steps_fullsize_test.jl
  - [ ] ヘルプメッセージ更新

- [ ] **タスク1-3**: SlidingWindowSolver.jl の修正
  - [ ] par パラメータの伝播確認・修正

### Phase 2: デフォルト設定の変更

- [ ] CGMSolver.jl line 260
- [ ] DHCPSolver.jl line 205
- [ ] AdjointSolver.jl line 227
- [ ] SensitivitySolver.jl 該当箇所
- [ ] SlidingWindowSolver.jl 該当箇所

### テスト（Phase 1&2完了後）

- [ ] **テスト1**: 動作確認
- [ ] **テスト2**: 数値精度検証（相対誤差 < 1e-10）
- [ ] **テスト3**: 性能ベンチマーク（8コアで5-7倍目標）
- [ ] 結果レポート作成

---

## ⚠️ 重要な注意事項

### スレッド数は明示的に指定（必須）

```bash
# 正しい
JULIA_NUM_THREADS=8 julia --project=julia julia/scripts/run_sliding_window.jl --par thread

# 間違い（1スレッドになる）
julia --project=julia julia/scripts/run_sliding_window.jl --par thread
```

### テスト順序

1. Phase 1実装 → 動作確認
2. Phase 2実装 → 動作確認
3. Phase 1&2完了後 → 全テスト実行（精度・性能）

---

## 📊 期待される成果

- **Phase 1完了**: 全スクリプトで並列化モード選択可能
- **Phase 2完了**: デフォルトで並列化が有効
- **テスト完了**: 8コアで5-7倍高速化達成

---

## 🔗 次のセッション開始手順

### 1. ブランチ確認

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git branch
# parallelization ブランチにいることを確認
```

### 2. このファイルと計画書を確認

```bash
cat TODO_NEXT_SESSION.md
cat docs/plans/parallelization_implementation_plan.md
```

### 3. Phase 1実装開始

計画書の手順に従って、全スクリプトにスレッド数表示と --par オプションを追加していく。

---

**次のセッション開始時**: 必ずこのファイルを読んでから作業を開始してください。
