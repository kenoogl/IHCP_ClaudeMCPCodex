# Julia ディレクトリ整理計画

**作成日**: 2025年10月11日
**対象**: julia/ディレクトリ配下の整理

## 整理方針

### 保持するもの
- src/: ソースコード全体
- scripts/: 現在使用中のスクリプト（2ファイル）
- test/: 主要テストファイル（runtests.jl, test_*.jl）
- config/: 設定ファイル
- data/: 参照データ（テスト用）
- plots/: プロット画像（必要に応じて）
- Project.toml, Manifest.toml: Julia環境設定
- README.md: ドキュメント

### アーカイブするもの（archive/julia/へ移動）

#### 1. 古いベンチマークファイル（14ファイル）
```
benchmark_baseline_old_layout.jl
benchmark_comparison.md
benchmark_large_comparison.md
benchmark_large_new.txt
benchmark_large_old.jl
benchmark_large_old.txt
benchmark_large.jl
benchmark_layouts.jl
benchmark_new.txt
benchmark_old_layout.jl
benchmark_old.txt
benchmark_phase22.jl
compare_layouts.jl
verify_layout_new.jl
verify_layout_old.jl
```

**理由**: Phase A-D完了、これらのベンチマークは古い

#### 2. 古い分析・計画ドキュメント（6ファイル）
```
grid_structure_comparison.md
guard_cell_analysis.md
guard_cell_implementation_summary.md
matrix_free_analysis.md
phase_b_implementation_plan.md
tuning1_vs_tuning2_analysis.md
```

**理由**: Phase A-D完了、分析は完了済み

#### 3. プロファイル関連（古い、3ファイル）
```
profile_new_mf.jl
profile_new_mf.txt
debug_performance.jl
```

**理由**: 古いプロファイリング結果

#### 4. examples/ディレクトリ（旧版、4ファイル）
```
examples/load_test_data.jl
examples/run_exact_match.jl
examples/run_test_1min.jl
examples/run_test_small.jl
```

**理由**: scripts/に統合済み、またはアーカイブ済みのスクリプトの旧版

#### 5. test/内の旧ファイル（9ファイル）
```
test/compare_medium_problem.jl
test/compare_phase5_detailed.jl
test/compare_sliding_window.jl
test/debug_sliding_window.jl
test/generate_reference_4c87fcf.jl
test/medium_problem_4c87fcf.txt
test/medium_problem_current.txt
test/run_dhcp_test_simple.jl
test/test_grid_transform.jl
test/test_rhs_core.jl
```

**理由**: 古いテスト・デバッグスクリプト、または統合済み

### 削除対象
- なし（すべてアーカイブへ移動）

---

## 整理後のjuliaディレクトリ構造

```
julia/
├── config/              # 設定ファイル
├── data/                # 参照データ（テスト用）
├── plots/               # プロット画像（必要に応じて）
├── scripts/             # 現在使用中のスクリプト
│   ├── generate_test_data.jl
│   └── run_10steps_fullsize_test.jl
├── src/                 # ソースコード
│   ├── IHCP_CGM.jl
│   ├── ThermalProperties.jl
│   ├── DataLoaders.jl
│   ├── main.jl
│   ├── solvers/
│   └── utils/
├── test/                # テストファイル
│   ├── runtests.jl
│   ├── test_thermal_properties.jl
│   ├── test_data_loaders.jl
│   ├── test_cgm_solver.jl
│   ├── test_sliding_window.jl
│   ├── test_validators.jl
│   ├── data/
│   └── deprecated/
├── Project.toml
├── Manifest.toml
└── README.md
```

---

## 実行コマンド

### Step 1: アーカイブディレクトリ作成
```bash
mkdir -p archive/julia/{benchmarks,analysis,profiling,examples,tests}
```

### Step 2: ベンチマークファイルの移動
```bash
mv julia/benchmark_*.{jl,md,txt} archive/julia/benchmarks/ 2>/dev/null || true
mv julia/compare_layouts.jl archive/julia/benchmarks/
mv julia/verify_layout_*.jl archive/julia/benchmarks/
```

### Step 3: 分析ドキュメントの移動
```bash
mv julia/grid_structure_comparison.md archive/julia/analysis/
mv julia/guard_cell_*.md archive/julia/analysis/
mv julia/matrix_free_analysis.md archive/julia/analysis/
mv julia/phase_b_implementation_plan.md archive/julia/analysis/
mv julia/tuning1_vs_tuning2_analysis.md archive/julia/analysis/
```

### Step 4: プロファイル関連の移動
```bash
mv julia/profile_*.{jl,txt} archive/julia/profiling/ 2>/dev/null || true
mv julia/debug_performance.jl archive/julia/profiling/
```

### Step 5: examples/ディレクトリの移動
```bash
mv julia/examples archive/julia/examples_old
```

### Step 6: test/内の旧ファイルの移動
```bash
mv julia/test/compare_*.jl archive/julia/tests/ 2>/dev/null || true
mv julia/test/debug_*.jl archive/julia/tests/ 2>/dev/null || true
mv julia/test/generate_reference_*.jl archive/julia/tests/ 2>/dev/null || true
mv julia/test/*_problem_*.txt archive/julia/tests/ 2>/dev/null || true
mv julia/test/run_dhcp_test_simple.jl archive/julia/tests/
mv julia/test/test_grid_transform.jl archive/julia/tests/
mv julia/test/test_rhs_core.jl archive/julia/tests/
```

---

## 統計

### ファイル数
- アーカイブ対象: 約36ファイル
- 保持: 主要ファイル（src/, scripts/, test/主要テスト、設定）

### 期待される効果
- ディレクトリ構造の明確化
- 現在使用中のファイルのみ残る
- 古いベンチマーク・分析ファイルはアーカイブで参照可能

---

**実施タイミング**: ユーザー確認後
