# ディレクトリ整理計画

**作成日**: 2025年10月11日
**状況**: Python→Julia移植完了、性能改善フェーズへ移行

## 整理方針

- **保存**: アーカイブディレクトリへ移動（将来参照用）
- **保持**: 現在の場所に維持（現在使用中 or スライディングウィンドウ検証用）
- **削除**: 完全に削除（不要な中間ファイル）

---

## Python Scripts

### 保持すべきもの（現在の場所に維持）

#### オリジナルコード
- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py` - Python版オリジナル（参照用）

#### スライディングウィンドウ検証用（未検証）
- `python/benchmark/run_test_1min_simple.py` - スライディングウィンドウテスト
- `python/benchmark/run_test_1min.py` - スライディングウィンドウテスト
- `python/generate_phase5_reference.py` - Phase 5参照データ生成

#### 現在使用中の検証スクリプト
- `python/validation/compare_python_julia_10steps_fullsize.py` - 10ステップ比較（現在使用中）
- `python/validation/run_10steps_fullsize_test.py` - 10ステップテスト（現在使用中）

### アーカイブすべきもの（archive/python/へ移動）

#### 完了済みの検証スクリプト
- `python/validation/compare_exact_match.py` - exact match検証完了
- `python/validation/run_exact_match.py` - exact match検証完了
- `python/validation/extract_params.py` - 一時的なユーティリティ
- `python/validation/check_previous_params.py` - 一時的なチェック
- `python/validation/compare_python_julia_10steps_fullsize_before_phaseA.py` - Phase A前の旧版

#### 古いベンチマーク
- `python/benchmark/analyze_python_results.py` - 古い分析
- `python/benchmark/estimate_timesteps_simple.py` - 見積もり用
- `python/benchmark/estimate_timesteps.py` - 見積もり用

#### 古いテスト
- `python/validation/compare_python_julia_100steps.py` - 100ステップ検証完了
- `python/validation/run_100steps_test.py` - 100ステップ検証完了

---

## Julia Scripts

### 保持すべきもの（現在の場所に維持）

#### 現在使用中
- `julia/scripts/run_10steps_fullsize_test.jl` - 性能測定用（現在使用中）

#### スライディングウィンドウ検証用（未検証）
- スライディングウィンドウ専用のスクリプトは未作成のため、現時点ではなし

#### データ生成用
- `julia/scripts/generate_test_data.jl` - テストデータ生成（必要に応じて使用）

### アーカイブすべきもの（archive/julia/scripts/へ移動）

#### Phase A前後の旧版
- `julia/scripts/run_10steps_fullsize_test_before_phaseA.jl` - Phase A前の旧版
- `julia/scripts/run_10steps_fullsize_test_phaseA.jl` - Phase A版（もう不要）

#### 古いテスト
- `julia/scripts/run_100steps_cgm1_test.jl` - 100ステップCGM1テスト完了
- `julia/scripts/run_100steps_test.jl` - 100ステップテスト完了
- `julia/scripts/run_10steps_mediumsize_test.jl` - 中間サイズテスト（不要）
- `julia/scripts/test_minimal.jl` - 最小テスト（デバッグ用、不要）

#### ベンチマーク
- `julia/scripts/benchmark_fill.jl` - fill!ベンチマーク（調査完了）

---

## Results Files

### 保持すべきもの（現在の場所に維持）

#### 最新の性能測定結果
- `shared/results/performance_22fde2d.md` - 性能ベースライン（重要）
- `shared/results/run_10steps_fullsize_22fde2d_detailed.log` - 最新ログ（重要）

#### 最新の比較結果
- `shared/results/comparison_10steps_fullsize_error_hist.png` - 最新比較プロット
- `shared/results/comparison_10steps_fullsize_heat_flux.png` - 最新比較プロット
- `shared/results/comparison_10steps_fullsize_temperature.png` - 最新比較プロット
- `shared/results/julia_10steps_fullsize.npz` - 最新Julia結果
- `shared/results/python_10steps_fullsize.npz` - 最新Python結果

#### 最終検証結果
- `shared/results/validation/` - 最終検証結果全体（保持）

#### ドキュメント
- `shared/results/README.md` - 結果ディレクトリの説明

#### スライディングウィンドウ検証用（未検証）
- `shared/results/python_300steps.npz` - 300ステップ結果（308MB、スライディングウィンドウ検証用）
- `shared/results/python_300steps_q.npy` - 熱流束データ
- `shared/results/python_300steps_T.npy` - 温度データ
- `shared/results/validation/python_test_1min.npz` - 1分テスト結果（スライディングウィンドウ検証用）

### アーカイブすべきもの（archive/results/へ移動）

#### 古いログ
- `shared/results/run_10steps_fullsize_4cef5c5_detailed.log` - 古いログ
- `shared/results/run_10steps_fullsize_4cef5c5.log` - 古いログ

#### 古い100ステップ結果
- `shared/results/comparison_100steps_error_hist.png` - 100ステップ比較プロット
- `shared/results/comparison_100steps_heat_flux.png` - 100ステップ比較プロット
- `shared/results/comparison_100steps_temperature.png` - 100ステップ比較プロット
- `shared/results/julia_100steps_q.npy` - Julia 100ステップ結果
- `shared/results/julia_100steps_T.npy` - Julia 100ステップ結果
- `shared/results/julia_100steps.npz` - Julia 100ステップ結果（135MB）
- `shared/results/python_100steps_cgm1_q.npy` - Python 100ステップCGM1結果
- `shared/results/python_100steps_cgm1_T.npy` - Python 100ステップCGM1結果
- `shared/results/python_100steps_cgm1.npz` - Python 100ステップCGM1結果（99MB）
- `shared/results/python_100steps_q.npy` - Python 100ステップ結果
- `shared/results/python_100steps_T.npy` - Python 100ステップ結果
- `shared/results/python_100steps.npz` - Python 100ステップ結果（95MB）

### 削除すべきもの

#### 完了済みタスク
- `shared/results/NEXT_SESSION_TASK.md` - 完了済みタスク（不要）

---

## アーカイブディレクトリ構造

```
archive/
├── python/
│   ├── benchmark/          # 古いベンチマークスクリプト
│   └── validation/         # 完了済み検証スクリプト
├── julia/
│   └── scripts/           # 旧版スクリプト
└── results/
    ├── logs/              # 古いログファイル
    ├── 100steps/          # 100ステップ検証結果
    └── plots/             # 古い比較プロット
```

---

## 実行計画

### Step 1: アーカイブディレクトリ作成
```bash
mkdir -p archive/{python/{benchmark,validation},julia/scripts,results/{logs,100steps,plots}}
```

### Step 2: ファイル移動（Python）
```bash
# 完了済み検証スクリプト
mv python/validation/compare_exact_match.py archive/python/validation/
mv python/validation/run_exact_match.py archive/python/validation/
mv python/validation/extract_params.py archive/python/validation/
mv python/validation/check_previous_params.py archive/python/validation/
mv python/validation/compare_python_julia_10steps_fullsize_before_phaseA.py archive/python/validation/
mv python/validation/compare_python_julia_100steps.py archive/python/validation/
mv python/validation/run_100steps_test.py archive/python/validation/

# 古いベンチマーク
mv python/benchmark/analyze_python_results.py archive/python/benchmark/
mv python/benchmark/estimate_timesteps_simple.py archive/python/benchmark/
mv python/benchmark/estimate_timesteps.py archive/python/benchmark/
```

### Step 3: ファイル移動（Julia）
```bash
# 旧版スクリプト
mv julia/scripts/run_10steps_fullsize_test_before_phaseA.jl archive/julia/scripts/
mv julia/scripts/run_10steps_fullsize_test_phaseA.jl archive/julia/scripts/
mv julia/scripts/run_100steps_cgm1_test.jl archive/julia/scripts/
mv julia/scripts/run_100steps_test.jl archive/julia/scripts/
mv julia/scripts/run_10steps_mediumsize_test.jl archive/julia/scripts/
mv julia/scripts/test_minimal.jl archive/julia/scripts/
mv julia/scripts/benchmark_fill.jl archive/julia/scripts/
```

### Step 4: ファイル移動（Results）
```bash
# 古いログ
mv shared/results/run_10steps_fullsize_4cef5c5_detailed.log archive/results/logs/
mv shared/results/run_10steps_fullsize_4cef5c5.log archive/results/logs/

# 100ステップ結果
mv shared/results/comparison_100steps_*.png archive/results/plots/
mv shared/results/*100steps*.npy archive/results/100steps/
mv shared/results/*100steps*.npz archive/results/100steps/
```

### Step 5: ファイル削除
```bash
# 完了済みタスク
rm shared/results/NEXT_SESSION_TASK.md
```

### Step 6: gitignoreへのアーカイブディレクトリ追加
```bash
echo "" >> .gitignore
echo "# Archive directory" >> .gitignore
echo "archive/" >> .gitignore
```

---

## 整理後のディレクトリ構成

### Python Scripts（残存）
```
python/
├── benchmark/
│   ├── run_test_1min_simple.py      # スライディングウィンドウ検証用
│   └── run_test_1min.py             # スライディングウィンドウ検証用
├── original/
│   └── IHCP_CGM_Sliding_Window_Calculation_ver2.py  # オリジナル
├── validation/
│   ├── compare_python_julia_10steps_fullsize.py  # 現在使用中
│   └── run_10steps_fullsize_test.py              # 現在使用中
└── generate_phase5_reference.py     # Phase 5参照データ生成
```

### Julia Scripts（残存）
```
julia/scripts/
├── generate_test_data.jl              # データ生成
└── run_10steps_fullsize_test.jl       # 性能測定用
```

### Results（残存）
```
shared/results/
├── comparison_10steps_fullsize_*.png      # 最新比較プロット（3ファイル）
├── julia_10steps_fullsize.npz            # 最新Julia結果
├── python_10steps_fullsize.npz           # 最新Python結果
├── python_300steps.npz                   # スライディングウィンドウ検証用
├── python_300steps_q.npy                 # スライディングウィンドウ検証用
├── python_300steps_T.npy                 # スライディングウィンドウ検証用
├── performance_22fde2d.md                # 性能ベースライン
├── run_10steps_fullsize_22fde2d_detailed.log  # 最新ログ
├── README.md
└── validation/                           # 最終検証結果（全体保持）
```

---

## ディスク使用量の削減見込み

### 移動されるファイル
- 100ステップ結果: 約329MB (135MB + 99MB + 95MB)
- 古いログ: 約1MB

### 削除されるファイル
- NEXT_SESSION_TASK.md: 約5KB

### 合計削減量: 約330MB（アーカイブへ移動、gitignore追加により）

---

## 注意事項

1. **スライディングウィンドウ検証前**のため、以下は保持：
   - `python_300steps.npz`（308MB）
   - `python/benchmark/run_test_1min*.py`
   - `python/generate_phase5_reference.py`
   - `shared/results/validation/python_test_1min.npz`

2. **アーカイブディレクトリはgitignore**に追加：
   - 移動したファイルはgit管理対象外
   - 必要に応じてローカルで参照可能

3. **削除前の確認**：
   - アーカイブ移動後、動作確認してから削除実行

---

## 実行タイミング

このクリーンアップは以下のタイミングで実行予定：
- ユーザー確認後、即座に実行
