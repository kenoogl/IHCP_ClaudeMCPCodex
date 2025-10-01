# IHCP CGM Solver - Julia実装

逆熱伝導問題（IHCP）を共役勾配法（CGM）で解くJuliaソルバー

## 概要

IRカメラからの温度測定データを使用して、SUS304材の表面熱流束を逆解析します。

- **Phase 1-5**: 完全実装済み（383テスト合格）
- **バージョン**: v0.5.0
- **Julia要件**: 1.10以上

## インストール

```bash
# プロジェクトディレクトリに移動
cd julia

# 依存パッケージのインストール
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## 使用方法

### 基本的な実行

```bash
julia --project=. src/main.jl --config <設定ファイル> --output <出力ファイル> [OPTIONS]
```

### コマンドライン引数

| 引数 | 説明 | 必須 | デフォルト |
|------|------|------|-----------|
| `--config`, `-c` | 設定ファイル（TOML形式） | ✓ | - |
| `--output`, `-o` | 出力ファイル（JLD2形式） | ✓ | - |
| `--mode`, `-m` | 実行モード: full/test/benchmark | - | full |
| `--threads`, `-t` | BLASスレッド数 | - | 自動 |
| `--verbose`, `-v` | 詳細ログ出力 | - | false |
| `--quiet`, `-q` | 最小限のログのみ | - | false |
| `--dry-run` | 設定検証のみ（計算しない） | - | false |
| `--help`, `-h` | ヘルプ表示 | - | - |
| `--version` | バージョン表示 | - | - |

### 実行モード

#### 1. fullモード（デフォルト）

全時間領域のスライディングウィンドウCGM計算を実行

```bash
# テスト用小規模計算
julia --project=. src/main.jl \
  --config config/test.toml \
  --output results/test_results.jld2 \
  --mode full

# 本番計算（詳細ログ付き）
julia --project=. src/main.jl \
  --config config/example.toml \
  --output results/full_results.jld2 \
  --mode full \
  --verbose
```

#### 2. testモード

Phase 1-5の統合テストを実行

```bash
julia --project=. src/main.jl \
  --config config/test.toml \
  --output dummy.jld2 \
  --mode test
```

#### 3. benchmarkモード

性能測定を実行

```bash
julia --project=. src/main.jl \
  --config config/test.toml \
  --output dummy.jld2 \
  --mode benchmark
```

### 設定検証のみ（dry-run）

計算を実行せず、設定ファイルの妥当性のみチェック

```bash
julia --project=. src/main.jl \
  --config config/example.toml \
  --output dummy.jld2 \
  --dry-run
```

## 設定ファイル（TOML形式）

### 基本構造

```toml
[problem]
ni = 10
nj = 10
nk = 20
dx = 0.12e-3
dy = 0.12e-3
dz = [...]  # nk要素の配列
nt = 1000
dt = 1.0e-3

[material]
rho = 7900.0
cp_coeffs = [4.184e-7, -0.0009168, 1.34436, 233.604]
k_coeffs = [9.2486e-6, -0.01356, 8.21814, 7821.93]

[cgm]
max_iter = 50
sigma = 0.1

[sliding_window]
window_size = 100
overlap = 20
q_init_value = 0.0

[input]
T_init_file = "data/T_init.npy"
Y_obs_file = "data/Y_obs.npy"

[output]
save_format = "jld2"
```

詳細な設定例は `config/example.toml` を参照してください。

### 重要な設定項目

#### 1. z方向格子（dz配列）

- **要素数**: nk個（セル数と同じ）
- **内容**: 各層の厚さ [m]
- **注意**: 表面に細かい格子を配置することを推奨

```toml
# 非均等格子の例（表面側5層: 0.05mm、底面側5層: 0.20mm）
dz = [
  0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4, 0.5e-4,  # 表面側
  1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4,  # 中間
  1.5e-4, 1.5e-4, 1.5e-4, 1.5e-4, 1.5e-4,  # 中間
  2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4, 2.0e-4   # 底面側
]
```

#### 2. 入力データ形状

- **初期温度場 (T_init)**: (ni, nj, nk)
- **観測温度 (Y_obs)**: (ni, nj, nt) ※読み込み時に自動変換

サポートファイル形式: `.npy`, `.txt`, `.mat`

## 出力ファイル

### JLD2形式の内容

```julia
using JLD2

# 結果読み込み
data = load("results.jld2")

# 内容
# - q_global: 推定熱流束 (nt-1, ni, nj) [W/m²]
# - windows_info: スライディングウィンドウ詳細情報
# - config: 計算設定
# - metadata: メタデータ（バージョン、日時等）
```

### 結果の読み込み例

```julia
using JLD2
using Statistics

# 結果ファイル読み込み
results = load("results.jld2")

q_global = results["q_global"]
windows_info = results["windows_info"]
config = results["config"]

# 統計情報
println("熱流束形状: ", size(q_global))
println("平均熱流束: ", mean(q_global), " W/m²")
println("最大熱流束: ", maximum(q_global), " W/m²")
println("最小熱流束: ", minimum(q_global), " W/m²")
println("ウィンドウ数: ", length(windows_info))
```

## テストデータ生成

小規模な合成データを生成してテスト

```bash
julia --project=. scripts/generate_test_data.jl
```

生成されるファイル:
- `data/test_T_init.npy`: 初期温度場（3×3×5）
- `data/test_Y_obs.npy`: 観測温度（3×3×50）

## テスト実行

### Phase 1-5統合テスト

```bash
# テストスイート実行
julia --project=. -e 'using Pkg; Pkg.test("IHCP_CGM")'

# または
julia --project=. src/main.jl \
  --config config/test.toml \
  --output dummy.jld2 \
  --mode test
```

### 個別Phaseのテスト

```bash
# Phase 1: 熱物性値計算
julia --project=. test/test_thermal_properties.jl

# Phase 2: DHCP（直接解法）
julia --project=. test/test_dhcp_solver.jl

# Phase 3: Adjoint（随伴解法）
julia --project=. test/test_adjoint_solver.jl

# Phase 4: CGM（共役勾配法）
julia --project=. test/test_cgm_solver.jl

# Phase 5: スライディングウィンドウ
julia --project=. test/test_sliding_window.jl
```

## 実行例

### 例1: 小規模テスト計算

```bash
# 1. テストデータ生成
julia --project=. scripts/generate_test_data.jl

# 2. dry-runで設定確認
julia --project=. src/main.jl \
  --config config/test.toml \
  --output test_output.jld2 \
  --dry-run

# 3. 計算実行
julia --project=. src/main.jl \
  --config config/test.toml \
  --output test_output.jld2 \
  --mode full \
  --verbose

# 出力例:
# 計算時間: 0.47秒
# 熱流束範囲: [-1.10e+05, -1.08e+04] W/m²
```

### 例2: 本番規模計算

```bash
# 大規模計算（1000時間ステップ）
julia --project=. src/main.jl \
  --config config/example.toml \
  --output results/production_results.jld2 \
  --mode full \
  --threads 8
```

## プロジェクト構造

```
julia/
├── src/
│   ├── main.jl                    # メインエントリポイント
│   ├── IHCP_CGM.jl                # メインモジュール
│   ├── ThermalProperties.jl       # Phase 1: 熱物性値計算
│   ├── DataLoaders.jl             # Phase 1: データ読み込み
│   ├── solvers/
│   │   ├── DHCPSolver.jl          # Phase 2: 直接解法
│   │   ├── AdjointSolver.jl       # Phase 3: 随伴解法
│   │   ├── StoppingCriteria.jl    # Phase 4: 停止判定
│   │   ├── CGMSolver.jl           # Phase 4: 共役勾配法
│   │   └── SlidingWindowSolver.jl # Phase 5: スライディングウィンドウ
│   └── utils/
│       ├── config.jl              # 設定ファイル処理
│       └── io.jl                  # 入出力ユーティリティ
├── config/
│   ├── example.toml               # 設定ファイル例
│   └── test.toml                  # テスト用設定
├── test/                          # テストスイート
├── scripts/                       # ユーティリティスクリプト
├── data/                          # データディレクトリ
├── Project.toml                   # 依存関係定義
└── README.md                      # このファイル
```

## トラブルシューティング

### エラー: "dz配列の長さがnkと一致しません"

dz配列の要素数をnkと同じにしてください（各層の厚さを指定）

### エラー: "観測温度の空間サイズが不正です"

Y_obs配列の形状を確認してください。期待形状: (ni, nj, nt)

### 警告: "ディスク空き容量不足"

出力ファイル用に十分な空き容量を確保してください

### 計算が遅い

- `--threads`オプションでスレッド数を増やす
- ウィンドウサイズを調整してメモリ使用量を最適化
- 問題サイズ（格子点数）を見直す

## 性能情報

### 小規模問題（3×3×5格子、50時間ステップ）

- 計算時間: 約0.5秒
- メモリ使用量: 約1MB
- ウィンドウ数: 8

### 中規模問題（10×10×20格子、1000時間ステップ）

- 計算時間: 数分〜数十分
- メモリ使用量: 約8MB
- ウィンドウ数: 約50

## ライセンス

このプロジェクトはIHCP研究グループのものです。

## 貢献者

- SHI ZHENGQI
- Claude Code (AI Assistant)

## 関連ファイル

- 元のPython実装: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`
- プロジェクト設定: `.claude/CLAUDE.md`

## 次のステップ

- A-2: 結果保存機能の拡張（NPZ形式、プロット機能）
- A-3: 並列計算の最適化
- A-4: GUI開発
