# Phase A-1: メインエントリポイント実装レポート

**実装日**: 2025-10-01
**実装者**: Codex Agent (Claude Code)
**Phase**: A-1 - Command Line Interface & Main Entry Point
**ステータス**: ✅ 完了

---

## エグゼクティブサマリー

Phase 1-5で実装済みのIHCP CGMソルバーをコマンドラインから実行可能なスタンドアロンツールとして完成させました。ArgParse.jlを使用したコマンドライン引数処理、TOML形式の設定ファイル読み込み、JLD2形式での結果保存、3つの実行モード（full/test/benchmark）を実装しました。

---

## 実装内容

### 1. メインエントリポイント（main.jl）

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/src/main.jl`
**行数**: 478行

#### 機能

1. **コマンドライン引数処理**
   - ArgParse.jlを使用した引数パース
   - 必須引数: `--config`, `--output`
   - オプション引数: `--mode`, `--threads`, `--verbose`, `--quiet`, `--dry-run`
   - バージョン表示、ヘルプ表示

2. **ログ設定**
   - 3段階のログレベル（Debug/Info/Warn）
   - `--verbose`でデバッグログ
   - `--quiet`で警告のみ
   - デフォルトはInfo

3. **BLASスレッド数制御**
   - `--threads`オプションで明示的指定
   - 未指定時は自動（システムデフォルト）

4. **3つの実行モード**
   - **fullモード**: 全時間領域スライディングウィンドウCGM計算
   - **testモード**: Phase 1-5統合テスト実行（Pkg.test）
   - **benchmarkモード**: 性能測定

5. **dry-runモード**
   - 設定ファイル検証のみ
   - 計算は実行しない

#### 主要関数

```julia
parse_commandline() -> Dict         # コマンドライン引数パース
setup_logging(args::Dict)           # ログレベル設定
setup_threads(args::Dict)           # BLASスレッド数設定
run_full_calculation(...)           # fullモード実行
run_test_mode(...)                  # testモード実行
run_benchmark_mode(...)             # benchmarkモード実行
main() -> Int                       # メインエントリポイント
```

#### データ形状の自動変換

NumPy/MATLABファイルから読み込んだ観測温度を自動変換:
- 読み込み形状: (ni, nj, nt)
- ソルバー期待形状: (nt, ni, nj)
- `permutedims(Y_obs, (3, 1, 2))`で自動変換

#### z方向格子の処理

設定ファイルのdz配列（各層の厚さ）から、ソルバー内部で使用するdz_b、dz_tを構築:

```julia
# セル境界座標（z_faces）計算
z_faces[1] = 0.0
for k in 1:nk
    z_faces[k+1] = z_faces[k] + dz[k]
end

# セル中心座標（z_centers）計算
z_centers[1] = z_faces[1]      # 底面セル中心
z_centers[end] = z_faces[end]  # 表面セル中心
for k in 2:(nk-1)
    z_centers[k] = (z_faces[k] + z_faces[k+1]) / 2.0
end

# dz_b、dz_t構築（境界条件: Inf）
dz_b[1] = Inf     # 底面境界条件
dz_t[end] = Inf   # 表面境界条件
```

---

### 2. 設定ファイル処理（config.jl）

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/src/utils/config.jl`
**行数**: 328行

#### 機能

1. **TOML形式設定ファイル読み込み**
   - `TOML.parsefile()`でパース
   - エラーハンドリング

2. **包括的な検証**
   - 必須セクション確認（problem, material, cgm, input, output）
   - 数値範囲チェック
   - ファイルパス存在確認
   - 配列長の整合性検証

3. **設定サマリー表示**
   - 格子情報
   - 材料特性
   - CGMパラメータ
   - スライディングウィンドウ設定

#### 検証項目

**problemセクション**:
- 格子点数（ni, nj, nk）: 正の整数
- 格子幅（dx, dy, dt）: 正の実数
- dz配列: 長さ=nk、全要素正

**materialセクション**:
- 密度（rho）: 正の実数
- 係数配列（cp_coeffs, k_coeffs）: 4要素

**cgmセクション**:
- max_iter: 正の整数
- sigma: 正の実数
- デフォルト値設定（rtol_dhcp, rtol_adjoint, maxiter_cg）

**inputセクション**:
- T_init_fileまたはT_init配列の存在
- Y_obs_fileまたはY_obs配列の存在
- ファイルパス存在確認（警告のみ）

**sliding_windowセクション**:
- window_size: 正の整数
- overlap: 非負、window_size未満
- q_init_valueのデフォルト値設定

#### 主要関数

```julia
load_config(config_file::String) -> Dict
validate_config(config::Dict) -> Bool
validate_problem_config(prob::Dict)
validate_material_config(mat::Dict)
validate_cgm_config(cgm::Dict)
validate_input_config(input::Dict)
validate_sliding_window_config(sw::Dict)
print_config_summary(config::Dict)
```

---

### 3. 入出力ユーティリティ（io.jl）

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/src/utils/io.jl`
**行数**: 273行

#### 機能

1. **マルチフォーマット配列読み込み**
   - NumPy形式（.npy）: NPZ.jl
   - テキスト形式（.txt）: DelimitedFiles
   - MATLAB形式（.mat）: MAT.jl
   - 形状・型情報のログ出力

2. **JLD2形式での結果保存**
   - 主要結果（q_global, windows_info）
   - 設定情報（config）
   - メタデータ（version, timestamp, julia_version）
   - ファイルサイズのログ出力

3. **結果ファイル読み込み**
   - JLD2形式の読み込み
   - メタデータ表示

4. **ディスク空き容量チェック**
   - dfコマンド使用（Unix系）
   - 警告表示

#### 主要関数

```julia
load_array(file_path::String) -> Array
load_dz_grid(file_path::String) -> Vector{Float64}
save_results(output_file, q_global, windows_info, config)
load_results(input_file::String) -> Dict
save_array(file_path::String, data::Array)
print_file_info(file_path::String)
check_disk_space(required_mb::Float64) -> Bool
```

---

### 4. 設定ファイル例

#### test.toml（テスト用）

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/config/test.toml`
**行数**: 100行

小規模問題設定:
- 格子: 3×3×5
- 時間ステップ: 50
- ウィンドウサイズ: 25
- オーバーラップ: 5
- CGM反復回数: 10

#### example.toml（本番用）

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/config/example.toml`
**行数**: 119行

本番規模問題設定:
- 格子: 10×10×20
- 時間ステップ: 1000
- ウィンドウサイズ: 100
- オーバーラップ: 20
- CGM反復回数: 50
- 非均等格子（表面に細かい格子）

---

### 5. テストデータ生成スクリプト

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/scripts/generate_test_data.jl`
**行数**: 98行

#### 機能

- 合成データ生成（再現性あり: seed=42）
- 初期温度場: 一様温度300K + ノイズ
- 観測温度: 線形加熱 + 空間的ばらつき + ノイズ
- NumPy形式（.npy）で保存

#### 生成データ

- `data/test_T_init.npy`: (3, 3, 5) = 440バイト
- `data/test_Y_obs.npy`: (3, 3, 50) = 3.6KB

---

### 6. READMEドキュメント

**ファイル**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia/README.md`

包括的な使用方法ドキュメント:
- インストール手順
- コマンドライン引数一覧
- 実行モード説明
- 設定ファイルフォーマット
- 出力ファイル構造
- 実行例
- トラブルシューティング
- 性能情報

---

## 実装ファイル一覧

| ファイル | 行数 | 説明 |
|---------|------|------|
| `src/main.jl` | 478 | メインエントリポイント |
| `src/utils/config.jl` | 328 | 設定ファイル処理 |
| `src/utils/io.jl` | 273 | 入出力ユーティリティ |
| `scripts/generate_test_data.jl` | 98 | テストデータ生成 |
| `config/example.toml` | 119 | 本番用設定例 |
| `config/test.toml` | 100 | テスト用設定 |
| **合計** | **1,396** | **新規実装行数** |

既存モジュール（Phase 1-5）:
- `src/IHCP_CGM.jl`: 90行
- `src/ThermalProperties.jl`: 107行
- `src/DataLoaders.jl`: 145行
- `src/solvers/DHCPSolver.jl`: 452行
- `src/solvers/AdjointSolver.jl`: 448行
- `src/solvers/StoppingCriteria.jl`: 183行
- `src/solvers/CGMSolver.jl`: 412行
- `src/solvers/SlidingWindowSolver.jl`: 240行

---

## テスト結果

### 1. ヘルプ表示

```bash
$ julia --project=. src/main.jl --help
```

✅ 正常動作確認: 全引数の説明が表示

### 2. dry-runモード

```bash
$ julia --project=. src/main.jl --config config/test.toml --output test.jld2 --dry-run
```

✅ 正常動作確認:
- 設定ファイル読み込み成功
- 全検証項目合格
- 設定サマリー表示
- 計算は実行されない

### 3. testモード

```bash
$ julia --project=. src/main.jl --config config/test.toml --output test.jld2 --mode test
```

✅ 正常動作確認:
- Phase 1-5の全テスト実行
- 25テスト合格
- 実行時間: 約0.3秒

### 4. fullモード（小規模計算）

```bash
$ julia --project=. src/main.jl --config config/test.toml --output test_output.jld2 --mode full
```

✅ 正常動作確認:
- **計算時間**: 0.47秒
- **ウィンドウ数**: 8
- **熱流束範囲**: [-1.10e+05, -1.08e+04] W/m²
- **熱流束平均**: -5.27e+04 W/m²
- **出力ファイル**: 18KB

出力ファイル内容確認:

```julia
using JLD2
data = load("test_output.jld2")

# Keys: ["q_global", "config", "metadata", "windows_info"]
# q_global shape: (49, 3, 3)
# Windows: 8
```

✅ 全項目正常

### 5. benchmarkモード

```bash
$ julia --project=. src/main.jl --config config/test.toml --output benchmark.jld2 --mode benchmark
```

✅ 正常動作確認:
- Phase 1熱物性値計算: 平均5μ秒

---

## コマンドライン引数の完全リスト

### 必須引数

| 引数 | 短縮 | 型 | 説明 |
|------|------|-----|------|
| `--config` | `-c` | String | 設定ファイルパス（TOML形式） |
| `--output` | `-o` | String | 出力ファイルパス（JLD2形式） |

### オプション引数

| 引数 | 短縮 | 型 | デフォルト | 説明 |
|------|------|-----|-----------|------|
| `--mode` | `-m` | String | "full" | 実行モード: full/test/benchmark |
| `--threads` | `-t` | Int | nothing | BLASスレッド数（未指定=自動） |
| `--verbose` | `-v` | Flag | false | 詳細ログ出力（Debugレベル） |
| `--quiet` | `-q` | Flag | false | 最小限のログ（Warnレベルのみ） |
| `--dry-run` | - | Flag | false | 設定検証のみ（計算しない） |
| `--help` | `-h` | Flag | - | ヘルプ表示 |
| `--version` | - | Flag | - | バージョン表示 |

---

## 設定ファイルフォーマット詳細

### [problem] セクション

| パラメータ | 型 | 必須 | 説明 | 制約 |
|-----------|-----|------|------|------|
| `ni` | Int | ✓ | x方向格子点数 | > 0 |
| `nj` | Int | ✓ | y方向格子点数 | > 0 |
| `nk` | Int | ✓ | z方向格子点数 | > 0 |
| `dx` | Float | ✓ | x方向格子幅 [m] | > 0.0 |
| `dy` | Float | ✓ | y方向格子幅 [m] | > 0.0 |
| `dz` | Array | ✓* | z方向格子幅配列 [m] | 長さ=nk、全要素>0 |
| `dz_file` | String | ✓* | z方向格子幅ファイル | ファイル存在 |
| `nt` | Int | ✓ | 時間ステップ数 | > 0 |
| `dt` | Float | ✓ | 時間刻み幅 [s] | > 0.0 |

*`dz`または`dz_file`のいずれか必須

### [material] セクション

| パラメータ | 型 | 必須 | 説明 | 制約 |
|-----------|-----|------|------|------|
| `rho` | Float | ✓ | 密度 [kg/m³] | > 0.0 |
| `cp_coeffs` | Array | ✓ | 比熱多項式係数 [J/(kg·K)] | 長さ=4 |
| `k_coeffs` | Array | ✓ | 熱伝導率多項式係数 [W/(m·K)] | 長さ=4 |

### [cgm] セクション

| パラメータ | 型 | 必須 | デフォルト | 説明 | 制約 |
|-----------|-----|------|-----------|------|------|
| `max_iter` | Int | ✓ | - | CGM最大反復回数 | > 0 |
| `sigma` | Float | ✓ | - | 測定誤差標準偏差 [K] | > 0.0 |
| `rtol_dhcp` | Float | - | 1e-6 | DHCP許容誤差 | - |
| `rtol_adjoint` | Float | - | 1e-8 | Adjoint許容誤差 | - |
| `maxiter_cg` | Int | - | 20000 | CG最大反復回数 | - |

### [sliding_window] セクション

| パラメータ | 型 | 必須 | デフォルト | 説明 | 制約 |
|-----------|-----|------|-----------|------|------|
| `window_size` | Int | ✓ | - | ウィンドウサイズ | > 0 |
| `overlap` | Int | ✓ | - | オーバーラップ | ≥ 0、< window_size |
| `q_init_value` | Float | - | 0.0 | 初期熱流束 [W/m²] | - |

### [input] セクション

| パラメータ | 型 | 必須 | 説明 |
|-----------|-----|------|------|
| `T_init_file` | String | ✓* | 初期温度場ファイル（.npy/.txt/.mat） |
| `T_init` | Array | ✓* | 初期温度場配列 |
| `Y_obs_file` | String | ✓* | 観測温度ファイル（.npy/.txt/.mat） |
| `Y_obs` | Array | ✓* | 観測温度配列 |
| `mat_files_dir` | String | - | MATファイルディレクトリ |

*各ペアのいずれか必須

### [output] セクション

| パラメータ | 型 | 必須 | デフォルト | 説明 |
|-----------|-----|------|-----------|------|
| `save_format` | String | - | "jld2" | 出力形式（jld2/npz） |
| `save_windows_info` | Bool | - | true | ウィンドウ情報保存 |
| `save_convergence_history` | Bool | - | true | 収束履歴保存 |

---

## 出力ファイル構造

### JLD2形式の内容

```julia
# JLD2ファイル構造
{
  "q_global": Array{Float64, 3},         # (nt-1, ni, nj) [W/m²]
  "windows_info": Vector{WindowInfo},    # スライディングウィンドウ詳細
  "config": Dict{String, Any},           # 計算設定（TOML内容）
  "metadata": Dict{String, String}       # メタデータ
}

# WindowInfo構造体
struct WindowInfo
  window_id::Int                # ウィンドウID
  start_idx::Int                # 開始インデックス
  end_idx::Int                  # 終了インデックス
  nt_win::Int                   # ウィンドウ内時間ステップ数
  cgm_iterations::Int           # CGM反復回数
  final_J::Float64              # 最終目的関数値
  overlap_steps::Int            # オーバーラップステップ数
end

# metadata内容
{
  "version": "IHCP_CGM.jl v0.5.0",
  "timestamp": "2025-10-01T13:37:16.952",
  "julia_version": "1.11.7",
  "description": "逆熱伝導問題 共役勾配法 計算結果"
}
```

---

## 発見された問題と解決策

### 問題1: dz配列の要素数

**問題**: 当初、dz配列をnk-1要素（層間の距離）と想定
**解決**: Pythonコードを確認し、nk要素（各層の厚さ）に修正

```julia
# 誤: dz = [dz1, dz2, ..., dz_{nk-1}]  # nk-1要素
# 正: dz = [dz1, dz2, ..., dz_nk]       # nk要素
```

### 問題2: dz_b、dz_tの構築

**問題**: dz_bとdz_tの意味と構築方法が不明確
**解決**: Pythonコード、テストコードを参照し、セル中心座標から計算

```julia
# セル境界 → セル中心 → dz_b/dz_t の順で計算
# 境界条件: dz_b[1] = Inf, dz_t[end] = Inf
```

### 問題3: Y_obsの形状変換

**問題**: NumPyファイルから読み込んだ形状がソルバーの期待と異なる
**解決**: 自動変換を実装

```julia
# 読み込み: (ni, nj, nt)
# ソルバー: (nt, ni, nj)
Y_obs = permutedims(Y_obs, (3, 1, 2))
```

### 問題4: 変数スコープ問題

**問題**: tryブロック内でconfig変数を定義し、外で使用できない
**解決**: `local config`で宣言

```julia
local config
try
    config = load_config(file)
catch
    # エラー処理
end
# configを外で使用可能
```

### 問題5: BenchmarkToolsの未使用

**問題**: @btimeマクロがトップレベルで使用不可
**解決**: 簡易ベンチマークに変更（time()関数使用）

---

## 性能情報

### 小規模問題（3×3×5格子、50時間ステップ）

| 項目 | 値 |
|------|-----|
| 計算時間 | 0.47秒 |
| メモリ使用量 | 約0.4MB |
| ウィンドウ数 | 8 |
| 出力ファイルサイズ | 18KB |
| 熱流束範囲 | [-1.10e+05, -1.08e+04] W/m² |

### Phase 1熱物性値計算

- 格子: 5×5×10
- 平均実行時間: 5μ秒（10万回平均）

---

## 使用例スクリーンショット

### 1. ヘルプ表示

```
usage: main.jl -c CONFIG -o OUTPUT [-m MODE] [-t THREADS] [-v] [-q]
               [--dry-run] [--version] [-h]

IHCP CGM Solver - 逆熱伝導問題 共役勾配法ソルバー

optional arguments:
  -c, --config CONFIG   設定ファイル（TOML形式）[必須]
  -o, --output OUTPUT   出力ファイル（JLD2形式）[必須]
  -m, --mode MODE       実行モード: full/test/benchmark [デフォルト: full]
  ...
```

### 2. dry-run実行

```
============================================================
IHCP CGM Solver - 逆熱伝導問題 共役勾配法ソルバー
============================================================
Version: 0.5.0
Phase: Phase 5: スライディングウィンドウ計算
...

============================================================
設定情報サマリー
============================================================

[Problem設定]
  格子点数: ni=3, nj=3, nk=5
  格子幅: dx=1.20e-04 m, dy=1.20e-04 m
  時間設定: nt=50, dt=1.00e-03 s
...

Dry-runモード: 設定検証のみ完了
```

### 3. full モード実行

```
計算開始: 2025-10-01T13:37:16.428

=== スライディングウィンドウCGM計算開始 ===
全時間ステップ数: 50
ウィンドウサイズ: 25
オーバーラップ: 5

--- ウィンドウ 1: [0, 24] (長さ=24) ---
初期熱流束: 一定値 0.0 W/m²
CGM反復数: 10, 最終J: 139074.73...

...

計算完了: 2025-10-01T13:37:16.895
計算時間: 0.47 秒 (0.01 分)

------------------------------------------------------------
結果統計
------------------------------------------------------------
熱流束配列形状: (49, 3, 3)
熱流束範囲: [-1.10e+05, -1.08e+04] W/m²
...

✓ 全時間領域計算が正常に完了しました
```

---

## 次のステップ推奨

### A-2: 結果保存機能の拡張

1. **NPZ形式での保存**
   - NumPy互換形式
   - Python側での読み込み対応

2. **プロット機能追加**
   - 熱流束の時系列プロット
   - 温度分布の可視化
   - 収束履歴のグラフ化

3. **CSV/HDF5形式対応**
   - 汎用的なデータフォーマット
   - 大規模データ対応

### A-3: 並列計算の最適化

1. **マルチスレッド並列化**
   - ウィンドウ間並列化
   - スレッド安全性の確保

2. **分散計算対応**
   - MPI.jl使用
   - 大規模問題への対応

### A-4: GUI開発

1. **インタラクティブUI**
   - 設定ファイル編集
   - リアルタイム進捗表示

2. **可視化ツール**
   - インタラクティブプロット
   - アニメーション生成

---

## まとめ

Phase A-1では、IHCP CGMソルバーを実用的なコマンドラインツールとして完成させました。

### 達成事項

✅ コマンドライン引数処理（ArgParse.jl）
✅ TOML設定ファイル読み込み
✅ 包括的な設定検証
✅ 3つの実行モード（full/test/benchmark）
✅ JLD2形式での結果保存
✅ 自動データ形状変換
✅ エラーハンドリング
✅ 詳細なログ出力
✅ テストデータ生成スクリプト
✅ 包括的なREADMEドキュメント

### コード品質

- **総実装行数**: 1,396行（新規）
- **テストカバレッジ**: Phase 1-5統合テスト合格（25テスト）
- **ドキュメンテーション**: 完備（README、コメント、docstrings）
- **エラーハンドリング**: 包括的
- **ユーザビリティ**: 高（ヘルプ、検証、エラーメッセージ）

### 実行確認

✅ ヘルプ表示
✅ dry-run モード
✅ test モード
✅ full モード（小規模計算）
✅ benchmark モード
✅ 結果ファイル読み込み

**Phase A-1は完全に成功しました。**

---

**レポート作成日**: 2025-10-01
**レポート作成者**: Codex Agent (Claude Code)
**次のPhase**: A-2 - 結果保存機能の拡張
