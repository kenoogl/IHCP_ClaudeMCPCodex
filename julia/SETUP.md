# Julia版DHCPソルバー 環境設定ガイド

このドキュメントでは、`test_dhcp_solver.jl`を他の環境で実行するための環境設定手順を説明します。

---

## 📋 前提条件

### 必須
- **Julia**: バージョン 1.10 以上
- **OS**: Linux, macOS, Windows（いずれも対応）
- **メモリ**: 最低 8GB 推奨（フルスケール計算の場合）
- **CPU**: マルチコア推奨（並列実行のため）

### オプション（測定データモードの場合のみ）
- **測定データファイル**: `T_measure_700um_1ms.npy` (1.1GB)

---

## 🚀 セットアップ手順

### ステップ1: Juliaのインストール

#### 方法A: juliaup（推奨）
```bash
# macOS/Linux
curl -fsSL https://install.julialang.org | sh

# Windows (PowerShell)
winget install julia -s msstore
```

#### 方法B: 公式サイトからダウンロード
https://julialang.org/downloads/ から対応するバイナリをダウンロード

#### インストール確認
```bash
julia --version
# Julia Version 1.10.0 以上であることを確認
```

---

### ステップ2: プロジェクトファイルの配置

以下のディレクトリ構成でファイルを配置します：

```
TrialClaudeMCPCodex/
├── julia/
│   ├── Project.toml          # 依存関係定義（必須）
│   ├── Manifest.toml         # 依存関係ロックファイル（推奨）
│   ├── src/                  # ソースコード（必須）
│   │   ├── IHCP_CGM.jl
│   │   ├── ThermalProperties.jl
│   │   ├── DataLoaders.jl
│   │   ├── solvers/
│   │   │   ├── DHCPSolver.jl
│   │   │   ├── CommonSolver.jl
│   │   │   └── ...
│   │   └── utils/
│   │       ├── validators.jl
│   │       ├── boundary_conditions.jl
│   │       └── ...
│   ├── scripts/              # 実行スクリプト
│   │   └── test_dhcp_solver.jl  # メインスクリプト
│   └── test/                 # テストファイル（オプション）
├── shared/
│   └── data/                 # データファイル
│       ├── metal_thermal_properties.csv  # 熱物性値（必須）
│       └── T_measure_700um_1ms.npy      # 測定データ（オプション）
└── README.md
```

**最小構成（syntheticモード）**:
- `julia/Project.toml`
- `julia/src/` 以下の全ファイル
- `julia/scripts/test_dhcp_solver.jl`
- `shared/data/metal_thermal_properties.csv`

**完全構成（測定データモード）**:
- 上記に加えて
- `shared/data/T_measure_700um_1ms.npy`

---

### ステップ3: 依存パッケージのインストール

#### 方法A: Manifest.tomlを使用（推奨、再現性保証）

```bash
# プロジェクトルートに移動
cd /path/to/TrialClaudeMCPCodex

# Julia環境に入る
julia --project=julia

# Juliaプロンプトで
julia> using Pkg
julia> Pkg.instantiate()  # Manifest.tomlから完全に復元
julia> exit()
```

`Pkg.instantiate()`は`Manifest.toml`に記録された**完全に同じバージョン**の依存関係をインストールします。

#### 方法B: Project.tomlのみを使用（Manifest.tomlがない場合）

```bash
cd /path/to/TrialClaudeMCPCodex
julia --project=julia

julia> using Pkg
julia> Pkg.resolve()      # 依存関係を解決
julia> Pkg.instantiate()  # パッケージをインストール
julia> exit()
```

この場合、最新の互換バージョンがインストールされます。

#### 方法C: コマンドライン（ワンライナー）

```bash
cd /path/to/TrialClaudeMCPCodex
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

#### インストール確認

```bash
julia --project=julia -e 'using Pkg; Pkg.status()'
```

以下のような出力が表示されれば成功：
```
Status `~/TrialClaudeMCPCodex/julia/Project.toml`
  [c7e460c6] ArgParse v1.2.0
  [6e4b80f9] BenchmarkTools v1.6.0
  [336ed68f] CSV v0.10.14
  ...
```

---

### ステップ4: データファイルの準備

#### 熱物性値ファイル（必須）

`shared/data/metal_thermal_properties.csv`を配置します。このファイルはSUS304の温度依存熱物性値を含みます（540バイト）。

#### 測定データファイル（測定データモードのみ）

`shared/data/T_measure_700um_1ms.npy`を配置します（1.1GB）。

**注意**: このファイルは大容量のため、gitignoreされています。以下の方法で取得してください：
- プロジェクト管理者から直接受け取る
- 共有ストレージからダウンロード
- 実験装置から生成

測定データがない場合は、**syntheticモード**で実行可能です（後述）。

---

## ✅ 動作確認

### 基本動作確認（syntheticモード）

測定データなしで動作確認できます：

```bash
cd /path/to/TrialClaudeMCPCodex

julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --synthetic --ni 40 --nj 50 --nk 10 --nt 5
```

**期待される出力**:
```
================================================================================
DHCP（直接熱伝導問題）ソルバー単体テスト
================================================================================
...
[Configuration]
  Test mode: Synthetic
  Grid size: 40 × 50 × 10 (N=20000)
...
Summary
  Total runtime: X.XX s
  DHCP share: XX.X%
================================================================================
```

### 測定データモード実行

測定データがある場合：

```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10
```

---

## 🎯 実行方法

### 基本実行

```bash
cd /path/to/TrialClaudeMCPCodex

# デフォルト設定（測定データモード、10ステップ）
julia --project=julia julia/scripts/test_dhcp_solver.jl

# GS前処理、PCGソルバー
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10
```

### 並列実行

```bash
# 4スレッドで実行
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10
```

**注意**: スレッド数は物理コア数または論理コア数に設定してください。

### Syntheticモード（測定データ不要）

```bash
# 小規模テスト
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --synthetic --ni 40 --nj 50 --nk 10 --nt 5

# フルサイズテスト
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --synthetic --ni 80 --nj 100 --nk 20 --nt 10
```

---

## 🔧 オプション一覧

| オプション | 説明 | デフォルト | 例 |
|-----------|------|----------|-----|
| `--solver` | ソルバータイプ | `pbicgstab` | `pcg`, `pbicgstab` |
| `--precond` | 前処理 | `diagonal` | `none`, `diagonal`, `gs` |
| `--nt` | 時間ステップ数 | `10` | `5`, `20`, `100` |
| `--basesize` | FLoops basesize | `600` | `300`, `1200` |
| `--ni` | X方向格子数 | `80` | `40`, `160` |
| `--nj` | Y方向格子数 | `100` | `50`, `200` |
| `--nk` | Z方向格子数 | `20` | `10`, `40` |
| `--synthetic` | 合成データモード | なし | （フラグのみ） |

### 推奨設定

#### 最速実行（Julia版最速設定）
```bash
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10
```

#### 高精度実行
```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pbicgstab --precond diagonal --nt 10
```

---

## 📊 ベンチマーク結果（参考値）

**環境**: Apple Silicon (aarch64), 4コア, Julia 1.10

| 前処理 | スレッド数 | 実行時間 | 反復回数 | 推奨度 |
|--------|-----------|---------|---------|--------|
| gs | 1 | 13.89秒 | 177回 | ⭐⭐⭐⭐⭐ |
| gs | 4 | 4.29秒 | 177回 | ⭐⭐⭐⭐⭐ |
| none | 1 | 30.21秒 | 884回 | ⭐⭐⭐ |
| none | 4 | 5.01秒 | 884回 | ⭐⭐⭐⭐ |
| diagonal | 1 | 32.62秒 | 737回 | ⭐⭐ |
| diagonal | 4 | 5.68秒 | 737回 | ⭐⭐⭐ |

**推奨設定**: `--solver pcg --precond gs` + マルチスレッド

---

## 🐛 トラブルシューティング

### エラー: "Package X not found"

**原因**: 依存パッケージがインストールされていない

**解決**:
```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

### エラー: "SystemError: opening file shared/data/T_measure_700um_1ms.npy"

**原因**: 測定データファイルがない

**解決**: syntheticモードで実行
```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl --synthetic
```

### エラー: "MethodError" または "LoadError"

**原因**: Juliaのバージョンが古い、または依存関係の不整合

**解決**:
```bash
# Juliaバージョン確認
julia --version  # 1.10以上であることを確認

# 依存関係の再インストール
julia --project=julia -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
```

### エラー: "WARNING: conflicting import"

**原因**: `include()`で直接読み込もうとしている

**解決**: コマンドラインから実行（`include()`は使わない）
```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl
```

### 実行が遅い

**原因**: 初回実行時はJITコンパイルのため遅い

**解決**:
- 2回目以降は高速化します
- マルチスレッドを使用: `JULIA_NUM_THREADS=4`
- GS前処理を使用: `--precond gs`

---

## 📚 関連ドキュメント

- **プロジェクトREADME**: `julia/README.md`
- **CLAUDE.md**: `.claude/CLAUDE.md`（開発ガイドライン）
- **性能比較レポート**: `docs/reports/julia_python_dhcp_comparison_report_20251030.md`

---

## 🆘 サポート

問題が解決しない場合：
1. `julia/README.md`を確認
2. `docs/INDEX.md`で関連ドキュメントを検索
3. プロジェクト管理者に連絡

---

**最終更新**: 2025年10月30日
**バージョン**: 1.0
