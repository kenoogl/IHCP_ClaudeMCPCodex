# クイックスタートガイド

Julia版DHCPソルバーを最速で実行するためのガイドです。

---

## ⚡ 3ステップで実行

### 1️⃣ Juliaのインストール

```bash
# macOS/Linux
curl -fsSL https://install.julialang.org | sh

# バージョン確認（1.10以上であること）
julia --version
```

### 2️⃣ 依存パッケージのインストール

```bash
cd /path/to/TrialClaudeMCPCodex
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

### 3️⃣ 実行

```bash
# 測定データなしで動作確認（5秒程度）
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --synthetic --ni 40 --nj 50 --nk 10 --nt 5

# 測定データありの場合（推奨設定、4スレッド）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10
```

---

## 📦 必要なファイル

### 最小構成（測定データなし）
```
TrialClaudeMCPCodex/
├── julia/
│   ├── Project.toml         # 必須
│   ├── src/                 # 必須（全ファイル）
│   └── scripts/
│       └── test_dhcp_solver.jl
└── shared/data/
    └── metal_thermal_properties.csv  # 必須（540B）
```

### 完全構成（測定データあり）
上記 + `shared/data/T_measure_700um_1ms.npy`（1.1GB）

---

## 🎯 推奨実行例

```bash
# 最速設定（Julia版最速、4.29秒 @4スレッド）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10

# 1スレッド実行（13.89秒）
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pcg --precond gs --nt 10

# 小規模テスト（測定データ不要、数秒）
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --synthetic --ni 40 --nj 50 --nk 10 --nt 5
```

---

## 🐛 よくあるエラーと解決

| エラー | 解決方法 |
|--------|---------|
| "Package X not found" | `julia --project=julia -e 'using Pkg; Pkg.instantiate()'` |
| "opening file T_measure_700um_1ms.npy" | `--synthetic`フラグを追加 |
| "MethodError" | Julia 1.10以上か確認、依存関係を再インストール |
| "conflicting import" | `include()`は使わず、コマンドラインから実行 |

---

## 📖 詳細情報

詳しい設定方法は `julia/SETUP.md` を参照してください。

---

**最終更新**: 2025年10月30日
