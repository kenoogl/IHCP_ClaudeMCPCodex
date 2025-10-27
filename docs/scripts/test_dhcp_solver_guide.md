# test_dhcp_solver.jl 説明書

**スクリプトパス**: `julia/scripts/test_dhcp_solver.jl`
**作成日**: 2025年10月23日
**最終更新**: 2025年10月27日
**目的**: DHCP（直接熱伝導問題）ソルバーの性能と精度検証

---

## 概要

このスクリプトはDHCPソルバー単体の性能を測定するための軽量テストツールです。実際のIHCP問題の全体（逆解析）ではなく、**順問題（DHCP）のみ**を実行して以下を検証します：

- ソルバーの計算性能（実行時間、スケーラビリティ）
- 収束性（反復回数、収束挙動）
- 数値精度（basesize変更時の再現性）

---

## 使用方法

### 基本実行

```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl
```

### オプション指定

```bash
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pbicgstab \
  --precond gs \
  --nt 10 \
  --basesize 600
```

### 利用可能なオプション

| オプション | 説明 | デフォルト値 | 選択肢 |
|-----------|------|------------|--------|
| `--solver` | ソルバータイプ | `pbicgstab` | `pbicgstab`, `pcg` |
| `--precond` | 前処理タイプ | `diagonal` | `diagonal`, `gs`, `none` |
| `--nt` | タイムステップ数 | `10` | 任意の正整数 |
| `--basesize` | FLoops並列化粒度 | `600` | 任意の正整数 |
| `-h`, `--help` | ヘルプ表示 | - | - |

---

## スクリプトの構成

### 1. コマンドライン引数パース

```julia
solver_type, precond_type, nt, basesize = parse_command_line_args()
```

引数を解析して実行パラメータを取得します。

### 2. Backend設定

```julia
set_backend_config(basesize=basesize)
```

FLoopsの並列化パラメータ（basesize）を設定します。

### 3. 問題定義

```julia
dt = 1.0e-3              # 時間ステップ幅 [s]
ni, nj, nk = 80, 100, 20 # 格子数
dx = 0.12e-3             # x方向格子幅 [m]
dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))  # y方向格子幅 [m]
Lz = 0.5e-3              # z方向全長 [m]
stretch_factor = 3.0      # z方向格子の集中度
```

80×100×20格子（N=160,000要素）の3D問題を定義します。

### 4. 格子生成

```julia
z_centers, ΔZ, z_faces, dz = build_z_grid(nk, Lz, stretch_factor)
```

z方向に非均等格子を生成（表面側に集中）。

### 5. 熱物性値読み込み

```julia
rho, cp_coeffs, k_coeffs = load_material_properties()
```

SUS304の熱物性値（密度、比熱係数、熱伝導率係数）を読み込みます。

### 6. 測定データ読み込み

```julia
Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)
Y_obs = permutedims(Y_obs_python, (2, 3, 1))  # (nt,ni,nj) → (ni,nj,nt)
```

**データファイル**: `shared/data/T_measure_700um_1ms.npy`（1.1GB）

**⚠️ 重要**: このデータは以下の用途にのみ使用されます：
- **初期温度場の生成**（最初のフレームをz方向に複製）
- **結果の検証**（計算結果との残差計算）

**DHCP計算自体には使われません**（熱流束はゼロに設定）。

### 7. 初期温度場の生成

```julia
T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
```

測定データの最初のフレーム（t=1）を20層にコピーして3D初期温度場を生成します。

### 8. 熱流束の設定

```julia
q_test = zeros(Float64, ni, nj, nt - 1)
```

**⚠️ 重要**: 熱流束は**ゼロ**に設定されています。

これは以下を意味します：
- 表面から熱の出入りがない（断熱に近い）
- 温度はほぼ変化しない（わずかな熱拡散のみ）
- **物理的な妥当性テストではない**
- **ソルバーの性能テストとして適切**

### 9. DHCP solve

```julia
T_result, iter_counts, residuals = solve_dhcp!(
  T_init,           # 初期温度場 (ni, nj, nk)
  q_test,           # 表面熱流束 (ni, nj, nt-1) ← ゼロ！
  work,             # ワークバッファ
  nt,               # タイムステップ数
  rho,              # 密度
  cp_coeffs,        # 比熱係数
  k_coeffs,         # 熱伝導率係数
  dx,               # x方向格子幅
  dy,               # y方向格子幅
  z_centers,        # セル中心座標（ガイドセル込み）
  ΔZ,               # セル幅（ガイドセル込み）
  dt;               # 時間ステップ幅
  rtol = 1.0e-6,
  maxiter = 20000,
  verbose = true,
  solver = solver_type,
  smoother = precond_type
)
```

**境界条件**:
- X方向: 両端断熱
- Y方向: 両端断熱
- Z下面: 断熱
- Z上面: 熱流束指定（q_test）← **ゼロなので実質断熱**

### 10. 残差診断

```julia
T_bottom = @view T_result[:, :, 1, :]  # 底面温度 (ni, nj, nt)
residual = T_bottom .- Y_obs
rms_error = sqrt(mean(residual .^ 2))
max_error = maximum(abs.(residual))
```

計算結果の底面温度と測定データを比較して残差を計算します。

**残差の意味**:
- q=0で計算しているため、物理的な意味は限定的
- 主に「初期温度場を全層で同一にした影響」を反映
- 数値的一貫性の確認には有効

---

## 測定データの利用状況

### データの流れ

```
T_measure_700um_1ms.npy (1.1GB)
  ↓
Y_obs (nt, ni, nj)  ← 測定温度データ
  ↓
├─→ T_init (ni, nj, nk)  ← 初期条件（最初のフレームをz方向複製）
│   └─→ solve_dhcp!() に渡す
│
└─→ residual = T_result - Y_obs  ← 結果検証
```

### データが使われる場所

1. **初期条件の生成** (行275):
   ```julia
   T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
   ```
   - 測定データの最初のフレーム（t=1）を使用
   - z方向に20層複製して3D温度場を生成

2. **結果検証** (行326):
   ```julia
   residual = T_bottom .- Y_obs
   ```
   - 計算結果の底面温度と測定データを比較
   - RMS誤差、最大誤差を計算

### データが使われない場所

- **DHCP計算の駆動力（熱流束）**: `q_test = zeros(...)` ← **ゼロ！**
- **境界条件**: 全面断熱（熱流束ゼロ）

---

## 出力結果

### コンソール出力

```
================================================================================
DHCP（直接熱伝導問題）ソルバー単体テスト
================================================================================
Project root: /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/

[Configuration]
  Solver: pbicgstab
  Preconditioner: gs
  Time steps: 10
  FLoops basesize: 600
  Julia threads: 4

[1/5] Grid and material parameters
  grid: ni=80, nj=100, nk=20
  spacing: dx=1.200e-04 m, dy=1.671e-04 m
  time steps: nt=10, dt=1.000e-03 s
  rho (reference): 7823.493963 kg/m^3

[2/5] Loading measurement data
  loading measurement data from: .../T_measure_700um_1ms.npy
  measurement subset shape (Python order): (10, 80, 100)
  load time: 1.56 s
  initial temperature range: 550.11 ~ 587.98 K

[3/5] 設定された熱流束パターン
  q_test shape: (80, 100, 9)
  q_test range: 0.0000e+00 ~ 0.0000e+00 W/m^2  ← ゼロ！

[4/5] Running DHCP solver
============================================================
Start DHCP direct solver
============================================================
格子: 80×100×20 (N=160000)
時間ステップ: 10, dt=0.001s
CG許容誤差: rtol=1.0e-6, maxiter=20000
============================================================
=== Boundary Conditions ===
X-minus: Adiabatic
X-plus : Adiabatic
Y-minus: Adiabatic
Y-plus : Adiabatic
Z-minus: Adiabatic
Z-plus : Distribution  ← q=0なので実質断熱
===================
[t=2/10] converged: Iteration= 11 : Res_0= 0.3654 : time=0.978s
[t=3/10] converged: Iteration= 12 : Res_0= 0.1703 : time=0.321s
...
[t=10/10] converged: Iteration= 11 : Res_0= 0.0135 : time=0.294s
============================================================
DHCP直接ソルバー完了
  最終温度範囲: 550.11 - 587.98 K
============================================================
  DHCP elapsed time: 4.86 s

[5/5] Residual analysis
  RMS residual: 2.9516e-01 K
  Max residual: 2.1465e+00 K
  Mean temperature: 572.20 K
  Temperature range: 550.11 ~ 587.98 K

Summary
  Total runtime: 6.83 s
  DHCP share: 71.2%
================================================================================
```

### 主要な出力指標

| 指標 | 説明 | 典型値 |
|------|------|--------|
| DHCP elapsed time | DHCP計算のみの実行時間 | 4.86s (nt=10) |
| Total runtime | 全体の実行時間 | 6.83s (nt=10) |
| DHCP share | DHCP計算の割合 | 71.2% (nt=10) |
| RMS residual | 残差のRMS値 | 0.30K (q=0) |
| Max residual | 残差の最大値 | 2.15K (q=0) |
| Iteration count | 各ステップの反復回数 | 11-13反復 |
| Convergence rate | 収束速度（Res_0） | 0.3654→0.0135 |

---

## 典型的な使用例

### 例1: デフォルト設定でテスト

```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl
```

### 例2: basesize最適化テスト

```bash
# basesize=1（ベースライン）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pbicgstab --precond gs --nt 10 --basesize 1

# basesize=1000（最適値候補）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pbicgstab --precond gs --nt 10 --basesize 1000
```

**期待される結果**: basesize=1000で16倍高速化

### 例3: スケーラビリティテスト

```bash
# nt=10
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl --nt 10

# nt=50
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl --nt 50

# nt=100
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl --nt 100
```

**期待される結果**:
- ステップあたり計算時間が減少（0.486s → 0.250s）
- 反復回数が減少（12.0 → 8.65）

---

## 制限事項と注意点

### 1. 物理的妥当性の限界

**⚠️ このテストは物理的に妥当な問題を解いていません**

- 熱流束がゼロ（q=0）
- 温度はほぼ変化しない
- 実際のIHCP問題とは異なる挙動

**用途**: ソルバーの性能テストとしてのみ有効

### 2. 測定データの利用が限定的

測定データは以下にのみ使用：
- 初期条件の生成（最初のフレームのみ）
- 結果検証（残差計算）

**DHCP計算の駆動力（熱流束）としては使われない**

### 3. 残差の解釈

```
RMS residual: 2.9516e-01 K
Max residual: 2.1465e+00 K
```

この残差は以下を反映：
- 初期温度場を全層で同一にした影響
- わずかな熱拡散の影響
- **物理的な意味は限定的**（q=0のため）

### 4. 実際のIHCP問題との違い

| 項目 | test_dhcp_solver.jl | 実際のIHCP |
|------|-------------------|-----------|
| 熱流束 | q=0（ゼロ） | 非ゼロ（逆解析対象） |
| 温度変化 | ほぼなし | 大きく変化 |
| 収束性 | 容易 | より困難 |
| 計算時間 | 比較的短い | より長い |
| 物理的意味 | 限定的 | 重要 |

---

## 実際のIHCP問題でテストするには

### 修正が必要な箇所

```julia
# 現在（q=0）
q_test = zeros(Float64, ni, nj, nt - 1)

# 修正例1: 一定熱流束
q_test = fill(1000.0, ni, nj, nt - 1)  # 1000 W/m²

# 修正例2: 測定温度から逆算した熱流束（要実装）
q_test = compute_heat_flux_from_measurement(Y_obs, ...)
```

---

## 性能ベンチマーク

### basesize最適化（nt=10、4スレッド）

| basesize | DHCP時間 | Total時間 | 高速化率 |
|----------|----------|-----------|---------|
| 1 | 108.93 s | 111.21 s | 基準 |
| 600 | 4.86 s | 6.83 s | **22.4倍** |
| 1000 | 6.52 s | 9.30 s | **16.7倍** |
| 10000 | 9.91 s | 11.93 s | 11.0倍 |

**推奨値**: basesize=600（デフォルト）

### スケーラビリティ（basesize=600、4スレッド）

| nt | DHCP時間 | ステップあたり時間 | 平均反復回数 |
|----|----------|------------------|-------------|
| 10 | 4.86 s | 0.486 s/step | 12.0反復 |
| 50 | 14.05 s | 0.281 s/step | 9.31反復 |
| 100 | 25.03 s | 0.250 s/step | 8.65反復 |

**観察**:
- ステップあたり時間が48.6%短縮（nt=10→100）
- 反復回数が28%減少（ホットスタート効果）

---

## 関連ドキュメント

- **レポート**: `docs/reports/phase5_2_step1_dhcp_basesize_validation.md`
- **詳細分析**: DHCP単体テスト スケーラビリティレポート
- **ソースコード**: `julia/src/solvers/DHCPSolver.jl`
- **テストスクリプト**: `julia/test/test_dhcp_solver.jl`（別物）

---

## まとめ

### このスクリプトの目的

✅ **DHCPソルバーの性能検証**
- 計算時間の測定
- スケーラビリティの評価
- basesizeパラメータの最適化
- 数値的一貫性の確認

❌ **物理的妥当性の検証ではない**
- q=0（ゼロ熱流束）で計算
- 温度変化がほぼない
- 実際のIHCP問題とは異なる

### 使い分け

- **性能テスト**: このスクリプトが最適 ✅
- **物理的検証**: `run_10steps_fullsize_test.jl` を使用
- **逆解析**: `run_sliding_window.jl` を使用

---

**作成者**: Claude Code
**最終更新**: 2025年10月27日
