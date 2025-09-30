# Julia移植計画書

**プロジェクト**: 逆熱伝導問題（IHCP）CGMソルバー Python→Julia移植
**策定日**: 2025-09-30
**対象ファイル**: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`（1701行、32関数）
**目標**: 高性能・保守性・数値精度の三位一体改善

---

## 1. 技術対応表（Python → Julia）

### 1.1 コアライブラリ対応

| Python | Julia | 移植戦略 | 注意点 |
|--------|-------|---------|--------|
| **numpy** | `Base`/`LinearAlgebra` | 組み込み多次元配列 | 1始まりインデックス、列優先（Column-major） |
| **pandas** | `DataFrames.jl` + `CSV.jl` | CSV読込・データフレーム | `DataFrame`の操作APIは類似 |
| **scipy.sparse** | `SparseArrays.jl` | CSR/CSC疎行列 | Juliaは標準でCSC優先、CSRも利用可 |
| **scipy.sparse.linalg.cg** | `IterativeSolvers.jl` (`cg`) | 共役勾配法 | 前処理器オプションが豊富 |
| **scipy.io.loadmat** | `MAT.jl` | MATLABファイル読込 | v5～v7.3形式対応 |
| **numba @njit** | **不要**（Julia自体がJIT） | 型推論最適化 | `@inbounds`、`@simd`で高速化 |
| **numba prange** | `Threads.@threads` | マルチスレッド並列 | 環境変数 `JULIA_NUM_THREADS` で制御 |
| **psutil** | `Sys` モジュール | システムメモリ情報 | `Sys.total_memory()`、`Sys.free_memory()` |
| **gc.collect()** | `GC.gc()` | ガベージコレクション | Juliaでは通常不要（自動管理が優秀） |

### 1.2 数値計算パターン対応

| Python手法 | Julia等価実装 | 性能改善見込み |
|-----------|--------------|--------------|
| `np.empty((n,m,k))` | `Array{Float64}(undef, n, m, k)` | 同等 |
| `np.zeros((n,m))` | `zeros(Float64, n, m)` | 同等 |
| `field.ravel(order="F")` | `vec(field)` ※列優先がデフォルト | より自然 |
| `for i in prange(ni):` | `Threads.@threads for i in 1:ni` | 同等～向上 |
| `@njit(fastmath=False)` | デフォルト（`@fastmath`で明示高速化） | 精度優先がデフォルト |
| `scipy.sparse.diags()` | `spdiagm()` | APIが簡潔 |
| LinearOperator | 関数定義 or `LinearMaps.jl` | より直感的 |

### 1.3 インデックス変換ルール

**重要**: Pythonは0始まり、Juliaは1始まり

```julia
# Python: T[0, :, -1]  →  Julia: T[1, :, end]
# Python: range(nt)     →  Julia: 1:nt
# Python: k_ijk in [0, nk-1]  →  Julia: k_ijk in [1, nk]
```

**配列レイアウト**:
- Python: 行優先（C順序、`order="C"`）
- Julia: **列優先（Fortran順序）** ← 既存コードが `order="F"` 使用で相性良好

---

## 2. 移植優先順位（5フェーズ戦略）

### Phase 1: 基盤構築（Week 1-2）

**目標**: Julia環境整備 + ユーティリティ関数移植

#### 対象関数（10個）
1. `polyval_numba()` → `polyval()`
2. `thermal_properties_calculator()` → `thermal_properties!()`
3. `check_field_finite()` → `isfinite_field()`
4. `check_temperature_range()` → `check_temperature_range()`
5. `check_flux_range()` → `check_flux_range()`
6. `check_gradient_magnitude()` → `check_gradient_magnitude()`
7. `extract_sorted_mat_files()` → `extract_sorted_mat_files()`
8. `load_region_temperature()` → `load_region_temperature()`
9. `get_memory_info()` → `get_memory_info()`
10. `safe_load_large_data()` → `safe_load_large_data()`

#### 成果物
- `julia/src/utils/thermal_properties.jl`
- `julia/src/utils/validators.jl`
- `julia/src/utils/data_loader.jl`
- `julia/src/utils/memory_monitor.jl`
- `julia/test/test_thermal_properties.jl`（TDD）
- `julia/test/test_validators.jl`

#### 成功基準
- [ ] 熱物性値計算がPython結果と1e-12以内で一致
- [ ] MATファイル読込が正常動作
- [ ] メモリ監視が1.1GBデータで警告発火
- [ ] 全単体テスト通過

#### リスク
- **低**: 数学関数のみ、外部依存少

---

### Phase 2: 直接ソルバー（DHCP）移植（Week 3-5）

**目標**: 前進時間差分 + 疎行列CG法の実装

#### 対象関数（4個）
1. `coeffs_and_rhs_building_DHCP()` → `build_dhcp_system!()`
2. `assemble_A_DHCP()` → `assemble_dhcp_matrix()`
3. `multiple_time_step_solver_DHCP()` → `solve_dhcp!()`
4. `check_diffusion_stability()` → `check_diffusion_stability()`

#### 技術課題
- **疎行列構築**: `SparseArrays.sparse(I, J, V, n, n)` でCSC形式
- **CG法**: `IterativeSolvers.cg!(x, A, b; Pl=前処理器)` 使用
- **前処理器**: 対角スケーリング → `Diagonal(inv_diag)`
- **ホットスタート**: 初期値 `x0` 引数で前時間ステップ解渡し

#### 成果物
- `julia/src/solvers/dhcp_solver.jl`
- `julia/test/test_dhcp_solver.jl`
- `julia/benchmarks/bench_dhcp.jl`（性能比較）

#### 成功基準
- [ ] 既知解問題（製作解）で相対誤差 < 1e-6
- [ ] Python結果との温度場差分 < 1e-8
- [ ] CG収束回数が同等以下
- [ ] 計算時間がPython比 **2倍以上高速化**

#### リスク
- **中**: 疎行列のメモリレイアウト差異、CG収束特性

---

### Phase 3: 随伴ソルバー（Adjoint）移植（Week 6-8）

**目標**: 後退時間積分 + 残差注入機構

#### 対象関数（4個）
1. `coeffs_and_rhs_building_Adjoint()` → `build_adjoint_system!()`
2. `assemble_A_Adjoint()` → `assemble_adjoint_matrix()`
3. `multiple_time_step_solver_Adjoint()` → `solve_adjoint!()`
4. 異常検出ラッパー群（5個） → `validators/anomaly_handlers.jl`

#### 技術課題
- **時間反転**: `for t in (nt-1):-1:1` （後退ループ）
- **残差注入**: 底面（`k=1`）での `2*(T_cal - Y_obs)` 項追加
- **初期条件**: `λ[:,:,:,nt] .= 0.0`
- **収束判定**: rtol=1e-8（DHCPより厳しい）

#### 成果物
- `julia/src/solvers/adjoint_solver.jl`
- `julia/src/utils/anomaly_detection.jl`
- `julia/test/test_adjoint_solver.jl`
- `julia/test/test_anomaly_detection.jl`

#### 成功基準
- [ ] 随伴場 `λ` がPython結果と1e-7以内
- [ ] 勾配場 `λ[:,:,end,:]` の一致（CGM勾配計算用）
- [ ] 異常検出が同一条件で同じ警告発火
- [ ] 計算時間がDHCPと同オーダー

#### リスク
- **中**: 時間反転ループの正確性、境界条件実装

---

### Phase 4: CGM最適化（Week 9-11）

**目標**: 共役勾配法本体 + 停止判定機構

#### 対象関数（1個 ← 大規模）
1. `global_CGM_time()` → `solve_cgm!()`
   - 勾配計算
   - ポラック・リビエール探索方向
   - 感度問題（`dT`計算）
   - ステップサイズ計算
   - 熱流束更新
   - Discrepancy + プラトー停止判定

#### 技術課題
- **感度問題**: DHCPソルバーの流用（`q`を摂動）
- **内積計算**: `dot(r, Sp)` で分子・分母
- **停止判定**:
  ```julia
  J < ε && maximum(abs.(residual)) ≤ σ  # Discrepancy
  mean(J_window[end-P:end]) / J_window[end-P] < η  # プラトー
  ```
- **探索方向リセット**: `mod(iter, 5) == 0` で勾配方向に戻す

#### 成果物
- `julia/src/solvers/cgm_solver.jl`
- `julia/src/solvers/stopping_criteria.jl`
- `julia/test/test_cgm_solver.jl`
- `julia/examples/single_window_cgm.jl`

#### 成功基準
- [ ] 逆解析熱流束 `q` がPython結果と相関係数 > 0.99
- [ ] 収束履歴 `J_hist` の一致（同一反復数）
- [ ] 総計算時間がPython比 **3倍以上高速化**
- [ ] メモリ使用量が同等以下

#### リスク
- **高**: CGMの数値安定性、ステップサイズ調整の微妙な差異

---

### Phase 5: スライディングウィンドウ統合（Week 12-14）

**目標**: 全時間領域計算 + 結果保存

#### 対象関数（2個）
1. `sliding_window_CGM_q_saving()` → `sliding_window_cgm!()`
2. `main()` → `main()`（エントリポイント）

#### 技術課題
- **ウィンドウ分割**: 可変長配列スライス `Y_obs[t_start:t_end, :, :]`
- **オーバーラップ平均**: `q_combined = vcat(q1[1:end-o], 0.5*(q1[end-o+1:end] .+ q2[1:o]), q2[o+1:end])`
- **初期値継承**: 前窓の `q[end-overlap:end]` を次窓の初期推定に
- **結果保存**: `NPZ.jl` or `JLD2.jl` でバイナリ保存

#### 成果物
- `julia/src/sliding_window.jl`
- `julia/src/main.jl`
- `julia/examples/full_calculation.jl`
- `julia/benchmarks/bench_full_pipeline.jl`

#### 成功基準
- [ ] 全時間領域での熱流束 `q_all` がPython結果と相関 > 0.95
- [ ] 計算時間がPython比 **5倍以上高速化**（総合）
- [ ] メモリピークが8GB以下（Python同等）
- [ ] 結果ファイルが正常保存・読込可能

#### リスク
- **中**: 大規模データのメモリ管理、ウィンドウ境界の数値処理

---

## 3. テスト戦略（TDD原則）

### 3.1 単体テスト構造

**ファイル配置**: `julia/test/`
```
test/
├── runtests.jl                  # テストエントリポイント
├── test_thermal_properties.jl   # Phase 1
├── test_validators.jl           # Phase 1
├── test_data_loader.jl          # Phase 1
├── test_dhcp_solver.jl          # Phase 2
├── test_adjoint_solver.jl       # Phase 3
├── test_cgm_solver.jl           # Phase 4
└── test_sliding_window.jl       # Phase 5
```

### 3.2 テストレベル定義

| レベル | 目的 | Python参照データ | 許容誤差 |
|--------|------|-----------------|---------|
| **L1: 数学関数** | 多項式、勾配計算 | 手計算値 | 1e-14（丸め誤差） |
| **L2: 単体関数** | 熱物性値、検証器 | Python出力保存 | 1e-12 |
| **L3: サブシステム** | DHCPソルバー単体 | 製作解 or Python | 1e-8 |
| **L4: 統合** | CGM1ウィンドウ | Python完全実行 | 相関0.99以上 |
| **L5: 性能** | ベンチマーク | Python時間比較 | 2倍以上高速 |

### 3.3 TDDワークフロー

各Phase開始時：
1. **テスト先行作成**: 期待入出力を定義
2. **Python参照データ生成**: 同一入力でPython実行し結果保存（`.npy`形式）
3. **Julia実装**: テスト失敗を確認しながら実装
4. **テスト通過**: 許容誤差内で一致確認
5. **コミット**: Phase単位でgit commit

**例**: Phase 2（DHCP）
```julia
# test/test_dhcp_solver.jl
@testset "DHCP Solver" begin
    # Python参照データ読込
    T_ref = npzread("test/data/dhcp_T_reference.npz")["T"]

    # Julia実装実行
    T_julia = solve_dhcp!(q, T0, Y_obs, params)

    # 数値比較
    @test maximum(abs.(T_julia .- T_ref)) < 1e-8
    @test cor(vec(T_julia), vec(T_ref)) > 0.9999
end
```

---

## 4. 予想課題とリスク管理

### 4.1 インデックス変換（リスク: 高）

**課題**: Pythonの0始まり → Juliaの1始まり

**対策**:
- **統一ルール**: 全コードで `1:n` 形式、`end` キーワード活用
- **境界条件注意**: `k=0`（Python底面） → `k=1`（Julia）
- **検証**: 小規模問題で全格子点の値を比較

**チェックリスト**:
- [ ] ループ範囲: `for k in 1:nk`
- [ ] 配列アクセス: `T[i, j, end]`（表面）、`T[i, j, 1]`（底面）
- [ ] スライス: `T[1:end-1, :, :]` ≠ `T[:-1, :, :]`（Python）

### 4.2 配列レイアウト（リスク: 中）

**課題**: Juliaは列優先（Fortran順序）、NumPyはデフォルト行優先

**対策**:
- **既存コード確認**: Python側が `order="F"` 多用 → 移植有利
- **明示的指定**: `Array{Float64, 3}(undef, ni, nj, nk)` でコンストラクタ
- **検証**: `vec(T)` の要素順がPython `T.ravel(order="F")` と一致

### 4.3 疎行列実装差異（リスク: 中）

**課題**: SciPyのCSR vs JuliaのCSC

**対策**:
- **CSR使用可**: `SparseMatricesCSR.jl` パッケージ
- **または転置**: Pythonで行構築 → Julia列構築に読み替え
- **性能測定**: CSC/CSR両方で実装し速度比較

**選択肢**:
```julia
# Option 1: CSC（Julia標準）
using SparseArrays
A = sparse(I, J, V, n, n)  # CSC形式

# Option 2: CSR（Python互換）
using SparseMatricesCSR
A = sparsecsr(I, J, V, n, n)  # CSR形式
```

### 4.4 CG収束特性（リスク: 高）

**課題**: IterativeSolvers.cgとscipy.sparse.linalg.cgの内部実装差異

**対策**:
- **前処理器統一**: 対角スケーリングを明示的に実装
- **収束判定厳密化**: `reltol=1e-8`, `maxiter=20000` 同値設定
- **ホットスタート**: 初期値渡し `cg!(x, A, b; initially_zero=false)`
- **デバッグ**: 各反復の残差ノルム比較

### 4.5 Numba並列化の移行（リスク: 低）

**課題**: `@njit(parallel=True)` → `Threads.@threads`

**対策**:
- **スレッド数設定**: 環境変数 `export JULIA_NUM_THREADS=8`
- **並列安全性**: ループ内で共有変数への書き込み回避
- **性能検証**: `@benchmark` でシングル/マルチコア比較

**実装例**:
```julia
using Base.Threads

function thermal_properties!(cp, k, T, cp_coeffs, k_coeffs)
    ni, nj, nk = size(T)
    @threads for i in 1:ni
        for j in 1:nj, k_ijk in 1:nk
                T_current = T[i, j, k_ijk]
                cp[i, j, k_ijk] = polyval(cp_coeffs, T_current)
                k[i, j, k_ijk] = polyval(k_coeffs, T_current)
            end
        end
    end
end
```

### 4.6 大容量データ読込（リスク: 低）

**課題**: 1.1GB `.npy` ファイルの読込

**対策**:
- **NPZ.jl使用**: `npzread("T_measure_700um_1ms.npy")`
- **メモリマップ**: 必要なら `Mmap.mmap()` で部分読込
- **検証**: ファイルサイズとメモリ使用量監視

---

## 5. 性能最適化戦略

### 5.1 Julia固有の高速化技術

| 技術 | 用途 | 適用箇所 |
|------|------|---------|
| `@inbounds` | 境界チェック省略 | 内側3重ループ |
| `@simd` | SIMD命令 | 1Dループ（係数計算） |
| `@views` | コピー回避 | 配列スライス渡し |
| `@turbo` (LoopVectorization.jl) | 自動最適化 | 密行列演算 |
| `StaticArrays.jl` | 小配列高速化 | 3x3勾配計算 |

### 5.2 メモリ最適化

**インプレース演算**: `!` 付き関数で上書き
```julia
# 悪い例（メモリコピー）
T_new = solve_dhcp(q, T0, Y_obs, params)

# 良い例（インプレース）
solve_dhcp!(T, q, T0, Y_obs, params)  # Tを直接更新
```

**事前確保**: ループ外で配列確保
```julia
cp = Array{Float64}(undef, ni, nj, nk)
k = similar(cp)
thermal_properties!(cp, k, T, cp_coeffs, k_coeffs)
```

### 5.3 型安定性の保証

**型注釈**: 関数引数・戻り値に型指定
```julia
function polyval(coeffs::Vector{Float64}, x::Float64)::Float64
    result = 0.0
    for i in eachindex(coeffs)
        result += coeffs[i] * x^(length(coeffs) - i)
    end
    return result
end
```

**型推論確認**: `@code_warntype` でチェック
```julia
@code_warntype solve_dhcp!(T, q, T0, Y_obs, params)
# 赤色の"Any"が出たら型不安定 → 修正必要
```

---

## 6. スケジュールとマイルストーン

### 6.1 全体スケジュール（14週間）

| Phase | 期間 | 主要成果物 | 累積達成率 |
|-------|------|-----------|-----------|
| **Phase 1** | Week 1-2 | 基盤関数10個 | 20% |
| **Phase 2** | Week 3-5 | DHCPソルバー | 45% |
| **Phase 3** | Week 6-8 | Adjointソルバー | 70% |
| **Phase 4** | Week 9-11 | CGM最適化 | 90% |
| **Phase 5** | Week 12-14 | スライディングウィンドウ | 100% |

### 6.2 マイルストーン定義

**M1（Week 2）**: Phase 1完了
- [ ] 熱物性値計算がPython一致
- [ ] MATファイル読込成功
- [ ] 単体テスト全通過

**M2（Week 5）**: Phase 2完了
- [ ] DHCP製作解テスト通過（誤差<1e-6）
- [ ] Python結果との温度場一致（<1e-8）
- [ ] 計算時間2倍以上高速化

**M3（Week 8）**: Phase 3完了
- [ ] Adjoint随伴場がPython一致（<1e-7）
- [ ] 異常検出システム動作確認
- [ ] DHCPと同等の計算時間

**M4（Week 11）**: Phase 4完了
- [ ] CGM収束履歴がPython一致
- [ ] 逆解析熱流束の相関0.99以上
- [ ] 総計算時間3倍以上高速化

**M5（Week 14）**: Phase 5完了
- [ ] 全時間領域計算成功
- [ ] 最終性能目標達成（5倍高速化）
- [ ] ドキュメント・ベンチマーク完備

### 6.3 週次進捗管理

**毎週金曜**:
- [ ] 単体テスト実行（`julia test/runtests.jl`）
- [ ] Python比較レポート作成
- [ ] 次週タスク洗い出し
- [ ] リスク再評価

---

## 7. 成功基準（最終検証）

### 7.1 機能要件

| 項目 | 基準 | 検証方法 |
|------|------|---------|
| **数値精度** | 温度場誤差 < 1e-8 | Python結果とのdiff |
| **逆解析精度** | 熱流束相関 > 0.99 | `cor(q_julia, q_python)` |
| **収束性** | CGM反復数が同等±10% | 収束履歴比較 |
| **データ処理** | 1.1GB読込成功 | メモリ監視ログ |
| **異常検出** | Python同一条件で同じ警告 | テストケース実行 |

### 7.2 非機能要件

| 項目 | 基準 | 測定方法 |
|------|------|---------|
| **計算速度** | Python比 **5倍以上** | `@benchmark`（BenchmarkTools.jl） |
| **メモリ使用** | ピーク8GB以下 | `@time`のメモリ表示 |
| **スレッド効率** | 8コアで6倍以上 | 並列スケーラビリティ測定 |
| **コード品質** | テストカバレッジ > 80% | `Coverage.jl` |

### 7.3 保守性要件

- [ ] 全関数にdocstring（Google Style）
- [ ] コード分割（1ファイル<500行）
- [ ] 設定外部化（`TOML`ファイル）
- [ ] CI/CD設定（GitHub Actions）

---

## 8. 依存パッケージ一覧

### 8.1 必須パッケージ

```julia
# Project.toml
[deps]
LinearAlgebra = "標準ライブラリ"
SparseArrays = "標準ライブラリ"
Statistics = "標準ライブラリ"
Sys = "標準ライブラリ"

DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
MAT = "23992714-dd62-5051-b70f-ba57cb901cac"
IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
NPZ = "15e1cf62-19b3-5cfa-8e77-841668bca605"  # .npy/.npz読書き

[extras]
Test = "標準ライブラリ"
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
```

### 8.2 オプション（高速化）

```julia
LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"  # @turbo
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"       # 固定サイズ配列
SparseMatricesCSR = "a0a7dd2c-ebf4-11e9-1f05-cf50bc540ca1"  # CSR形式
```

---

## 9. リスクマトリックス

| リスク項目 | 発生確率 | 影響度 | 優先度 | 対策 |
|-----------|---------|-------|-------|------|
| **CG収束特性差異** | 中 | 高 | ★★★ | 前処理器統一、デバッグ詳細化 |
| **インデックス変換ミス** | 高 | 高 | ★★★ | 小規模問題で全点検証 |
| **疎行列性能劣化** | 低 | 中 | ★★☆ | CSR/CSC両方実装し選択 |
| **並列化効率低下** | 低 | 中 | ★★☆ | プロファイリング、false sharing回避 |
| **メモリリーク** | 低 | 低 | ★☆☆ | Julia自動GC（通常問題なし） |

---

## 10. 参考資料

### 10.1 Julia公式ドキュメント
- [Performance Tips](https://docs.julialang.org/en/v1/manual/performance-tips/)
- [Multi-Threading](https://docs.julialang.org/en/v1/manual/multi-threading/)
- [Calling C/Fortran](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/) ※必要に応じて

### 10.2 主要パッケージドキュメント
- [SparseArrays.jl](https://docs.julialang.org/en/v1/stdlib/SparseArrays/)
- [IterativeSolvers.jl](https://iterativesolvers.julialang.org/stable/)
- [DataFrames.jl](https://dataframes.juliadata.org/stable/)
- [BenchmarkTools.jl](https://juliaci.github.io/BenchmarkTools.jl/stable/)

### 10.3 数値計算ベストプラクティス
- [SciML Style Guide](https://github.com/SciML/SciMLStyle)
- [Julia for HPC](https://github.com/johnfgibson/whyjulia)

---

## 11. まとめ

### 11.1 移植の期待効果

| 指標 | 現状（Python） | 目標（Julia） | 改善率 |
|------|--------------|-------------|-------|
| **計算時間** | 基準 | **1/5以下** | 5倍以上高速化 |
| **メモリ効率** | 8GB | 同等以下 | GC自動管理 |
| **コード行数** | 1701行単一 | 分割構造 | 保守性向上 |
| **並列効率** | Numba依存 | ネイティブ | スケーラビリティ向上 |

### 11.2 重点管理項目

1. **数値精度の保証**: 各Phaseで参照データ比較
2. **性能目標の達成**: 2週ごとにベンチマーク
3. **テストカバレッジ**: 常時80%以上維持
4. **リスク早期検出**: 週次レビューで軌道修正

### 11.3 次ステップ

**Phase 1開始準備**:
1. Julia環境構築（v1.10以上推奨）
2. 必須パッケージインストール（`Pkg.add()`）
3. Python参照データ生成スクリプト作成
4. ディレクトリ構造整備（`julia/src/`, `julia/test/`）
5. TDDテンプレート作成

---

**計画書終了**
*本計画は14週間でPython→Julia完全移植を達成し、5倍以上の性能向上を実現することを目標としています。各Phaseのマイルストーン達成を確実にし、リスク管理を徹底することで、高品質な数値計算コードを構築します。*

**総文字数**: 約11,500トークン（制限内）