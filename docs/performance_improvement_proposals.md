# 性能改善提案書

**作成日**: 2025年10月11日
**最終更新**: 2025年10月12日 21:00
**対象コミット**: 22fde2d（ベースライン）
**現在のブランチ**: tuning5
**ベースライン性能**: 1072.43秒（約17.9分、80×100×20格子、10ステップ）

## 📋 実装状況サマリー

| フェーズ | 提案 | 状態 | 実装日 | ブランチ |
|---------|------|------|-------|---------|
| Phase 0 | **配列アロケーション削減** | ✅ **完了** | 2025-10-12 | tuning4 |
| Phase 1-A | **型安定化（提案2）** | ✅ **完了** | 2025-10-12 | tuning5 |
| Phase 1-B | 初期推定値改善（提案4） | 🔜 予定 | - | - |
| Phase 1-C | 前処理改善（提案1） | 🔜 予定 | - | - |
| Phase 1-D | BLASスレッド最適化（提案7） | 🔜 予定 | - | - |
| Phase 2-A | WorkBuffers共有化（提案6） | 🔜 予定 | - | - |
| Phase 2-B | メモリアクセス最適化（提案3） | 🔄 一部完了 | 2025-10-12 | tuning4 |
| Phase 2-C | 適応的収束判定（提案5） | 🔜 予定 | - | - |
| Phase 3 | 並列化（提案8） | 📅 将来対応 | - | - |

## エグゼクティブサマリー

現在の実装における主要なボトルネックは以下の通り：

1. **Adjoint（Gradient）ソルバー**: 447.2秒（41.7%）- 最大ボトルネック
2. **CG反復回数**: 平均530-715回/ステップ - 前処理改善の余地大
3. **型不安定性**: 型推論失敗による動的ディスパッチの可能性
4. **メモリアクセスパターン**: キャッシュミス、アロケーション最適化の余地

**目標**: 総実行時間を50-70%削減（500-600秒台）

---

## 現状分析

### 実行時間内訳

| ソルバー | 時間(秒) | 割合 | 平均反復数/ステップ | 総反復数 |
|---------|---------|------|-------------------|---------|
| DHCP | 321.4 | 30.0% | 560.3 | 5,043 |
| **Gradient** | **447.2** | **41.7%** | **715.1** | **6,436** |
| Sensitivity | 300.7 | 28.0% | 530.3 | 4,773 |
| その他 | 3.0 | 0.3% | - | - |
| **合計** | **1072.4** | **100%** | - | **16,252** |

### ボトルネック特定

#### 1. Adjointソルバーの高反復回数
- ステップ9（後退時間積分の最初）: 1348回 ← 突出
- 初期残差: 6.02e-05（比較的大きい）
- その後のステップ: 586-952回（依然として高い）

#### 2. 初期残差の不均一性
- DHCP: 3.54e+06 → 4.23e-04（6桁以上の変動）
- Sensitivity: 1.40e-12 → 9.45e-12（過剰精度要求の可能性）

#### 3. CG反復効率
- 1反復あたり約0.063秒（DHCP基準）
- 16,252反復 × 0.063秒 ≈ 1024秒（総実行時間の95%）

---

## 性能改善提案

### 【提案1】前処理手法の改善（最優先・高効果）

**優先度**: ★★★★★
**期待効果**: CG反復回数30-50%削減 → 総実行時間20-30%削減
**実装難易度**: 中
**実装期間**: 2-3週間

#### 現状の問題点
- Gauss-Seidel (GS) smootherは収束が遅い
- 特にAdjoint問題で効果が不十分

#### 改善案

##### 1-A: ILU(0)前処理の導入
```julia
# julia/src/solvers/Preconditioner.jl (新規作成)
using SparseArrays, LinearAlgebra

"""
ILU(0)前処理行列の構築
"""
struct ILU0Preconditioner
  L::SparseMatrixCSC{Float64, Int64}
  U::SparseMatrixCSC{Float64, Int64}
end

function build_ilu0(A::SparseMatrixCSC)
  n = size(A, 1)
  LU = copy(A)

  for k in 1:n-1
    for i in k+1:n
      if LU[i,k] != 0
        LU[i,k] /= LU[k,k]
        for j in k+1:n
          if LU[k,j] != 0
            LU[i,j] -= LU[i,k] * LU[k,j]
          end
        end
      end
    end
  end

  L = tril(LU, -1) + I
  U = triu(LU)
  return ILU0Preconditioner(L, U)
end

function apply_preconditioner!(y, P::ILU0Preconditioner, x)
  # Lz = x を解く（前進代入）
  z = P.L \ x
  # Uy = z を解く（後退代入）
  y .= P.U \ z
end
```

##### 1-B: Jacobi前処理 + 粗格子補正（マルチグリッド的アプローチ）
```julia
"""
2レベルマルチグリッド前処理
"""
struct TwoLevelMGPreconditioner
  jacobi_diag::Vector{Float64}
  coarse_grid_solver::Function
  restriction::SparseMatrixCSC
  prolongation::SparseMatrixCSC
end

function apply_mg_preconditioner!(y, P::TwoLevelMGPreconditioner, r)
  # 1. Fine grid smoothing (weighted Jacobi)
  @. y = 0.7 * r / P.jacobi_diag

  # 2. Restriction to coarse grid
  r_coarse = P.restriction * r

  # 3. Coarse grid solve
  e_coarse = P.coarse_grid_solver(r_coarse)

  # 4. Prolongation and correction
  y .+= P.prolongation * e_coarse

  # 5. Post-smoothing
  @. y = 0.7 * (r - A * y) / P.jacobi_diag
end
```

#### 実装ステップ
1. **Week 1**: ILU(0)実装とベンチマーク
2. **Week 2**: マルチグリッド実装（オプション）
3. **Week 3**: 3ソルバー統合とテスト

#### 検証基準
- CG反復回数: 560 → 300回/ステップ以下（DHCP）
- 総実行時間: 1072秒 → 750秒以下

---

### 【提案2】型安定化とコンパイラ最適化（高効果・中実装）

**優先度**: ★★★★★
**期待効果**: 実行時間10-20%削減
**実装難易度**: 中
**実装期間**: 1-2週間

#### 現状の問題点
- 型推論失敗による動的ディスパッチの可能性
- 不要な型変換とボクシング
- コンパイラ最適化の阻害

#### 改善案

##### 2-A: 型安定性の検証と修正

```julia
# 型安定性チェックスクリプト
using JET

# julia/scripts/check_type_stability.jl
function check_solver_type_stability()
  # DHCPソルバーの型推論チェック
  @report_opt solve_dhcp!(T_init, q, work, nt, rho, cp_coeffs, k_coeffs,
                          dx, dy, Z, ΔZ, dt; rtol=1e-6, maxiter=20000)

  # PBiCGSTAB!の型推論チェック
  @report_opt PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ;
                         tol=1e-6, maxItr=20000, smoother="gs")
end
```

##### 2-B: 型アノテーションの追加

```julia
# julia/src/solvers/CommonSolver.jl
function PBiCGSTAB!(
    wk::WorkBuffers,
    Δh::AbstractArray{T,3},  # 型パラメータT追加
    dt::T,
    Z::AbstractVector{T},
    ΔZ::AbstractVector{T},
    z_range::UnitRange{Int},
    HT::AbstractArray{T,3},
    ρ::AbstractArray{T,3};
    tol::T = T(1e-6),  # 型安定な初期値
    maxItr::Int = 20000,
    smoother::String = "none",
    par::Int = 1
) where T <: AbstractFloat
    # 内部変数も型注釈
    rho::T = zero(T)
    alpha::T = zero(T)
    omega::T = zero(T)

    # ...
end
```

##### 2-C: 関数特殊化の活用

```julia
# smoother引数のコンパイル時特殊化
@inline function apply_smoother!(
    z::AbstractArray{T,3},
    r::AbstractArray{T,3},
    ::Val{:none},
    args...
) where T
    mycopy!(z, r)
end

@inline function apply_smoother!(
    z::AbstractArray{T,3},
    r::AbstractArray{T,3},
    ::Val{:gs},
    args...
) where T
    GS_smoother_in_place!(z, r, args...)
end

# 呼び出し側
function PBiCGSTAB!(...; smoother::Symbol = :none, ...)
    smoother_val = Val(smoother)
    # ...
    apply_smoother!(wk.pcg_z, wk.pcg_r, smoother_val, ...)
end
```

#### 実装ステップ
1. **Week 1**: JETによる型不安定箇所の特定
2. **Week 1-2**: 型アノテーション追加と検証
3. **Week 2**: Val型による関数特殊化

#### 検証基準
- `@code_warntype`でwarning無し
- 動的ディスパッチ箇所ゼロ（`@report_opt`）
- 実行時間5-10%改善

---

### 【提案3】メモリアクセス最適化（高効果・高難易度）

**優先度**: ★★★★☆
**期待効果**: 実行時間10-20%削減
**実装難易度**: 高
**実装期間**: 2-3週間
**実装状況**: 🔄 **一部完了**（2025-10-12）

#### 実装済み項目（Phase 0）

✅ **3-B: @simd最適化** - `tensor_dot`関数を`@simd`ループに変更完了
✅ **3-D: プリアロケーションの徹底** - CGMループ内の全バッファを事前確保完了

詳細は「実装完了レポート > Phase 0」参照。

#### 現状の問題点（残課題）
- キャッシュミスの発生（ループ順序最適化未実施）
- SIMD化の余地（一部のみ実施、他関数への展開可能）

#### 改善案

##### 3-A: ループ順序の最適化（キャッシュフレンドリー）

```julia
# 現状（列優先だが、最内ループがi）
for k in z_range
  for j in 1:nj
    for i in 1:ni
      # T[i,j,k]へのアクセス
    end
  end
end

# 改善案（最内ループをkに）
for j in 1:nj
  for i in 1:ni
    @simd ivdep for k in z_range
      # T[i,j,k]へのアクセス - キャッシュライン効率的
    end
  end
end
```

##### 3-B: @tturboマクロによるSIMD最適化

```julia
using LoopVectorization

# julia/src/utils/vector_operations.jl
function Fdot2_optimized(
    x::AbstractArray{T,3},
    y::AbstractArray{T,3},
    z_range::UnitRange{Int}
) where T
    s = zero(T)
    ni, nj, nk = size(x)

    @tturbo for k in z_range
        for j in 1:nj
            for i in 1:ni
                s += x[i,j,k] * y[i,j,k]
            end
        end
    end

    return s
end
```

##### 3-C: メモリプールとアロケーション削減

```julia
# julia/src/solvers/MemoryPool.jl (新規作成)
mutable struct MemoryPool{T}
    buffers::Dict{Symbol, Array{T,3}}
    available::Set{Symbol}

    function MemoryPool{T}(ni::Int, nj::Int, nk::Int, n_buffers::Int) where T
        buffers = Dict{Symbol, Array{T,3}}()
        available = Set{Symbol}()

        for i in 1:n_buffers
            key = Symbol("buffer_", i)
            buffers[key] = zeros(T, ni, nj, nk)
            push!(available, key)
        end

        new{T}(buffers, available)
    end
end

function acquire!(pool::MemoryPool, key::Symbol)
    if key in pool.available
        delete!(pool.available, key)
        return pool.buffers[key]
    else
        error("Buffer $key is not available")
    end
end

function release!(pool::MemoryPool, key::Symbol)
    push!(pool.available, key)
end
```

##### 3-D: プリアロケーションの徹底

```julia
# 現状: 関数内でアロケーション
function compute_residual(T_cal, Y_obs, nt)
    residual = zeros(nt, ni, nj)  # アロケーション発生
    for t in 1:nt
        residual[t, :, :] .= T_cal[:, :, nk, t] .- Y_obs[:, :, t]
    end
    return residual
end

# 改善案: プリアロケーション済みバッファ使用
function compute_residual!(residual, T_cal, Y_obs, nt, ni, nj, nk)
    @tturbo for t in 1:nt
        for j in 1:nj
            for i in 1:ni
                residual[t, i, j] = T_cal[i, j, nk, t] - Y_obs[i, j, t]
            end
        end
    end
end
```

#### 実装ステップ
1. **Week 1**: プロファイリングでホットスポット特定
2. **Week 1-2**: ループ順序最適化とSIMD化
3. **Week 2-3**: メモリプール実装とアロケーション削減

#### 検証基準
- アロケーション回数: 現状の50%以下
- キャッシュミス率: 20%改善（`perf stat`使用）
- 実行時間10%以上改善

---

### 【提案4】初期推定値の改善（中効果・低実装）

**優先度**: ★★★★☆
**期待効果**: CG反復回数15-25%削減
**実装難易度**: 低
**実装期間**: 1週間

#### 改善案

##### 4-A: 前ステップ解の利用

```julia
# julia/src/solvers/DHCPSolver.jl
function solve_dhcp!(T_init, q, work, nt, ...; use_previous_solution=true)
    T_all = zeros(ni, nj, nk, nt)
    T_all[:, :, :, 1] .= T_init

    # 前ステップの解を次ステップの初期値として使用
    for t in 2:nt
        if use_previous_solution
            # T(t)の初期推定値 = T(t-1)
            work.Δh .= T_all[:, :, :, t-1]
        else
            work.Δh .= 0.0
        end

        # PBiCGSTAB!呼び出し（work.Δhが初期値兼解ベクトル）
        isconverged, itr, res0 = PBiCGSTAB!(work, ...)

        T_all[:, :, :, t] .= work.Δh
    end

    return T_all, iter_counts
end
```

##### 4-B: 線形外挿による予測初期値

```julia
# 2次精度の線形外挿
if t >= 3 && use_extrapolation
    # T(t)の初期推定 = 2*T(t-1) - T(t-2)
    @. work.Δh = 2 * T_all[:, :, :, t-1] - T_all[:, :, :, t-2]
elseif t >= 2
    work.Δh .= T_all[:, :, :, t-1]
else
    work.Δh .= 0.0
end
```

##### 4-C: Adjointソルバーの特別対応

```julia
# Adjoint問題の初期値（後退時間積分）
function compute_gradient!(T_cal, Y_obs, work, ...)
    λ = zeros(ni, nj, nk, nt)

    # t=ntから逆順に解く
    for t in nt:-1:2
        if t == nt
            # 最終ステップ（後退の最初）: 残差の負を初期値
            residual = Y_obs[:, :, t] .- T_cal[:, :, nk, t]
            work.Δh[:, :, nk] .= -residual
            work.Δh[:, :, 1:nk-1] .= 0.0
        else
            # 次ステップ（時間的に後）の解を初期値
            work.Δh .= λ[:, :, :, t+1]
        end

        # Adjoint方程式を解く
        isconverged, itr, res0 = PBiCGSTAB!(work, ...)

        λ[:, :, :, t] .= work.Δh
    end

    return grad, grad_iters
end
```

#### 実装ステップ
1. **Day 1-2**: DHCPソルバーに前ステップ解初期値実装
2. **Day 3-4**: 線形外挿実装とベンチマーク
3. **Day 5-7**: Adjointソルバー特別対応とテスト

#### 検証基準
- DHCP反復回数: 560 → 450回/ステップ以下
- Adjoint反復回数: 715 → 550回/ステップ以下
- ステップ9（Adjoint）: 1348 → 800回以下

---

### 【提案5】適応的収束判定（小効果・低実装）

**優先度**: ★★★☆☆
**期待効果**: 実行時間5-10%削減
**実装難易度**: 低
**実装期間**: 3日

#### 現状の問題点
- Sensitivity問題で初期残差1e-12に対してrtol=1e-8は過剰
- 一律の許容誤差では効率が悪い

#### 改善案

```julia
# julia/src/solvers/CommonSolver.jl
function adaptive_tolerance(
    res0::T,
    base_rtol::T,
    min_rtol::T = T(1e-10),
    max_rtol::T = T(1e-4)
) where T
    # 初期残差が小さい場合は許容誤差を緩和
    if res0 < min_rtol
        # 絶対残差基準に切り替え
        return max(min_rtol, base_rtol * sqrt(res0))
    else
        return clamp(base_rtol, min_rtol, max_rtol)
    end
end

function PBiCGSTAB!(...; tol=1e-6, adaptive=true, ...)
    res0 = sqrt(Fdot2(wk.pcg_r, wk.pcg_r, z_range, par))

    effective_tol = adaptive ? adaptive_tolerance(res0, tol) : tol

    for itr in 1:maxItr
        res = sqrt(Fdot2(wk.pcg_r, wk.pcg_r, z_range, par))

        if res < effective_tol * res0
            isconverged = true
            break
        end

        # ...
    end
end
```

#### 実装ステップ
1. **Day 1**: adaptive_tolerance実装
2. **Day 2**: 3ソルバー統合
3. **Day 3**: ベンチマークと調整

---

### 【提案6】WorkBuffers内部配列の共有化（中効果・中実装）

**優先度**: ★★★☆☆
**期待効果**: メモリ使用量10-20%削減、キャッシュ効率向上
**実装難易度**: 中
**実装期間**: 1-2週間
**提案元**: Codex（cgm_allocation_reduction_report.md）

#### 現状の問題点
- 各ソルバー（DHCP/Adjoint/Sensitivity）が独立した`WorkBuffers`を使用
- 同時に複数のWorkBuffersがメモリに存在
- 内部配列（`pcg_r`, `pcg_z`, `pcg_p`など）の重複

#### 改善案

```julia
# julia/src/solvers/SharedWorkBuffers.jl (新規作成)
"""
複数ソルバー間で共有可能なWorkBuffers
"""
mutable struct SharedWorkBuffers
  # 共有可能なバッファ
  shared_temp1::Array{Float64,3}
  shared_temp2::Array{Float64,3}
  shared_temp3::Array{Float64,3}

  # ソルバー別ビュー
  dhcp_buffers::WorkBuffers
  adjoint_buffers::WorkBuffers
  sensitivity_buffers::WorkBuffers

  function SharedWorkBuffers(ni::Int, nj::Int, nk::Int)
    # 共有バッファを1回だけ確保
    temp1 = zeros(Float64, ni, nj, nk)
    temp2 = zeros(Float64, ni, nj, nk)
    temp3 = zeros(Float64, ni, nj, nk)

    # 各ソルバーは共有バッファへのビューを使用
    dhcp = WorkBuffers(ni, nj, nk, shared_arrays=(temp1, temp2, temp3))
    adjoint = WorkBuffers(ni, nj, nk, shared_arrays=(temp1, temp2, temp3))
    sensitivity = WorkBuffers(ni, nj, nk, shared_arrays=(temp1, temp2, temp3))

    new(temp1, temp2, temp3, dhcp, adjoint, sensitivity)
  end
end
```

#### 実装ステップ
1. **Week 1**: SharedWorkBuffers設計と実装
2. **Week 1-2**: 各ソルバーの統合とテスト
3. **Week 2**: ベンチマークと効果測定

#### 検証基準
- メモリ使用量: 現状の70-80%以下
- 実行時間: 現状維持または改善
- 全テスト通過

---

### 【提案7】BLASスレッド設定の最適化（小効果・低実装）

**優先度**: ★★☆☆☆
**期待効果**: 実行時間5-10%削減（環境依存）
**実装難易度**: 低
**実装期間**: 2-3日
**提案元**: Codex（cgm_allocation_reduction_report.md）

#### 現状の問題点
- BLASスレッド数がデフォルト設定（環境依存）
- 過剰なスレッド数でオーバーヘッド増加の可能性
- ユーザーが手動で`BLAS.set_num_threads()`を設定する必要

#### 改善案

```julia
# julia/src/utils/config.jl
using LinearAlgebra

"""
実行環境に応じた最適なBLASスレッド数を自動設定
"""
function configure_blas_threads()
  # 利用可能なCPUコア数を取得
  available_cores = Sys.CPU_THREADS

  # 最適なスレッド数を計算
  # - 小規模問題: 2-4スレッド（オーバーヘッド削減）
  # - 大規模問題: コア数の50-75%（HT考慮）
  optimal_threads = if available_cores <= 4
    min(2, available_cores)
  else
    div(available_cores * 3, 4)  # 75%使用
  end

  BLAS.set_num_threads(optimal_threads)

  @info "BLAS threads configured" threads=optimal_threads available_cores

  return optimal_threads
end

# モジュールロード時に自動実行
__init__() = configure_blas_threads()
```

#### 実装ステップ
1. **Day 1**: 自動設定関数実装
2. **Day 2**: 異なるスレッド数でベンチマーク
3. **Day 3**: 最適値の決定とドキュメント化

#### 検証基準
- 小規模問題（10×10×10）: 2-4スレッドで最速
- 大規模問題（80×100×20）: 6-8スレッドで最速
- 環境変数`JULIA_NUM_THREADS`との整合性確保

---

### 【提案8】並列化（超高効果・高実装）

**優先度**: ★★★☆☆（将来対応）
**期待効果**: 実行時間50-90%削減（環境依存）
**実装難易度**: 高
**実装期間**: 4-8週間

#### 改善案

##### 6-A: CPUマルチスレッド並列化

```julia
using ThreadsX

# ベクトル内積の並列化
function Fdot2_parallel(x, y, z_range, par)
    ni, nj, nk = size(x)

    result = ThreadsX.sum(z_range) do k
        s = 0.0
        @simd for j in 1:nj
            for i in 1:ni
                s += x[i,j,k] * y[i,j,k]
            end
        end
        s
    end

    return result
end
```

**期待効果**:
- 8コアで3-5倍高速化
- 実行時間: 1072秒 → 250-350秒

##### 6-B: GPU実装（CUDA.jl）

```julia
using CUDA

# GPU版ソルバー
function solve_dhcp_gpu!(T_init_gpu, q_gpu, work_gpu, ...)
    # データはCuArray
    T_all = CUDA.zeros(Float32, ni, nj, nk, nt)

    for t in 2:nt
        # GPUカーネル実行
        @cuda threads=(16,16,4) blocks=... heat_conduction_kernel!(...)
    end

    return T_all
end
```

**期待効果**:
- RTX 3090で20-40倍高速化
- 実行時間: 1072秒 → 25-50秒

---

## 実装ロードマップ

### ✅ Phase 0（配列アロケーション削減、完了）
- **実装日**: 2025年10月12日
- **実装者**: Codex
- **内容**: CGMループのインプレース化、バッファ再利用、`tensor_dot`最適化
- **結果**: テスト505件全通過、後方互換性維持

**次のアクション**: ベンチマーク実施（BenchmarkToolsによる定量評価）

---

### Phase 1（優先実装、1ヶ月）
1. **Week 1**: 型安定化（提案2）
   - 型推論チェック
   - 型アノテーション追加
   - 検証とベンチマーク

2. **Week 2**: 初期推定値改善（提案4）
   - 前ステップ解利用
   - Adjoint特別対応
   - 性能測定

3. **Week 3-4**: 前処理改善（提案1）
   - ILU(0)実装
   - 3ソルバー統合
   - 包括的テスト

**目標**: 総実行時間600-700秒（30-40%削減）

### Phase 2（メモリ最適化、1ヶ月）
1. **Week 5-6**: メモリアクセス最適化（提案3）
   - プロファイリング
   - ループ順序最適化
   - SIMD化

2. **Week 7**: アロケーション削減
   - メモリプール実装
   - プリアロケーション徹底

3. **Week 8**: 適応的収束判定（提案5）

**目標**: 総実行時間400-500秒（50-60%削減）

### Phase 3（並列化、2-3ヶ月、将来対応）
1. マルチスレッド並列化
2. GPU実装検討
3. 分散並列計算（MPI）

**目標**: 総実行時間100秒以下

---

## 性能測定計画

### ベンチマーク環境
- **ハードウェア**: 現在の実行環境を継続
- **測定条件**: 80×100×20格子、10ステップ、CGM反復1回
- **測定項目**:
  - 総実行時間
  - 各ソルバー実行時間
  - CG反復回数
  - メモリアロケーション回数・量
  - キャッシュミス率（可能であれば）

### 比較基準
- **ベースライン**: 22fde2d（1072.43秒）
- **各Phase完了時**: 前Phase比較 + ベースライン比較
- **最終目標**: 500-600秒（50-60%削減）

### 検証基準
- **正確性**: Python版との精度維持（相対誤差2%以内）
- **安定性**: 全テスト505項目合格維持
- **再現性**: 3回測定の標準偏差5%以内

---

## リスクと対策

### リスク1: 前処理改善が効かない
- **対策**: 複数の前処理手法を並行検証
- **代替案**: マルチグリッド法の導入

### リスク2: 型安定化で性能向上が見られない
- **対策**: プロファイリングで真のボトルネック再特定
- **代替案**: 低レベル最適化（SIMD、キャッシュ）を優先

### リスク3: メモリ最適化で複雑化
- **対策**: 段階的実装、各ステップでテスト
- **代替案**: 効果の高い部分のみ適用

### リスク4: 並列化でデバッグ困難化
- **対策**: 逐次版を保持、並列版は別モジュール
- **代替案**: 並列化は将来対応として保留

---

## まとめ

### 推奨実装順序

#### ✅ 完了済み
0. **配列アロケーション削減**（Phase 0）: 2025-10-12完了

#### 次のフェーズ
1. **型安定化**（提案2）: 低リスク・中効果、基盤整備
2. **初期推定値改善**（提案4）: 低リスク・中効果、早期成果
3. **前処理改善**（提案1）: 中リスク・高効果、最大の改善
4. **BLASスレッド最適化**（提案7）: 低リスク・小効果、簡単実装
5. **WorkBuffers共有化**（提案6）: 中リスク・中効果、メモリ効率
6. **メモリ最適化**（提案3残り）: 高リスク・中効果、専門性要
7. **適応的収束判定**（提案5）: 低リスク・小効果、補助的
8. **並列化**（提案8）: 高リスク・超高効果、将来対応

### 期待される最終性能
- **Phase 1完了後**: 600-700秒（30-40%削減）
- **Phase 2完了後**: 400-500秒（50-60%削減）
- **Phase 3完了後**: 100秒以下（90%以上削減、環境依存）

### 次のアクション
1. 提案内容のレビューと優先順位確認
2. Phase 1実装ブランチ作成（`performance-opt-phase1`）
3. 型安定性チェックスクリプト作成
4. 実装開始

---

## 🎉 実装完了レポート

### ✅ Phase 0: 配列アロケーション削減（2025-10-12実装完了）

**実装者**: Codex (GPT-5ベース)
**実装ブランチ**: tuning4
**レポート**: `docs/reports/benchmarks/cgm_allocation_reduction_report.md`

#### 実装内容

1. **CGMループのインプレース化**
   - 残差・勾配・探索方向用バッファを事前確保（`CGMSolver.jl:272`）
   - ビュー演算と`@.`マクロによるインプレース演算に変更
   - 一時配列生成を大幅に削減

2. **ソルバーAPIのバッファ再利用対応**
   - DHCP/Adjoint/Sensitivityソルバーに外部バッファ受け取り機能追加
   - キーワード引数: `T_buffer`, `lambda_buffer`, `iter_buffer`
   - バッファサイズバリデーション（`ArgumentError`による安全性確保）

3. **テンソル内積の最適化**
   ```julia
   # 変更前
   function tensor_dot(a::Array{T}, b::Array{T}) where T <: Real
     return Float64(dot(vec(a), vec(b)))  # vec()でコピー発生
   end

   # 変更後
   function tensor_dot(a::AbstractArray{T1}, b::AbstractArray{T2}) where {T1 <: Real, T2 <: Real}
     acc = zero(promote_type(T1, T2))
     @inbounds @simd for idx in eachindex(a, b)
       acc += a[idx] * b[idx]  # コピーなし、SIMD最適化
     end
     return Float64(acc)
   end
   ```

4. **勾配計算のバッファ再活用**
   - `compute_gradient!`で外部バッファを受け取り可能に
   - CGMループ内での繰り返しアロケーションを回避

#### 変更ファイル

- `julia/src/solvers/CGMSolver.jl`: CGMループ最適化
- `julia/src/solvers/DHCPSolver.jl`: バッファ再利用対応
- `julia/src/solvers/AdjointSolver.jl`: バッファ再利用対応
- `julia/src/solvers/SensitivitySolver.jl`: バッファ再利用対応

#### テスト結果

- ✅ **505件のテスト全て通過**
- ✅ 既存機能の後方互換性維持（デフォルトでは新規バッファ確保）
- ✅ Julia 1.12.0で検証済み

#### 期待効果

- **GC負荷削減**: CGM外層ループでのアロケーションが1回のみ
- **メモリ効率向上**: スライディングウィンドウ計算での大幅なメモリ削減
- **キャッシュ効率向上**: バッファ再利用によるキャッシュヒット率向上

#### 次のアクション

1. ✅ テスト実行（完了）
2. ⏳ ベンチマーク実施（BenchmarkTools使用）
   - アロケーション数の定量化
   - 実行時間の測定
   - ベースライン（22fde2d）との比較
3. ⏳ 性能改善レポート作成

---

**付録**: 詳細な実装例とベンチマークコードは各提案の実装時に作成
