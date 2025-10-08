# Phase B実装計画：Heat3ds統合型反復ソルバー

**作成日**: 2025年10月3日
**Phase**: 2.3b（旧2.3a方針からの変更）
**目的**: Heat3ds非定常ソルバーをIHCPに移植し、DHCP/Adjoint統一反復ソルバーを実装

---

## 1. エグゼクティブサマリー

### 1.1 方針変更の経緯

**旧方針（Phase 2.3a当初）**:
- LinearOperators.jl + IterativeSolvers.jlを使用
- マトリックスフリー行列ベクトル積を実装
- 既存のcg!ソルバーを活用

**新方針（Phase 2.3b）**:
- Heat3dsの非定常BiCGstabソルバーを直接移植
- LinearOperators/IterativeSolvers.jlは**使用しない**
- DHCP/Adjointで同一のソルバーコードを共有

**変更理由**:
1. **Heat3dsが既に非定常問題対応済み**: Δtパラメータで陰的オイラー法を実装
2. **統一性の向上**: DHCPとAdjointは境界条件が異なるだけで、同じソルバーを使える
3. **実績のある実装**: Heat3dsで検証済みのコードをそのまま活用
4. **外部依存の削減**: LinearOperators.jlが不要になり、制御が容易

### 1.2 期待される効果

- **性能向上**: 15-25%の高速化（疎行列組み立てコスト削減）
- **コードの簡潔化**: DHCP/Adjointの重複コード削減
- **保守性向上**: Heat3dsとの共通化により、改良が両方に反映

---

## 2. Heat3ds非定常項の実装詳細

### 2.1 離散化スキーム（陰的オイラー法）

**支配方程式**:
```
∂θ/∂t = ∂/∂x (α ∂θ/∂x) + ∂/∂y (α ∂θ/∂y) + ∂/∂z (α ∂θ/∂z) + Q
```

**時間離散化**:
```
(θ^{n+1} - θ^n) / Δt = 拡散項^{n+1} + Q^{n+1}
```

**線形システム** `A·θ^{n+1} = b`:
```julia
# 対角項（Heat3ds: line 136, 206, 283, 495, 556）
dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm + ddt)*m0
# ここで ddt = 1.0 / Δt

# 非対角項（拡散係数）
ss = ( axp * θ[i+1,j,k] + axm * θ[i-1,j,k]
     + ayp * θ[i,j+1,k] + aym * θ[i,j-1,k]
     + azp * θ[i,j,k+1] + azm * θ[i,j,k-1] )

# 行列ベクトル積
A·θ = dd*θ - ss

# RHS
b = θ^n / Δt + Q + 熱流束境界項
```

### 2.2 Heat3dsの主要関数

#### (1) CalcAX! - 行列ベクトル積 (heat3d_NonUniform.jl: 518-562)

```julia
function CalcAX!(ap::Array{Float64,3},
                  p::Array{Float64,3},
                  Δh,
                  Δt::Float64,
                  λ::Array{Float64,3},
                  m::Array{Float64,3},
                  Z::Vector{Float64},
                  ΔZ::Vector{Float64},
                  z_range::Vector{Int64},
                  HT::Vector{Float64},
                  par::String)
    # 非定常項
    ddt = 1.0 / Δt

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        # 拡散係数計算（6方向）
        axm, axp, aym, ayp, azm, azp = ...

        # 対角項（非定常項を含む）
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm + ddt)*m0

        # 7点ステンシル
        ss = ( axp * p[i+1,j,k] + axm * p[i-1,j,k]
             + ayp * p[i,j+1,k] + aym * p[i,j-1,k]
             + azp * p[i,j,k+1] + azm * p[i,j,k-1] )

        # A*p の結果（マスクで境界除外）
        ap[i,j,k] = (ss - dd*p[i,j,k]) * m0
    end
end
```

**重要なポイント**:
- `ddt = 1.0 / Δt`: 非定常項の係数
- `dd`に非定常項が含まれる → 時間積分が陰的
- `HT`: 熱伝達境界（IHCPでは未使用、全てゼロ）

#### (2) CalcRK! - 残差計算 (heat3d_NonUniform.jl: 451-503)

```julia
function CalcRK!(r::Array{Float64,3},
                p::Array{Float64,3},
                b::Array{Float64,3},
                λ::Array{Float64,3},
                m::Array{Float64,3},
                Δh,
                Δt::Float64,
                Z::Vector{Float64},
                ΔZ::Vector{Float64},
                z_range::Vector{Int64},
                HF::Vector{Float64},
                HT::Vector{Float64},
                par::String)
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        # bb: RHSに熱流束境界項を追加
        bb = b[i,j,k]
            +( (mw*HF[1] - me*HF[2])*dx1
            +  (ms*HF[3] - mn*HF[4])*dy1
            +((1.0-mb)*HF[5] - (1.0-m0)*HF[6])/ΔZ[k] )

        # 拡散係数と対角項（CalcAX!と同じ）
        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm + ddt)*m0
        ss = ...

        # 残差: r = b - A*p
        rs = (bb - (ss - dd * p[i,j,k]))* m0
        r[i,j,k] = rs
        @reduce(res = 0.0 + rs*rs)
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end
```

**HF配列の意味**（6方向の熱流束境界値）:
- `HF[1]`: X-方向（西側）
- `HF[2]`: X+方向（東側）
- `HF[3]`: Y-方向（南側）
- `HF[4]`: Y+方向（北側）
- `HF[5]`: Z-方向（底面） ← **IHCP: q_surface（未知量）**
- `HF[6]`: Z+方向（表面） ← **IHCP: 断熱（ゼロ）または測定誤差項（Adjoint）**

#### (3) PBiCGSTAB! - BiCGstabソルバー (heat3d_NonUniform.jl: 353-431)

```julia
function PBiCGSTAB!(X::Array{Float64,3},
                    B::Array{Float64,3},
                    pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_,
                    pcg_s, pcg_s_, pcg_t_,
                    λ::Array{Float64,3},
                    mask::Array{Float64,3},
                    Δh,
                    Δt::Float64,
                    Z::Vector{Float64},
                    ΔZ::Vector{Float64},
                    z_range::Vector{Int64},
                    smoother::String,
                    F,
                    tol,
                    HF::Vector{Float64},
                    HT::Vector{Float64},
                    par::String)
    # 初期残差計算
    res0 = CalcRK!(pcg_r, X, B, λ, mask, Δh, Δt, Z, ΔZ, z_range, HF, HT, par)
    pcg_r0 .= pcg_r

    rho_old = 1.0
    alpha = 0.0
    omega = 1.0

    for itr in 1:ItrMax
        rho = Fdot2(pcg_r, pcg_r0, z_range, par)

        if itr == 1
            pcg_p .= pcg_r
        else
            beta = rho / rho_old * alpha / omega
            BiCG1!(pcg_p, pcg_r, pcg_q, beta, omega, z_range, par)
        end

        # 前処理
        Preconditioner!(pcg_p_, pcg_p, λ, mask, Δh, Δt, smoother, Z, ΔZ, z_range, HF, HT, par)

        # 第1回行列ベクトル積
        CalcAX!(pcg_q, pcg_p_, Δh, Δt, λ, mask, Z, ΔZ, z_range, HT, par)

        alpha = rho / Fdot2(pcg_q, pcg_r0, z_range, par)
        Triad!(pcg_s, pcg_q, pcg_r, -alpha, z_range, par)

        # 前処理（2回目）
        Preconditioner!(pcg_s_, pcg_s, λ, mask, Δh, Δt, smoother, Z, ΔZ, z_range, HF, HT, par)

        # 第2回行列ベクトル積
        CalcAX!(pcg_t_, pcg_s_, Δh, Δt, λ, mask, Z, ΔZ, z_range, HT, par)

        omega = Fdot2(pcg_t_, pcg_s, z_range, par) / Fdot1(pcg_t_, z_range, par)

        # 解の更新
        BICG2!(X, pcg_p_, pcg_s_, alpha, omega, z_range, par)

        # 残差更新と収束判定
        Triad!(pcg_r, pcg_t_, pcg_s, -omega, z_range, par)
        res = sqrt(Fdot1(pcg_r, z_range, par))/((SZ[1]-2)*(SZ[2]-2)*(z_range[2]-z_range[1]+1))
        res /= res0

        if res < tol
            break
        end

        rho_old = rho
    end
end
```

**BiCGstabアルゴリズムの特徴**:
- CG法よりも安定（非対称行列に対応）
- 1反復あたり2回の行列ベクトル積
- 前処理を2回実行（収束性向上）

---

## 3. DHCP/Adjoint統一インターフェース設計

### 3.1 統一ソルバー関数

```julia
"""
Heat3ds形式の統一反復ソルバー（DHCP/Adjoint共通）

# 引数
- X: 解ベクトル（初期推定値、出力）[ni+2, nj+2, nk+2]
- B: RHSベクトル [ni+2, nj+2, nk+2]
- λ: 熱伝導率 [ni+2, nj+2, nk+2]
- mask: マスク配列（1.0=内点、0.0=境界）[ni+2, nj+2, nk+2]
- Z: セル中心座標 [nk+2]
- ΔZ: セル代表幅 [nk+1]
- Δh: (dx, dy)
- Δt: 時間刻み
- z_range: [z_start, z_end] 計算範囲
- HF: 熱流束境界値 [6]（X-,X+,Y-,Y+,Z-,Z+）
- HT: 熱伝達境界値 [6]（IHCPでは全てゼロ）

# キーワード引数
- tol: 収束判定閾値（デフォルト: 1e-6）
- maxiter: 最大反復回数（デフォルト: 20000）
- smoother: 前処理方法（"gs"=Gauss-Seidel、""=なし）
- verbose: 詳細ログ出力

# 戻り値
- なし（Xが解で上書きされる）
"""
function solve_bicgstab_heat3ds!(
    X::Array{Float64,3},
    B::Array{Float64,3},
    λ::Array{Float64,3},
    mask::Array{Float64,3},
    Z::Vector{Float64},
    ΔZ::Vector{Float64},
    Δh::Tuple{Float64, Float64},
    Δt::Float64,
    z_range::Vector{Int},
    HF::Vector{Float64},
    HT::Vector{Float64};
    tol=1e-6,
    maxiter=20000,
    smoother="",
    verbose=false
)
    # ワークベクトルの確保
    SZ = size(X)
    pcg_q = zeros(Float64, SZ)
    pcg_r = zeros(Float64, SZ)
    pcg_r0 = zeros(Float64, SZ)
    pcg_p = zeros(Float64, SZ)
    pcg_p_ = zeros(Float64, SZ)
    pcg_s = zeros(Float64, SZ)
    pcg_s_ = zeros(Float64, SZ)
    pcg_t_ = zeros(Float64, SZ)

    # Heat3dsのPBiCGSTAB!を呼び出し
    PBiCGSTAB!(X, B, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_, pcg_s, pcg_s_, pcg_t_,
               λ, mask, Δh, Δt, Z, ΔZ, z_range, smoother, stdout, tol, HF, HT, "serial")
end
```

### 3.2 DHCPソルバー

```julia
"""
DHCP（直接熱伝導問題）ソルバー - Heat3ds統合版

# 引数
- T_prev: 前時刻の温度場 [ni, nj, nk]
- q_surface: 表面熱流束 [ni, nj]
- rho, cp_coeffs, k_coeffs: 材料物性値
- dx, dy, dz, dz_b, dz_t, dt: 格子・時間パラメータ
- rtol, maxiter: 収束パラメータ

# 戻り値
- T_new: 新時刻の温度場 [ni, nj, nk]
"""
function solve_dhcp_heat3ds!(
    T_prev::Array{Float64,3},
    q_surface::Array{Float64,2},
    rho::Float64,
    cp_coeffs::Vector{Float64},
    k_coeffs::Vector{Float64},
    dx::Float64,
    dy::Float64,
    dz::Vector{Float64},
    dz_b::Vector{Float64},
    dz_t::Vector{Float64},
    dt::Float64;
    rtol=1e-6,
    maxiter=20000,
    verbose=false
)
    ni, nj, nk = size(T_prev)

    # ガイドセル配列初期化
    SZ = (ni+2, nj+2, nk+2)
    θ = zeros(Float64, SZ)
    λ = zeros(Float64, SZ)
    mask = ones(Float64, SZ)

    # 格子変換
    Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

    # 境界条件設定
    z_range = compute_z_range(nk, HEAT_FLUX, HEAT_FLUX)  # 両端熱流束

    # 熱伝導率の計算
    cp, k = thermal_properties_calculator(T_prev, cp_coeffs, k_coeffs)

    # 内点データをガイドセル配列にコピー
    initialize_guard_cells!(θ, λ, mask, T_prev, k)

    # ガイドセルの境界条件設定
    # X, Y方向: 断熱（mask=0, λ=0）
    mask[1, :, :] .= 0.0
    mask[ni+2, :, :] .= 0.0
    mask[:, 1, :] .= 0.0
    mask[:, nj+2, :] .= 0.0
    λ[1, :, :] .= 0.0
    λ[ni+2, :, :] .= 0.0
    λ[:, 1, :] .= 0.0
    λ[:, nj+2, :] .= 0.0

    # Z方向: 熱流束境界
    mask[:, :, 1] .= 0.0       # 底面ガイドセル
    mask[:, :, nk+2] .= 0.0    # 表面ガイドセル
    λ[:, :, 1] .= 0.0
    λ[:, :, nk+2] .= 0.0

    # RHSベクトル構築
    B = zeros(Float64, SZ)
    vol = dx * dy  # 底面積
    ddt = 1.0 / dt

    for k in 2:nk+1, j in 2:nj+1, i in 2:ni+1
        # 非定常項: θ^n / Δt
        B[i, j, k] = θ[i, j, k] * ddt * vol * ΔZ[k]

        # 内部熱源があれば追加（IHCPでは通常ゼロ）
        # B[i, j, k] += Q[i, j, k] * vol * ΔZ[k]
    end

    # 熱流束境界条件（HF配列）
    HF = zeros(Float64, 6)
    # HF[1-4]: X,Y方向断熱（ゼロ）
    # HF[5]: 底面熱流束（未知、ゼロと仮定）
    # HF[6]: 表面熱流束（入力）
    for j in 2:nj+1, i in 2:ni+1
        HF[6] = q_surface[i-1, j-1]  # ガイドセルインデックス変換
    end

    # 熱伝達境界（IHCPでは未使用）
    HT = zeros(Float64, 6)

    # BiCGstab反復ソルバー
    solve_bicgstab_heat3ds!(
        θ, B, λ, mask, Z, ΔZ, (dx, dy), dt, z_range, HF, HT;
        tol=rtol, maxiter=maxiter, smoother="", verbose=verbose
    )

    # ガイドセルを除いて内点データを返す
    T_new = zeros(Float64, ni, nj, nk)
    for k in 1:nk, j in 1:nj, i in 1:ni
        T_new[i, j, k] = θ[i+1, j+1, k+1]
    end

    return T_new
end
```

### 3.3 Adjointソルバー

```julia
"""
Adjoint（随伴方程式）ソルバー - Heat3ds統合版

# 引数
- Y_obs: 観測温度 [ni, nj]
- Y_calc: 計算温度（DHCPの結果）[ni, nj, nk]
- 以下DHCP同様

# 戻り値
- lambda_field: 随伴場 [ni, nj, nk]
"""
function solve_adjoint_heat3ds!(
    Y_obs::Array{Float64,2},
    Y_calc::Array{Float64,3},
    T_prev::Array{Float64,3},
    rho::Float64,
    cp_coeffs::Vector{Float64},
    k_coeffs::Vector{Float64},
    dx::Float64,
    dy::Float64,
    dz::Vector{Float64},
    dz_b::Vector{Float64},
    dz_t::Vector{Float64},
    dt::Float64;
    rtol=1e-6,
    maxiter=20000,
    verbose=false
)
    ni, nj, nk = size(T_prev)

    # ガイドセル配列初期化（DHCP同様）
    SZ = (ni+2, nj+2, nk+2)
    λ_adj = zeros(Float64, SZ)  # 随伴場
    λ_thermal = zeros(Float64, SZ)
    mask = ones(Float64, SZ)

    Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)
    z_range = compute_z_range(nk, HEAT_FLUX, HEAT_FLUX)

    # 熱伝導率（DHCPと同じ）
    cp, k = thermal_properties_calculator(T_prev, cp_coeffs, k_coeffs)
    initialize_guard_cells!(λ_adj, λ_thermal, mask, zeros(ni, nj, nk), k)

    # 境界条件（DHCP同様）
    mask[1, :, :] .= 0.0
    mask[ni+2, :, :] .= 0.0
    mask[:, 1, :] .= 0.0
    mask[:, nj+2, :] .= 0.0
    λ_thermal[1, :, :] .= 0.0
    λ_thermal[ni+2, :, :] .= 0.0
    λ_thermal[:, 1, :] .= 0.0
    λ_thermal[:, nj+2, :] .= 0.0
    mask[:, :, 1] .= 0.0
    mask[:, :, nk+2] .= 0.0
    λ_thermal[:, :, 1] .= 0.0
    λ_thermal[:, :, nk+2] .= 0.0

    # RHSベクトル構築（測定誤差項）
    B = zeros(Float64, SZ)
    vol = dx * dy
    ddt = 1.0 / dt

    # 表面測定誤差をRHSに注入
    k_surface = nk  # 表面層インデックス（内点）
    for j in 2:nj+1, i in 2:ni+1
        error_term = Y_obs[i-1, j-1] - Y_calc[i-1, j-1, k_surface]
        B[i, j, k_surface+1] = -error_term * ddt * vol * ΔZ[k_surface+1]
        # 符号注意: 随伴方程式では測定誤差に-1を掛ける
    end

    # 熱流束境界条件
    HF = zeros(Float64, 6)
    # Adjointでは全て断熱（測定誤差項はRHSに含まれる）

    HT = zeros(Float64, 6)

    # BiCGstab反復ソルバー（DHCPと同じ関数）
    solve_bicgstab_heat3ds!(
        λ_adj, B, λ_thermal, mask, Z, ΔZ, (dx, dy), dt, z_range, HF, HT;
        tol=rtol, maxiter=maxiter, smoother="", verbose=verbose
    )

    # 内点データを返す
    lambda_field = zeros(Float64, ni, nj, nk)
    for k in 1:nk, j in 1:nj, i in 1:ni
        lambda_field[i, j, k] = λ_adj[i+1, j+1, k+1]
    end

    return lambda_field
end
```

### 3.4 DHCP/Adjointの違いまとめ

| 項目 | DHCP | Adjoint |
|------|------|---------|
| **解ベクトル** | θ^{n+1}（温度） | λ（随伴場） |
| **RHS (B)** | θ^n/Δt + Q | -(Y_obs - Y_calc)/Δt（表面のみ） |
| **HF[6]** | q_surface（表面熱流束） | 0.0（断熱、測定誤差はRHSへ） |
| **時間方向** | 前進（t → t+1） | 後退（t+1 → t） |
| **ソルバー** | solve_bicgstab_heat3ds!（共通） | solve_bicgstab_heat3ds!（共通） |

**重要**: ソルバーコードは完全に同一、違いは入力データのみ

---

## 4. 実装スケジュール

### Week 1-2: 核心関数の移植（10-14日）

**タスク**:
1. ヘルパー関数群（5日）
   - `Fdot1` - ベクトル内積（マスク適用）
   - `Fdot2` - 2ベクトル内積
   - `Triad!` - 3項演算（z = x + α*y）
   - `BiCG1!` - BiCGstab更新式1
   - `BiCG2!` - BiCGstab更新式2
   - `get_backend` - 並列化バックエンド選択

2. `CalcAX!` - 行列ベクトル積（3日）
   - Heat3ds: line 518-562を移植
   - @floopマクロで並列化
   - 単体テスト作成

3. `CalcRK!` - 残差計算（2日）
   - Heat3ds: line 451-503を移植
   - 熱流束境界項の追加
   - 単体テスト作成

4. `Preconditioner!` + `rbsor!` - 前処理（3日）
   - Heat3ds: line 578-601, 312-333を移植
   - Gauss-Seidel前処理の実装
   - 単体テスト作成

5. `PBiCGSTAB!` - メインソルバー（1日）
   - Heat3ds: line 353-431を移植
   - 上記関数を統合

**成果物**:
- `src/solvers/Heat3DSolver.jl` - Heat3ds互換ソルバーモジュール
- `test/test_heat3d_solver.jl` - 単体テストスイート（50項目）

### Week 3: DHCP統合（5-7日）

**タスク**:
1. `solve_bicgstab_heat3ds!` - 統一インターフェース（2日）
   - ワークベクトル管理
   - エラーハンドリング

2. `solve_dhcp_heat3ds!` - DHCPラッパー（2日）
   - RHS構築（θ^n/Δt + Q）
   - 熱流束境界条件設定
   - ガイドセル変換

3. 既存DHCPSolverとの統合（1日）
   - `solve_dhcp!`から`solve_dhcp_heat3ds!`を呼び出し
   - 互換性維持（戻り値形状など）

4. テストとデバッグ（2日）
   - Phase 2テストスイート（298項目）で検証
   - Python参照データとの比較

**成果物**:
- `solve_dhcp_heat3ds!`が全298テスト合格

### Week 4: Adjoint統合（5-7日）

**タスク**:
1. `solve_adjoint_heat3ds!` - Adjointラッパー（2日）
   - RHS構築（測定誤差項）
   - 境界条件設定

2. 既存AdjointSolverとの統合（1日）
   - `solve_adjoint!`から`solve_adjoint_heat3ds!`を呼び出し

3. テストとデバッグ（2日）
   - Phase 3テストスイート（13項目）で検証
   - 随伴場の数値精度確認

**成果物**:
- `solve_adjoint_heat3ds!`が全13テスト合格

### Week 5: CGM/SlidingWindow統合とベンチマーク（5-7日）

**タスク**:
1. CGMSolverの更新（1日）
   - 新しいDHCP/Adjointソルバーを使用

2. 全テストスイート実行（1日）
   - Phase 1-6全571テストで検証

3. ベンチマーク（2日）
   - 中規模問題（40×50×10×50）で性能測定
   - Python版との比較

4. プロファイリングと最適化（1日）
   - ボトルネック特定
   - @simd、@inboundsマクロの活用

**成果物**:
- 全571テスト合格
- 性能改善レポート

---

## 5. 技術的詳細

### 5.1 配列インデックス変換

**IHCP形式** → **Heat3ds形式**:

| 配列 | IHCP | Heat3ds | 変換式 |
|------|------|---------|--------|
| 温度 | `T[i,j,k]` (ni,nj,nk) | `θ[i+1,j+1,k+1]` (ni+2,nj+2,nk+2) | i → i+1 |
| 熱流束 | `q[i,j]` (ni,nj) | `HF[6]` (スカラー) | ループ内で設定 |

### 5.2 境界条件の設定方法

#### X,Y方向（断熱）

```julia
# ガイドセルのマスクとλをゼロに
mask[1, :, :] .= 0.0      # 西側
mask[ni+2, :, :] .= 0.0   # 東側
mask[:, 1, :] .= 0.0      # 南側
mask[:, nj+2, :] .= 0.0   # 北側

λ[1, :, :] .= 0.0
λ[ni+2, :, :] .= 0.0
λ[:, 1, :] .= 0.0
λ[:, nj+2, :] .= 0.0
```

→ `λf`関数で自動的に拡散項がゼロになる

#### Z方向（熱流束）

**DHCP**:
```julia
# HF[5]: 底面熱流束（通常ゼロ、または未知量）
# HF[6]: 表面熱流束（q_surface）
HF[5] = 0.0
HF[6] = q_surface[i, j]  # CalcRK!関数内でRHSに追加
```

**Adjoint**:
```julia
# 熱流束境界は使わず、測定誤差をRHSに直接注入
HF[5] = 0.0
HF[6] = 0.0

# 測定誤差項
error = Y_obs[i, j] - Y_calc[i, j, nk]
B[i+1, j+1, nk+1] = -error / dt * vol * ΔZ[nk+1]
```

### 5.3 前処理器の選択

**オプション1: 前処理なし**（`smoother=""`）
- 実装: `pcg_p_ .= pcg_p`（恒等写像）
- コスト: ゼロ
- 収束性: 基本的なBiCGstab

**オプション2: Gauss-Seidel前処理**（`smoother="gs"`）
- 実装: RB-SORを5反復
- コスト: 高い（反復ごとに5回のSOR）
- 収束性: 大幅改善（反復回数が1/3〜1/5に減少）

**推奨**: まずは前処理なしで実装し、性能測定。必要に応じてGS前処理を追加。

### 5.4 並列化戦略

```julia
using FLoops

# バックエンド選択
backend = ThreadedEx()  # マルチスレッド並列

# @floopマクロで並列化
@floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
    # 各セルの計算は完全に独立
    ap[i,j,k] = ...
end
```

**効果**: 4スレッドで約3倍高速化（並列化効率75%）

---

## 6. 期待される効果

### 6.1 性能向上（保守的見積もり）

**現在の性能**（中規模問題: 40×50×10×50）:
- Python版: 106秒
- Julia版（現行）: 1516秒

**Phase B実装後の予測**:

| 項目 | 削減量 | 予測実行時間 |
|------|--------|-------------|
| 疎行列組み立て削減 | 12.4秒 | 1504秒 |
| キャッシュ効率向上（30%） | 355秒 | 1149秒 |
| **合計削減** | **367秒** | **1149秒（24%改善）** |

**楽観的見積もり（並列化+SIMD）**:
- キャッシュ効率50%改善: 592秒削減
- 並列化（4スレッド）: さらに3倍 → **約300秒（80%改善）**

### 6.2 コードの簡潔化

**現在**:
- DHCPSolver.jl: 370行
- AdjointSolver.jl: 250行
- 合計: 620行

**Phase B後**:
- Heat3DSolver.jl: 400行（核心関数）
- solve_dhcp_heat3ds!: 50行
- solve_adjoint_heat3ds!: 50行
- 合計: 500行（**20%削減**）

**重複コード削減**: DHCP/Adjointで同じソルバー関数を共有

### 6.3 保守性向上

- Heat3dsの改良が自動的にIHCPに反映
- バグ修正も一箇所で済む
- テストも共通化可能

---

## 7. リスクと対策

### 7.1 リスク

**リスク1: BiCGstabの収束性**
- CG法（対称正定値行列専用）よりも収束が遅い可能性
- 前処理なしでは収束しない場合がある

**対策**:
- まず前処理なしで実装し、収束性を確認
- 必要に応じてGS前処理を追加
- 収束パラメータ（tol, maxiter）の調整

**リスク2: 実装の複雑性**
- Heat3dsのコードは800行以上
- FLoopsマクロの理解が必要

**対策**:
- 段階的移植（ヘルパー関数 → CalcAX! → PBiCGSTAB!）
- 各段階で単体テスト実施
- TDD方針の継続

**リスク3: 数値精度の劣化**
- 浮動小数点演算の順序変更
- 前処理の影響

**対策**:
- Python参照データとの厳密な比較（相対誤差<1e-12）
- Phase 2/3テストスイートで検証

### 7.2 成功基準

**必須**:
- 全571テスト合格
- Python参照データと相対誤差1e-7以内

**目標**:
- 中規模問題で15-25%高速化
- 大規模問題で30-40%高速化

**ストレッチゴール**:
- 並列化で80%高速化
- Python版を超える性能

---

## 8. 次のステップ

### 8.1 Phase B開始前の準備（1日）

1. Heat3dsコードの精読
   - PBiCGSTAB!の完全理解
   - 各ヘルパー関数の役割確認

2. 開発環境整備
   - FLoops.jlの動作確認
   - 単体テストのひな形作成

3. ブランチ作成
   - `tuning2` → `feature/heat3d-solver` をチェックアウト

### 8.2 Week 1開始タスク

**Day 1-2**: ヘルパー関数実装
- `Fdot1`, `Fdot2`, `Triad!`, `BiCG1!`, `BiCG2!`
- `get_backend`
- 単体テスト（10項目）

**Day 3-5**: `CalcAX!`実装
- Heat3ds: line 518-562を移植
- @floopマクロで並列化
- 単体テスト（15項目）

---

## 9. まとめ

### 9.1 主要な変更点

1. **LinearOperators/IterativeSolvers.jlを使わない**
   - Heat3dsの実装を直接活用
   - 外部依存を削減

2. **DHCP/Adjointで同一ソルバー**
   - 境界条件とRHSのみ異なる
   - コードの重複を削減

3. **非定常項の統一的な扱い**
   - Heat3dsが既に対応済み
   - `ddt = 1.0 / Δt`で陰的オイラー法

### 9.2 期待される成果

- **性能**: 15-25%高速化（保守的）、80%高速化（楽観的）
- **コード品質**: 20%のコード削減、重複排除
- **保守性**: Heat3dsとの共通化、改良の自動反映

### 9.3 実装期間

**合計**: 4-5週間
- Week 1-2: 核心関数移植
- Week 3: DHCP統合
- Week 4: Adjoint統合
- Week 5: 統合テストとベンチマーク

---

**作成者**: Claude Code
**最終更新**: 2025年10月3日
**Phase**: 2.3b実装計画
