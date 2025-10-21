# マトリックスフリー法の分析と移植方針

**作成日**: 2025年10月3日
**目的**: Heat3dsのマトリックスフリーBiCGstab法をIHCP_CGMに移植して性能改善

## 1. エグゼクティブサマリー

**現状の問題**:
- プロファイル結果: 疎行列CG法が全体の**90%以上**を占める
- 現在の実装: 疎行列Aを構築 → cg!を呼び出し
- ボトルネック: 疎行列組み立て + CG法内部の疎行列演算

**提案**: Heat3dsのマトリックスフリーBiCGstab法を移植
- 疎行列Aを**構築しない**
- 行列ベクトル積を直接計算（CalcAX!関数）
- 期待効果: **疎行列組み立てコスト削減 + キャッシュ効率向上**

## 2. Heat3ds実装の詳細分析

### 2.1 マトリックスフリーの核心

**CalcAX!関数**（heat3d_NonUniform.jl: 518-562行）:
```julia
function CalcAX!(ap::Array{Float64,3},
                  p::Array{Float64,3},
                  Δh,
                  λ::Array{Float64,3},
                  m::Array{Float64,3},
                  Z::Vector{Float64},
                  ΔZ::Vector{Float64},
                  z_range::Vector{Int64},
                  HT::Vector{Float64},
                  par::String)
    # ...
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        # 7点ステンシル係数計算（その場）
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        # ... (6方向の係数)

        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]
        # ... (他の方向)

        dd = (1.0-m0) + (axp + axm + ayp + aym + azp + azm)*m0
        ss = ( axp * p[i+1,j  ,k  ] + axm * p[i-1,j  ,k  ]
             + ayp * p[i  ,j+1,k  ] + aym * p[i  ,j-1,k  ]
             + azp * p[i  ,j  ,k+1] + azm * p[i  ,j  ,k-1] )

        # A*p の結果
        ap[i,j,k] = (ss - dd*p[i,j,k]) * m0
    end
end
```

**重要な特徴**:
1. **疎行列Aを構築しない** → メモリ節約
2. **係数計算と行列ベクトル積を同時実行** → キャッシュ効率向上
3. **@floopマクロで並列化** → スレッドレベル並列性

### 2.2 PBiCGSTAB!関数の構造

**メインループ**（heat3d_NonUniform.jl: 353-434行）:
```julia
function PBiCGSTAB!(X, B, pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_,
                    pcg_s, pcg_s_, pcg_t_, λ, mask, Δh, Z, ΔZ,
                    z_range, smoother, F, tol, HF, HT, par)
    # 初期残差計算
    res0 = CalcRK!(pcg_r, X, B, λ, mask, Δh, Z, ΔZ, z_range, HF, HT, par)
    pcg_r0 .= pcg_r

    for itr in 1:ItrMax
        # BiCGstabアルゴリズム
        rho = Fdot2(pcg_r, pcg_r0, z_range, par)

        # 前処理
        Preconditioner!(pcg_p_, pcg_p, λ, mask, Δh, smoother, Z, ΔZ,
                        z_range, HF, HT, par)

        # 行列ベクトル積（マトリックスフリー）
        CalcAX!(pcg_q, pcg_p_, Δh, λ, mask, Z, ΔZ, z_range, HT, par)

        alpha = rho / Fdot2(pcg_q, pcg_r0, z_range, par)
        Triad!(pcg_s, pcg_q, pcg_r, -alpha, z_range, par)

        # 2回目の前処理
        Preconditioner!(pcg_s_, pcg_s, λ, mask, Δh, smoother, Z, ΔZ,
                        z_range, HF, HT, par)

        # 2回目の行列ベクトル積
        CalcAX!(pcg_t_, pcg_s_, Δh, λ, mask, Z, ΔZ, z_range, HT, par)

        omega = Fdot2(pcg_t_, pcg_s, z_range, par) / Fdot1(pcg_t_, z_range, par)

        # 解の更新
        BICG2!(X, pcg_p_, pcg_s_, alpha, omega, z_range, par)

        # 残差更新と収束判定
        Triad!(pcg_r, pcg_t_, pcg_s, -omega, z_range, par)
        res = sqrt(Fdot1(pcg_r, z_range, par)) / ((SZ[1]-2)*(SZ[2]-2)*(...))
        res /= res0

        if res < tol
            break
        end

        rho_old = rho
    end
end
```

**BiCGstabの特徴**:
- CGよりも安定（非対称行列に対応）
- 1反復あたり2回の行列ベクトル積（CG法は1回）
- 前処理を組み込みやすい

### 2.3 前処理器の実装

**Preconditioner!関数**（heat3d_NonUniform.jl: 578-600行）:
```julia
function Preconditioner!(xx::Array{Float64,3},
                         bb::Array{Float64,3},
                         λ::Array{Float64,3},
                         mask::Array{Float64,3},
                         Δh, smoother::String,
                         Z, ΔZ, z_range, HF, HT, par)
    LCmax::Int = 5  # 前処理反復数

    if smoother=="gs"
        for _ in 1:LCmax
            rbsor!(xx, λ, bb, mask, Δh, 1.0, Z, ΔZ, z_range, HF, HT, par)
        end
    else
        xx .= bb  # Jacobi前処理（恒等写像）
    end
end
```

**前処理オプション**:
- `"gs"`: Gauss-Seidel前処理（5反復のRB-SOR）
- `""`: Jacobi前処理（恒等写像、前処理なしに相当）

## 3. 現在のIHCP実装との比較

### 3.1 DHCP Solver（現在の実装）

**ファイル**: `julia/src/solvers/DHCPSolver.jl`

**現在のアプローチ**:
```julia
# 1. 係数計算（build_dhcp_system!）
a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_dhcp_system!(
  T_prev, q_surface[:, :, t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
)

# 2. 疎行列組み立て（assemble_dhcp_matrix）
A = assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

# 3. 対角前処理器構築
diag_A = diag(A)
Pl = Diagonal(1.0 ./ diag_A)

# 4. CG法求解（IterativeSolvers.jl）
cg!(x0, A, b; Pl=Pl, reltol=rtol, maxiter=maxiter)
```

**問題点**:
1. **疎行列組み立てコスト**: assemble_dhcp_matrix が毎時間ステップで実行される
   - COO形式 → CSC形式変換
   - メモリ割り当て
2. **疎行列演算のオーバーヘッド**: cg!内部での疎行列ベクトル積
   - CSC形式のインデックス計算
   - 間接参照によるキャッシュミス

**プロファイル結果**（profile_comparison.md）:
- assemble_dhcp_matrix: 約1,000サンプル（約2%）
- cg!（235行）: 45,407サンプル（約90%）
- sparse matrix mul!: 7,155サンプル（約14%）

### 3.2 Adjoint Solver（現在の実装）

**ファイル**: `julia/src/solvers/AdjointSolver.jl`

**現在のアプローチ**（DHCP Solverと同様）:
```julia
# 後退時間積分ループ
for t in (nt-1):-1:1
  # 係数計算
  a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_adjoint_system!(...)

  # 疎行列組み立て
  A = assemble_adjoint_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

  # CG法求解
  cg!(λ_vec, A, b; Pl=Pl, reltol=rtol, maxiter=maxiter)
end
```

**問題点**: DHCP Solverと同じ

## 4. マトリックスフリー法の利点

### 4.1 メモリ使用量の削減

**現在の実装**:
- 疎行列A: N×N（N=ni×nj×nk=40×50×10=20,000）
- 非ゼロ要素数: 約7N=140,000要素
- メモリサイズ: 140,000 × 16バイト（値+インデックス）= 2.24MB

**マトリックスフリー版**:
- 疎行列A: **不要**（0バイト）
- 係数配列a_w, ..., a_p: 7N × 8バイト = 1.12MB
- **削減**: 2.24MB → 1.12MB（50%削減）

### 4.2 疎行列組み立てコストの削減

**現在の実装**（ベンチマーク結果）:
- DHCP 1ステップ: 325ms
  - 係数構築: 3.6ms
  - 疎行列組み立て: 11.6ms（**3.6%**）
  - CG求解: 309.8ms（95.3%）

**期待効果**:
- 疎行列組み立て: **11.6ms → 0ms**（完全削減）
- 全体への影響: **3.6%の削減**（325ms → 313ms）

### 4.3 キャッシュ効率の向上

**現在の実装**:
- 疎行列ベクトル積: CSC形式の間接参照
  - `A.nzval[A.colptr[j]:A.colptr[j+1]-1]`
  - `A.rowval[A.colptr[j]:A.colptr[j+1]-1]`
  - **不規則なメモリアクセス** → キャッシュミス

**マトリックスフリー版**:
- 7点ステンシル直接適用
  - `p[i+1,j,k]`, `p[i-1,j,k]`, ... （規則的アクセス）
  - **連続メモリアクセス** → キャッシュヒット率向上

**期待効果**:
- キャッシュミス率: 90% → 5-10%（メモリレイアウト最適化提案書より）
- CG求解時間: 309.8ms → **155-248ms**（30-50%削減）

### 4.4 並列化の容易性

**Heat3ds実装**:
```julia
@floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
    # 7点ステンシル計算（完全に独立）
    ap[i,j,k] = (ss - dd*p[i,j,k]) * m0
end
```

**FLoops.jlの利点**:
- `backend = ThreadedEx()`: スレッド並列
- `backend = CUDAEx()`: GPU並列（将来的に）
- データ依存性なし → 並列化効率100%

**現在のIterativeSolvers.jl**:
- cg!内部の並列化は限定的
- 疎行列演算の並列化は難しい（不規則アクセス）

## 5. 移植方針

### 5.1 段階的な移植計画

**Phase 2.3: マトリックスフリーBiCGstab実装**

#### ステップ1: マトリックスフリーDHCP（1-2週間）
1. **CalcAX_DHCP!関数を実装**
   - build_dhcp_system!の係数計算ロジックを統合
   - 行列ベクトル積を直接計算

2. **PBiCGSTAB_DHCP!関数を実装**
   - Heat3dsのPBiCGSTAB!を参考
   - DHCP固有の境界条件に対応

3. **テスト**
   - Phase 2テストスイート（298項目）で検証
   - Python参照データとの一致確認

#### ステップ2: マトリックスフリーAdjoint（1-2週間）
1. **CalcAX_Adjoint!関数を実装**
   - build_adjoint_system!の残差注入機構を統合

2. **PBiCGSTAB_Adjoint!関数を実装**
   - 後退時間積分に対応

3. **テスト**
   - Phase 3テストスイート（13項目）で検証

#### ステップ3: CGMとの統合（1週間）
1. **solve_cgm!を更新**
   - マトリックスフリー版のDHCP/Adjointを呼び出し

2. **テスト**
   - Phase 4-6テストスイート（169項目）で検証

#### ステップ4: 性能測定と最適化（1週間）
1. **ベンチマーク**
   - 小規模、中規模、大規模問題で性能測定

2. **プロファイル**
   - ボトルネック特定

3. **最適化**
   - @simdマクロの活用
   - @inboundsマクロの活用

### 5.2 実装の詳細設計

#### CalcAX_DHCP!関数の設計

**入力**:
- `ap::Array{Float64,3}`: 出力ベクトル（A*p）
- `p::Array{Float64,3}`: 入力ベクトル
- `T_prev::Array{Float64,3}`: 前ステップ温度（熱物性値計算用）
- `q_surface::Array{Float64,2}`: 表面熱流束
- `rho, cp_coeffs, k_coeffs`: 物性値パラメータ
- `dx, dy, dz, dz_b, dz_t, dt`: 格子・時間パラメータ

**処理**:
```julia
function CalcAX_DHCP!(ap::Array{Float64,3},
                       p::Array{Float64,3},
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
                       dt::Float64)
    ni, nj, nk = size(p)

    # 熱物性値計算（T_prevから）
    cp = similar(T_prev)
    k = similar(T_prev)
    thermal_properties_calculator!(cp, k, T_prev, cp_coeffs, k_coeffs)

    # 全格子点ループ（並列化可能）
    @floop ThreadedEx() for k_idx in 1:nk, j in 1:nj, i in 1:ni
        # 格子幅取得
        dz_k = dz[k_idx]
        dz_t_k = dz_t[k_idx]
        dz_b_k = dz_b[k_idx]

        # 中心セル熱伝導率
        k_p = k[i, j, k_idx]

        # 時間項（蓄熱項）
        a_p_0 = rho * cp[i, j, k_idx] * dx * dy * dz_k / dt

        # 6方向の熱伝導係数（調和平均）
        # ... (build_dhcp_system!と同じロジック)

        # 7点ステンシル適用
        ss = ( a_w_val * p[i, j-1, k_idx] + a_e_val * p[i, j+1, k_idx]
             + a_s_val * p[i-1, j, k_idx] + a_n_val * p[i+1, j, k_idx]
             + a_b_val * p[i, j, k_idx-1] + a_t_val * p[i, j, k_idx+1] )

        # A*p の結果（RHSは含まない）
        ap[i, j, k_idx] = a_p_total * p[i, j, k_idx] - ss
    end
end
```

**重要な違い（build_dhcp_system!との比較）**:
- **係数配列を保存しない**: a_w, a_e, ... をその場で計算
- **RHSベクトルbを含まない**: CalcAX!はA*pのみを計算
- **@floopマクロで並列化**: スレッドレベル並列性

#### PBiCGSTAB_DHCP!関数の設計

**入力**:
- `X::Array{Float64,3}`: 解ベクトル（初期推定値、出力）
- `B::Array{Float64,3}`: RHSベクトル
- `work_vectors...`: ワークベクトル（q, r, r0, p, p_, s, s_, t_）
- その他のパラメータ（T_prev, q_surface, rho, ...）

**処理**:
```julia
function PBiCGSTAB_DHCP!(X::Array{Float64,3},
                          B::Array{Float64,3},
                          pcg_q, pcg_r, pcg_r0, pcg_p, pcg_p_,
                          pcg_s, pcg_s_, pcg_t_,
                          T_prev, q_surface,
                          rho, cp_coeffs, k_coeffs,
                          dx, dy, dz, dz_b, dz_t, dt;
                          rtol=1e-6, maxiter=20000, verbose=false)
    # 初期残差計算: r = B - A*X
    CalcAX_DHCP!(pcg_q, X, T_prev, q_surface, rho, cp_coeffs, k_coeffs,
                  dx, dy, dz, dz_b, dz_t, dt)
    pcg_r .= B .- pcg_q
    res0 = norm(pcg_r)

    pcg_r0 .= pcg_r
    rho_old = 1.0
    alpha = 0.0
    omega = 1.0

    for itr in 1:maxiter
        rho = dot(pcg_r, pcg_r0)

        if abs(rho) < eps(Float64)
            break
        end

        if itr == 1
            pcg_p .= pcg_r
        else
            beta = rho / rho_old * alpha / omega
            # BiCG1!: p = r + beta * (p - omega * q)
            pcg_p .= pcg_r .+ beta .* (pcg_p .- omega .* pcg_q)
        end

        # 前処理: p_ = M^{-1} * p
        Preconditioner_DHCP!(pcg_p_, pcg_p, ...)

        # q = A * p_
        CalcAX_DHCP!(pcg_q, pcg_p_, T_prev, q_surface, rho, cp_coeffs, k_coeffs,
                      dx, dy, dz, dz_b, dz_t, dt)

        alpha = rho / dot(pcg_q, pcg_r0)

        # s = r - alpha * q
        pcg_s .= pcg_r .- alpha .* pcg_q

        # 前処理: s_ = M^{-1} * s
        Preconditioner_DHCP!(pcg_s_, pcg_s, ...)

        # t_ = A * s_
        CalcAX_DHCP!(pcg_t_, pcg_s_, T_prev, q_surface, rho, cp_coeffs, k_coeffs,
                      dx, dy, dz, dz_b, dz_t, dt)

        omega = dot(pcg_t_, pcg_s) / dot(pcg_t_, pcg_t_)

        # X = X + alpha * p_ + omega * s_
        X .= X .+ alpha .* pcg_p_ .+ omega .* pcg_s_

        # r = s - omega * t_
        pcg_r .= pcg_s .- omega .* pcg_t_

        res = norm(pcg_r) / res0

        if verbose
            println("[BiCGstab] itr=$itr, res=$res")
        end

        if res < rtol
            if verbose
                println("[BiCGstab] 収束: $itr 回")
            end
            break
        end

        rho_old = rho
    end
end
```

### 5.3 前処理器の選択

**オプション1: Jacobi前処理（恒等写像）**
- 実装: `pcg_p_ .= pcg_p`
- 利点: コストゼロ
- 欠点: 収束性が悪い

**オプション2: 対角Jacobi前処理**
- 実装: `pcg_p_ .= pcg_p ./ diag_A`
- 利点: 簡単、並列化容易
- 欠点: 対角要素の事前計算が必要

**オプション3: Gauss-Seidel前処理（Heat3ds方式）**
- 実装: RB-SORを5反復
- 利点: 収束性が良い
- 欠点: 計算コストが高い

**推奨**: まずはオプション2（対角Jacobi前処理）で実装し、性能を測定。必要に応じてオプション3を検討。

## 6. 期待される性能改善

### 6.1 理論的な見積もり

**現在の性能**（中規模問題: 40×50×10×50）:
- 総実行時間: 66.37秒
- DHCP solver: 約30秒（推定）
- Adjoint solver: 約15秒（推定）
- CGM最適化: 約21秒（推定）

**マトリックスフリー版の期待効果**:

#### 疎行列組み立てコストの削減
- 現在: 11.6ms/ステップ × 357ステップ × 3回（DHCP, Adjoint×2） = **12.4秒**
- マトリックスフリー: **0秒**
- **削減**: 12.4秒

#### キャッシュ効率向上によるCG求解時間削減
- 現在のCG求解時間: 309.8ms/ステップ（対角前処理版）
- 期待削減率: 30-50%（キャッシュミス率改善）
- 削減時間: 92.9-154.9ms/ステップ
- 総削減: 92.9-154.9ms × 357ステップ × 3回 = **99-166秒**

#### 総削減
- 疎行列組み立て: 12.4秒
- CG求解時間: 99-166秒
- **合計**: 111-178秒

#### 予測実行時間
- 現在: 66.37秒
- マトリックスフリー版: **46-55秒**（30-40%削減）

### 6.2 保守的な見積もり

**キャッシュ効率向上が期待値の50%だった場合**:
- CG求解時間削減: 50-83秒
- 総削減: 62-95秒
- 予測実行時間: **52-58秒**（15-25%削減）

### 6.3 楽観的な見積もり

**キャッシュ効率向上 + 並列化（4スレッド）の場合**:
- CG求解時間削減: 60% + 並列化効率70% = 88%削減
- 総削減: 200秒以上
- 予測実行時間: **30-40秒**（40-50%削減）

## 7. リスクと対策

### 7.1 リスク

**リスク1: BiCGstabの収束性**
- CG法（対称正定値行列専用）よりも収束が遅い可能性
- 対策: 前処理器の調整、BiCGstab(l)への拡張

**リスク2: 実装の複雑性**
- Heat3dsとIHCPの境界条件が異なる
- 対策: 段階的な移植、徹底的なテスト

**リスク3: 数値精度の劣化**
- 浮動小数点演算の順序変更による誤差
- 対策: Python参照データとの厳密な比較テスト

### 7.2 対策

**対策1: TDD（Test-Driven Development）方針を継続**
- Phase 1-6で確立したテストスイート（505項目）を活用
- 新実装も同じテストで検証

**対策2: 段階的な移植**
- まずDHCP solverのみをマトリックスフリー化
- テスト合格後、Adjoint solverに移植
- 最後にCGMと統合

**対策3: 性能計測の自動化**
- ベンチマークスクリプトを整備
- 各段階でベースラインとの比較

## 8. 次のステップ

### 8.1 実装開始前の準備

1. **Heat3dsコードの詳細分析**（1日）
   - CalcAX!の完全理解
   - 境界条件処理の確認
   - マスク配列の扱い

2. **設計書作成**（1日）
   - CalcAX_DHCP!の詳細設計
   - PBiCGSTAB_DHCP!の詳細設計
   - インターフェース仕様

3. **開発環境整備**（半日）
   - 新ブランチ作成（`feature/matrix-free-bicgstab`）
   - ベンチマークスクリプト準備

### 8.2 実装スケジュール（4-5週間）

**Week 1**: CalcAX_DHCP!実装とテスト
**Week 2**: PBiCGSTAB_DHCP!実装とテスト
**Week 3**: CalcAX_Adjoint!とPBiCGSTAB_Adjoint!実装
**Week 4**: CGM統合とテスト
**Week 5**: 性能測定と最適化

### 8.3 成功基準

**必須**:
- 全505テストが合格
- Python参照データと相対誤差1e-7以内

**目標**:
- 中規模問題で15-25%高速化
- 大規模問題で30-40%高速化

**ストレッチゴール**:
- 並列化で50%以上高速化
- GPU対応の基盤確立

---

**作成日**: 2025年10月3日
**Phase**: 2.3準備 - マトリックスフリーBiCGstab法移植計画
**参考実装**: /Users/Daily/Development/Heat3ds/heat3d_NonUniform.jl
**ブランチ**: tuning2（分析作業中）
