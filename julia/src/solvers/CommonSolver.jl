"""
CommonSolver.jl


主要関数:
- solve_dhcp!: 複数時間ステップCG求解（ホットスタート対応）
"""

module CommonSolver

using Printf
using FLoops
using ThreadsX

import ..Commons
using ..Commons: WorkBuffers, FloatMin, λf, get_backend, myfill!, mycopy!


export PBiCGSTAB!, PCG!, solve_linear_system!


"""
@brief PCG反復（前処理付き共役勾配法）

@param [in,out] wk       ワークベクトル
@param [in]     Δh       セル幅
@param [in]     Δt       時間積分幅
@param [in]     ZC       CVセンター座標
@param [in]     ΔZ       CV幅
@param [in]     ρ        SUS密度

# キーワード引数
@param [in]     tol      反復閾値
@param [in]     maxItr   最大反復数
@param [in]     smoother [:none, :gs, :jacobi]
@param [in]     par      バックエンド（"sequential", "thread"）
@param [in]     verbose  詳細出力フラグ

収束判定： 相対残差ノルム (||r_k|| / ||r_0||) < tol
        IterativeSolvers.jlのcg!関数 : 相対残差ノルム
@ret            収束/未収束、反復回数、初期残差
"""
function PCG!(wk::WorkBuffers,
              Δh::NTuple{3,T},
              Δt::T,
              ZC::AbstractVector{T},
              ΔZ::AbstractVector{T},
              ρ::T,
              HC::AbstractVector{T};
              tol::T = T(1e-6),
              maxItr::Int = 20_000,
              smoother::Symbol = :none,
              par::String = "sequential",
              verbose::Bool=false) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(wk.θ)

    # 初期残差を計算: r = b - Ax
    myfill!(wk.pcg_q, zero(T), par)
    res0 = CalcRK!(wk.pcg_r, wk.θ, wk.b, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, ZC, ΔZ, HC, par)
    if verbose
        println("Inital residual = ", res0)
    end

    # 初期残差がゼロの場合は収束済み（数値安定性対策）
    if res0 ≈ zero(T)
        return true, 0, res0
    end

    # Smoother選択: Symbol → Val型変換（コンパイル時分岐解決）
    smoother_val = smoother_selector(smoother)

    # 前処理: z = M^-1 * r (pcg_sをzとして使用)
    myfill!(wk.pcg_s, zero(T), par)
    Preconditioner!(wk.pcg_s, wk.pcg_r, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother_val, ZC, ΔZ, HC, par, wk.tmp)

    # p = z
    mycopy!(wk.pcg_p, wk.pcg_s, par)

    # rho_old = (r, z)
    rho_old::T = Fdot2(wk.pcg_r, wk.pcg_s, par)
    isconverged::Bool = false
    itr::Int = 0
    float_min_T = T(FloatMin)

    for k in 1:maxItr
        itr = k

        # q = A * p
        CalcAX!(wk.pcg_q, wk.pcg_p, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, ZC, ΔZ, HC, par)

        # 分母ゼロ対策（数値安定性）
        denom = Fdot2(wk.pcg_p, wk.pcg_q, par)
        if abs(denom) < float_min_T
            # 分母がゼロに近い場合は数値的に不安定（未収束として扱う）
            isconverged = false
            break
        end

        # alpha = rho_old / (p, q)
        alpha::T = rho_old / denom

        # x = x + alpha * p, r = r - alpha * q
        @floop backend for kk in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
            wk.θ[i,j,kk] += alpha * wk.pcg_p[i,j,kk]
            wk.pcg_r[i,j,kk] -= alpha * wk.pcg_q[i,j,kk]
        end

        # 残差チェック
        res = sqrt(Fdot1(wk.pcg_r, par))
        res /= res0

        if res < tol
            isconverged = true
            break
        end

        # 前処理: z = M^-1 * r
        myfill!(wk.pcg_s, zero(T), par)
        Preconditioner!(wk.pcg_s, wk.pcg_r, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother_val, ZC, ΔZ, HC, par, wk.tmp)

        # rho_new = (r, z)
        rho_new = Fdot2(wk.pcg_r, wk.pcg_s, par)

        # rhoがゼロに近い場合は数値的に不安定（未収束として扱う）
        if abs(rho_new) < float_min_T
            isconverged = false
            break
        end

        # beta = rho_new / rho_old
        beta::T = rho_new / rho_old

        # p = z + beta * p
        @floop backend for kk in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
            wk.pcg_p[i,j,kk] = wk.pcg_s[i,j,kk] + beta * wk.pcg_p[i,j,kk]
        end

        rho_old = rho_new
    end # itr

    return isconverged, itr, res0
end

"""
@brief Smoother選択器（Symbol → Val型変換）

Symbol型から型パラメータVal{S}に変換することで、コンパイル時に分岐を解決。
これによりループ内でも高速な前処理ディスパッチが可能になる。

処理フロー:
  PBiCGSTAB! → smoother_selector() → Preconditioner! → _Preconditioner!(Val{S})

対応smoother:
  - :none   → 前処理なし（恒等変換）
  - :gs     → Gauss-Seidel法（RB-SOR、5回反復）
  - :jacobi → Jacobi法（加重Jacobi、ω=0.8、5回反復）

@param [in] s  Smoother種別（:none, :gs, :jacobi）
@ret           Val型（コンパイル時ディスパッチ用）
"""
@inline function smoother_selector(s::Symbol)
  if s === :none || s === :gs || s === :jacobi
    return Val(s)
  end
  throw(ArgumentError("Unsupported smoother: $s"))
end

# 前処理反復回数（Gauss-Seidel、Jacobiで使用）
const PRECONDITIONER_SWEEPS = 5

# Jacobi法の緩和係数（加重Jacobi法: x_new = (1-ω)x_old + ω*D^{-1}(b-Rx)）
const JACOBI_RELAXATION = 0.8


"""
@brief PBiCGSTAB反復
@param [in]     wk   ワークベクトル
@param [in]     Δh   セル幅
@param [in]     Δt   時間積分幅
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ   CV幅
@param [in]     ρ    SUS密度

# キーワード引数
@param [in]     tol    反復閾値
@param [in]     maxItr 最大反復数
@param [in]     smoother [:none, :gs, :jacobi]
@param [in]     par    バックエンド（"sequential", "thread"）

収束判定： 相対残差ノルム (||r_k|| / ||r_0||) < tol
        IterativeSolvers.jlのcg!関数 : 相対残差ノルム
@ret            収束/未収束、反復回数、初期残差
"""
function PBiCGSTAB!(wk::WorkBuffers,
                    Δh::NTuple{3,T},
                    Δt::T,
                    ZC::AbstractVector{T},
                    ΔZ::AbstractVector{T},
                    ρ::T,
                    HC::AbstractVector{T};
                    tol::T = T(1e-6),
                    maxItr::Int = 20_000,
                    smoother::Symbol = :none,
                    par::String = "sequential",
                    verbose::Bool=false) where {T <: AbstractFloat}
    SZ = size(wk.θ)
    myfill!(wk.pcg_q, zero(T), par)
    res0 = CalcRK!(wk.pcg_r, wk.θ, wk.b, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, ZC, ΔZ, HC, par)
    if verbose
        println("Inital residual = ", res0)
    end

    # 初期残差がゼロの場合は収束済み（数値安定性対策）
    if res0 ≈ zero(T)
        return true, 0, res0
    end

    mycopy!(wk.pcg_r0, wk.pcg_r, par) # wk.pcg_r0 .= wk.pcg_r  #copy!(pcg_r0, pcg_r)

    rho_old::T = one(T)
    alpha::T = zero(T)
    omega::T  = one(T)
    r_omega::T = -omega
    beta::T = zero(T)
    isconverged::Bool = false
    itr::Int = 0

    # Smoother選択: Symbol → Val型変換（コンパイル時分岐解決）
    smoother_val = smoother_selector(smoother)

    #inv_cell_count = inv(T((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2)))
    float_min_T = T(FloatMin)

    for k in 1:maxItr
        itr = k
        rho = Fdot2(wk.pcg_r, wk.pcg_r0, par) # 非計算部分はゼロのこと

        if abs(rho) < float_min_T
            # rhoがゼロに近い場合は数値的に不安定（未収束として扱う）
            isconverged = false
            break
        end

        if k == 1
            mycopy!(wk.pcg_p, wk.pcg_r, par)  #copy!(pcg_p, pcg_r)
        else
            beta = rho / rho_old * alpha / omega
            BiCG1!(wk.pcg_p, wk.pcg_r, wk.pcg_q, beta, omega, par)
        end

        myfill!(wk.pcg_p_, zero(T), par)  #fill!(pcg_p_, 0.0)
        # 前処理: pcg_p_ ≈ M^{-1} pcg_p （smoother_valで分岐、wk.tmpはJacobi用）
        Preconditioner!(wk.pcg_p_, wk.pcg_p, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother_val, ZC, ΔZ, HC, par, wk.tmp)

        CalcAX!(wk.pcg_q, wk.pcg_p_, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, ZC, ΔZ, HC, par)
        alpha = rho / Fdot2(wk.pcg_q, wk.pcg_r0, par)
        r_alpha = -alpha
        Triad!(wk.pcg_s, wk.pcg_q, wk.pcg_r, r_alpha, par)

        myfill!(wk.pcg_s_, zero(T), par)  #fill!(pcg_s_, 0.0)
        # 前処理: pcg_s_ ≈ M^{-1} pcg_s （2回目の前処理適用）
        Preconditioner!(wk.pcg_s_, wk.pcg_s, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother_val, ZC, ΔZ, HC, par, wk.tmp);

        CalcAX!(wk.pcg_t_, wk.pcg_s_, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, ZC, ΔZ, HC, par)

        # 分母ゼロ対策（数値安定性）
        denom = Fdot1(wk.pcg_t_, par)
        if abs(denom) < float_min_T
            # 分母がゼロに近い場合は数値的に不安定（未収束として扱う）
            isconverged = false
            break
        end
        omega = Fdot2(wk.pcg_t_, wk.pcg_s, par) / denom
        r_omega = -omega

        BICG2!(wk.θ, wk.pcg_p_, wk.pcg_s_, alpha , omega, par)

        Triad!(wk.pcg_r, wk.pcg_t_, wk.pcg_s, r_omega, par)
        res = sqrt(Fdot1(wk.pcg_r, par)) #* inv_cell_count
        res /= res0

        if res < tol
            isconverged = true
            break
        end

        rho_old = rho
    end # itr
    
    return  isconverged, itr, res0
end


"""
@brief 残差ベクトルの計算
@param [out]    r    残差ベクトル
@param [in]     θ    解ベクトル
@param [in]     b    右辺ベクトル
@param [in]     λ    熱伝導率
@param [in]     cp   比熱
@param [in]     m    マスク配列
@param [in]     ρ    密度
@param [in]     Δh   セル幅
@param [in]     Δt   時間積分幅
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ   CV幅
@param [in]     HC   熱伝達係数 [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]

@ret                 セルあたりの残差RMS
"""
function CalcRK!(
                r::AbstractArray{T,3},
                θ::AbstractArray{T,3},
                b::AbstractArray{T,3},
                λ::AbstractArray{T,3},
                cp::AbstractArray{T,3},
                m::AbstractArray{T,3},
                ρ::T,
                Δh::NTuple{3,T},
                Δt::T,
                ZC::AbstractVector{T},
                ΔZ::AbstractVector{T},
                HC::AbstractVector{T},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    ddt = inv(Δt)
    ddx = inv(dx0)
    ddy = inv(dy0)
    oneT = one(T)  # ループ外で事前計算

    @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        dz_k = ΔZ[k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]

        # 熱伝導項（面積×コンダクタンス）
        cond_xm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dy0 * dz_k * ddx
        cond_xp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dy0 * dz_k * ddx
        cond_ym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dx0 * dz_k * ddy
        cond_yp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dx0 * dz_k * ddy
        cond_zm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dx0 * dy0 / (ZC[k]-ZC[k-1])
        cond_zp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dx0 * dy0 / (ZC[k+1]-ZC[k])

        # 対流項（対角項のみに追加）
        conv_xm = HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
        conv_xp = HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
        conv_ym = HC[3] * dx0 * dz_k * (oneT - m[i,j-1,k])
        conv_yp = HC[4] * dx0 * dz_k * (oneT - m[i,j+1,k])
        conv_zm = HC[5] * dx0 * dy0 * (oneT - m[i,j,k-1])
        conv_zp = HC[6] * dx0 * dy0 * (oneT - m[i,j,k+1])

        # 時間項
        a_p_0 = ρ * cp[i,j,k] * dx0 * dy0 * dz_k * ddt

        # 対角項: 熱伝導項 + 対流項 + 時間項
        dd = (oneT-m0) + (cond_xp + cond_xm + cond_yp + cond_ym + cond_zp + cond_zm +
                          conv_xp + conv_xm + conv_yp + conv_ym + conv_zp + conv_zm + a_p_0)*m0

        # 隣接項: 熱伝導項のみ（対流項は除外）
        ss = ( cond_xp * θ[i+1,j  ,k  ] + cond_xm * θ[i-1,j  ,k  ]
             + cond_yp * θ[i  ,j+1,k  ] + cond_ym * θ[i  ,j-1,k  ]
             + cond_zp * θ[i  ,j  ,k+1] + cond_zm * θ[i  ,j  ,k-1] )

        rs = (b[i,j,k] - (ss - dd * θ[i,j,k]))* m0
        r[i,j,k] = rs
        @reduce(res = zero(T) + rs*rs)
    end
    return sqrt(res)
end


"""
@brief AXの計算
@param [out] ax   AX
@param [in]  θ    解ベクトル
@param [in]  Δh   セル幅
@param [in]  Δt   時間積分幅
@param [in]  λ    熱伝導率
@param [in]  cp   比熱
@param [in]  m    マスク配列
@param [in]  ρ    密度
@param [in]  ZC   CVセンター座標
@param [in]  ΔZ   CV幅
@param [in]  HC   熱伝達係数 [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]

"""
function CalcAX!(ax::AbstractArray{T,3},
                  θ::AbstractArray{T,3},
                  Δh::NTuple{3,T},
                  Δt::T,
                  λ::AbstractArray{T,3},
                  cp::AbstractArray{T,3},
                  m::AbstractArray{T,3},
                  ρ::T,
                  ZC::AbstractVector{T},
                  ΔZ::AbstractVector{T},
                  HC::AbstractVector{T},
                  par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    ddt = inv(Δt)
    ddx = inv(dx0)
    ddy = inv(dy0)
    oneT = one(T)  # ループ外で事前計算

    @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        dz_k = ΔZ[k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]

        # 熱伝導項（面積×コンダクタンス）
        cond_xm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dy0 * dz_k * ddx
        cond_xp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dy0 * dz_k * ddx
        cond_ym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dx0 * dz_k * ddy
        cond_yp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dx0 * dz_k * ddy
        cond_zm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dx0 * dy0 / (ZC[k]-ZC[k-1])
        cond_zp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dx0 * dy0 / (ZC[k+1]-ZC[k])

        # 対流項（対角項のみに追加）
        conv_xm = HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
        conv_xp = HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
        conv_ym = HC[3] * dx0 * dz_k * (oneT - m[i,j-1,k])
        conv_yp = HC[4] * dx0 * dz_k * (oneT - m[i,j+1,k])
        conv_zm = HC[5] * dx0 * dy0 * (oneT - m[i,j,k-1])
        conv_zp = HC[6] * dx0 * dy0 * (oneT - m[i,j,k+1])

        # 時間項
        a_p_0 = ρ * cp[i,j,k] * dx0 * dy0 * dz_k * ddt

        # 対角項: 熱伝導項 + 対流項 + 時間項
        dd = (oneT-m0) + (cond_xp + cond_xm + cond_yp + cond_ym + cond_zp + cond_zm +
                          conv_xp + conv_xm + conv_yp + conv_ym + conv_zp + conv_zm + a_p_0)*m0

        # 隣接項: 熱伝導項のみ（対流項は除外）
        ss = ( cond_xp * θ[i+1,j  ,k  ] + cond_xm * θ[i-1,j  ,k  ]
             + cond_yp * θ[i  ,j+1,k  ] + cond_ym * θ[i  ,j-1,k  ]
             + cond_zp * θ[i  ,j  ,k+1] + cond_zm * θ[i  ,j  ,k-1] )

        ax[i,j,k] = (ss - dd*θ[i,j,k]) * m0
    end
end


"""
@brief 前処理（公開API）

BiCGSTAB法の収束を加速するための前処理を適用する。

【設計パターン】
この関数は公開APIとして以下の役割を持つ：
  1. 型チェックとドキュメント提供
  2. 内部実装 _Preconditioner! へのディスパッチ

【処理フロー】
  Preconditioner!(smoother::Val{S})
    → _Preconditioner!(smoother::Val{S})  ← Juliaの多重ディスパッチ
       ├─ Val{:none}   → mycopy!（恒等変換）
       ├─ Val{:gs}     → rbsor!（RB-SOR法、5回反復）
       └─ Val{:jacobi} → jacobi_preconditioner!（加重Jacobi法、5回反復）

【Val型の利点】
  - コンパイル時に分岐解決（if文なし）
  - ループ内でもゼロオーバーヘッド
  - 型安定性による高速化

@param [in,out] xx       解ベクトル（前処理後の値で上書き）
@param [in]     bb       RHSベクトル
@param [in]     λ        熱伝導率
@param [in]     cp       比熱
@param [in]     mask     マスク配列
@param [in]     ρ        密度
@param [in]     Δh       セル幅
@param [in]     Δt       時間積分幅
@param [in]     smoother Val{:none}, Val{:gs}, Val{:jacobi}
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ       CV幅
@param [in]     par      バックエンド（"sequential", "thread"）
@param [in]     scratch  ワーク配列（Jacobi前処理で使用、WorkBuffers.tmp）
"""
function Preconditioner!(xx::AbstractArray{T,3},
                         bb::AbstractArray{T,3},
                         λ::AbstractArray{T,3},
                         cp::AbstractArray{T,3},
                         mask::AbstractArray{T,3},
                         ρ::T,
                         Δh::NTuple{3,T},
                         Δt::T,
                         smoother::Val,
                         ZC::AbstractVector{T},
                         ΔZ::AbstractVector{T},
                         HC::AbstractVector{T},
                         par::String,
                         scratch::AbstractArray{T,3}) where {T <: AbstractFloat}
    # 内部実装にディスパッチ（Val型による分岐はコンパイル時に解決）
    _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, smoother, ZC, ΔZ, HC, par, scratch)
end

"""
@brief 前処理実装: :none（恒等変換）

前処理なし。右辺ベクトルをそのまま解ベクトルにコピーする。
xx = bb

【用途】
  - 前処理のオーバーヘッドを避けたい場合
  - 行列の条件数が良好な場合（収束が速い）
  - ベースライン性能測定

【計算量】 O(N)（単純なコピー）
"""
@inline function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{:none}, ZC, ΔZ, HC, par, _)
    mycopy!(xx, bb, par)
    return nothing
end

"""
@brief 前処理実装: :gs（Gauss-Seidel法）

Red-Black SOR法を5回反復して前処理を行う。
xx ≈ M^{-1} bb （M: 前処理行列、ここではGS反復5回による近似逆行列）

【アルゴリズム】
  - Red-Black順序付きSOR法（並列化対応）
  - 緩和係数 ω = 1.0（純粋なGauss-Seidel）
  - 反復回数: PRECONDITIONER_SWEEPS = 5

【用途】
  - 従来から使用されている標準的な前処理
  - 滑らかな誤差成分を効果的に除去

【計算量】 O(5N)（5回反復 × N点更新）
"""
function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{:gs}, ZC, ΔZ, HC, par, _)
    for _ in 1:PRECONDITIONER_SWEEPS
        rbsor!(xx, λ, cp, bb, mask, ρ, Δh, Δt, one(ρ), ZC, ΔZ, HC, par)
    end
    return nothing
end

"""
@brief 前処理実装: :jacobi（Jacobi法）

加重Jacobi法を5回反復して前処理を行う。
xx ≈ M^{-1} bb （M: 対角行列Dによる前処理）

【アルゴリズム】
  - 加重Jacobi法: x_{k+1} = (1-ω)x_k + ω D^{-1}(bb - R x_k)
  - 緩和係数: ω = JACOBI_RELAXATION = 0.8
  - 反復回数: PRECONDITIONER_SWEEPS = 5
  - ワーク配列: scratch（WorkBuffers.tmp）を使用

【特徴】
  - Gauss-Seidelより並列性が高い（理論上）
  - 対角要素のみ計算（マトリクスフリー）
  - メモリアクセスパターンが単純

【用途】
  - Phase 1-Cで導入（GS法との性能比較）
  - 並列計算での性能改善を期待

【計算量】 O(5N)（5回反復 × N点更新）
"""
function _Preconditioner!(xx::AbstractArray{T,3},
                          bb::AbstractArray{T,3},
                          λ::AbstractArray{T,3},
                          cp::AbstractArray{T,3},
                          mask::AbstractArray{T,3},
                          ρ::T,
                          Δh::NTuple{3,T},
                          Δt::T,
                          ::Val{:jacobi},
                          ZC::AbstractVector{T},
                          ΔZ::AbstractVector{T},
                          HC::AbstractVector{T},
                          par::String,
                          scratch::AbstractArray{T,3}) where {T <: AbstractFloat}
    jacobi_preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ZC, ΔZ, HC, par, scratch)
    return nothing
end

"""
@brief 前処理実装: エラーハンドラ

サポートされていないsmoother種別に対するエラー処理。
この関数は通常実行されない（smoother_selector()で事前チェック済み）。

型システムの完全性のために定義。
"""
@inline function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{s}, ZC, ΔZ, HC, par, scratch) where {s}
    throw(ArgumentError("Unsupported smoother: $(s)"))
end

"""
@brief Jacobi前処理の実装（マトリクスフリー）

加重Jacobi法による前処理を5回反復で実行する。
xx ≈ (I - ωD^{-1}(A-D))^5 bb

【アルゴリズム詳細】
1. 反復: PRECONDITIONER_SWEEPS = 5回
2. 各反復で加重Jacobi更新:
   - 現在の解xxをscratchにコピー
   - 各格子点で対角要素ddを計算（マトリクスフリー）
   - 残差: rhs - Ax を計算
   - 更新: xx[i,j,k] = xx[i,j,k] + ω * (rhs - Ax) / dd
3. 境界セルはbbの値をそのままコピー（mask==0の場合）

【パラメータ】
- 緩和係数: ω = JACOBI_RELAXATION = 0.8
- ワーク配列: scratch（WorkBuffers.tmp、サイズはxxと同じ）

【数値安定性】
- 対角要素ゼロ回避: dd = max(dd, FloatMin)
- mask処理: 境界セル（m0==0）はbbの値を維持

【計算コスト】
- 5回反復 × N点 × 7点ステンシル
- GS法と同等のFLOPS、但し並列性が高い

@param [in,out] xx       解ベクトル（前処理後の値で上書き）
@param [in]     bb       RHSベクトル
@param [in]     λ        熱伝導率
@param [in]     cp       比熱
@param [in]     m     マスク配列（1.0=計算点、0.0=境界点）
@param [in]     ρ        密度
@param [in]     Δh       セル幅 (dx, dy)
@param [in]     Δt       時間積分幅
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ       CV幅
@param [in]     par      バックエンド（"sequential", "thread"）
@param [in]     scratch  ワーク配列（WorkBuffers.tmp）
"""
function jacobi_preconditioner!(xx::AbstractArray{T,3},
                                bb::AbstractArray{T,3},
                                λ::AbstractArray{T,3},
                                cp::AbstractArray{T,3},
                                m::AbstractArray{T,3},
                                ρ::T,
                                Δh::NTuple{3,T},
                                Δt::T,
                                ZC::AbstractVector{T},
                                ΔZ::AbstractVector{T},
                                HC::AbstractVector{T},
                                par::String,
                                scratch::AbstractArray{T,3}) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(xx)
    dx0 = Δh[1]
    dy0 = Δh[2]
    ddt = inv(Δt)
    ddx = inv(dx0)
    ddy = inv(dy0)
    float_min_T = T(FloatMin)
    ω = T(JACOBI_RELAXATION)
    zeroT = zero(T)  # ループ外で事前計算
    oneT = one(T)    # ループ外で事前計算

    # 加重Jacobi法を5回反復
    for _ in 1:PRECONDITIONER_SWEEPS
      mycopy!(scratch, xx, par)

      @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        dz_k = ΔZ[k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]

        # Dirichlet/ghost cells are copied through
        if m0 == zeroT
          scratch[i,j,k] = bb[i,j,k]
          continue
        end

        # 体積積分形式：面積×コンダクタンス + 対流項
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dy0 * dz_k * ddx + HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dy0 * dz_k * ddx + HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dx0 * dz_k * ddy + HC[3] * dx0 * dz_k * (oneT - m[i,j-1,k])
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dx0 * dz_k * ddy + HC[4] * dx0 * dz_k * (oneT - m[i,j+1,k])
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dx0 * dy0 / (ZC[k]-ZC[k-1]) + HC[5] * dx0 * dy0 * (oneT - m[i,j,k-1])
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dx0 * dy0 / (ZC[k+1]-ZC[k]) + HC[6] * dx0 * dy0 * (oneT - m[i,j,k+1])

        # 時間項（体積積分形式） - 定常解析の場合は0
        a_p_0 = ρ * cp[i,j,k] * dx0 * dy0 * dz_k * ddt

        dd = (oneT-m0) + (axp + axm + ayp + aym + azp + azm + a_p_0)*m0
        ss = ( axp * scratch[i+1,j  ,k  ] + axm * scratch[i-1,j  ,k  ]
             + ayp * scratch[i  ,j+1,k  ] + aym * scratch[i  ,j-1,k  ]
             + azp * scratch[i  ,j  ,k+1] + azm * scratch[i  ,j  ,k-1] )

        Ax = ss - dd * scratch[i,j,k]
        rhs = bb[i,j,k]
        diag = max(dd, float_min_T)
        scratch[i,j,k] = xx[i,j,k] + ω * (rhs - Ax) / diag
      end

      mycopy!(xx, scratch, par)
    end

  return nothing
end


"""
@brief ベクトルの内積
@param [in]     x    ベクトル
@ret            内積
"""
function Fdot1(x::AbstractArray{T,3}, par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(x)
    z_st = 2
    z_ed = SZ[3] - 1
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        @reduce(sum = zero(T) + x[i,j,k] * x[i,j,k])
    end
    return sum
end


"""
@brief 2ベクトルの内積
@param [in]     x    ベクトル
@param [in]     y    ベクトル
@ret            内積
"""
function Fdot2(x::AbstractArray{T,3}, y::AbstractArray{T,3}, par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(x)
    z_st = 2
    z_ed = SZ[3] - 1
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        @reduce(sum = zero(T) + x[i,j,k] * y[i,j,k])
    end
    return sum
end


"""
@brief BiCGstabの部分演算1
@param [in,out] p    ベクトル
@param [in]     r    ベクトル
@param [in]     q    ベクトル
@param [in]     beta 係数
@param [in]     omg  係数
"""
function BiCG1!(p::AbstractArray{T,3},
                r::AbstractArray{T,3},
                q::AbstractArray{T,3},
                beta::T,
                omg::T,
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(p)
    z_st = 2
    z_ed = SZ[3] - 1
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        p[i,j,k] = r[i,j,k] + beta * (p[i,j,k] - omg * q[i,j,k])
    end
end


"""
@brief AXPYZ
@param [out]    z    ベクトル
@param [in]     y    ベクトル
@param [in]     x    ベクトル
@param [in]     a    係数
"""
function Triad!(z::AbstractArray{T,3},
                x::AbstractArray{T,3},
                y::AbstractArray{T,3},
                a::T,
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(z)
    z_st = 2
    z_ed = SZ[3] - 1
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] = a * x[i,j,k] + y[i,j,k]
    end
end


"""
@brief BiCGstab 2
@param [in,out] z    ベクトル
@param [in]     y    ベクトル
@param [in]     x    ベクトル
@param [in]     a    係数
@param [in]     b    係数
"""
function BICG2!(z::AbstractArray{T,3},
                x::AbstractArray{T,3},
                y::AbstractArray{T,3},
                a::T,
                b::T,
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(z)
    z_st = 2
    z_ed = SZ[3] - 1
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        z[i,j,k] += a * x[i,j,k] + b * y[i,j,k]
    end
end


"""
@brief SOR法の残差
@param [in,out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     cp   比熱
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     ρ    密度
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ   格子幅
@param [in]     HC   熱伝達係数 [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]
@ret                 残差RMS
"""
function resSOR(θ::AbstractArray{T,3},
                λ::AbstractArray{T,3},
                cp::AbstractArray{T,3},
                b::AbstractArray{T,3},
                m::AbstractArray{T,3},
                ρ::T,
                Δh::NTuple{3,T},
                Δt::T,
                ω::T,
                ZC::AbstractVector{T},
                ΔZ::AbstractVector{T},
                HC::AbstractVector{T},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    ddt = inv(Δt)
    ddx = inv(dx0)
    ddy = inv(dy0)
    oneT = one(T)  # ループ外で事前計算

    @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        dz_k = ΔZ[k]
        pp = θ[i,j,k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]

        # 体積積分形式：面積×コンダクタンス + 対流項
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dy0 * dz_k * ddx + HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dy0 * dz_k * ddx + HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dx0 * dz_k * ddy + HC[3] * dx0 * dz_k * (oneT - m[i,j-1,k])
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dx0 * dz_k * ddy + HC[4] * dx0 * dz_k * (oneT - m[i,j+1,k])
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dx0 * dy0 / (ZC[k]-ZC[k-1]) + HC[5] * dx0 * dy0 * (oneT - m[i,j,k-1])
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dx0 * dy0 / (ZC[k+1]-ZC[k]) + HC[6] * dx0 * dy0 * (oneT - m[i,j,k+1])

        # 時間項（体積積分形式） - 定常解析の場合は0
        a_p_0 = ρ * cp[i,j,k] * dx0 * dy0 * dz_k * ddt

        dd = (oneT-m0) + (axp + axm + ayp + aym + azp + azm + a_p_0)*m0
        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )

        dp = (((ss-b[i,j,k])/dd - pp)) * m0
        θ[i,j,k] = pp + ω * dp
        r = (dd + ω*(axm+aym+azm))*dp / ω
        @reduce(res = zero(T) + r*r)
    end

    return res
end

"""
@brief RB-SOR法のカーネル
@param [in,out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     cp   比熱
@param [in]     b    右辺ベクトル
@param [in]     m    マスク配列
@param [in]     ρ    密度
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ   格子幅
@param [in]     HC   熱伝達係数 [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]
@param [in]     color R or B
@ret                 残差2乗和
"""
function rbsor_core!(θ::AbstractArray{T,3},
                     λ::AbstractArray{T,3},
                     cp::AbstractArray{T,3},
                     b::AbstractArray{T,3},
                     m::AbstractArray{T,3},
                     ρ::T,
                     Δh::NTuple{3,T},
                     Δt::T,
                     ω::T,
                     ZC::AbstractVector{T},
                     ΔZ::AbstractVector{T},
                     HC::AbstractVector{T},
                     color::Int,
                     par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    ddt = inv(Δt)
    ddx = inv(dx0)
    ddy = inv(dy0)
    oneT = one(T)  # ループ外で事前計算

    @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            dz_k = ΔZ[k]
            pp = θ[i,j,k]
            λ0 = λ[i,j,k]
            m0 = m[i,j,k]

            # 体積積分形式：面積×コンダクタンス + 対流項
            axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dy0 * dz_k * ddx + HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
            axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dy0 * dz_k * ddx + HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
            aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dx0 * dz_k * ddy + HC[3] * dx0 * dz_k * (oneT - m[i,j-1,k])
            ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dx0 * dz_k * ddy + HC[4] * dx0 * dz_k * (oneT - m[i,j+1,k])
            azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dx0 * dy0 / (ZC[k]-ZC[k-1]) + HC[5] * dx0 * dy0 * (oneT - m[i,j,k-1])
            azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dx0 * dy0 / (ZC[k+1]-ZC[k]) + HC[6] * dx0 * dy0 * (oneT - m[i,j,k+1])

            # 時間項（体積積分形式） - 定常解析の場合は0
            a_p_0 = ρ * cp[i,j,k] * dx0 * dy0 * dz_k * ddt

            dd = (oneT-m0) + (axp + axm + ayp + aym + azp + azm + a_p_0)*m0
            ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
                 + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
                 + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )

            dp = (((ss-b[i,j,k])/dd - pp)) * m0
            θ[i,j,k] = pp + ω * dp
            r = (dd + ω*(axm+aym+azm))*dp / ω
            @reduce(res = zero(T) + r*r)
        end
    end

    return res
end


"""
@brief RB-SOR法
@param [in,out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     cp   比熱
@param [in]     b    右辺ベクトル
@param [in]     mask マスク配列
@param [in]     ρ    密度
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     ZC   CVセンター座標
@param [in]     ΔZ   格子幅
@param [in]     HC   熱伝達係数 [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]
@ret                 残差RMS
"""
function rbsor!(θ::AbstractArray{T,3},
                λ::AbstractArray{T,3},
                cp::AbstractArray{T,3},
                b::AbstractArray{T,3},
                mask::AbstractArray{T,3},
                ρ::T,
                Δh::NTuple{3,T},
                Δt::T,
                ω::T,
                ZC::AbstractVector{T},
                ΔZ::AbstractVector{T},
                HC::AbstractVector{T},
                par::String) where {T <: AbstractFloat}
    SZ = size(b)
    res = zero(T)

    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        res += rbsor_core!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, ZC, ΔZ, HC, c, par)
    end
    #norm_factor = T((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
    return sqrt(res) #/ norm_factor
end



"""
@brief SOR法による求解
@param [in/out] θ    解ベクトル
@param [in]     λ    熱伝導率
@param [in]     cp   比熱
@param [in]     b    RHSベクトル
@param [in]     mask マスク配列
@param [in]     ρ    密度
@param [in]     Δh   セル幅
@param [in]     ω    加速係数
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     F    ファイルディスクリプタ
@param [in]     tol  収束判定基準
@param [in]     ItrMax 最大反復回数（デフォルト: 1000）

"""
function solveSOR!(θ::AbstractArray{T,3},
                    λ::AbstractArray{T,3},
                    cp::AbstractArray{T,3},
                    b::AbstractArray{T,3},
                    mask::AbstractArray{T,3},
                    ρ::T,
                    Δh::NTuple{3,T},
                    Δt::T,
                    ω::T,
                    ZC::AbstractVector{T},
                    ΔZ::AbstractVector{T},
                    F,
                    tol::T,
                    par::String;
                    ItrMax::Int=1000) where {T <: AbstractFloat}

    res0 = resSOR(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, ZC, ΔZ, par)
    if res0 == zero(T)
        res0 = one(T)
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        #res = sor!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, ZC, ΔZ) / res0
        res = rbsor!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, ZC, ΔZ, par) / res0
        @printf(F, "%10d %24.14E\n", n, Float64(res)) # 時間計測の場合にはコメントアウト
        @printf(stdout, "%10d %24.14E\n", n, Float64(res)) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end


"""
@brief ソルバー選択ヘルパー関数

線形方程式系 Ax = b を解くソルバーを選択して実行する。

【サポートされているソルバー】
  - :pbicgstab → PBiCGSTAB!（前処理付き双共役勾配安定化法）
  - :pcg → PCG!（前処理付き共役勾配法）

【使用例】
```julia
isconverged, itr, res0 = solve_linear_system!(
    wk, Δh, dt, ZC, dz, ρ,
    solver=:pcg,
    tol=1e-6, maxItr=20000, smoother=:gs, par="sequential"
)
```

@param [in,out] wk       ワークベクトル
@param [in]     Δh       セル幅
@param [in]     Δt       時間積分幅
@param [in]     ZC       CVセンター座標
@param [in]     ΔZ       CV幅
@param [in]     ρ        SUS密度

# キーワード引数
@param [in]     solver   ソルバー種別 [:pbicgstab, :pcg]
@param [in]     tol      反復閾値
@param [in]     maxItr   最大反復数
@param [in]     smoother [:none, :gs, :jacobi]
@param [in]     par      バックエンド（"sequential", "thread"）
@param [in]     verbose  詳細出力フラグ

@ret            収束/未収束、反復回数、初期残差
"""
function solve_linear_system!(
    wk::WorkBuffers,
    Δh::NTuple{3,T},
    Δt::T,
    ZC::AbstractVector{T},
    ΔZ::AbstractVector{T},
    ρ::T,
    HC::AbstractVector{T};
    solver::Symbol = :pbicgstab,
    tol::T = T(1e-6),
    maxItr::Int = 20_000,
    smoother::Symbol = :none,
    par::String = "sequential",
    verbose::Bool=false
) where {T <: AbstractFloat}
    if solver === :pbicgstab
        return PBiCGSTAB!(wk, Δh, Δt, ZC, ΔZ, ρ, HC,
            tol=tol, maxItr=maxItr, smoother=smoother, par=par, verbose=verbose)
    elseif solver === :pcg
        return PCG!(wk, Δh, Δt, ZC, ΔZ, ρ, HC,
            tol=tol, maxItr=maxItr, smoother=smoother, par=par, verbose=verbose)
    else
        throw(ArgumentError("Unsupported solver: $(solver). Use :pbicgstab or :pcg."))
    end
end


end  # module CommonSolver
