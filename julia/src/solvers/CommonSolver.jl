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
using ..Commons: WorkBuffers, FloatMin, λf, get_backend


export PBiCGSTAB!

"""
@brief 配列のコピー（並列対応版）
@param [out]    a    コピー先配列
@param [in]     b    コピー元配列
@param [in]     par  バックエンド

シングルスレッド時は組み込みのcopyto!を使用（SIMD最適化）
マルチスレッド時は並列ループで処理
"""
function mycopy!(a::AbstractArray{T,3}, b::AbstractArray{T,3}, par::String) where {T <: AbstractFloat}
  if par == "sequential" || Threads.nthreads() == 1
    # シングルスレッド: 高速な組み込み関数を使用
    copyto!(a, b)
  else
    # マルチスレッド: 並列ループ
    backend = get_backend(par)
    SZ = size(a)
    @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
      a[i,j,k] = b[i,j,k]
    end
  end
end

"""
@brief 並列用のfill（任意値の代入）
@param [out]    a    配列
@param [in]     val  代入値
@param [in]     par  バックエンド

シングルスレッド時は組み込みのfill!を使用（SIMD最適化）
マルチスレッド時は並列ループで処理
"""
function myfill!(a::AbstractArray{T,3}, val::T, par::String) where {T <: AbstractFloat}
  if par == "sequential" || Threads.nthreads() == 1
    # シングルスレッド: 高速な組み込み関数を使用
    fill!(a, val)
  else
    # マルチスレッド: 並列ループ
    backend = get_backend(par)
    SZ = size(a)
    @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
      a[i,j,k] = val
    end
  end
end

@inline function smoother_selector(s::Symbol)
  if s === :none || s === :gs || s === :jacobi
    return Val(s)
  end
  throw(ArgumentError("Unsupported smoother: $s"))
end

const PRECONDITIONER_SWEEPS = 5
const JACOBI_RELAXATION = 0.8


"""
@brief PBiCGSTAB反復
@param [in]     wk   ワークベクトル
@param [in]     Δh     セル幅
@param [in]     Δt     時間積分幅
@param [in]     Z      Z座標
@param [in]     ΔZ     格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値
@param [in]     ρ    SUS密度

# キーワード引数
@param [in]     tol    反復閾値
@param [in]     maxItr 最大反復数
@param [in]     smoother [:none, :gs, :jacobi]
@param [in]     par    バックエンド（"sequential", "thread"）

収束判定： 相対残差ノルム (||r_k|| / ||r_0||) < tol
一方、IterativeSolvers.jlのcg!関数は相対残差ノルム
@ret            収束/未収束、反復回数、初期残差
"""
function PBiCGSTAB!(wk::WorkBuffers,
                    Δh::NTuple{3,T},
                    Δt::T,
                    Z::AbstractVector{T},
                    ΔZ::AbstractVector{T},
                    z_range::AbstractVector{<:Integer},
                    HT::AbstractVector{T},
                    ρ::T;
                    tol::T = T(1e-6),
                    maxItr::Int = 20_000,
                    smoother::Symbol = :none,
                    par::String = "sequential",
                    verbose::Bool=false) where {T <: AbstractFloat}
    SZ = size(wk.θ)
    myfill!(wk.pcg_q, zero(T), par)
    res0 = CalcRK!(wk.pcg_r, wk.θ, wk.b, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, Z, ΔZ, z_range, HT, par)
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
    smoother_val = smoother_selector(smoother)
    inv_cell_count = inv(T((SZ[1]-2)*(SZ[2]-2)*(Int(last(z_range)) - Int(first(z_range)) + 1)))
    float_min_T = T(FloatMin)

    for k in 1:maxItr
        itr = k
        rho = Fdot2(wk.pcg_r, wk.pcg_r0, z_range, par) # 非計算部分はゼロのこと

        if abs(rho) < float_min_T
            # rhoがゼロに近い場合は数値的に不安定（未収束として扱う）
            isconverged = false
            break
        end

        if k == 1
            mycopy!(wk.pcg_p, wk.pcg_r, par)  #copy!(pcg_p, pcg_r)
        else
            beta = rho / rho_old * alpha / omega
            BiCG1!(wk.pcg_p, wk.pcg_r, wk.pcg_q, beta, omega, z_range, par)
        end

        myfill!(wk.pcg_p_, zero(T), par)  #fill!(pcg_p_, 0.0)
        Preconditioner!(wk.pcg_p_, wk.pcg_p, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother_val, Z, ΔZ, z_range, HT, par; scratch=wk.pcg_q)

        CalcAX!(wk.pcg_q, wk.pcg_p_, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, Z, ΔZ, z_range, HT, par)
        alpha = rho / Fdot2(wk.pcg_q, wk.pcg_r0, z_range, par)
        r_alpha = -alpha
        Triad!(wk.pcg_s, wk.pcg_q, wk.pcg_r, r_alpha, z_range, par)

        myfill!(wk.pcg_s_, zero(T), par)  #fill!(pcg_s_, 0.0)
        Preconditioner!(wk.pcg_s_, wk.pcg_s, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother_val, Z, ΔZ, z_range, HT, par; scratch=wk.pcg_q);

        CalcAX!(wk.pcg_t_, wk.pcg_s_, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, Z, ΔZ, z_range, HT, par)

        # 分母ゼロ対策（数値安定性）
        denom = Fdot1(wk.pcg_t_, z_range, par)
        if abs(denom) < float_min_T
            # 分母がゼロに近い場合は数値的に不安定（未収束として扱う）
            isconverged = false
            break
        end
        omega = Fdot2(wk.pcg_t_, wk.pcg_s, z_range, par) / denom
        r_omega = -omega

        BICG2!(wk.θ, wk.pcg_p_, wk.pcg_s_, alpha , omega, z_range, par)

        Triad!(wk.pcg_r, wk.pcg_t_, wk.pcg_s, r_omega, z_range, par)
        res = sqrt(Fdot1(wk.pcg_r, z_range, par)) * inv_cell_count
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
@param [in]     Z    CV境界座標
@param [in]     ΔZ   CV幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値

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
                Z::AbstractVector{T},
                ΔZ::AbstractVector{T},
                z_range::AbstractVector{<:Integer},
                HT::AbstractVector{T},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = inv(dx0 * dx0)
    dy2 = inv(dy0 * dy0)
    dx1 = inv(dx0)
    dy1 = inv(dy0)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
    ddt = inv(Δt)
    cell_count_inv = inv(T((SZ[1]-2)*(SZ[2]-2)*(z_ed - z_st + 1)))

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        r_ρc = inv(ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = one(T) - m[i-1,j  ,k  ]
        me = one(T) - m[i+1,j  ,k  ]
        ms = one(T) - m[i  ,j-1,k  ]
        mn = one(T) - m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc 
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc 
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc 
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc 
        zb = (Z[k]-Z[k-1])*mb + (one(T)-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (one(T)-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (one(T)-mb)/ΔZ[k]*HT[5]*r_ρc 
        azp = (λ0        / (ΔZ[k]*zt))*mt + (one(T)-mt)/ΔZ[k]*HT[6]*r_ρc 
        dd = (one(T)-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
        rs = (b[i,j,k] - (ss - dd * θ[i,j,k])) * m0
        r[i,j,k] = rs
        @reduce(res = zero(T) + rs*rs)
    end
    return sqrt(res) * cell_count_inv
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
@param [in]  Z    CV境界座標
@param [in]  ΔZ   CV幅
@param [in]  z_range Zループ開始/終了インデクス
@param [in]  HT   熱伝達境界の値

"""
function CalcAX!(ax::AbstractArray{T,3},
                  θ::AbstractArray{T,3},
                  Δh::NTuple{3,T},
                  Δt::T,
                  λ::AbstractArray{T,3},
                  cp::AbstractArray{T,3},
                  m::AbstractArray{T,3},
                  ρ::T,
                  Z::AbstractVector{T},
                  ΔZ::AbstractVector{T},
                  z_range::AbstractVector{<:Integer},
                  HT::AbstractVector{T},
                  par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = inv(dx0*dx0)
    dy2 = inv(dy0*dy0)
    dx1 = inv(dx0)
    dy1 = inv(dy0)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
    ddt = inv(Δt)

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        r_ρc = inv(ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = one(T)-m[i-1,j  ,k  ]
        me = one(T)-m[i+1,j  ,k  ]
        ms = one(T)-m[i  ,j-1,k  ]
        mn = one(T)-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
        zb = (Z[k]-Z[k-1])*mb + (one(T)-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (one(T)-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (one(T)-mb)/ΔZ[k]*HT[5]*r_ρc
        azp = (λ0        / (ΔZ[k]*zt))*mt + (one(T)-mt)/ΔZ[k]*HT[6]*r_ρc
        dd = (one(T)-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
        ax[i,j,k] = (ss - dd*θ[i,j,k]) * m0
    end
end


"""
@brief 前処理
@param [in,out] xx   解ベクトル
@param [in]     bb   RHSベクトル
@param [in]     λ    熱伝導率
@param [in]     cp   比熱
@param [in]     mask マスク配列
@param [in]     ρ    密度
@param [in]     Δh   セル幅
@param [in]     Δt   時間積分幅
@param [in]     smoother  [:none, :gs, :jacobi]
@param [in]     Z    CV境界座標
@param [in]     ΔZ   CV幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値
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
                         Z::AbstractVector{T},
                         ΔZ::AbstractVector{T},
                         z_range::AbstractVector{<:Integer},
                         HT::AbstractVector{T},
                         par::String;
                         scratch::Union{Nothing,AbstractArray{T,3}}=nothing) where {T <: AbstractFloat}
    _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, smoother, Z, ΔZ, z_range, HT, par, scratch)
end

@inline function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{:none}, Z, ΔZ, z_range, HT, par, _)
    mycopy!(xx, bb, par)
    return nothing
end

function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{:gs}, Z, ΔZ, z_range, HT, par, _)
    for _ in 1:PRECONDITIONER_SWEEPS
        rbsor!(xx, λ, cp, bb, mask, ρ, Δh, Δt, one(ρ), Z, ΔZ, z_range, HT, par)
    end
    return nothing
end

function _Preconditioner!(xx::AbstractArray{T,3},
                          bb::AbstractArray{T,3},
                          λ::AbstractArray{T,3},
                          cp::AbstractArray{T,3},
                          mask::AbstractArray{T,3},
                          ρ::T,
                          Δh::NTuple{3,T},
                          Δt::T,
                          ::Val{:jacobi},
                          Z::AbstractVector{T},
                          ΔZ::AbstractVector{T},
                          z_range::AbstractVector{<:Integer},
                          HT::AbstractVector{T},
                          par::String,
                          scratch::Union{Nothing,AbstractArray{T,3}}) where {T <: AbstractFloat}
    scratch === nothing && throw(ArgumentError("Jacobi preconditioner requires scratch workspace"))
    jacobi_preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, Z, ΔZ, z_range, HT, par, scratch)
    return nothing
end

@inline function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{s}, Z, ΔZ, z_range, HT, par, scratch) where {s}
    throw(ArgumentError("Unsupported smoother: $(s)"))
end

function jacobi_preconditioner!(xx::AbstractArray{T,3},
                                bb::AbstractArray{T,3},
                                λ::AbstractArray{T,3},
                                cp::AbstractArray{T,3},
                                mask::AbstractArray{T,3},
                                ρ::T,
                                Δh::NTuple{3,T},
                                Δt::T,
                                Z::AbstractVector{T},
                                ΔZ::AbstractVector{T},
                                z_range::AbstractVector{<:Integer},
                                HT::AbstractVector{T},
                                par::String,
                                scratch::AbstractArray{T,3}) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(xx)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = inv(dx0 * dx0)
    dy2 = inv(dy0 * dy0)
    dx1 = inv(dx0)
    dy1 = inv(dy0)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
    ddt = inv(Δt)
    float_min_T = T(FloatMin)
    ω = T(JACOBI_RELAXATION)
    oneT = one(T)
    zeroT = zero(T)

    for _ in 1:PRECONDITIONER_SWEEPS
        mycopy!(scratch, xx, par)
        @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
            λ0 = λ[i,j,k]
            m0 = mask[i,j,k]
            # Dirichlet/ghost cells are copied through
            if m0 == zeroT
                scratch[i,j,k] = bb[i,j,k]
                continue
            end
            r_ρc = inv(ρ * cp[i,j,k])
            mw = oneT - mask[i-1,j  ,k  ]
            me = oneT - mask[i+1,j  ,k  ]
            ms = oneT - mask[i  ,j-1,k  ]
            mn = oneT - mask[i  ,j+1,k  ]
            mb = mask[i  ,j  ,k-1]
            mt = mask[i  ,j  ,k+1]
            axm = λf(λ[i-1,j,k], λ0, mask[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
            axp = λf(λ[i+1,j,k], λ0, mask[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
            aym = λf(λ[i,j-1,k], λ0, mask[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
            ayp = λf(λ[i,j+1,k], λ0, mask[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
            zb = (Z[k]-Z[k-1])*mb + (oneT-mb)*ΔZ[k]
            zt = (Z[k+1]-Z[k])*m0 + (oneT-m0)*ΔZ[k]
            azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (oneT-mb)/ΔZ[k]*HT[5]*r_ρc
            azp = (λ0        / (ΔZ[k]*zt))*mt + (oneT-mt)/ΔZ[k]*HT[6]*r_ρc
            dd = axp + axm + ayp + aym + azp + azm + ddt
            ss = ( axp * xx[i+1,j  ,k  ] + axm * xx[i-1,j  ,k  ]
                 + ayp * xx[i  ,j+1,k  ] + aym * xx[i  ,j-1,k  ]
                 + azp * xx[i  ,j  ,k+1] + azm * xx[i  ,j  ,k-1] )
            Ax = ss - dd * xx[i,j,k]
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
@param [in]     z_range Zループ開始/終了インデクス
@ret            内積
"""
function Fdot1(x::AbstractArray{T,3}, z_range::AbstractVector{<:Integer}, par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(x)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        @reduce(sum = zero(T) + x[i,j,k] * x[i,j,k])
    end
    return sum
end


"""
@brief 2ベクトルの内積
@param [in]     x    ベクトル
@param [in]     y    ベクトル
@param [in]     z_range Zループ開始/終了インデクス
@ret            内積
"""
function Fdot2(x::AbstractArray{T,3}, y::AbstractArray{T,3}, z_range::AbstractVector{<:Integer}, par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(x)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
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
@param [in]     z_range Zループ開始/終了インデクス
"""
function BiCG1!(p::AbstractArray{T,3},
                r::AbstractArray{T,3},
                q::AbstractArray{T,3},
                beta::T,
                omg::T,
                z_range::AbstractVector{<:Integer},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(p)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
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
@param [in]     z_range Zループ開始/終了インデクス
"""
function Triad!(z::AbstractArray{T,3},
                x::AbstractArray{T,3},
                y::AbstractArray{T,3},
                a::T,
                z_range::AbstractVector{<:Integer},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(z)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
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
@param [in]     z_range Zループ開始/終了インデクス
"""
function BICG2!(z::AbstractArray{T,3},
                x::AbstractArray{T,3},
                y::AbstractArray{T,3},
                a::T,
                b::T,
                z_range::AbstractVector{<:Integer},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(z)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値
@ret                 1セルあたりの残差RMS
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
                Z::AbstractVector{T},
                ΔZ::AbstractVector{T},
                z_range::AbstractVector{<:Integer},
                HT::AbstractVector{T},
                par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = inv(dx0*dx0)
    dy2 = inv(dy0*dy0)
    dx1 = inv(dx0)
    dy1 = inv(dy0)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
    ddt = inv(Δt)

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = θ[i,j,k]
        λ0 = λ[i,j,k]
        r_ρc = inv(ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = one(T) - m[i-1,j  ,k  ]
        me = one(T) - m[i+1,j  ,k  ]
        ms = one(T) - m[i  ,j-1,k  ]
        mn = one(T) - m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
        zb = (Z[k]-Z[k-1])*mb + (one(T)-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (one(T)-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (one(T)-mb)/ΔZ[k]*HT[5]*r_ρc
        azp = (λ0        / (ΔZ[k]*zt))*mt + (one(T)-mt)/ΔZ[k]*HT[6]*r_ρc
        dd = (one(T)-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値
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
                     Z::AbstractVector{T},
                     ΔZ::AbstractVector{T},
                     z_range::AbstractVector{<:Integer},
                     HT::AbstractVector{T},
                     color::Int,
                     par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = inv(dx0*dx0)
    dy2 = inv(dy0*dy0)
    dx1 = inv(dx0)
    dy1 = inv(dy0)
    z_st = Int(first(z_range))
    z_ed = Int(last(z_range))
    ddt = inv(Δt)

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            pp = θ[i,j,k]
            λ0 = λ[i,j,k]
            r_ρc = inv(ρ * cp[i,j,k])
            m0 = m[i,j,k]
            me = one(T) - m[i+1,j  ,k  ]
            mw = one(T) - m[i-1,j  ,k  ]
            mn = one(T) - m[i  ,j+1,k  ]
            ms = one(T) - m[i  ,j-1,k  ]
            mt = m[i  ,j  ,k+1]
            mb = m[i  ,j  ,k-1]
            axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
            axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
            aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
            ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
            zb = (Z[k]-Z[k-1])*mb + (one(T)-mb)*ΔZ[k]
            zt = (Z[k+1]-Z[k])*m0 + (one(T)-m0)*ΔZ[k]
            azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (one(T)-mb)/ΔZ[k]*HT[5]*r_ρc
            azp = (λ0        / (ΔZ[k]*zt))*mt + (one(T)-mt)/ΔZ[k]*HT[6]*r_ρc
            dd = (one(T)-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
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
@param [in]     Z    Z座標
@param [in]     ΔZ   格子幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値
@ret                 セルあたりの残差RMS
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
                Z::AbstractVector{T},
                ΔZ::AbstractVector{T},
                z_range::AbstractVector{<:Integer},
                HT::AbstractVector{T},
                par::String) where {T <: AbstractFloat}
    SZ = size(b)
    res = zero(T)

    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        res += rbsor_core!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, c, par)
    end
    norm_factor = T((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
    return sqrt(res) / norm_factor
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
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     F    ファイルディスクリプタ
@param [in]     tol  収束判定基準
@param [in]     HT   熱伝達境界の値
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
                    Z::AbstractVector{T},
                    ΔZ::AbstractVector{T},
                    z_range::AbstractVector{<:Integer},
                    F,
                    tol::T,
                    HT::AbstractVector{T},
                    par::String;
                    ItrMax::Int=1000) where {T <: AbstractFloat}

    res0 = resSOR(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, par)
    if res0 == zero(T)
        res0 = one(T)
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        #res = sor!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, par) / res0
        res = rbsor!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, par) / res0
        @printf(F, "%10d %24.14E\n", n, Float64(res)) # 時間計測の場合にはコメントアウト
        @printf(stdout, "%10d %24.14E\n", n, Float64(res)) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end


end  # module CommonSolver
