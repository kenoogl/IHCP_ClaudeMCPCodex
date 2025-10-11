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
function mycopy!(a::Array{Float64,3}, b::Array{Float64,3}, par::String)
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
function myfill!(a::Array{Float64,3}, val::Float64, par::String)
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
@param [in]     smoother ["gs", ""]
@param [in]     par    バックエンド（"sequential", "thread"）

収束判定： 有効セル数(Nf)あたりの相対残差ノルム (||r_k|| / ||r_0||)/Nf < tol
一方、IterativeSolvers.jlのcg!関数は相対残差ノルム
@ret            収束/未収束、反復回数、初期残差
"""
function PBiCGSTAB!(wk::WorkBuffers,
                    Δh::Tuple{Float64, Float64, Float64},
                    Δt::Float64, 
                    Z::Vector{Float64},
                    ΔZ::Vector{Float64},
                    z_range::Vector{Int64},
                    HT::Vector{Float64},
                    ρ::Float64;
                    tol::Float64,
                    maxItr::Int64,
                    smoother::String,
                    par::String,
                    verbose::Bool=false)
    SZ = size(wk.θ)
    myfill!(wk.pcg_q, 0.0, par)
    res0 = CalcRK!(wk.pcg_r, wk.θ, wk.b, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, Z, ΔZ, z_range, HT, par)
    if verbose
        println("Inital residual = ", res0)
    end

    # 初期残差がゼロの場合は収束済み（数値安定性対策）
    if res0 ≈ 0.0
        return true, 0, res0
    end

    mycopy!(wk.pcg_r0, wk.pcg_r, par) # wk.pcg_r0 .= wk.pcg_r  #copy!(pcg_r0, pcg_r)

    rho_old::Float64 = 1.0
    alpha::Float64 = 0.0
    omega::Float64  = 1.0
    r_omega::Float64 = -omega
    beta::Float64 = 0.0
    isconverged::Bool = false
    itr::Int = 0

    for k in 1:maxItr
        itr = k
        rho = Fdot2(wk.pcg_r, wk.pcg_r0, z_range, par) # 非計算部分はゼロのこと

        if abs(rho) < FloatMin
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

        myfill!(wk.pcg_p_, 0.0, par)  #fill!(pcg_p_, 0.0)
        Preconditioner!(wk.pcg_p_, wk.pcg_p, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother, Z, ΔZ, z_range, HT, par)

        CalcAX!(wk.pcg_q, wk.pcg_p_, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, Z, ΔZ, z_range, HT, par)
        alpha = rho / Fdot2(wk.pcg_q, wk.pcg_r0, z_range, par)
        r_alpha = -alpha
        Triad!(wk.pcg_s, wk.pcg_q, wk.pcg_r, r_alpha, z_range, par)

        myfill!(wk.pcg_s_, 0.0, par)  #fill!(pcg_s_, 0.0)
        Preconditioner!(wk.pcg_s_, wk.pcg_s, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, smoother, Z, ΔZ, z_range, HT, par);

        CalcAX!(wk.pcg_t_, wk.pcg_s_, Δh, Δt, wk.λ, wk.cp, wk.mask, ρ, Z, ΔZ, z_range, HT, par)

        # 分母ゼロ対策（数値安定性）
        denom = Fdot1(wk.pcg_t_, z_range, par)
        if abs(denom) < FloatMin
            # 分母がゼロに近い場合は数値的に不安定（未収束として扱う）
            isconverged = false
            break
        end
        omega = Fdot2(wk.pcg_t_, wk.pcg_s, z_range, par) / denom
        r_omega = -omega

        BICG2!(wk.θ, wk.pcg_p_, wk.pcg_s_, alpha , omega, z_range, par)

        Triad!(wk.pcg_r, wk.pcg_t_, wk.pcg_s, r_omega, z_range, par)
        res = sqrt(Fdot1(wk.pcg_r, z_range, par))/((SZ[1]-2)*(SZ[2]-2)*(z_range[2]-z_range[1]+1))
        res /= res0

        if res<tol
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
                r::Array{Float64,3},
                θ::Array{Float64,3},
                b::Array{Float64,3},
                λ::Array{Float64,3},
                cp::Array{Float64,3},
                m::Array{Float64,3},
                ρ::Float64,
                Δh::Tuple{Float64, Float64, Float64},
                Δt::Float64,
                Z::Vector{Float64},
                ΔZ::Vector{Float64},
                z_range::Vector{Int64},
                HT::Vector{Float64},
                par::String)
    backend = get_backend(par)
    SZ = size(θ)
    dx0::Float64 = Δh[1]
    dy0::Float64 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]
    ddt = 1.0 / Δt

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        r_ρc = 1.0 / (ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc 
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc 
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc 
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc 
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]*r_ρc 
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]*r_ρc 
        dd = (1.0-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
        rs = (b[i,j,k] - (ss - dd * θ[i,j,k]))* m0
        r[i,j,k] = rs
        @reduce(res = 0.0 + rs*rs)
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
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
function CalcAX!(ax::Array{Float64,3},
                  θ::Array{Float64,3},
                  Δh::Tuple{Float64, Float64, Float64},
                  Δt::Float64,
                  λ::Array{Float64,3},
                  cp::Array{Float64,3},
                  m::Array{Float64,3},
                  ρ::Float64,
                  Z::Vector{Float64},
                  ΔZ::Vector{Float64},
                  z_range::Vector{Int64},
                  HT::Vector{Float64},
                  par::String)
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]
    ddt = 1.0 / Δt

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        λ0 = λ[i,j,k]
        r_ρc = 1.0 / (ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]*r_ρc
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]*r_ρc
        dd = (1.0-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
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
@param [in]     smoother  ["gs", ""]
@param [in]     Z    CV境界座標
@param [in]     ΔZ   CV幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     HT   熱伝達境界の値
"""
function Preconditioner!(xx::Array{Float64,3},
                         bb::Array{Float64,3},
                          λ::Array{Float64,3},
                          cp::Array{Float64,3},
                       mask::Array{Float64,3},
                        ρ::Float64,
                         Δh::Tuple{Float64, Float64, Float64},
                         Δt::Float64,
                   smoother::String,
                   Z, ΔZ,
                   z_range::Vector{Int64},
                   HT::Vector{Float64},
                   par::String)

    LCmax::Int = 5

    if smoother=="gs"
        for _ in 1:LCmax
            #sor!(xx, λ, cp, bb, mask, ρ, Δh, 1.0, Z, ΔZ, z_range, HT, par)
            rbsor!(xx, λ, cp, bb, mask, ρ, Δh, Δt, 1.0, Z, ΔZ, z_range, HT, par)
        end
    else
        xx .= bb  #copy!(xx, bb)
    end
end


"""
@brief ベクトルの内積
@param [in]     x    ベクトル
@param [in]     z_range Zループ開始/終了インデクス
@ret            内積
"""
function Fdot1(x::Array{Float64,3}, z_range::Vector{Int64}, par::String)
    backend = get_backend(par)
    SZ = size(x)
    z_st = z_range[1]
    z_ed = z_range[2]
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        @reduce(sum = 0.0 + x[i,j,k] * x[i,j,k])
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
function Fdot2(x::Array{Float64,3}, y::Array{Float64,3}, z_range::Vector{Int64}, par::String)
    backend = get_backend(par)
    SZ = size(x)
    z_st = z_range[1]
    z_ed = z_range[2]
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        @reduce(sum = 0.0 + x[i,j,k] * y[i,j,k])
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
function BiCG1!(p::Array{Float64,3},
                r::Array{Float64,3},
                q::Array{Float64,3},
                beta::Float64,
                omg::Float64,
                z_range::Vector{Int64},
                par::String)
    backend = get_backend(par)
    SZ = size(p)
    z_st = z_range[1]
    z_ed = z_range[2]
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
function Triad!(z::Array{Float64,3},
                x::Array{Float64,3},
                y::Array{Float64,3},
                a::Float64,
                z_range::Vector{Int64},
                par::String)
    backend = get_backend(par)
    SZ = size(z)
    z_st = z_range[1]
    z_ed = z_range[2]
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
function BICG2!(z::Array{Float64,3},
                x::Array{Float64,3},
                y::Array{Float64,3},
                a::Float64,
                b::Float64,
                z_range::Vector{Int64},
                par::String)
    backend = get_backend(par)
    SZ = size(z)
    z_st = z_range[1]
    z_ed = z_range[2]
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
function resSOR(θ::Array{Float64,3},
                λ::Array{Float64,3},
                cp::Array{Float64,3},
                b::Array{Float64,3},
                m::Array{Float64,3},
                ρ::Float64,
                Δh::Tuple{Float64, Float64, Float64},
                Δt::Float64, 
                ω::Float64,
                Z::Vector{Float64},
                ΔZ::Vector{Float64},
                z_range::Vector{Int64},
                HT::Vector{Float64},
                par::String)
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]
    ddt = 1.0 / Δt

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = θ[i,j,k]
        λ0 = λ[i,j,k]
        r_ρc = 1.0 / (ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]*r_ρc
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]*r_ρc

        dd = (1.0-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
        dp = (((ss-b[i,j,k])/dd - pp)) * m0
        r = (dd + ω*(axm+aym+azm))*dp / ω
        @reduce(res = 0.0 + r*r)
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
end

"""
@brief SOR法
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
@ret                 セルあたりの残差RMS
"""
function sor!(θ::Array{Float64,3},
              λ::Array{Float64,3},
              cp::Array{Float64,3},
              b::Array{Float64,3},
              m::Array{Float64,3},
              ρ::Float64,
              Δh::Tuple{Float64, Float64, Float64},
              Δt::Float64, 
              ω::Float64,
              Z::Vector{Float64},
              ΔZ::Vector{Float64},
              z_range::Vector{Int64},
              HT::Vector{Float64},
              par::String)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]
    ddt = 1.0 / Δt

    res::Float64 = 0.0
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        pp = θ[i,j,k]
        bb = b[i,j,k]
        λ0 = λ[i,j,k]
        r_ρc = 1.0 / (ρ * cp[i,j,k])
        m0 = m[i,j,k]
        mw = 1.0-m[i-1,j  ,k  ]
        me = 1.0-m[i+1,j  ,k  ]
        ms = 1.0-m[i  ,j-1,k  ]
        mn = 1.0-m[i  ,j+1,k  ]
        mb = m[i  ,j  ,k-1]
        mt = m[i  ,j  ,k+1]
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
        zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
        zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
        azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]*r_ρc
        azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]*r_ρc
        dd = (1.0-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
        dp = (((ss-b[i,j,k])/dd - pp)) * m0
        pn = pp + ω * dp
        θ[i,j,k] = pn
        r = (dd + ω*(axm+aym+azm))*dp / ω
        res += r*r
    end

    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(z_ed-z_st+1))
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
function rbsor_core!(θ::Array{Float64,3},
                     λ::Array{Float64,3},
                     cp::Array{Float64,3},
                     b::Array{Float64,3},
                     m::Array{Float64,3},
                     ρ::Float64,
                     Δh::Tuple{Float64, Float64, Float64},
                     Δt::Float64, 
                     ω::Float64,
                     Z::Vector{Float64},
                     ΔZ::Vector{Float64},
                     z_range::Vector{Int64},
                     HT::Vector{Float64},
                     color::Int,
                     par::String)
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    dx2 = 1.0 / (dx0*dx0)
    dy2 = 1.0 / (dy0*dy0)
    dx1 = 1.0 / dx0
    dy1 = 1.0 / dy0
    z_st = z_range[1]
    z_ed = z_range[2]
    ddt = 1.0 / Δt

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1
        @simd for i in 2+mod(k+j+color,2):2:SZ[1]-1
            pp = θ[i,j,k]
            λ0 = λ[i,j,k]
            r_ρc = 1.0 / (ρ * cp[i,j,k])
            m0 = m[i,j,k]
            me = 1.0-m[i+1,j  ,k  ]
            mw = 1.0-m[i-1,j  ,k  ]
            mn = 1.0-m[i  ,j+1,k  ]
            ms = 1.0-m[i  ,j-1,k  ]
            mt = m[i  ,j  ,k+1]
            mb = m[i  ,j  ,k-1]
            axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dx2 + mw*dx1*HT[1]*r_ρc
            axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dx2 + me*dx1*HT[2]*r_ρc
            aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dy2 + ms*dy1*HT[3]*r_ρc
            ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dy2 + mn*dy1*HT[4]*r_ρc
            zb = (Z[k]-Z[k-1])*mb + (1.0-mb)*ΔZ[k] # 境界の半セル処理
            zt = (Z[k+1]-Z[k])*m0 + (1.0-m0)*ΔZ[k]
            azm = (λ[i,j,k-1]/ (ΔZ[k]*zb))*mb + (1.0-mb)/ΔZ[k]*HT[5]*r_ρc
            azp = (λ0        / (ΔZ[k]*zt))*mt + (1.0-mt)/ΔZ[k]*HT[6]*r_ρc
            dd = (1.0-m0) + ( axp + axm + ayp + aym + azp + azm + ddt)*m0
            ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
            dp = (((ss-b[i,j,k])/dd - pp)) * m0
            θ[i,j,k] = pp + ω * dp
            r = (dd + ω*(axm+aym+azm))*dp / ω
            @reduce(res = 0.0 + r*r)
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
function rbsor!(θ::Array{Float64,3},
                λ::Array{Float64,3},
                cp::Array{Float64,3},
                b::Array{Float64,3},
                mask::Array{Float64,3},
                ρ::Float64,
                Δh::Tuple{Float64, Float64, Float64},
                Δt::Float64, 
                ω::Float64,
                Z::Vector{Float64},
                ΔZ::Vector{Float64},
                z_range::Vector{Int64},
                HT::Vector{Float64},
                par::String)
    SZ = size(b)
    res::Float64 = 0.0

    # 2色のマルチカラー(Red&Black)のセットアップ
    for c in 0:1
        res += rbsor_core!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, c, par)
    end
    return sqrt(res)/((SZ[1]-2)*(SZ[2]-2)*(SZ[3]-2))
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
function solveSOR!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ,
                    z_range::Vector{Int64},
                    F, tol,
                    HT::Vector{Float64},
                    par::String;
                    ItrMax::Int=1000)

    res0 = resSOR(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, par)
    if res0==0.0
        res0 = 1.0
    end
    println("Inital residual = ", res0)

    n = 0
    for n in 1:ItrMax
        #res = sor!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, par) / res0
        res = rbsor!(θ, λ, cp, b, mask, ρ, Δh, Δt, ω, Z, ΔZ, z_range, HT, par) / res0
        @printf(F, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        @printf(stdout, "%10d %24.14E\n", n, res) # 時間計測の場合にはコメントアウト
        if res < tol
            println("Converged at ", n)
            return
        end
    end
end


end  # module CommonSolver