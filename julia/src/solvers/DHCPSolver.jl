"""
DHCPSolver.jl

直接熱伝導問題（DHCP）ソルバー

陰解法による熱伝導方程式の数値解法：
    (I - α∇²)T^(n+1) = T^n + q_boundary
    α = (k·dt) / (ρ·cp·dx²)


主要関数:
- solve_dhcp!: 複数時間ステップCG求解（ホットスタート対応）
"""

module DHCPSolver

using FLoops
using Printf

import ..Commons
using ..Commons: WorkBuffers, get_backend

import ..ThermalProperties
using ..ThermalProperties: set_properties!

import ..BoundaryConditions
using ..BoundaryConditions: BoundaryCondition, BoundaryConditionSet,
            isothermal_bc, heat_flux_bc, adiabatic_bc, convection_bc,
            create_boundary_conditions, apply_boundary_conditions!,
            print_boundary_conditions, ISOTHERMAL, HEAT_FLUX,
            apply_face_boundary!, set_BC_coef

import ..RHSCore
using ..RHSCore: calRHS_core!

import ..CommonSolver
using ..CommonSolver: PBiCGSTAB!, PCG!, solve_linear_system!

import ..AdaptiveTolerance
using ..AdaptiveTolerance: AdaptiveToleranceParams, compute_adaptive_tol

export solve_dhcp!

"""
順問題の境界条件  
Z下面: 断熱、Z上面: 熱流束、側面: 断熱
"""

function set_dhcp_bc_parameters()
    x_minus_bc = adiabatic_bc()
    x_plus_bc  = adiabatic_bc()
    y_minus_bc = adiabatic_bc()
    y_plus_bc  = adiabatic_bc()
    z_minus_bc = adiabatic_bc()
    z_plus_bc  = heat_flux_bc(0.0, true) # 分布を与える > true, 0.0はダミー

    # 境界条件セットを作成
    return create_boundary_conditions(
                              x_minus_bc, x_plus_bc,
                              y_minus_bc, y_plus_bc,
                              z_minus_bc, z_plus_bc)
end


"""
@brief 右辺項b
@param [in,out] wk.b   RHS
@param [in]     dx, dy セル幅
@param [in]     Δt   時間積分幅
@param [in]     dz   CV幅
@param [in]     qsrf 熱流束分布
@param [in]     q_dist 熱流束分布を与える場合true

"""
function calRHS!(
    wk::WorkBuffers,
    dx::T,
    dy::T,
    Δt::T,
    dz::AbstractVector{T},
    qsrf::AbstractArray{T,2},
    bc_set::BoundaryConditionSet,
    distribution::Bool,
    ρ::T,
    par::String
    ) where {T <: AbstractFloat}

    # コア処理（共通部分）: 初期化 + 6面一様境界条件
    calRHS_core!(wk.b, wk.θ, dx, dy, dz, bc_set, par)
    backend = get_backend(par)
    SZ = size(wk.b)

    # DHCP固有: Z上面の熱流束分布
    if distribution == true
        let k = SZ[3]-1, area = dx * dy
            @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
                wk.b[i,j,k] -= qsrf[i-1,j-1] * area
            end
        end
    end

    # DHCP固有: 最終RHS計算（内部熱源項含む）
    ddt = inv(Δt)
    @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1, i in 2:SZ[1]-1
      dz_k = dz[k]
      a_p_0 = ρ * wk.cp[i,j,k] * dx * dy * dz_k * ddt
      wk.b[i,j,k] = -( a_p_0 * wk.θ[i,j,k]
                  + wk.hsrc[i,j,k] * dx * dy * dz_k + wk.b[i,j,k] )
    end

end


"""
DHCP（直接熱伝導問題）ソルバー

各時間ステップで、前ステップ温度から熱物性値計算と反復求解

# 引数
T_prev: 前時刻の温度場 (ni, nj, nk) [K]
q_surface: 表面熱流束 (ni, nj, nt-1) [W/m²] 
work: ワーク配列群 (ni+2, nj+2, nk+2)
nt: 時間ステップ数
ρ: 密度 [kg/m³]
cp_coeffs: 比熱多項式係数 [c0, c1, c2, c3]
k_coeffs: 熱伝導率多項式係数 [k0, k1, k2, k3]
dx, dy, ZC, dz, Δt: 格子・時間パラメータ
rtol, maxiter: 収束パラメータ
verbose: 進捗表示フラグ（デフォルト: false）
par: バックエンド

# 戻り値
T_all: 新時刻の温度場 (ni, nj, nk, nt) [K]
iter_counts: 各時間ステップの反復回数 (nt) [回]
"""

const VALID_EXTRAPOLATION_METHODS = (:none, :linear, :quadratic)

@inline function apply_temporal_initial_guess!(
  θ::AbstractArray{T,3},
  history::AbstractArray{T,4},
  t::Int,
  use_previous_solution::Bool,
  extrapolation::Symbol
) where {T}
  if !use_previous_solution
    fill!(θ, zero(T))
    return
  end

  ni, nj, nk, _ = size(history)
  interior = @view θ[2:ni+1, 2:nj+1, 2:nk+1]

  if extrapolation === :quadratic && t >= 4
    prev1 = @view history[:, :, :, t-1]
    prev2 = @view history[:, :, :, t-2]
    prev3 = @view history[:, :, :, t-3]
    @. interior = 3 * prev1 - 3 * prev2 + prev3
  elseif (extrapolation === :linear || extrapolation === :quadratic) && t >= 3
    prev1 = @view history[:, :, :, t-1]
    prev2 = @view history[:, :, :, t-2]
    @. interior = 2 * prev1 - prev2
  else
    prev1 = @view history[:, :, :, t-1]
    interior .= prev1
  end
end

function solve_dhcp!(
  T_initial::AbstractArray{T,3},
  q_surface::AbstractArray{T,3},
  wk::WorkBuffers,
  nt::Int,
  ρ::T,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  dx::T,
  dy::T,
  ZC::Vector{Float64},
  dz::Vector{Float64},
  dt::T;
  rtol::T=T(1e-6),
  maxiter::Int=20000,
  verbose::Bool=false,
  par::String="sequential",
  T_buffer::Union{Nothing,Array{T,4}}=nothing,
  iter_buffer::Union{Nothing,Vector{Int}}=nothing,
  use_previous_solution::Bool=true,
  extrapolation::Symbol=:none,
  smoother::Symbol=:gs,
  solver::Symbol=:pbicgstab,
  adaptive_tol::Bool=false,
  adaptive_tol_params::Union{Nothing,AdaptiveToleranceParams{T}}=nothing
) where {T <: AbstractFloat}
  ni, nj, nk = size(T_initial)
  N = ni * nj * nk
  Δh = (dx, dy, one(T)) # 1.0はダミー

  T_all = isnothing(T_buffer) ? zeros(T, ni, nj, nk, nt) : T_buffer
  expected_shape = (ni, nj, nk, nt)
  if size(T_all) != expected_shape
    throw(ArgumentError("T_buffer size mismatch: expected $(expected_shape), got $(size(T_all))"))
  end
  fill!(T_all, zero(T))
  @views T_all[:, :, :, 1] .= T_initial

  # 反復回数を記録
  iter_counts = isnothing(iter_buffer) ? zeros(Int, nt) : iter_buffer
  expected_len = nt
  if length(iter_counts) != expected_len
    throw(ArgumentError("iter_buffer length mismatch: expected $(expected_len), got $(length(iter_counts))"))
  end
  fill!(iter_counts, 0)

  if verbose
    println("="^60)
    println("Start DHCP direct solver")
    println("="^60)
    println("格子: $(ni)×$(nj)×$(nk) (N=$(N))")
    println("時間ステップ: $(nt), dt=$(dt)s")
    println("CG許容誤差: rtol=$(rtol), maxiter=$(maxiter)")
    println("="^60)
  end
  
  if extrapolation ∉ VALID_EXTRAPOLATION_METHODS
    throw(ArgumentError("Unsupported extrapolation=:$(extrapolation). Use one of $(VALID_EXTRAPOLATION_METHODS)."))
  end


   # PBICGSTAB 初期値設定
  for k in 1:nk, j in 1:nj, i in 1:ni
    wk.θ[i+1, j+1, k+1] = T_initial[i, j, k]
    wk.mask[i+1, j+1, k+1] = one(T)
  end

  # 温度場から初期値の熱物性値計算
  set_properties!(T_initial, wk.cp, wk.λ, cp_coeffs, k_coeffs)

  # Boundary condition
  bc_set = set_dhcp_bc_parameters()
  print_boundary_conditions(bc_set)
  apply_boundary_conditions!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set)

  # HC配列を生成
  HC = set_BC_coef(bc_set)

# 時間積分ループ
  for t in 2:nt
    step_start = time()

    # 前ステップ温度から熱物性値計算
    set_properties!(@view(T_all[:, :, :, t-1]), wk.cp, wk.λ, cp_coeffs, k_coeffs)

    apply_temporal_initial_guess!(wk.θ, T_all, t, use_previous_solution, extrapolation)

    # Z+方向のみ時間とともに更新
    apply_face_boundary!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set.z_plus, :z_plus)

    # work.b (RHS)の計算
    calRHS!(wk, dx, dy, dt, dz,
      @view(q_surface[:, :, t-1]),
      bc_set, true, ρ, par)

    # 適応的収束基準の計算
    current_tol = rtol
    if adaptive_tol && !isnothing(adaptive_tol_params)
      current_tol = compute_adaptive_tol(
        @view(T_all[:, :, :, t-1]),  # 前ステップ
        t >= 3 ? @view(T_all[:, :, :, t-2]) : @view(T_all[:, :, :, t-1]),  # 2ステップ前（初期はt-1を再利用）
        t,
        rtol,
        adaptive_tol_params
      )
      if verbose
        println("  [t=$(t)] Adaptive tol: $(current_tol) (default: $(rtol))")
      end
    end

    isconverged, itr, res0 = solve_linear_system!(wk, Δh, dt, ZC, dz, ρ, HC,
        solver=solver, tol=current_tol, maxItr=maxiter, smoother=smoother, par=par)

    iter_counts[t] = itr
    step_time = time() - step_start

    if verbose
      if isconverged
        println("[t=$(t)/$(nt)] converged: Iteration= $(itr) : Res_0= $(res0) : time=$(step_time)")
      else
        @warn "[t=$(t)/$(nt)] not converged: Iteration=$(itr) : Res_0= $(res0) : time=$(step_time)"
      end
    end

    # 数値異常チェック
    if any(isnan.(wk.θ)) || any(isinf.(wk.θ))
      error("[t=$(t)] 数値異常が発生しました（NaN/Inf検出）")
    end

    # ガイドセルを除いて内点データを返す
    backend = get_backend(par)
    @floop backend for k in 1:nk, j in 1:nj, i in 1:ni
      T_all[i, j, k, t] = wk.θ[i+1, j+1, k+1]
    end

  end

  if verbose
    println("="^60)
    println("DHCP直接ソルバー完了")
    println("  最終温度範囲: $(minimum(T_all)) - $(maximum(T_all)) K")
    println("="^60)
  end

  return T_all, iter_counts
end



end  # module DHCPSolver
