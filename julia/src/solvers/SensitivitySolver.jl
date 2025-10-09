"""
SensitivitySolver.jl

感度問題ソルバー


主要関数:
- solve_sensitivity!: 複数時間ステップCG求解（ホットスタート対応）
"""

module SensitivitySolver

using FLoops

import ..Commons
using ..Commons: WorkBuffers, λf, get_backend

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
using ..CommonSolver: PBiCGSTAB!

export solve_sensitivity!

"""
順問題の境界条件  
Z下面: 断熱、Z上面: 熱流束、側面: 断熱
"""

function set_sensitivity_bc_parameters(nk::Int)
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
                              z_minus_bc, z_plus_bc, nk)
end


"""
@brief 右辺項b
@param [in,out] wk.b   RHS
@param [in]     HF  熱流束境界の値
@param [in]     dx, dy セル幅
@param [in]     Δt   時間積分幅
@param [in]     ΔZ   CV幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     qsrf 熱流束分布
@param [in]     q_dist 熱流束分布を与える場合true

"""
function calRHS!(wk::WorkBuffers,
    HF::Vector{Float64},
    dx::Float64,
    dy::Float64,
    Δt::Float64,
    ΔZ::Vector{Float64},
    z_range::Vector{Int64},
    qsrf::AbstractArray{Float64,2},
    distribution::Bool,
    ρ::Float64,
    par::String
    )

    # コア処理（共通部分）: 初期化 + 6面一様境界条件
    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)
    backend = get_backend(par)

    # Sensitivity固有: Z上面の熱流束分布（DHCPと同じ処理）
    if distribution == true
        let k = z_ed, a = 1.0 / ΔZ[z_ed]
            @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
                wk.b[i,j,k] -= qsrf[i-1,j-1] * a
            end
        end
    end

    # Sensitivity固有: 最終RHS計算（内部熱源項なし）
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        wk.b[i,j,k] = -( ddt * wk.θ[i,j,k]
                        + wk.b[i,j,k] / (ρ * wk.cp[i,j,k])  )
    end

end


"""
Sensitivity（感度問題）ソルバー

微小な熱流束の変化が裏面S2温度に及ぼす影響(感度)を計算

# 引数
T_initial: 初期温度場 (ni, nj, nk) [K]
q_surface: 表面熱流束 (ni, nj, nt-1) [W/m²] 
work: ワーク配列群 (ni+2, nj+2, nk+2)
nt: 時間ステップ数
ρ: 密度 [kg/m³]
cp_coeffs: 比熱多項式係数 [c0, c1, c2, c3]
k_coeffs: 熱伝導率多項式係数 [k0, k1, k2, k3]
dx, dy, Z, ΔZ, Δt: 格子・時間パラメータ
rtol, maxiter: 収束パラメータ
verbose: 進捗表示フラグ（デフォルト: false）
par: バックエンド

# 戻り値
T_all: 新時刻の温度場 (ni, nj, nk, nt) [K] 
"""

function solve_sensitivity!(
  T_initial::Array{Float64,3},
  q_surf::Array{Float64,3},
  wk::WorkBuffers,
  nt::Int,
  ρ::Float64,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  dx::Float64,
  dy::Float64,
  Z::Vector{Float64},
  ΔZ::Vector{Float64},
  dt::Float64;
  rtol=1e-6,
  maxiter=20000,
  verbose=false,
  par::String="sequential"
)
  ni, nj, nk = size(T_initial)
  N = ni * nj * nk
  Δh = (dx, dy, 1.0) # 1.0はダミー

  T_all = zeros(Float64, ni, nj, nk, nt)
  T_all[:, :, :, 1] = T_initial


  if verbose
    println("="^60)
    println("Start Sensitivity solver")
    println("="^60)
    println("Grid: $(ni)×$(nj)×$(nk) (N=$(N))")
    println("Time steps: $(nt), dt=$(dt)s")
    println("CG tolerance: rtol=$(rtol), maxiter=$(maxiter)")
    println("="^60)
  end

   # PBICGSTAB 初期値設定
  for k in 1:nk, j in 1:nj, i in 1:ni
    wk.θ[i+1, j+1, k+1] = T_initial[i, j, k]
    wk.mask[i+1, j+1, k+1] = 1.0
  end

  # 温度場から初期値の熱物性値計算
  set_properties!(T_initial, wk.cp, wk.λ, cp_coeffs, k_coeffs)

  # Boundary condition
  z_range, bc_set = set_sensitivity_bc_parameters(nk)
  print_boundary_conditions(bc_set)
  apply_boundary_conditions!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set)
  HF, HT = set_BC_coef(bc_set) # 時間変化なし

  SZ = (ni+2, nj+2, nk+2)
  Nf = (SZ[1]-2)*(SZ[2]-2)*(z_range[2]-z_range[1]+1) # cg!と判定基準を合わせるため

# 時間積分ループ
  for t in 2:nt
    # 前ステップ温度から熱物性値計算
    set_properties!(@view(T_all[:, :, :, t-1]), wk.cp, wk.λ, cp_coeffs, k_coeffs)

    # Z+方向のみ時間とともに更新
    apply_face_boundary!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set.z_plus, :z_plus)

    # work.b (RHS)の計算
    calRHS!(wk, HF, dx, dy, dt, ΔZ, z_range,
      @view(q_surf[:, :, t-1]),
      true, ρ, par)
    

    if verbose
      isconverged, itr, res0 = PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol*Nf, maxItr=maxiter, smoother="", par=par)
      if isconverged
        println("[t=$(t)/$(nt)] CG converged: $(itr) iterations, initial residual: $(res0)")
      else
        @warn "[t=$(t)/$(nt)] CG not converged: $(itr) iterations, initial residual: $(res0)"
      end
    else
      PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol*Nf, maxItr=maxiter, smoother="", par=par)
    end

    # Numerical anomaly check
    if any(isnan.(wk.θ)) || any(isinf.(wk.θ))
      error("[t=$(t)] Numerical anomaly detected (NaN/Inf)")
    end

    # ガイドセルを除いて内点データを返す
    for k in 1:nk, j in 1:nj, i in 1:ni
      T_all[i, j, k, t] = wk.θ[i+1, j+1, k+1]
    end

  end

  if verbose
    println("="^60)
    println("Sensitivity solver completed")
    println("  Final temperature range: $(minimum(T_all)) - $(maximum(T_all)) K")
    println("="^60)
  end

  return T_all
end



end  # module DHCPSolver