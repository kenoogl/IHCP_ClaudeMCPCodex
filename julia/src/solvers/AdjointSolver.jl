"""
AdjointSolver.jl

Phase 3: 随伴問題（Adjoint）ソルバー（マトリクスフリー版）

随伴方程式による感度解析（後退時間積分）：
    (I - α∇²)λ^(n-1) = λ^n + 2·(T_cal - Y_obs)
    α = (k·dt) / (ρ·cp·dx²)

時間反転ループ:
    for t in (nt-1):-1:1

残差注入（底面 k=1）:
    rhs[k=1] += 2.0 * (T_cal[t+1, :, :, 1] - Y_obs[t+1, :, :]) * dx * dy

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- multiple_time_step_solver_Adjoint() (1297-1364行)

主要関数:
- solve_adjoint_mf!: マトリクスフリー版PBICGSTAB!による後退時間積分
- calRHS!: RHSベクトル直接計算（残差注入含む）
"""

module AdjointSolver

using FLoops

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
using ..CommonSolver: PBiCGSTAB!

export solve_adjoint!


"""
逆問題の境界条件  
Z下面: 値指定、Z上面: ノイマン、側面: ノイマン
"""

function set_adjoint_bc_parameters(nk::Int)
    x_minus_bc = adiabatic_bc()
    x_plus_bc  = adiabatic_bc()
    y_minus_bc = adiabatic_bc()
    y_plus_bc  = adiabatic_bc()
    z_minus_bc = isothermal_bc(0.0, true) # 分布を与える, 0.0はダミー
    z_plus_bc  = adiabatic_bc()

    # 境界条件セットを作成
    return create_boundary_conditions(
                              x_minus_bc, x_plus_bc,
                              y_minus_bc, y_plus_bc,
                              z_minus_bc, z_plus_bc, nk)
end


"""
    calRHS!(wk, Tsrf, Yobs, HF, dx, dy, Δt, ΔZ, z_range, distribution, ρ, par)

随伴問題用右辺項ベクトルbを計算（残差注入含む）

# 引数
- `wk::WorkBuffers`: ワークバッファ（wk.bを更新）
- `Tsrf::AbstractArray{Float64,2}`: 底面の計算温度 (ni, nj) [K]
- `Yobs::AbstractArray{Float64,2}`: 底面の観測温度 (ni, nj) [K]
- `HF::Vector{Float64}`: 熱流束境界条件係数（6面分）
- `dx::Float64`: x方向セル幅 [m]
- `dy::Float64`: y方向セル幅 [m]
- `Δt::Float64`: 時間刻み幅 [s]
- `ΔZ::Vector{Float64}`: z方向CV幅配列 (nk+1,) [m]
- `z_range::Vector{Int64}`: Zループ範囲 [開始, 終了]
- `distribution::Bool`: 分布を与える場合true
- `ρ::Float64`: 密度 [kg/m³]
- `par::String`: 並列化バックエンド（"sequential" / "thread"）
"""
function calRHS!(
    wk::WorkBuffers,
    Tsrf::AbstractArray{Float64,2},
    Yobs::AbstractArray{Float64,2},
    HF::Vector{Float64},
    dx::Float64,
    dy::Float64,
    Δt::Float64,
    ΔZ::Vector{Float64},
    z_range::Vector{Int64},
    distribution::Bool,
    ρ::Float64,
    par::String
    )

    # コア処理（共通部分）: 初期化 + 6面一様境界条件
    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)
    backend = get_backend(par)

    # Adjoint固有: Z下面の残差注入
    if distribution == true
        let k = z_st, a = 2.0 / ΔZ[z_st]
            @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
                wk.b[i,j,k] += (Tsrf[i-1,j-1]-Yobs[i-1,j-1]) * a
            end
        end
    end

    # Adjoint固有: 最終RHS計算（内部熱源項なし）
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        wk.b[i,j,k] = -( ddt * wk.θ[i,j,k]
                        + wk.b[i,j,k] / (ρ * wk.cp[i,j,k])  )
    end

end




"""
    solve_adjoint!(T_cal, Y_obs, wk, nt, ρ, cp_coeffs, k_coeffs, dx, dy, Z, ΔZ, dt;
                   rtol=1e-8, maxiter=1000, verbose=false, par="sequential") -> (λa_all, cg_iters)

随伴ソルバー（後退時間積分、マトリクスフリー版）

時間反転ループで随伴方程式を解く:
    for t in (nt-1):-1:1

対応Pythonコード: multiple_time_step_solver_Adjoint() (1297-1364行)

アルゴリズム:
  1. 終端条件: λ[nt] = 0.0
  2. 後退時間ループ（nt-1 → 1）:
     a. 前ステップ（時間的に後）の随伴場を初期値とする
     b. T_cal[t]から熱物性値を計算
     c. RHS構築（残差注入: T_cal[t,:,:,1] vs Y_obs[t]）
     d. マトリクスフリーPBICGSTAB!で求解
     e. 結果保存

# 引数
- `T_cal::Array{Float64,4}`: DHCP計算温度場 (ni, nj, nk, nt) [K]
- `Y_obs::Array{Float64,3}`: 観測温度（底面） (ni, nj, nt) [K]
- `wk::WorkBuffers`: ワーク配列群 (ni+2, nj+2, nk+2)
- `nt::Int`: 時間ステップ数
- `ρ::Float64`: 密度 [kg/m³]
- `cp_coeffs::Vector{Float64}`: 比熱多項式係数 [c0, c1, c2, c3]
- `k_coeffs::Vector{Float64}`: 熱伝導率多項式係数 [k0, k1, k2, k3]
- `dx::Float64`: x方向格子幅 [m]
- `dy::Float64`: y方向格子幅 [m]
- `Z::Vector{Float64}`: z方向座標配列 (nk+2,) [m]
- `ΔZ::Vector{Float64}`: z方向格子幅配列 (nk+1,) [m]
- `dt::Float64`: 時間刻み [s]
- `rtol::Float64`: CG相対許容誤差（デフォルト: 1e-8）
- `maxiter::Int`: CG最大反復回数（デフォルト: 1000）
- `verbose::Bool`: 詳細出力フラグ（デフォルト: false）
- `par::String`: 並列化バックエンド（デフォルト: "sequential"）

# 戻り値
- `λa_all::Array{Float64,4}`: 随伴場時系列 (ni, nj, nk, nt) ※熱伝導率のλと混同するのでλaとする
- `cg_iters::Vector{Int}`: CG反復回数履歴 (nt-1,)
"""
function solve_adjoint!(
  T_cal::Array{Float64,4},
  Y_obs::Array{Float64,3},
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
  rtol::Float64=1e-8,
  maxiter::Int=1000,
  verbose::Bool=false,
  par::String="sequential"
)
  ni, nj, nk = size(T_cal[:, :, :, 1])
  N = ni * nj * nk
  Δh = (dx, dy, 1.0) # 1.0はダミー

  # 随伴場の初期化
  λa_all = zeros(Float64, ni, nj, nk, nt)
  λa_all[:, :, :, nt] .= 0.0  # 終端条件（Pythonオリジナル1322行）

  cg_iters = zeros(Int, nt-1)

  if verbose
    println("Starting Adjoint solver (backward time integration)")
    println("  Grid: $ni×$nj×$nk, N=$N")
    println("  Time steps: $nt")
    println("  CG: rtol=$rtol, maxiter=$maxiter")
  end
  

  # PBICGSTAB 初期値設定
  for k in 1:nk, j in 1:nj, i in 1:ni
    wk.θ[i+1, j+1, k+1] = λa_all[i, j, k, nt]
    wk.mask[i+1, j+1, k+1] = 1.0
  end

  # 温度場から初期値の熱物性値計算
  set_properties!(@view(T_cal[:, :, :, nt]), wk.cp, wk.λ, cp_coeffs, k_coeffs)

  # Boundary condition
  z_range, bc_set = set_adjoint_bc_parameters(nk)
  print_boundary_conditions(bc_set)
  apply_boundary_conditions!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set)
  HF, HT = set_BC_coef(bc_set) # 時間変化なし


  # 後退時間ループ（Pythonオリジナル1328行: range(nt-2, -1, -1)）
  for t in (nt-1):-1:1
    # 次ステップ（時間的に後）の随伴場を初期値とする
    λa_initial = @view λa_all[:, :, :, t+1]

    # 温度場から熱物性値計算（Pythonオリジナル1332行: T_cal[t]）
    set_properties!(@view(T_cal[:, :, :, t]), wk.cp, wk.λ, cp_coeffs, k_coeffs)

    # Z-方向のみ時間とともに更新
    apply_face_boundary!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set.z_minus, :z_minus)

    # work.b (RHS)の計算
    # 底面（物理座標 k=1）の温度と観測値を使用
    calRHS!(wk, @view(T_cal[:, :, 1, t]), @view(Y_obs[:, :, t]), HF, dx, dy, dt, ΔZ, z_range,
      true, ρ, par)

    if verbose
      isconverged, itr, res0 = PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol, maxItr=maxiter, smoother="", par=par)
      if isconverged
        println("[t=$(t)/$(nt)] CG converged: $(itr) iterations, initial residual: $(res0)")
      else
        @warn "[t=$(t)/$(nt)] CG not converged: $(itr) iterations, initial residual: $(res0)"
      end
      # 反復回数記録
      cg_iters[t] = itr
    else
      PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol, maxItr=maxiter, smoother="", par=par)
    end

    # 結果保存 ガイドセルを除いて内点データを返す
    for k in 1:nk, j in 1:nj, i in 1:ni
      λa_all[i, j, k, t] = wk.θ[i+1, j+1, k+1]
    end

    if verbose && (t % 10 == 0 || t == nt-1)
      println("  t=$t: CG converged (iter=$(cg_iters[t]))")
    end
  end

  if verbose
    println("Adjoint solver completed")
    println("  Average CG iterations: $(mean(cg_iters))")
    println("  Adjoint field range: [$(minimum(λa_all)), $(maximum(λa_all))]")
  end

  return λa_all, cg_iters
end

end # module AdjointSolver