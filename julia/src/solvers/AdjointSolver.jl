"""
AdjointSolver.jl

Phase 3: 随伴問題（Adjoint）ソルバー

随伴方程式による感度解析（後退時間積分）：
    (I - α∇²)λ^(n-1) = λ^n + 2·(T_cal - Y_obs)
    α = (k·dt) / (ρ·cp·dx²)

時間反転ループ:
    for t in (nt-1):-1:1

残差注入（底面 k=1）:
    rhs[k=1] += 2.0 * (T_cal[t+1, :, :, 1] - Y_obs[t+1, :, :]) * dx * dy

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- coeffs_and_rhs_building_Adjoint() (1228-1270行)
- assemble_A_Adjoint() (1272-1294行)
- multiple_time_step_solver_Adjoint() (1297-1364行)

主要関数:
- solve_adjoint_mf!: 複数時間ステップPBICGSTAB!求解（マトリクスフリー版、後退時間積分、ホットスタート対応）
"""

module AdjointSolver

using LinearAlgebra
using FLoops
using Statistics

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

export solve_adjoint_mf!


"""
逆問題の境界条件  
Z下面: 値指定、Z上面: ノイマン、側面: ノイマン
"""

function set_adjoint_bc_parameters(nk::Int)
    x_minus_bc = adiabatic_bc()
    x_plus_bc  = adiabatic_bc()
    y_minus_bc = adiabatic_bc()
    y_plus_bc  = adiabatic_bc()
    z_minus_bc = heat_flux_bc(0.0, true) # 分布を与える, 0.0はダミー
    z_plus_bc  = adiabatic_bc()

    # 境界条件セットを作成
    return create_boundary_conditions(
                              x_minus_bc, x_plus_bc,
                              y_minus_bc, y_plus_bc,
                              z_minus_bc, z_plus_bc, nk)
end


"""
@brief 右辺項b
@param [in,out] wk.b   RHS
@param [in]     Tsrf   底面の計算温度
@param [in]     Yobs   底面の観測値
@param [in]     HF  熱流束境界の値
@param [in]     dx, dy セル幅
@param [in]     Δt   時間積分幅
@param [in]     ΔZ   CV幅
@param [in]     z_range Zループ開始/終了インデクス
@param [in]     distribution 分布を与える場合true

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
    solve_adjoint_mf!(T_cal, Y_obs, nt, ρ, cp_coeffs, k_coeffs,
                   dx, dy, dz, dz_b, dz_t, dt; rtol=1e-8, maxiter=1000, verbose=false)
      -> (λ_all, cg_iters)

複数時間ステップ随伴ソルバー（後退時間積分）

時間反転ループで随伴方程式を解く:
    for t in (nt-1):-1:1

対応Pythonコード: multiple_time_step_solver_Adjoint() (1297-1364行)

アルゴリズム:
  1. 終端条件: λ[nt] = 0.0
  2. 後退時間ループ（nt-1 → 1）:
     a. 前ステップ（時間的に後）の随伴場を初期値とする
     b. T_cal[t]から熱物性値を計算
     c. 係数とRHS構築（残差注入: T_cal[t,:,:,1] vs Y_obs[t]）
     d. 疎行列組み立て
     e. 対角前処理CG法で求解
     f. ホットスタート更新（x0 = x）

Args:
  T_cal: DHCP計算温度場 (ni, nj, nk, nt) [K] ※Phase 2.2: 時間次元を最後に配置
  Y_obs: 観測温度（底面） (ni, nj, nt) [K] ※Phase 2.2: 時間次元を最後に配置
  nt: 時間ステップ数
  ρ: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数 [c0, c1, c2, c3]
  k_coeffs: 熱伝導率多項式係数 [k0, k1, k2, k3]
  dx, dy: x, y方向格子幅 [m]
  dz: z方向格子幅配列 (nk,) [m]
  dz_b: 下側界面距離 (nk,) [m]
  dz_t: 上側界面距離 (nk,) [m]
  dt: 時間刻み [s]
  rtol: CG相対許容誤差（デフォルト: 1e-8）
  maxiter: CG最大反復回数（デフォルト: 1000）
  verbose: 詳細出力フラグ（デフォルト: false）

Returns:
  λa_all: 随伴場時系列 (ni, nj, nk, nt) ※熱伝導率のλと混同するのでλaとする
  cg_iters: CG反復回数履歴 (nt-1,)
"""
function solve_adjoint_mf!(
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
    println("Adjoint求解開始（後退時間積分）")
    println("  格子: $ni×$nj×$nk, N=$N")
    println("  時間ステップ: $nt")
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
    # 次ステップ（時間的に後）の随伴場を初期値とする（wk.θに保持されているためホットスタート）

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
        println("[t=$(t)/$(nt)] CG収束: $(itr)回 初期残差: $(res0)")
      else
        @warn "[t=$(t)/$(nt)] CG未収束: $(itr)回 初期残差: $(res0)"
      end
      # 反復回数記録
      cg_iters[t] = itr
    else
      PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol, maxItr=maxiter, smoother="", par=par)
    end

    # 結果保存 ガイドセルを除いて内点データを返す
    backend = get_backend(par)
    @floop backend for k in 1:nk, j in 1:nj, i in 1:ni
      λa_all[i, j, k, t] = wk.θ[i+1, j+1, k+1]
    end

    if verbose && (t % 10 == 0 || t == nt-1)
      println("  t=$t: CG収束 (iter=$(cg_iters[t]))")
    end
  end

  if verbose
    println("Adjoint求解完了")
    println("  平均CG反復回数: $(mean(cg_iters))")
    println("  随伴場範囲: [$(minimum(λa_all)), $(maximum(λa_all))]")
  end

  return λa_all, cg_iters
end

end # module AdjointSolver