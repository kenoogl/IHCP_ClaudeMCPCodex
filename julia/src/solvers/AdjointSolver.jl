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
- build_adjoint_system!: 係数とRHS構築（7点ステンシル + 残差注入）
- assemble_adjoint_matrix: CSC疎行列組み立て（DHCPと同じ構造）
- solve_adjoint!: 複数時間ステップCG求解（後退時間積分、ホットスタート対応）
"""

module AdjointSolver

using LinearAlgebra
using SparseArrays
using IterativeSolvers

# Phase 1モジュールの読み込み
include("../ThermalProperties.jl")
using .ThermalProperties

export build_adjoint_system!, assemble_adjoint_matrix, solve_adjoint!, adjoint_index


"""
    adjoint_index(i, j, k, ni, nj) -> Int

Fortran順序（列優先）でのグローバルインデックス計算

DHCPSolverのdhcp_index()と同じ実装（随伴方程式も同じインデックス）

Args:
  i: x方向インデックス (1:ni)
  j: y方向インデックス (1:nj)
  k: z方向インデックス (1:nk)
  ni, nj: 格子点数

Returns:
  p: グローバルインデックス (1:N)
"""
@inline function adjoint_index(i::Int, j::Int, k::Int, ni::Int, nj::Int)
  return i + (j - 1) * ni + (k - 1) * ni * nj
end


"""
    build_adjoint_system!(λ_initial, T_cal_bottom, Y_obs, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt)
      -> (a_w, a_e, a_s, a_n, a_b, a_t, a_p, b)

随伴方程式の係数とRHSベクトル構築

随伴方程式（DHCPと同じ行列構造、RHSが異なる）:
    A·λ^(n-1) = b
    b = (ρ·cp·V/dt)·λ^n + 残差注入項（底面のみ）

残差注入（底面 k=1）:
    b[i,j,1] += 2.0 * (T_cal[i,j] - Y_obs[i,j]) * dx * dy

対応Pythonコード: coeffs_and_rhs_building_Adjoint() (1228-1270行)

Args:
  λ_initial: 次ステップ随伴場 (ni, nj, nk) - 時間反転なので「次」=時間的に後
  T_cal_bottom: DHCP計算温度（底面） (ni, nj) [K]
  Y_obs: 観測温度（底面） (ni, nj) [K]
  rho: 密度 [kg/m³]
  cp: 比熱場 (ni, nj, nk) [J/(kg·K)]
  k: 熱伝導率場 (ni, nj, nk) [W/(m·K)]
  dx, dy: x, y方向格子幅 [m]
  dz: z方向格子幅配列 (nk,) [m]
  dz_b: 下側界面距離 (nk,) [m]
  dz_t: 上側界面距離 (nk,) [m]
  dt: 時間刻み [s]

Returns:
  a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列 (N,)
  b: RHSベクトル (N,)
"""
function build_adjoint_system!(
  λ_initial::Array{Float64,3},
  T_cal_bottom::Matrix{Float64},
  Y_obs::Matrix{Float64},
  rho::Float64,
  cp::Array{Float64,3},
  k::Array{Float64,3},
  dx::Float64,
  dy::Float64,
  dz::Vector{Float64},
  dz_b::Vector{Float64},
  dz_t::Vector{Float64},
  dt::Float64
)
  ni, nj, nk = size(λ_initial)
  N = ni * nj * nk

  # 係数配列の初期化
  a_w = zeros(Float64, N)
  a_e = zeros(Float64, N)
  a_s = zeros(Float64, N)
  a_n = zeros(Float64, N)
  a_b = zeros(Float64, N)
  a_t = zeros(Float64, N)
  a_p = zeros(Float64, N)
  b   = zeros(Float64, N)

  # 全格子点ループ
  for k_idx in 1:nk
    dz_k = dz[k_idx]
    dz_t_k = dz_t[k_idx]
    dz_b_k = dz_b[k_idx]

    for j in 1:nj
      for i in 1:ni
        # グローバルインデックス
        p = adjoint_index(i, j, k_idx, ni, nj)

        # 中心点熱伝導率
        k_p = k[i, j, k_idx]

        # 時間項（蓄熱項、DHCPと同じ）
        a_p_0 = rho * cp[i, j, k_idx] * dx * dy * dz_k / dt

        # 6方向の熱伝導係数（調和平均、DHCPと同じ）
        # 西（y方向負）
        if j > 1
          k_w = k[i, j-1, k_idx]
          a_w[p] = (2.0 * k_p * k_w / (k_p + k_w)) * dy * dz_k / dx
        else
          a_w[p] = 0.0
        end

        # 東（y方向正）
        if j < nj
          k_e = k[i, j+1, k_idx]
          a_e[p] = (2.0 * k_p * k_e / (k_p + k_e)) * dy * dz_k / dx
        else
          a_e[p] = 0.0
        end

        # 南（x方向負）
        if i > 1
          k_s = k[i-1, j, k_idx]
          a_s[p] = (2.0 * k_p * k_s / (k_p + k_s)) * dx * dz_k / dy
        else
          a_s[p] = 0.0
        end

        # 北（x方向正）
        if i < ni
          k_n = k[i+1, j, k_idx]
          a_n[p] = (2.0 * k_p * k_n / (k_p + k_n)) * dx * dz_k / dy
        else
          a_n[p] = 0.0
        end

        # 下（z方向負）
        if k_idx > 1
          k_b = k[i, j, k_idx-1]
          a_b[p] = (2.0 * k_p * k_b / (k_p + k_b)) * dx * dy / dz_b_k
        else
          a_b[p] = 0.0
        end

        # 上（z方向正）
        if k_idx < nk
          k_t = k[i, j, k_idx+1]
          a_t[p] = (2.0 * k_p * k_t / (k_p + k_t)) * dx * dy / dz_t_k
        else
          a_t[p] = 0.0
        end

        # 対角項（DHCPと同じ）
        a_p[p] = a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p] + a_p_0

        # RHS（時間項）
        rhs = a_p_0 * λ_initial[i, j, k_idx]

        # 残差注入（底面 k=1 のみ、Pythonオリジナル1266-1267行）
        if k_idx == 1
          residual = T_cal_bottom[i, j] - Y_obs[i, j]
          rhs += 2.0 * residual * dx * dy
        end

        b[p] = rhs
      end
    end
  end

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b
end


"""
    assemble_adjoint_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p) -> SparseMatrixCSC

随伴方程式の疎行列組み立て（CSC形式）

7点ステンシル（対角+6方向）のCSC疎行列を構築
DHCPと同じ行列構造なので、係数の符号のみ注意

対応Pythonコード: assemble_A_Adjoint() (1272-1294行)

疎行列構造:
  - 対角: a_p
  - オフセット-1: -a_s (南)
  - オフセット+1: -a_n (北)
  - オフセット-ni: -a_w (西)
  - オフセット+ni: -a_e (東)
  - オフセット-ni*nj: -a_b (下)
  - オフセット+ni*nj: -a_t (上)

Args:
  ni, nj, nk: 格子点数
  a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列 (N,)

Returns:
  A: CSC疎行列 (N, N)
"""
function assemble_adjoint_matrix(
  ni::Int,
  nj::Int,
  nk::Int,
  a_w::Vector{Float64},
  a_e::Vector{Float64},
  a_s::Vector{Float64},
  a_n::Vector{Float64},
  a_b::Vector{Float64},
  a_t::Vector{Float64},
  a_p::Vector{Float64}
)
  N = ni * nj * nk

  # オフセット計算（Fortran順序）
  off_i = 1           # i方向（南北）
  off_j = ni          # j方向（西東）
  off_k = ni * nj     # k方向（下上）

  # COO形式で係数を収集
  I_list = Int[]
  J_list = Int[]
  V_list = Float64[]

  sizehint!(I_list, 7 * N)  # メモリ事前確保
  sizehint!(J_list, 7 * N)
  sizehint!(V_list, 7 * N)

  for p in 1:N
    # 対角成分
    push!(I_list, p)
    push!(J_list, p)
    push!(V_list, a_p[p])

    # i-1 (南)
    if p > off_i
      push!(I_list, p)
      push!(J_list, p - off_i)
      push!(V_list, -a_s[p])
    end

    # i+1 (北)
    if p <= N - off_i
      push!(I_list, p)
      push!(J_list, p + off_i)
      push!(V_list, -a_n[p])
    end

    # j-1 (西)
    if p > off_j
      push!(I_list, p)
      push!(J_list, p - off_j)
      push!(V_list, -a_w[p])
    end

    # j+1 (東)
    if p <= N - off_j
      push!(I_list, p)
      push!(J_list, p + off_j)
      push!(V_list, -a_e[p])
    end

    # k-1 (下)
    if p > off_k
      push!(I_list, p)
      push!(J_list, p - off_k)
      push!(V_list, -a_b[p])
    end

    # k+1 (上)
    if p <= N - off_k
      push!(I_list, p)
      push!(J_list, p + off_k)
      push!(V_list, -a_t[p])
    end
  end

  # COO→CSC変換
  A = sparse(I_list, J_list, V_list, N, N)

  return A
end


"""
    solve_adjoint!(T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
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
  T_cal: DHCP計算温度場 (nt, ni, nj, nk) [K]
  Y_obs: 観測温度（底面） (nt, ni, nj) [K]
  nt: 時間ステップ数
  rho: 密度 [kg/m³]
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
  λ_all: 随伴場時系列 (nt, ni, nj, nk)
  cg_iters: CG反復回数履歴 (nt-1,)
"""
function solve_adjoint!(
  T_cal::Array{Float64,4},
  Y_obs::Array{Float64,3},
  nt::Int,
  rho::Float64,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  dx::Float64,
  dy::Float64,
  dz::Vector{Float64},
  dz_b::Vector{Float64},
  dz_t::Vector{Float64},
  dt::Float64;
  rtol::Float64=1e-8,
  maxiter::Int=1000,
  verbose::Bool=false
)
  ni, nj, nk = size(T_cal[1, :, :, :])
  N = ni * nj * nk

  # 随伴場の初期化
  λ_all = zeros(Float64, nt, ni, nj, nk)
  λ_all[nt, :, :, :] .= 0.0  # 終端条件（Pythonオリジナル1322行）

  cg_iters = zeros(Int, nt-1)

  # ホットスタート用初期推定値
  x0 = vec(λ_all[nt, :, :, :])

  if verbose
    println("Adjoint求解開始（後退時間積分）")
    println("  格子: $ni×$nj×$nk, N=$N")
    println("  時間ステップ: $nt")
    println("  CG: rtol=$rtol, maxiter=$maxiter")
  end

  # 後退時間ループ（Pythonオリジナル1328行: range(nt-2, -1, -1)）
  for t in (nt-1):-1:1
    # 次ステップ（時間的に後）の随伴場を初期値とする
    λ_initial = λ_all[t+1, :, :, :]

    # 温度場から熱物性値計算（Pythonオリジナル1332行: T_cal[t]）
    cp, k = thermal_properties_calculator(T_cal[t, :, :, :], cp_coeffs, k_coeffs)

    # 係数とRHS構築（Pythonオリジナル1334-1335行）
    # 残差注入: T_cal[t,:,:,1]（底面）と Y_obs[t] を使用
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_adjoint_system!(
      λ_initial, T_cal[t, :, :, 1], Y_obs[t, :, :],
      rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 疎行列組み立て
    A = assemble_adjoint_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    # 対角前処理器（Pythonオリジナル1339-1345行）
    diag_A = diag(A)
    inv_diag = [d != 0.0 ? 1.0/d : 0.0 for d in diag_A]
    M = Diagonal(inv_diag)

    # CG法求解（Pythonオリジナル1347行）
    x, cg_history = cg!(x0, A, b, Pl=M, reltol=rtol, maxiter=maxiter, log=true)

    # 収束チェック（Pythonオリジナル1350-1353行）
    if !cg_history.isconverged
      @warn "[t=$t] CG法が収束しませんでした (iter=$(length(cg_history.data[:resnorm])))"
    end

    # 結果保存（Pythonオリジナル1356-1357行）
    λ_all[t, :, :, :] = reshape(x, ni, nj, nk)
    x0 .= x  # ホットスタート更新

    # 反復回数記録
    cg_iters[t] = length(cg_history.data[:resnorm])

    if verbose && (t % 10 == 0 || t == nt-1)
      println("  t=$t: CG収束 (iter=$(cg_iters[t]))")
    end
  end

  if verbose
    println("Adjoint求解完了")
    println("  平均CG反復回数: $(mean(cg_iters))")
    println("  随伴場範囲: [$(minimum(λ_all)), $(maximum(λ_all))]")
  end

  return λ_all, cg_iters
end


end # module AdjointSolver