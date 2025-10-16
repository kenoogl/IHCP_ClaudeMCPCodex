"""
CGMSolver.jl

共役勾配法（CGM）による逆問題ソルバー

共役勾配法（ポラック・リビエール型）による表面熱流束の逆解析：
1. 勾配計算（随伴場λの表面値）
2. ポラック・リビエール探索方向
3. 感度問題（熱流束微小変化に対する温度応答）
4. ステップサイズ計算（ライン検索）
5. 熱流束更新
6. 停止判定（Discrepancy & プラトー）

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- global_CGM_time() (1371-1553行)

主要関数:
- compute_gradient!: 勾配計算（Adjointソルバーを呼び出し）
- compute_sensitivity!: 感度問題（DHCPソルバーを流用）
- compute_step_size: ステップサイズ計算（ライン検索）
- solve_cgm!: CGMメインループ
"""

module CGMSolver

using LinearAlgebra: dot
using Printf

# 共通モジュール
import ..Commons
using ..Commons: WorkBuffers, get_backend, reset_work_buffers!

# 親モジュールで既にinclude済み
using ..DHCPSolver
using ..AdjointSolver
using ..SensitivitySolver
using ..StoppingCriteria
using ..AdaptiveTolerance: AdaptiveToleranceParams

import ..CommonSolver
using ..CommonSolver: PBiCGSTAB!

export solve_cgm!, compute_gradient!, compute_sensitivity!, compute_step_size


"""
    tensor_dot(a, b) -> Float64

全要素のテンソル内積計算（Pythonオリジナルの_dot関数）

Pythonオリジナル:
  def _dot(a, b):
    return float(np.tensordot(a, b, axes=a.ndim))

Args:
  a, b: 同じサイズの配列

Returns:
  内積値（スカラー）
  CHECK: 型安定性のため、Array{Float64,3}

メモリ効率化:
  一時配列を生成しないdot(vec(), vec())を使用
"""
function tensor_dot(a::AbstractArray{T1}, b::AbstractArray{T2}) where {T1 <: Real, T2 <: Real}
  @assert size(a) == size(b)
  acc = zero(promote_type(T1, T2))
  @inbounds @simd for idx in eachindex(a, b)
    acc += a[idx] * b[idx]
  end
  return acc
end


"""
    compute_gradient!(T_cal, Y_obs, rho, cp_coeffs, k_coeffs, dx, dy, dz, dz_b, dz_t, dt,
                      rtol, maxiter) -> gradient

勾配計算（随伴場の表面値）

手順:
  1. Adjoint随伴場を求解（Phase 3）
  2. 表面（k=nk）での随伴場を抽出 → 勾配

勾配の物理的意味:
  gradient[n, i, j] = ∂J/∂q[n, i, j]
  ここでJは目的関数（残差二乗和）

Args:
  T_cal: DHCP計算温度場 (ni, nj, nk, nt) ※Phase 2.2: 時間次元を最後に配置
  Y_obs: 観測温度（底面） (ni, nj, nt) [K] ※Phase 2.2: 時間次元を最後に配置
  work: 
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数
  k_coeffs: 熱伝導率多項式係数
  dx, dy: x, y方向格子幅 [m]
  ZC, dz: Z方向格子情報
  dt: 時間刻み [s]
  rtol: CG相対許容誤差
  maxiter: CG最大反復数

Returns:
  gradient: 勾配場 (ni, nj, nt-1) ※Phase 2.2: 時間次元を最後に配置
  cg_iters: 各時間ステップの反復回数
"""
function compute_gradient!(
  T_cal::AbstractArray{T,4},
  Y_obs::AbstractArray{T,3},
  work::WorkBuffers,
  rho::T,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  dx::T, dy::T,
  ZC::Vector{Float64}, dz::Vector{Float64},
  dt::T;
  rtol::T=T(1e-8),
  maxiter::Int=20000,
  verbose::Bool=false,
  par::String="sequential",
  gradient_buffer::Union{Nothing,Array{T,3}}=nothing,
  adjoint_buffer::Union{Nothing,Array{T,4}}=nothing,
  iter_buffer::Union{Nothing,Vector{Int}}=nothing,
  initial_strategy::Symbol=:residual,
  residual_scale::T=T(1),
  smoother::Symbol=:gs,
  adaptive_tol::Bool=false,
  adaptive_tol_params::Union{Nothing,AdaptiveToleranceParams{T}}=nothing
) where {T <: AbstractFloat}
  ni, nj, nk, nt = size(T_cal)

  # 随伴場求解（Phase 3）
  lambda_field, cg_iters = solve_adjoint_mf!(
    T_cal, Y_obs, work, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, ZC, dz, dt;
    rtol=rtol, maxiter=maxiter, verbose=verbose, par=par,
    lambda_buffer=adjoint_buffer,
    iter_buffer=iter_buffer,
    initial_strategy=initial_strategy,
    residual_scale=residual_scale,
    smoother=smoother,
    adaptive_tol=adaptive_tol,
    adaptive_tol_params=adaptive_tol_params
  )

  # 勾配抽出（表面 k=nk での随伴場）
  gradient = isnothing(gradient_buffer) ? zeros(T, ni, nj, nt - 1) : gradient_buffer
  expected_shape = (ni, nj, nt - 1)
  if size(gradient) != expected_shape
    throw(ArgumentError("gradient_buffer size mismatch: expected $(expected_shape), got $(size(gradient))"))
  end
  for n in 1:(nt - 1)
    @views gradient[:, :, n] .= lambda_field[:, :, nk, n]  # 表面（上端）
  end

  return gradient, cg_iters
end


"""
    compute_step_size(res_T, Sp, eps) -> beta

ステップサイズ計算（ライン検索）

ステップサイズの導出:
  J(q - β*p) を最小化
  → β = <res_T, Sp> / <Sp, Sp>

ここで:
  res_T: 残差場（T_cal - Y_obs）
  Sp: 感度場の底面値（dT[:, :, :, 1]）

Args:
  res_T: 残差場 (ni, nj, nt-1) [K] ※Phase 2.2: 時間次元を最後に配置
  Sp: 感度場の底面値 (ni, nj, nt-1) [K] ※Phase 2.2: 時間次元を最後に配置
  eps: ゼロ除算防止の微小値（デフォルト1e-12）

Returns:
  beta: ステップサイズ
"""
function compute_step_size(res_T::AbstractArray{T1,3}, Sp::AbstractArray{T2,3}, eps::Real=1e-12) where {T1 <: Real, T2 <: Real}
  numerator = tensor_dot(res_T, Sp)
  denominator = tensor_dot(Sp, Sp)
  eps_T = convert(typeof(numerator), eps)
  beta = numerator / (denominator + eps_T)
  return beta
end


"""
    solve_cgm!(T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
               rho, cp_coeffs, k_coeffs; params)
      -> (q_final, T_cal_final, J_hist)

共役勾配法（CGM）による熱流束逆解析

CGMアルゴリズム（ポラック・リビエール型）:
  1. 初期化: q = q_init, gradient = ∇J(q), p = -gradient
  2. ループ (iter = 0, 1, 2, ...):
      a. DHCP求解: T_cal = DHCP(T_init, q)
      b. 停止判定: Discrepancy / プラトー / 最大反復数
      c. Adjoint求解: λ = Adjoint(T_cal, Y_obs)
      d. 勾配計算: gradient_new = λ[:, :, end, :]
      e. ポラック・リビエール係数:
         β = max(0, <gradient_new, gradient_new - gradient> / <gradient, gradient>)
      f. 探索方向更新: p = -gradient_new + β * p
      g. 感度問題: dT = DHCP(0, p)
      h. ステップサイズ: α = <res_T, Sp> / <Sp, Sp>
      i. 熱流束更新: q = q - α * p
      j. 勾配更新: gradient = gradient_new

対応Pythonコード: global_CGM_time() (1371-1553行)

Args:
  T_init: 初期温度場 (ni, nj, nk) [K]
  Y_obs: 観測温度（底面） (ni, nj, nt) [K] ※Phase 2.2: 時間次元を最後に配置
  q_init: 初期熱流束推定 (ni, nj, nt-1) [W/m²] ※Phase 2.2: 時間次元を最後に配置
  work: Heat3d用配列群
  dx, dy: x, y方向格子幅 [m]
  ZC: z方向格子配列 (nk,) [m]
  dz: CV幅 (nk,) [m]
  dt: 時間刻み [s]
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数
  k_coeffs: 熱伝導率多項式係数
  params: CGMパラメータ（NamedTuple）
    - max_iter: 最大反復数（デフォルト20000）
    - rtol_dhcp: DHCP CG許容誤差（デフォルト1e-6）
    - rtol_adjoint: Adjoint CG許容誤差（デフォルト1e-8）
    - maxiter_cg: CG最大反復数（デフォルト20000）
    - sigma: 測定誤差標準偏差 [K]（デフォルト1.8）
    - min_iter: 最小反復数（デフォルト10）
    - P: プラトー検出ウィンドウ（デフォルト10）
    - eta: プラトー相対下降閾値（デフォルト1e-4）
    - dire_reset_every: 方向リセット周期（デフォルト5）
    - eps: ゼロ除算防止値（デフォルト1e-12）
    - beta_max: ステップサイズ上限（デフォルト1e8）
    - verbose: 詳細出力フラグ（デフォルトtrue）

Returns:
  q_final: 最終逆解析熱流束 (ni, nj, nt-1) [W/m²] ※Phase 2.2: 時間次元を最後に配置
  T_cal_final: 最終温度場 (ni, nj, nk, nt) [K] ※Phase 2.2: 時間次元を最後に配置
  J_hist: 目的関数履歴 (vector)
"""
function solve_cgm!(
  T_init::AbstractArray{T,3},
  Y_obs::AbstractArray{T,3},
  q_init::AbstractArray{T,3},
  work::WorkBuffers,
  dx::T, dy::T,
  ZC::Vector{Float64}, dz::Vector{Float64},
  dt::T,
  rho::T,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64};
  params::NamedTuple=(;),
  par::String="sequential"
) where {T <: AbstractFloat}
  # パラメータ展開（デフォルト値付き）
  max_iter = get(params, :max_iter, 20000)
  rtol_dhcp = T(get(params, :rtol_dhcp, 1e-6))
  rtol_adjoint = T(get(params, :rtol_adjoint, 1e-8))
  maxiter_cg = get(params, :maxiter_cg, 20000)
  sigma = T(get(params, :sigma, 1.8))
  min_iter = get(params, :min_iter, 10)
  P = get(params, :P, 10)
  eta = T(get(params, :eta, 1e-4))
  dire_reset_every = get(params, :dire_reset_every, 5)
  eps = T(get(params, :eps, 1e-12))
  beta_max = T(get(params, :beta_max, 1e8))
  verbose = get(params, :verbose, true)
  dhcp_use_previous_solution = get(params, :dhcp_use_previous_solution, true)
  dhcp_extrapolation = get(params, :dhcp_extrapolation, :none)
  dhcp_smoother = get(params, :dhcp_smoother, :gs)
  sensitivity_use_previous_solution = get(params, :sensitivity_use_previous_solution, true)
  sensitivity_extrapolation = get(params, :sensitivity_extrapolation, :none)
  sensitivity_smoother = get(params, :sensitivity_smoother, :gs)
  adjoint_initial_strategy = get(params, :adjoint_initial_strategy, :residual)
  adjoint_residual_scale = T(get(params, :adjoint_residual_scale, 1.0))
  adjoint_smoother = get(params, :adjoint_smoother, :gs)

  # Phase 1-E: 適応的収束判定パラメータ
  adaptive_tol = get(params, :adaptive_tol, false)
  tol_min_dhcp = T(get(params, :tol_min_dhcp, 1e-7))
  tol_max_dhcp = T(get(params, :tol_max_dhcp, 1e-5))
  tol_min_adjoint = T(get(params, :tol_min_adjoint, 1e-9))
  tol_max_adjoint = T(get(params, :tol_max_adjoint, 1e-5))
  tol_min_sensitivity = T(get(params, :tol_min_sensitivity, 1e-7))
  tol_max_sensitivity = T(get(params, :tol_max_sensitivity, 1e-5))
  kappa_theta = T(get(params, :kappa_theta, 0.2))
  warmup_steps = get(params, :warmup_steps, 3)

  # 問題サイズ（メモリレイアウト最適化: Phase 2.2）
  ni, nj, nt = size(Y_obs)
  nk = size(T_init, 3)

  # 初期化
  q = copy(q_init)
  J_hist = T[]

  M = ni * nj
  epsilon = T(M) * (sigma^2) * T(nt - 1)  # Discrepancy基準値

  grad = zeros(T, ni, nj, nt - 1)
  grad_last = zeros(T, ni, nj, nt - 1)
  p_n = zeros(T, ni, nj, nt - 1)
  p_n_last = zeros(T, ni, nj, nt - 1)
  tmp_dir = zeros(T, ni, nj, nt - 1)
  res_T = zeros(T, ni, nj, nt - 1)

  dhcp_T_buffer = zeros(T, ni, nj, nk, nt)
  dhcp_iter_buffer = zeros(Int, nt)

  adjoint_buffer = zeros(T, ni, nj, nk, nt)
  adjoint_iter_buffer = zeros(Int, nt - 1)

  sensitivity_buffer = zeros(T, ni, nj, nk, nt)
  sensitivity_iter_buffer = zeros(Int, nt)
  dT_init = zeros(T, ni, nj, nk)


  bottom_idx = 1   # Julia 1-indexed（裏面 S2）

  # 停止判定パラメータ
  stop_params = (
    epsilon=epsilon, sigma=sigma, min_iter=min_iter,
    P=P, eta=eta, max_iter=max_iter, eps=eps
  )

  # Phase 1-E: 適応的収束判定パラメータオブジェクト作成
  adaptive_tol_params_dhcp = adaptive_tol ? AdaptiveToleranceParams{T}(
    tol_min=tol_min_dhcp, tol_max=tol_max_dhcp,
    κ_θ=kappa_theta, warmup_steps=warmup_steps, ε=eps
  ) : nothing

  adaptive_tol_params_adjoint = adaptive_tol ? AdaptiveToleranceParams{T}(
    tol_min=tol_min_adjoint, tol_max=tol_max_adjoint,
    κ_θ=kappa_theta, warmup_steps=warmup_steps, ε=eps
  ) : nothing

  adaptive_tol_params_sensitivity = adaptive_tol ? AdaptiveToleranceParams{T}(
    tol_min=tol_min_sensitivity, tol_max=tol_max_sensitivity,
    κ_θ=kappa_theta, warmup_steps=warmup_steps, ε=eps
  ) : nothing

  T_cal = nothing  # 最終温度場保存用

  # CGMループ
  for it in 0:(max_iter - 1)
    if verbose
      println("\n=== CGM反復 $(it) ===")
    end

    # Step 1: 順問題求解（DHCP）
    reset_work_buffers!(work)  # WorkBuffersをクリーンな状態にリセット
    dhcp_start = time()
    T_cal, iter_counts = solve_dhcp!(
      T_init, q, work,
      nt, rho, cp_coeffs, k_coeffs,
      dx, dy, ZC, dz, dt;
      rtol=rtol_dhcp, maxiter=maxiter_cg, verbose=verbose, par=par,
      T_buffer=dhcp_T_buffer,
      iter_buffer=dhcp_iter_buffer,
      use_previous_solution=dhcp_use_previous_solution,
      extrapolation=dhcp_extrapolation,
      smoother=dhcp_smoother,
      adaptive_tol=adaptive_tol,
      adaptive_tol_params=adaptive_tol_params_dhcp
    )
    dhcp_time = time() - dhcp_start

    total_iters = sum(iter_counts[2:end])
    avg_iters = total_iters / (nt - 1)
    @printf("  DHCP solve time: %.3f s (nt=%d steps, total_iters=%d, avg=%.1f)\n",
              dhcp_time, nt, total_iters, avg_iters)

    # Step 2: 目的関数計算
    @views @. res_T = T_cal[:, :, bottom_idx, 2:nt] - Y_obs[:, :, 2:nt]
    J = tensor_dot(res_T, res_T)
    push!(J_hist, J)

    if verbose
      @printf("J = %.5e\n", Float64(J))
    end

    # Step 3: 随伴問題求解（勾配計算）
    reset_work_buffers!(work)  # WorkBuffersをクリーンな状態にリセット
    gradient_start = time()
    grad, grad_iters = compute_gradient!(
      T_cal, Y_obs, work, rho, cp_coeffs, k_coeffs,
      dx, dy, ZC, dz, dt,
      rtol=rtol_adjoint, maxiter=maxiter_cg, verbose=verbose, par=par,
      gradient_buffer=grad,
      adjoint_buffer=adjoint_buffer,
      iter_buffer=adjoint_iter_buffer,
      initial_strategy=adjoint_initial_strategy,
      residual_scale=adjoint_residual_scale,
      smoother=adjoint_smoother,
      adaptive_tol=adaptive_tol,
      adaptive_tol_params=adaptive_tol_params_adjoint
    )
    gradient_time = time() - gradient_start

    total_iters = sum(grad_iters[2:end])
    avg_iters = total_iters / (nt - 1)
    @printf("  Gradient solve time: %.3f s (nt=%d steps, total_iters=%d, avg=%.1f)\n",
              gradient_time, nt, total_iters, avg_iters)

    # Step 4: 共役勾配方向計算（ポラック・リビエール）
    if it == 0 || tensor_dot(grad, p_n_last) <= 0 || it % dire_reset_every == 0
      # 最急降下方向にリセット
      @. p_n = grad
      gamma = zero(T)
    else
      # ポラック・リビエール係数
      @. tmp_dir = grad - grad_last
      denom = tensor_dot(grad_last, grad_last) + eps
      gamma = max(zero(T), tensor_dot(grad, tmp_dir) / denom)

      # 探索方向候補
      @. tmp_dir = grad + gamma * p_n_last

      # 下降方向チェック
      if tensor_dot(grad, tmp_dir) > zero(T)
        @. p_n = tmp_dir
      else
        # 下降方向でない場合はリセット
        @. p_n = grad
        gamma = zero(T)
      end
    end

    # Step 5: 感度問題求解（dT計算）
    reset_work_buffers!(work)  # WorkBuffersをクリーンな状態にリセット
    fill!(dT_init, zero(T))
    sensitivity_start = time()
    dT, sens_iters = solve_sensitivity!(
      dT_init, p_n, work,
      nt, rho, cp_coeffs, k_coeffs,
      dx, dy, ZC, dz, dt,
      rtol=rtol_adjoint,
      maxiter=maxiter_cg,
      verbose=verbose,
      par=par,
      T_buffer=sensitivity_buffer,
      iter_buffer=sensitivity_iter_buffer,
      use_previous_solution=sensitivity_use_previous_solution,
      extrapolation=sensitivity_extrapolation,
      smoother=sensitivity_smoother,
      adaptive_tol=adaptive_tol,
      adaptive_tol_params=adaptive_tol_params_sensitivity
    )
    sensitivity_time = time() - sensitivity_start

    total_iters = sum(sens_iters[2:end])
    avg_iters = total_iters / (nt - 1)
    @printf("  Sensitivity solve time: %.3f s (nt=%d steps, total_iters=%d, avg=%.1f)\n",
              sensitivity_time, nt, total_iters, avg_iters)

    # Step 6: ステップサイズ計算
    Sp = @view dT[:, :, bottom_idx, 2:nt]  # (ni, nj, nt-1)
    beta = compute_step_size(res_T, Sp, eps)

    # ステップサイズ制限（初回のみ）
    if it == 0 && abs(beta) > beta_max
      if verbose
        @printf("[警告] beta制限: %.2e => %.2e\n", Float64(beta), Float64(sign(beta) * beta_max))
      end
      beta = clamp(beta, -beta_max, beta_max)
    end

    if verbose
      @printf("beta = %.4e, gamma = %.4e\n", Float64(beta), Float64(gamma))
    end

    # Step 7: 熱流束更新
    @. q = q - beta * p_n

    # Step 8: 勾配更新
    @. grad_last = grad
    @. p_n_last = p_n

    # Step 9: 停止判定（更新後にチェック）
    status = check_stopping_criteria(J, J_hist, res_T, it, stop_params)
    if status.should_stop
      if verbose
        println(status.reason)
      end
      break
    end
  end

  return q, T_cal, J_hist
end

end # module CGMSolver
