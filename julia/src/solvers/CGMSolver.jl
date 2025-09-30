"""
CGMSolver.jl

Phase 4: 共役勾配法（CGM）による逆問題ソルバー

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

using LinearAlgebra
using SparseArrays
using Printf

# Phase 2, 3モジュールの読み込み
include("DHCPSolver.jl")
include("AdjointSolver.jl")
include("StoppingCriteria.jl")

using .DHCPSolver
using .AdjointSolver
using .StoppingCriteria

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
"""
function tensor_dot(a::Array{T}, b::Array{T}) where T <: Real
  return Float64(sum(a .* b))
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
  T_cal: DHCP計算温度場 (nt, ni, nj, nk)
  Y_obs: 観測温度（底面） (nt, ni, nj) [K]
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数
  k_coeffs: 熱伝導率多項式係数
  dx, dy: x, y方向格子幅 [m]
  dz: z方向格子幅配列 (nk,) [m]
  dz_b: 下側界面距離 (nk,) [m]
  dz_t: 上側界面距離 (nk,) [m]
  dt: 時間刻み [s]
  rtol: CG相対許容誤差
  maxiter: CG最大反復数

Returns:
  gradient: 勾配場 (nt-1, ni, nj)
"""
function compute_gradient!(
  T_cal::Array{Float64,4},
  Y_obs::Array{Float64,3},
  rho::Float64,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  dx::Float64, dy::Float64,
  dz::Vector{Float64}, dz_b::Vector{Float64}, dz_t::Vector{Float64},
  dt::Float64,
  rtol::Float64=1e-8,
  maxiter::Int=20000
)
  nt, ni, nj, nk = size(T_cal)

  # 随伴場求解（Phase 3）
  lambda_field, cg_iters = solve_adjoint!(
    T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt;
    rtol=rtol, maxiter=maxiter, verbose=false
  )

  # 勾配抽出（表面 k=nk での随伴場）
  gradient = zeros(nt - 1, ni, nj)
  for n in 1:(nt - 1)
    gradient[n, :, :] = lambda_field[n, :, :, nk]  # 表面（上端）
  end

  return gradient
end


"""
    compute_sensitivity!(T_init, p_n, rho, cp_coeffs, k_coeffs, dx, dy, dz, dz_b, dz_t, dt,
                          rtol, maxiter) -> dT

感度問題求解（熱流束微小変化に対する温度応答）

感度問題は、探索方向p_nを表面熱流束として与えたDHCP問題：
  dT = DHCP(T_init=0, q=p_n)

物理的意味:
  熱流束がp_n方向に微小変化した場合の温度場の変化

Args:
  T_init: 初期温度場（ゼロ） (ni, nj, nk)
  p_n: 探索方向（熱流束の方向） (nt-1, ni, nj)
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数
  k_coeffs: 熱伝導率多項式係数
  dx, dy: x, y方向格子幅 [m]
  dz: z方向格子幅配列 (nk,) [m]
  dz_b: 下側界面距離 (nk,) [m]
  dz_t: 上側界面距離 (nk,) [m]
  dt: 時間刻み [s]
  rtol: CG相対許容誤差
  maxiter: CG最大反復数

Returns:
  dT: 感度場（温度応答） (nt, ni, nj, nk)
"""
function compute_sensitivity!(
  T_init::Array{Float64,3},
  p_n::Array{Float64,3},
  nt::Int,
  rho::Float64,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64},
  dx::Float64, dy::Float64,
  dz::Vector{Float64}, dz_b::Vector{Float64}, dz_t::Vector{Float64},
  dt::Float64,
  rtol::Float64=1e-8,
  maxiter::Int=20000
)
  # DHCPソルバーを流用（初期温度ゼロ、熱流束=p_n）
  dT = solve_dhcp!(
    T_init, p_n, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt;
    rtol=rtol, maxiter=maxiter, verbose=false
  )

  return dT
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
  res_T: 残差場 (nt-1, ni, nj) [K]
  Sp: 感度場の底面値 (nt-1, ni, nj) [K]
  eps: ゼロ除算防止の微小値（デフォルト1e-12）

Returns:
  beta: ステップサイズ
"""
function compute_step_size(res_T::Array{Float64,3}, Sp::Array{Float64,3}, eps::Float64=1e-12)
  numerator = tensor_dot(res_T, Sp)
  denominator = tensor_dot(Sp, Sp)
  beta = numerator / (denominator + eps)
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
  Y_obs: 観測温度（底面） (nt, ni, nj) [K]
  q_init: 初期熱流束推定 (nt-1, ni, nj) [W/m²]
  dx, dy: x, y方向格子幅 [m]
  dz: z方向格子幅配列 (nk,) [m]
  dz_b: 下側界面距離 (nk,) [m]
  dz_t: 上側界面距離 (nk,) [m]
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
  q_final: 最終逆解析熱流束 (nt-1, ni, nj) [W/m²]
  T_cal_final: 最終温度場 (nt, ni, nj, nk) [K]
  J_hist: 目的関数履歴 (vector)
"""
function solve_cgm!(
  T_init::Array{Float64,3},
  Y_obs::Array{Float64,3},
  q_init::Array{Float64,3},
  dx::Float64, dy::Float64,
  dz::Vector{Float64}, dz_b::Vector{Float64}, dz_t::Vector{Float64},
  dt::Float64,
  rho::Float64,
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64};
  params::NamedTuple=(;)
)
  # パラメータ展開（デフォルト値付き）
  max_iter = get(params, :max_iter, 20000)
  rtol_dhcp = get(params, :rtol_dhcp, 1e-6)
  rtol_adjoint = get(params, :rtol_adjoint, 1e-8)
  maxiter_cg = get(params, :maxiter_cg, 20000)
  sigma = get(params, :sigma, 1.8)
  min_iter = get(params, :min_iter, 10)
  P = get(params, :P, 10)
  eta = get(params, :eta, 1e-4)
  dire_reset_every = get(params, :dire_reset_every, 5)
  eps = get(params, :eps, 1e-12)
  beta_max = get(params, :beta_max, 1e8)
  verbose = get(params, :verbose, true)

  # 問題サイズ
  nt, ni, nj = size(Y_obs)
  nk = size(T_init, 3)

  # 初期化
  q = copy(q_init)
  J_hist = Float64[]

  M = ni * nj
  epsilon = M * (sigma^2) * (nt - 1)  # Discrepancy基準値

  grad = zeros(nt - 1, ni, nj)
  grad_last = zeros(nt - 1, ni, nj)
  p_n_last = zeros(nt - 1, ni, nj)

  bottom_idx = 1   # Julia 1-indexed（底面）
  top_idx = nk     # Julia 1-indexed（表面）

  # 停止判定パラメータ
  stop_params = (
    epsilon=epsilon, sigma=sigma, min_iter=min_iter,
    P=P, eta=eta, max_iter=max_iter, eps=eps
  )

  T_cal = nothing  # 最終温度場保存用

  # CGMループ
  for it in 0:(max_iter - 1)
    if verbose
      println("\n=== CGM反復 $(it) ===")
    end

    # Step 1: 直接問題求解（DHCP）
    T_cal = solve_dhcp!(
      T_init, q, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt;
      rtol=rtol_dhcp, maxiter=maxiter_cg, verbose=false
    )

    # Step 2: 目的関数と停止判定
    res_T = T_cal[2:nt, :, :, bottom_idx] .- Y_obs[2:nt, :, :]  # (nt-1, ni, nj)
    J = tensor_dot(res_T, res_T)
    push!(J_hist, J)

    if verbose
      @printf("J = %.5e\n", J)
    end

    # 停止判定
    status = check_stopping_criteria(J, J_hist, res_T, it, stop_params)
    if status.should_stop
      if verbose
        println(status.reason)
      end
      break
    end

    # Step 3: 随伴問題求解（勾配計算）
    grad = compute_gradient!(
      T_cal, Y_obs, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol_adjoint, maxiter_cg
    )

    # Step 4: 共役勾配方向計算（ポラック・リビエール）
    if it == 0 || tensor_dot(grad, p_n_last) <= 0 || it % dire_reset_every == 0
      # 最急降下方向にリセット
      p_n = copy(grad)
      gamma = 0.0
    else
      # ポラック・リビエール係数
      y = grad .- grad_last
      denom = tensor_dot(grad_last, grad_last) + eps
      gamma = max(0.0, tensor_dot(grad, y) / denom)

      # 探索方向候補
      p_n_candidate = grad .+ gamma .* p_n_last

      # 下降方向チェック
      if tensor_dot(grad, p_n_candidate) > 0
        p_n = p_n_candidate
      else
        # 下降方向でない場合はリセット
        p_n = copy(grad)
        gamma = 0.0
      end
    end

    p_n_last = copy(p_n)

    # Step 5: 感度問題求解（dT計算）
    dT_init = zeros(ni, nj, nk)
    dT = compute_sensitivity!(
      dT_init, p_n, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol_adjoint, maxiter_cg
    )

    # Step 6: ステップサイズ計算
    Sp = dT[2:nt, :, :, bottom_idx]  # (nt-1, ni, nj)
    beta = compute_step_size(res_T, Sp, eps)

    # ステップサイズ制限（初回のみ）
    if it == 0 && abs(beta) > beta_max
      if verbose
        @printf("[警告] beta制限: %.2e => %.2e\n", beta, sign(beta) * beta_max)
      end
      beta = clamp(beta, -beta_max, beta_max)
    end

    if verbose
      @printf("beta = %.4e, gamma = %.4e\n", beta, gamma)
    end

    # Step 7: 熱流束更新
    q .= q .- beta .* p_n

    # Step 8: 勾配更新
    grad_last = copy(grad)
  end

  return q, T_cal, J_hist
end

end # module CGMSolver
