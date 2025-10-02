"""
DHCPSolver.jl

Phase 2: 直接熱伝導問題（DHCP）ソルバー

陰解法（後退差分）による熱伝導方程式の数値解法：
    (I - α∇²)T^(n+1) = T^n + q_boundary
    α = (k·dt) / (ρ·cp·dx²)

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- coeffs_and_rhs_building_DHCP() (1087-1129行)
- assemble_A_DHCP() (1131-1153行)
- multiple_time_step_solver_DHCP() (1156-1219行)

主要関数:
- build_dhcp_system!: 係数とRHS構築（7点ステンシル）
- assemble_dhcp_matrix: CSC疎行列組み立て
- solve_dhcp!: 複数時間ステップCG求解（ホットスタート対応）
"""

module DHCPSolver

using LinearAlgebra
using SparseArrays
using IterativeSolvers

# Phase 1モジュールの読み込み
include("../ThermalProperties.jl")
using .ThermalProperties

export build_dhcp_system!, assemble_dhcp_matrix, solve_dhcp!, dhcp_index, solve_dhcp_multiple_timesteps


"""
    dhcp_index(i, j, k, ni, nj) -> Int

Fortran順序（列優先）でのグローバルインデックス計算

Pythonオリジナルコード:
    p = i + j*ni + k*ni*nj  (0-indexed)

Julia変換（1-indexed）:
    p = (i-1) + (j-1)*ni + (k-1)*ni*nj + 1
      = i + (j-1)*ni + (k-1)*ni*nj

Args:
  i: x方向インデックス (1:ni)
  j: y方向インデックス (1:nj)
  k: z方向インデックス (1:nk)
  ni, nj: 格子点数

Returns:
  p: グローバルインデックス (1:N)
"""
@inline function dhcp_index(i::Int, j::Int, k::Int, ni::Int, nj::Int)
  return i + (j - 1) * ni + (k - 1) * ni * nj
end


"""
    build_dhcp_system!(T_initial, q_surface, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt)
      -> (a_w, a_e, a_s, a_n, a_b, a_t, a_p, b)

直接熱伝導問題（DHCP）の係数とRHSベクトル構築

有限体積法による3D熱伝導方程式の離散化:
    ρ·cp·(∂T/∂t) = ∇·(k∇T) + q

陰解法（後退差分）:
    (ρ·cp·V/dt)·T^(n+1) - Σ(a_f·T_f^(n+1)) = (ρ·cp·V/dt)·T^n + Q

7点ステンシル係数:
  - a_w, a_e: 西・東（y方向）
  - a_s, a_n: 南・北（x方向）
  - a_b, a_t: 下・上（z方向）
  - a_p: 対角項（中心）

境界条件:
  - x, y, z境界: 断熱（係数=0）
  - z上端: ノイマン条件（熱流束 q_surface）

対応Pythonコード: coeffs_and_rhs_building_DHCP() (1087-1129行)

Args:
  T_initial: 前ステップ温度場 (ni, nj, nk) [K]
  q_surface: 表面熱流束 (ni, nj) [W/m²]
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
function build_dhcp_system!(
  T_initial::AbstractArray{Float64,3},
  q_surface::AbstractArray{Float64,2},
  rho::Float64,
  cp::AbstractArray{Float64,3},
  k::AbstractArray{Float64,3},
  dx::Float64,
  dy::Float64,
  dz::Vector{Float64},
  dz_b::Vector{Float64},
  dz_t::Vector{Float64},
  dt::Float64
)
  ni, nj, nk = size(T_initial)
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

  # 全格子点ループ（並列化可能だが、まずは順次実行）
  for k_idx in 1:nk, j in 1:nj, i in 1:ni
    p = dhcp_index(i, j, k_idx, ni, nj)

    # 格子幅取得
    dz_k = dz[k_idx]
    dz_t_k = dz_t[k_idx]
    dz_b_k = dz_b[k_idx]

    # 中心セル熱伝導率
    k_p = k[i, j, k_idx]

    # 時間項（蓄熱項）
    a_p_0 = rho * cp[i, j, k_idx] * dx * dy * dz_k / dt

    # 6方向の熱伝導係数（調和平均）
    # 西（j-1）
    if j == 1
      a_w[p] = 0.0
    else
      k_w = k[i, j-1, k_idx]
      a_w[p] = (2.0 * k_p * k_w / (k_p + k_w)) * dy * dz_k / dx
    end

    # 東（j+1）
    if j == nj
      a_e[p] = 0.0
    else
      k_e = k[i, j+1, k_idx]
      a_e[p] = (2.0 * k_p * k_e / (k_p + k_e)) * dy * dz_k / dx
    end

    # 南（i-1）
    if i == 1
      a_s[p] = 0.0
    else
      k_s = k[i-1, j, k_idx]
      a_s[p] = (2.0 * k_p * k_s / (k_p + k_s)) * dx * dz_k / dy
    end

    # 北（i+1）
    if i == ni
      a_n[p] = 0.0
    else
      k_n = k[i+1, j, k_idx]
      a_n[p] = (2.0 * k_p * k_n / (k_p + k_n)) * dx * dz_k / dy
    end

    # 下（k-1）
    if k_idx == 1
      a_b[p] = 0.0
    else
      k_b = k[i, j, k_idx-1]
      a_b[p] = (2.0 * k_p * k_b / (k_p + k_b)) * dx * dy / dz_b_k
    end

    # 上（k+1）
    if k_idx == nk
      a_t[p] = 0.0
    else
      k_t = k[i, j, k_idx+1]
      a_t[p] = (2.0 * k_p * k_t / (k_p + k_t)) * dx * dy / dz_t_k
    end

    # 対角項（中心係数）
    a_p[p] = a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p] + a_p_0

    # RHS（初期温度項）
    rhs = a_p_0 * T_initial[i, j, k_idx]

    # 表面（上端、k=nk）での熱流束境界条件
    if k_idx == nk
      rhs += q_surface[i, j] * dx * dy
    end

    b[p] = rhs
  end

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b
end


"""
    assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p) -> SparseMatrixCSC

DHCPの疎行列（CSC形式）を組み立て

7点ステンシル（対角+6方向）の疎行列を構築:
    A = [
      [a_p[1]    -a_e[1]   0     ...  -a_t[1]    ...]
      [-a_w[2]   a_p[2]   -a_e[2] ...  0         ...]
      [...]
    ]

JuliaではCSC（Compressed Sparse Column）形式がデフォルト
Pythonと異なり、列優先でメモリ効率が良い

対応Pythonコード: assemble_A_DHCP() (1131-1153行)

Args:
  ni, nj, nk: 格子点数
  a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列 (N,)

Returns:
  A: CSC疎行列 (N, N)
"""
function assemble_dhcp_matrix(
  ni::Int, nj::Int, nk::Int,
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
  off_k = ni * nj     # k方向（上下）

  # COO形式（I, J, V）で係数を収集
  I_list = Int[]
  J_list = Int[]
  V_list = Float64[]

  # 予約サイズ（7点ステンシル: 最大7*N個の非ゼロ要素）
  sizehint!(I_list, 7 * N)
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

  # COO→CSC変換（Juliaのデフォルト疎行列形式）
  A = sparse(I_list, J_list, V_list, N, N)

  return A
end


"""
    solve_dhcp!(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
                dx, dy, dz, dz_b, dz_t, dt;
                rtol=1e-8, maxiter=1000, verbose=false) -> T_all

複数時間ステップDHCPソルバー（ホットスタート対応）

各時間ステップで:
1. 前ステップ温度から熱物性値計算
2. 係数とRHS構築
3. 疎行列組み立て
4. CG法求解（対角前処理、ホットスタート）
5. 次ステップへ

陰解法の特徴:
  - 無条件安定（dt制約なし）
  - 大規模問題でも安定
  - CG収束が高速（良条件数）

対応Pythonコード: multiple_time_step_solver_DHCP() (1156-1219行)

Args:
  T_initial: 初期温度場 (ni, nj, nk) [K]
  q_surface: 表面熱流束時系列 (ni, nj, nt-1) [W/m²] ※Phase 2.2: 時間次元を最後に配置
  nt: 時間ステップ数
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数 [c0, c1, c2, c3]
  k_coeffs: 熱伝導率多項式係数 [k0, k1, k2, k3]
  dx, dy, dz, dz_b, dz_t: 格子情報
  dt: 時間刻み [s]
  rtol: CG相対許容誤差（デフォルト: 1e-8）
  maxiter: CG最大反復回数（デフォルト: 1000）
  verbose: 進捗表示フラグ（デフォルト: false）

Returns:
  T_all: 温度場時系列 (ni, nj, nk, nt) [K] ※Phase 2.2: 時間次元を最後に配置
"""
function solve_dhcp!(
  T_initial::Array{Float64,3},
  q_surface::Array{Float64,3},
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
  ni, nj, nk = size(T_initial)
  N = ni * nj * nk

  # 結果配列の初期化（メモリレイアウト最適化: Phase 2.2）
  # 時間次元を最後に配置して空間方向のメモリ連続性を確保
  T_all = zeros(Float64, ni, nj, nk, nt)
  T_all[:, :, :, 1] = T_initial

  # ホットスタート用初期推定値
  x0 = zeros(Float64, N)
  for k_idx in 1:nk, j in 1:nj, i in 1:ni
    p = dhcp_index(i, j, k_idx, ni, nj)
    x0[p] = T_initial[i, j, k_idx]
  end

  if verbose
    println("="^60)
    println("DHCP直接ソルバー開始")
    println("="^60)
    println("格子: $(ni)×$(nj)×$(nk) (N=$(N))")
    println("時間ステップ: $(nt), dt=$(dt)s")
    println("CG許容誤差: rtol=$(rtol), maxiter=$(maxiter)")
    println("="^60)
  end

  # 時間積分ループ
  for t in 2:nt
    # 前ステップ温度から熱物性値計算（メモリビュー使用: Phase 2.2）
    T_prev = @view T_all[:, :, :, t-1]
    cp, k = thermal_properties_calculator(T_prev, cp_coeffs, k_coeffs)

    # 係数とRHS構築（熱流束も時間次元を最後に: Phase 2.2）
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_dhcp_system!(
      T_prev, @view(q_surface[:, :, t-1]), rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 疎行列組み立て
    A = assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    # 対角前処理器（Jacobi前処理）
    diag_A = diag(A)
    inv_diag = similar(diag_A)
    for i in 1:N
      inv_diag[i] = diag_A[i] != 0.0 ? 1.0 / diag_A[i] : 0.0
    end
    Pl = Diagonal(inv_diag)

    # CG法求解（IterativeSolvers.jl）
    # ホットスタート: 前ステップ解を初期推定値に使用
    if verbose
      result = cg!(x0, A, b; Pl=Pl, reltol=rtol, maxiter=maxiter, log=true)
      history = result[2]
      if history.isconverged
        println("[t=$(t)/$(nt)] CG収束: $(history.iters)回")
      else
        @warn "[t=$(t)/$(nt)] CG未収束: $(history.iters)回"
      end
    else
      cg!(x0, A, b; Pl=Pl, reltol=rtol, maxiter=maxiter)
    end

    # 結果を3D配列に復元（Fortran順序、メモリレイアウト最適化: Phase 2.2）
    for k_idx in 1:nk, j in 1:nj, i in 1:ni
      p = dhcp_index(i, j, k_idx, ni, nj)
      T_all[i, j, k_idx, t] = x0[p]
    end

    # 数値異常チェック
    if any(isnan.(x0)) || any(isinf.(x0))
      error("[t=$(t)] 数値異常が発生しました（NaN/Inf検出）")
    end
  end

  if verbose
    println("="^60)
    println("DHCP直接ソルバー完了")
    println("  最終温度範囲: $(minimum(T_all)) - $(maximum(T_all)) K")
    println("="^60)
  end

  return T_all
end


"""
    solve_dhcp_multiple_timesteps(...)

solve_dhcp!のエイリアス関数（完全一致検証スクリプト用）

Python版のmultiple_time_step_solver_DHCPに対応する名前で呼び出せるようにします。
"""
function solve_dhcp_multiple_timesteps(
  T_initial::Array{Float64,3},
  q_surface::Array{Float64,3},
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
  return solve_dhcp!(
    T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt;
    rtol=rtol, maxiter=maxiter, verbose=verbose
  )
end


end  # module DHCPSolver