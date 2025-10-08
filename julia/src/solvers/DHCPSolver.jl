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

using LinearAlgebra
using FLoops
using SparseArrays
using IterativeSolvers

import ..Commons
using ..Commons: WorkBuffers, λf, get_backend, compute_z_range

import ..ThermalProperties
using ..ThermalProperties: thermal_properties!, set_properties!

import ..BoundaryConditions
using ..BoundaryConditions: BoundaryCondition, BoundaryConditionSet,
            isothermal_bc, heat_flux_bc, adiabatic_bc, convection_bc,
            create_boundary_conditions, apply_boundary_conditions!,
            print_boundary_conditions, ISOTHERMAL, HEAT_FLUX,
            apply_face_boundary!, set_BC_coef

import ..CommonSolver
using ..CommonSolver: PBiCGSTAB!

export solve_dhcp!, solve_dhcp_cg!

"""
順問題の境界条件  
Z下面: 断熱、Z上面: 熱流束、側面: 断熱
"""

function set_dhcp_bc_parameters(nk::Int)
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
    q_dist::Bool,
    ρ::Float64,
    par::String
    )

    backend = get_backend(par)
    SZ = size(wk.b)
    dx1 = 1.0 / dx
    dy1 = 1.0 / dy
    z_st = z_range[1]
    z_ed = z_range[2]
    ddt = 1.0 / Δt
    
    @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
        wk.b[i,j,k] = 0.0
    end

    # 領域境界面に一様な熱流束の場合
    # X_minus
    if HF[1] != 0.0
        i=2
        @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1
            wk.b[i,j,k] += HF[1] * dx1
        end
    end

    # X_plus
    if HF[2] != 0.0
        i=SZ[1]-1
        @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1
            wk.b[i,j,k] -= HF[2] * dx1
        end
    end

    # Y_minus
    if HF[3] != 0.0
        j=2
        @floop backend for k in z_st:z_ed, i in 2:SZ[1]-1
            wk.b[i,j,k] += HF[3] * dy1
        end
    end

    # Y_plus
    if HF[4] != 0.0
        j=SZ[2]-1
        @floop backend for k in z_st:z_ed, i in 2:SZ[1]-1
            wk.b[i,j,k] -= HF[4] * dy1
        end
    end

    # Z_minus
    if HF[5] != 0.0
        k=z_st
        a = HF[5] / ΔZ[k]
        @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
            wk.b[i,j,k] += a
        end
    end

    # Z_plus
    if HF[6] != 0.0
        k=z_ed
        a = HF[6] / ΔZ[k]
        @floop for j in 2:SZ[2]-1, i in 2:SZ[1]-1
            wk.b[i,j,k] -= a
        end
    end

    # IHCPの場合、Z方向のみ分布を考慮した熱流束、必要があれば他の面も同様に処理
    # Z_plus
    if q_dist == true
        k=z_ed
        a = 1.0 / ΔZ[k]
        @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
            wk.b[i,j,k] -= qsrf[i-1,j-1] * a
        end
    end

    # RHS
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        wk.b[i,j,k] = -( ddt * wk.θ[i,j,k]
                        + (wk.hsrc[i,j,k] + wk.b[i,j,k]) / (ρ * wk.cp[i,j,k])  )
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
dx, dy, Z, ΔZ, Δt: 格子・時間パラメータ
rtol, maxiter: 収束パラメータ
verbose: 進捗表示フラグ（デフォルト: false）
par: バックエンド

# 戻り値
T_all: 新時刻の温度場 (ni, nj, nk, nt) [K] 
"""

function solve_dhcp!(
  T_initial::Array{Float64,3},
  q_surface::Array{Float64,3},
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
    println("Start DHCP direct solver")
    println("="^60)
    println("格子: $(ni)×$(nj)×$(nk) (N=$(N))")
    println("時間ステップ: $(nt), dt=$(dt)s")
    println("CG許容誤差: rtol=$(rtol), maxiter=$(maxiter)")
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
  z_range, bc_set = set_dhcp_bc_parameters(nk)
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
      @view(q_surface[:, :, t-1]),
      true, ρ, par)
    

    if verbose
      isconverged, itr, res0 = PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol*Nf, maxItr=maxiter, smoother="", par=par)
      if isconverged
        println("[t=$(t)/$(nt)] CG収束: $(itr)回 初期残差: $(res0)")
      else
        @warn "[t=$(t)/$(nt)] CG未収束: $(itr)回 初期残差: $(res0)"
      end
    else
      PBiCGSTAB!(wk, Δh, dt, Z, ΔZ, z_range, HT, ρ,
          tol=rtol*Nf, maxItr=maxiter, smoother="", par=par)
    end

    # 数値異常チェック
    if any(isnan.(wk.θ)) || any(isinf.(wk.θ))
      error("[t=$(t)] 数値異常が発生しました（NaN/Inf検出）")
    end

    # ガイドセルを除いて内点データを返す
    for k in 1:nk, j in 1:nj, i in 1:ni
      T_all[i, j, k, t] = wk.θ[i+1, j+1, k+1]
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


# ============================================================================
# 感度問題用の旧API関数群（疎行列cg!使用）
# ============================================================================

"""
    dhcp_index(i, j, k, ni, nj) -> Int

Fortran順序（列優先）でのグローバルインデックス計算

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

直接熱伝導問題（DHCP）の係数とRHSベクトル構築（旧版）

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
  T_initial::Array{Float64,3},
  q_surface::Array{Float64,2},
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

  # 全格子点ループ
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

DHCPの疎行列（CSC形式）を組み立て（旧版）

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

  # COO→CSC変換
  A = sparse(I_list, J_list, V_list, N, N)

  return A
end


"""
    solve_dhcp_cg!(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
                   dx, dy, dz, dz_b, dz_t, dt;
                   rtol=1e-8, maxiter=1000, verbose=false) -> T_all

感度問題用DHCPソルバー（旧API、疎行列cg!使用）

Args:
  T_initial: 初期温度場 (ni, nj, nk) [K]
  q_surface: 表面熱流束時系列 (nt-1, ni, nj) [W/m²]
  nt: 時間ステップ数
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数
  k_coeffs: 熱伝導率多項式係数
  dx, dy, dz, dz_b, dz_t: 格子情報
  dt: 時間刻み [s]
  rtol: CG相対許容誤差
  maxiter: CG最大反復回数
  verbose: 進捗表示フラグ

Returns:
  T_all: 温度場時系列 (nt, ni, nj, nk) [K]
"""
function solve_dhcp_cg!(
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

  # 結果配列の初期化
  T_all = zeros(Float64, nt, ni, nj, nk)
  T_all[1, :, :, :] = T_initial

  # ホットスタート用初期推定値
  x0 = zeros(Float64, N)
  for k_idx in 1:nk, j in 1:nj, i in 1:ni
    p = dhcp_index(i, j, k_idx, ni, nj)
    x0[p] = T_initial[i, j, k_idx]
  end

  if verbose
    println("="^60)
    println("DHCP-CG感度問題ソルバー開始")
    println("="^60)
    println("格子: $(ni)×$(nj)×$(nk) (N=$(N))")
    println("時間ステップ: $(nt), dt=$(dt)s")
    println("CG許容誤差: rtol=$(rtol), maxiter=$(maxiter)")
    println("="^60)
  end

  # 熱物性値配列（事前確保）
  cp = zeros(Float64, ni, nj, nk)
  k = zeros(Float64, ni, nj, nk)

  # 時間積分ループ
  for t in 2:nt
    # 前ステップ温度から熱物性値計算
    T_prev = T_all[t-1, :, :, :]
    thermal_properties!(T_prev, cp, k, cp_coeffs, k_coeffs)

    # 係数とRHS構築
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_dhcp_system!(
      T_prev, q_surface[t-1, :, :], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 疎行列組み立て
    A = assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    # 対角前処理器
    diag_A = diag(A)
    inv_diag = similar(diag_A)
    for i in 1:N
      inv_diag[i] = diag_A[i] != 0.0 ? 1.0 / diag_A[i] : 0.0
    end
    Pl = Diagonal(inv_diag)

    # CG法求解
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

    # 結果を3D配列に復元
    for k_idx in 1:nk, j in 1:nj, i in 1:ni
      p = dhcp_index(i, j, k_idx, ni, nj)
      T_all[t, i, j, k_idx] = x0[p]
    end

    # 数値異常チェック
    if any(isnan.(x0)) || any(isinf.(x0))
      error("[t=$(t)] 数値異常が発生しました（NaN/Inf検出）")
    end
  end

  if verbose
    println("="^60)
    println("DHCP-CG感度問題ソルバー完了")
    println("  最終温度範囲: $(minimum(T_all)) - $(maximum(T_all)) K")
    println("="^60)
  end

  return T_all
end


end  # module DHCPSolver