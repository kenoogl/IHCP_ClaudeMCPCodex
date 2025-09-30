#!/usr/bin/env python3
"""
Phase 3参照データ生成スクリプト（Adjoint随伴ソルバー）

このスクリプトは、Julia実装の検証用に小規模な随伴問題の参照データを生成します：
1. DHCP直接解を計算
2. 疑似観測データを生成（ノイズ付加）
3. Adjoint随伴場を計算

生成データ：
- DHCP温度場（T_cal）
- 疑似観測データ（Y_obs）
- 随伴場の時間発展（λ_all）
- 係数行列とRHS（後退時間ステップ）
- CG収束履歴
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg, LinearOperator
from scipy.sparse import diags
import json
import os

# ===== Pythonオリジナルコードからの関数抽出 =====

def polyval_numba_simple(coeffs, x):
  """
  多項式評価（Pythonオリジナル57-61行と同じ）

  Args:
    coeffs: 多項式係数 [a, b, c, d] for a*x^3 + b*x^2 + c*x + d
    x: 評価点

  Returns:
    多項式の値
  """
  result = 0.0
  for i in range(len(coeffs)):
    result += coeffs[i] * x ** (len(coeffs) - i - 1)
  return result


def thermal_properties_calculator_simple(T, cp_coeffs, k_coeffs):
  """
  簡易版熱物性値計算（Pythonオリジナル63-78行と同じ）

  Args:
    T: 温度場 (ni, nj, nk)
    cp_coeffs: 比熱多項式係数 [a, b, c, d] for a*T^3 + b*T^2 + c*T + d
    k_coeffs: 熱伝導率多項式係数 [a, b, c, d] for a*T^3 + b*T^2 + c*T + d

  Returns:
    cp: 比熱場 (ni, nj, nk)
    k: 熱伝導率場 (ni, nj, nk)
  """
  ni, nj, nk = T.shape
  cp = np.zeros_like(T)
  k = np.zeros_like(T)

  for i in range(ni):
    for j in range(nj):
      for k_idx in range(nk):
        T_val = T[i, j, k_idx]
        # polyval_numbaと同じ多項式評価
        cp[i, j, k_idx] = polyval_numba_simple(cp_coeffs, T_val)
        k[i, j, k_idx] = polyval_numba_simple(k_coeffs, T_val)

  return cp, k


def coeffs_and_rhs_building_DHCP(T_initial, q_surface, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt):
  """
  直接熱伝導問題（DHCP）の係数とRHS構築

  Args:
    T_initial: 初期温度場 (ni, nj, nk)
    q_surface: 表面熱流束 (ni, nj)
    rho: 密度 [kg/m³]
    cp: 比熱場 (ni, nj, nk)
    k: 熱伝導率場 (ni, nj, nk)
    dx, dy: x, y方向格子幅 [m]
    dz: z方向格子幅配列 (nk,) [m]
    dz_b: 下側界面距離 (nk,) [m]
    dz_t: 上側界面距離 (nk,) [m]
    dt: 時間刻み [s]

  Returns:
    a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列 (N,)
    b: RHSベクトル (N,)
  """
  ni, nj, nk = T_initial.shape
  N = ni * nj * nk

  a_w = np.zeros(N)
  a_e = np.zeros(N)
  a_s = np.zeros(N)
  a_n = np.zeros(N)
  a_b = np.zeros(N)
  a_t = np.zeros(N)
  a_p = np.zeros(N)
  b   = np.zeros(N)

  for p in range(N):
    # Fortran順序（列優先）のインデックス変換
    i = p % ni
    j = (p // ni) % nj
    k_idx = p // (ni * nj)

    dz_k = dz[k_idx]
    dz_t_k = dz_t[k_idx]
    dz_b_k = dz_b[k_idx]

    k_p = k[i, j, k_idx]

    # 時間項（蓄熱項）
    a_p_0 = rho * cp[i, j, k_idx] * dx * dy * dz_k / dt

    # 6方向の熱伝導係数（調和平均）
    a_w[p] = 0.0 if j == 0      else (2 * k_p * k[i, j-1, k_idx] / (k_p + k[i, j-1, k_idx])) * dy * dz_k / dx
    a_e[p] = 0.0 if j == nj-1   else (2 * k_p * k[i, j+1, k_idx] / (k_p + k[i, j+1, k_idx])) * dy * dz_k / dx
    a_s[p] = 0.0 if i == 0      else (2 * k_p * k[i-1, j, k_idx] / (k_p + k[i-1, j, k_idx])) * dx * dz_k / dy
    a_n[p] = 0.0 if i == ni-1   else (2 * k_p * k[i+1, j, k_idx] / (k_p + k[i+1, j, k_idx])) * dx * dz_k / dy
    a_b[p] = 0.0 if k_idx == 0  else (2 * k_p * k[i, j, k_idx-1] / (k_p + k[i, j, k_idx-1])) * dx * dy / dz_b_k
    a_t[p] = 0.0 if k_idx == nk-1 else (2 * k_p * k[i, j, k_idx+1] / (k_p + k[i, j, k_idx+1])) * dx * dy / dz_t_k

    # 対角項
    a_p[p] = a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p] + a_p_0

    # RHS（初期温度項 + 境界条件）
    rhs = a_p_0 * T_initial[i, j, k_idx]

    # 表面（上端）での熱流束境界条件
    if k_idx == nk - 1:
      rhs += q_surface[i, j] * dx * dy

    b[p] = rhs

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p):
  """DHCPの疎行列を組み立て"""
  N = ni * nj * nk

  # オフセット計算
  off_i = 1           # i方向（南北）
  off_j = ni          # j方向（西東）
  off_k = ni * nj     # k方向（上下）

  # COO形式で係数を収集
  I_list = []
  J_list = []
  V_list = []

  for p in range(N):
    # 対角成分
    I_list.append(p)
    J_list.append(p)
    V_list.append(a_p[p])

    # i-1 (南)
    if p >= off_i:
      I_list.append(p)
      J_list.append(p - off_i)
      V_list.append(-a_s[p])

    # i+1 (北)
    if p < N - off_i:
      I_list.append(p)
      J_list.append(p + off_i)
      V_list.append(-a_n[p])

    # j-1 (西)
    if p >= off_j:
      I_list.append(p)
      J_list.append(p - off_j)
      V_list.append(-a_w[p])

    # j+1 (東)
    if p < N - off_j:
      I_list.append(p)
      J_list.append(p + off_j)
      V_list.append(-a_e[p])

    # k-1 (下)
    if p >= off_k:
      I_list.append(p)
      J_list.append(p - off_k)
      V_list.append(-a_b[p])

    # k+1 (上)
    if p < N - off_k:
      I_list.append(p)
      J_list.append(p + off_k)
      V_list.append(-a_t[p])

  # COO→CSR変換
  A_csr = sp.coo_matrix((V_list, (I_list, J_list)), shape=(N, N)).tocsr()

  return A_csr


def multiple_time_step_solver_DHCP(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
                                   dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000):
  """複数時間ステップDHCPソルバー"""
  ni, nj, nk = T_initial.shape
  N = ni * nj * nk

  T_all = np.zeros((nt, ni, nj, nk))
  T_all[0] = T_initial.copy()

  cg_iters = []

  # ホットスタート用初期推定値
  x0 = T_initial.ravel(order='F').copy()

  for t in range(1, nt):
    # 前ステップ温度から熱物性値計算
    cp, k = thermal_properties_calculator_simple(T_all[t-1], cp_coeffs, k_coeffs)

    # 係数とRHS構築
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
      T_all[t-1], q_surface[t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 疎行列組み立て
    A_csr = assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    # 対角前処理器
    diag = A_csr.diagonal()
    inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

    def Mv(z):
      return inv_diag * z

    M = LinearOperator(A_csr.shape, matvec=Mv)

    # CG法求解
    iter_count = [0]
    def callback(xk):
      iter_count[0] += 1

    x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0, callback=callback)

    if info > 0:
      print(f"[警告][DHCP t={t}] CG法が収束しませんでした (info={info}, iter={iter_count[0]})")
    elif info < 0:
      print(f"[エラー][DHCP t={t}] CG法の入力が不正です (info={info})")

    # 結果保存とホットスタート更新
    T_all[t] = x.reshape((ni, nj, nk), order='F')
    x0 = x.copy()
    cg_iters.append(iter_count[0])

  return T_all, np.array(cg_iters)


# ===== Adjoint随伴ソルバー =====

def coeffs_and_rhs_building_Adjoint(lambda_initial, T_cal, Y_obs, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt):
  """
  随伴方程式の係数とRHS構築（Pythonオリジナル1228-1270行）

  Args:
    lambda_initial: 前ステップ随伴場 (ni, nj, nk) - 時間反転なので「前」=時間的に後
    T_cal: DHCP計算温度（底面のみ使用） (ni, nj)
    Y_obs: 観測温度（底面） (ni, nj)
    rho, cp, k, dx, dy, dz, dz_b, dz_t, dt: DHCP同様

  Returns:
    a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列 (N,) - DHCPと同じ構造
    b: RHSベクトル (N,) - 残差注入あり
  """
  ni, nj, nk = lambda_initial.shape
  N = ni * nj * nk

  a_w = np.zeros(N)
  a_e = np.zeros(N)
  a_s = np.zeros(N)
  a_n = np.zeros(N)
  a_b = np.zeros(N)
  a_t = np.zeros(N)
  a_p = np.zeros(N)
  b   = np.zeros(N)

  for p in range(N):
    # Fortran順序（列優先）のインデックス変換
    i = p % ni
    j = (p // ni) % nj
    k_idx = p // (ni * nj)

    dz_k = dz[k_idx]
    dz_t_k = dz_t[k_idx]
    dz_b_k = dz_b[k_idx]

    k_p = k[i, j, k_idx]

    # 時間項（DHCPと同じ）
    a_p_0 = rho * cp[i, j, k_idx] * dx * dy * dz_k / dt

    # 6方向の熱伝導係数（DHCPと同じ - 随伴方程式も同じ行列）
    a_w[p] = 0.0 if j == 0      else (2 * k_p * k[i, j-1, k_idx] / (k_p + k[i, j-1, k_idx])) * dy * dz_k / dx
    a_e[p] = 0.0 if j == nj-1   else (2 * k_p * k[i, j+1, k_idx] / (k_p + k[i, j+1, k_idx])) * dy * dz_k / dx
    a_s[p] = 0.0 if i == 0      else (2 * k_p * k[i-1, j, k_idx] / (k_p + k[i-1, j, k_idx])) * dx * dz_k / dy
    a_n[p] = 0.0 if i == ni-1   else (2 * k_p * k[i+1, j, k_idx] / (k_p + k[i+1, j, k_idx])) * dx * dz_k / dy
    a_b[p] = 0.0 if k_idx == 0  else (2 * k_p * k[i, j, k_idx-1] / (k_p + k[i, j, k_idx-1])) * dx * dy / dz_b_k
    a_t[p] = 0.0 if k_idx == nk-1 else (2 * k_p * k[i, j, k_idx+1] / (k_p + k[i, j, k_idx+1])) * dx * dy / dz_t_k

    # 対角項（DHCPと同じ）
    a_p[p] = a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p] + a_p_0

    # RHS（時間項）
    rhs = a_p_0 * lambda_initial[i, j, k_idx]

    # 残差注入（底面のみ、Python 1266-1267行）
    if k_idx == 0:
      rhs += 2.0 * (T_cal[i, j] - Y_obs[i, j]) * dx * dy

    b[p] = rhs

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def assemble_A_Adjoint(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p):
  """
  随伴方程式の疎行列組み立て（Pythonオリジナル1272-1294行）

  DHCPと同じ構造なので、assemble_A_DHCPを流用
  """
  return assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)


def multiple_time_step_solver_Adjoint(T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
                                       dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000):
  """
  複数時間ステップAdjointソルバー（Pythonオリジナル1297-1364行）

  重要: 時間反転ループ（後退時間積分）

  Args:
    T_cal: DHCP計算温度場 (nt, ni, nj, nk)
    Y_obs: 観測温度（底面） (nt, ni, nj)
    nt: 時間ステップ数
    rho, cp_coeffs, k_coeffs, dx, dy, dz, dz_b, dz_t, dt: DHCP同様
    rtol, maxiter: CG法パラメータ

  Returns:
    lambda_all: 随伴場時系列 (nt, ni, nj, nk)
    cg_iters: CG反復回数履歴 (nt-1,)
  """
  ni, nj, nk = T_cal[0].shape
  N = ni * nj * nk

  lambda_all = np.zeros_like(T_cal)
  lambda_all[-1] = 0.0  # 終端条件（Python 1322行）

  cg_iters = []

  # ホットスタート用初期推定値
  x0 = lambda_all[-1].ravel(order='F').copy()

  # 後退時間ループ（Python 1328行: range(nt-2, -1, -1)）
  for t in range(nt-2, -1, -1):
    # 次ステップ（時間的に後）の随伴場を初期値とする
    lambda_initial = lambda_all[t+1]

    # 温度場から熱物性値計算（Python 1332行: T_cal[t]）
    cp, k = thermal_properties_calculator_simple(T_cal[t], cp_coeffs, k_coeffs)

    # 係数とRHS構築（Python 1334-1335行）
    # 残差注入: T_cal[t][:,:,0]（底面）と Y_obs[t] を使用
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_Adjoint(
      lambda_initial, T_cal[t][:,:,0], Y_obs[t], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 疎行列組み立て
    A_csr = assemble_A_Adjoint(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    # 対角前処理器
    diag = A_csr.diagonal()
    inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

    def Mv(z):
      return inv_diag * z

    M = LinearOperator(A_csr.shape, matvec=Mv)

    # CG法求解
    iter_count = [0]
    def callback(xk):
      iter_count[0] += 1

    x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0, callback=callback)

    if info > 0:
      print(f"[警告][Adjoint t={t}] CG法が収束しませんでした (info={info}, iter={iter_count[0]})")
    elif info < 0:
      print(f"[エラー][Adjoint t={t}] CG法の入力が不正です (info={info})")

    # 結果保存とホットスタート更新（Python 1356-1357行）
    lambda_all[t] = x.reshape((ni, nj, nk), order='F')
    x0 = x.copy()
    cg_iters.append(iter_count[0])

  return lambda_all, np.array(cg_iters)


# ===== 参照データ生成関数 =====

def generate_simple_adjoint_problem():
  """
  Adjoint問題1: 1D小規模問題

  問題設定:
    - 1×1×5の1D問題
    - 定数熱物性値（ρ=8000, cp=500, k=15）
    - 初期温度: T = 300K
    - 表面熱流束: q = 1000 W/m²
    - 10ステップDHCP実行
    - 疑似観測データ: Y_obs = T_cal[:,:,:,0] + ノイズ
    - Adjoint計算
  """
  print("\n" + "="*60)
  print("Adjoint問題1: 1D小規模問題")
  print("="*60)

  # 格子設定
  ni, nj, nk = 1, 1, 5
  dx, dy = 1e-3, 1e-3
  Lz = 2.5e-3

  # 等間隔格子
  z_faces = np.linspace(0, Lz, nk+1)
  z_centers = (z_faces[:-1] + z_faces[1:]) / 2
  dz = np.diff(z_faces)

  # 界面距離
  dz_t = np.full(nk, dz[0])
  dz_t[-1] = np.inf
  dz_b = np.full(nk, dz[0])
  dz_b[0] = np.inf

  # 物性値（定数）
  rho = 8000.0
  cp_const = 500.0
  k_const = 15.0
  # 係数順序: [a, b, c, d] for a*T^3 + b*T^2 + c*T + d
  # 定数の場合は [0, 0, 0, cp_const]
  cp_coeffs = np.array([0.0, 0.0, 0.0, cp_const])
  k_coeffs = np.array([0.0, 0.0, 0.0, k_const])

  # 初期条件
  T0 = 300.0
  T_initial = np.full((ni, nj, nk), T0)

  # 境界条件
  q = 1000.0
  nt = 11
  dt = 0.1
  q_surface = np.full((nt-1, ni, nj), q)

  # DHCP求解
  print("Step 1: DHCP計算中...")
  T_cal, dhcp_iters = multiple_time_step_solver_DHCP(
    T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-10, maxiter=1000
  )
  print(f"  DHCP完了: 平均CG反復回数={np.mean(dhcp_iters):.1f}")

  # 疑似観測データ生成（底面温度 + 小ノイズ）
  print("Step 2: 疑似観測データ生成...")
  Y_obs = np.zeros((nt, ni, nj))
  noise_level = 0.1  # K
  np.random.seed(42)
  for t in range(nt):
    Y_obs[t] = T_cal[t, :, :, 0] + noise_level * np.random.randn(ni, nj)
  print(f"  ノイズレベル: {noise_level} K")

  # Adjoint求解
  print("Step 3: Adjoint計算中...")
  lambda_all, adjoint_iters = multiple_time_step_solver_Adjoint(
    T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000
  )
  print(f"  Adjoint完了: 平均CG反復回数={np.mean(adjoint_iters):.1f}")

  # 統計情報
  residual = T_cal[:, :, :, 0] - Y_obs
  print(f"\n結果統計:")
  print(f"  DHCP温度範囲: {np.min(T_cal):.2f} - {np.max(T_cal):.2f} K")
  print(f"  残差範囲: {np.min(residual):.2e} - {np.max(residual):.2e} K")
  print(f"  随伴場範囲: {np.min(lambda_all):.2e} - {np.max(lambda_all):.2e}")

  # データ保存
  data = {
    'problem_type': 'adjoint_1D_simple',
    'grid': {'ni': ni, 'nj': nj, 'nk': nk, 'dx': dx, 'dy': dy, 'dz': dz.tolist()},
    'z_coords': {'faces': z_faces.tolist(), 'centers': z_centers.tolist(), 'dz_b': dz_b.tolist(), 'dz_t': dz_t.tolist()},
    'properties': {'rho': rho, 'cp_coeffs': cp_coeffs.tolist(), 'k_coeffs': k_coeffs.tolist()},
    'boundary': {'q_surface': q, 'T0': T0, 'noise_level': noise_level},
    'time': {'nt': nt, 'dt': dt},
    'T_initial': T_initial.tolist(),
    'T_cal': T_cal.tolist(),
    'Y_obs': Y_obs.tolist(),
    'lambda_all': lambda_all.tolist(),
    'dhcp_iters': dhcp_iters.tolist(),
    'adjoint_iters': adjoint_iters.tolist(),
    'stats': {
      'T_min': float(np.min(T_cal)),
      'T_max': float(np.max(T_cal)),
      'residual_min': float(np.min(residual)),
      'residual_max': float(np.max(residual)),
      'lambda_min': float(np.min(lambda_all)),
      'lambda_max': float(np.max(lambda_all)),
      'dhcp_cg_mean': float(np.mean(dhcp_iters)),
      'adjoint_cg_mean': float(np.mean(adjoint_iters))
    }
  }

  return data


def generate_3d_adjoint_problem():
  """
  Adjoint問題2: 3D小規模問題（温度依存熱物性値）

  問題設定:
    - 3×3×8の3D問題
    - 温度依存熱物性値（SUS304近似）
    - 初期温度: T = 300K
    - 表面熱流束: q = 3000 * sin(π*x/Lx) * sin(π*y/Ly) W/m²
    - 15ステップDHCP実行
    - 疑似観測データ: Y_obs = T_cal[:,:,:,0] + ノイズ
    - Adjoint計算
  """
  print("\n" + "="*60)
  print("Adjoint問題2: 3D小規模問題（温度依存熱物性値）")
  print("="*60)

  # 格子設定
  ni, nj, nk = 3, 3, 8
  dx, dy = 0.12e-3, 0.12e-3
  Lz = 0.4e-3

  # 非等間隔格子（表面集中）
  stretch_factor = 2.5
  z_faces_norm = np.linspace(1, 0, nk+1)
  z_faces = Lz - (Lz / (np.exp(stretch_factor) - 1)) * (np.exp(stretch_factor * z_faces_norm) - 1)

  z_centers = np.zeros(nk)
  z_centers[0] = z_faces[0]
  z_centers[-1] = z_faces[-1]
  z_centers[1:-1] = (z_faces[1:-2] + z_faces[2:-1]) / 2

  dz = np.diff(z_faces)

  dz_t = np.zeros(nk)
  dz_t[-1] = np.inf
  dz_t[:-1] = z_centers[1:] - z_centers[:-1]

  dz_b = np.zeros(nk)
  dz_b[0] = np.inf
  dz_b[1:] = z_centers[1:] - z_centers[:-1]

  # 物性値（SUS304近似、温度依存）
  rho = 7900.0
  # 係数順序: [a, b, c, d] for a*T^3 + b*T^2 + c*T + d
  # SUS304近似: cp ≈ d + c*T, k ≈ d + c*T
  cp_coeffs = np.array([0.0, 0.0, 0.134, 462.0])
  k_coeffs = np.array([0.0, 0.0, 0.0127, 14.6])

  # 初期条件
  T0 = 300.0
  T_initial = np.full((ni, nj, nk), T0)

  # 境界条件（空間変動熱流束）
  nt = 16
  dt = 1e-3
  q_surface = np.zeros((nt-1, ni, nj))

  Lx, Ly = dx * ni, dy * nj
  x = np.linspace(0, Lx, ni)
  y = np.linspace(0, Ly, nj)
  X, Y = np.meshgrid(x, y, indexing='ij')

  q_base = 3000.0
  for t in range(nt-1):
    q_surface[t] = q_base * np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)

  # DHCP求解
  print("Step 1: DHCP計算中...")
  T_cal, dhcp_iters = multiple_time_step_solver_DHCP(
    T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000
  )
  print(f"  DHCP完了: 平均CG反復回数={np.mean(dhcp_iters):.1f} ± {np.std(dhcp_iters):.1f}")

  # 疑似観測データ生成（底面温度 + ノイズ）
  print("Step 2: 疑似観測データ生成...")
  Y_obs = np.zeros((nt, ni, nj))
  noise_level = 0.2  # K
  np.random.seed(123)
  for t in range(nt):
    Y_obs[t] = T_cal[t, :, :, 0] + noise_level * np.random.randn(ni, nj)
  print(f"  ノイズレベル: {noise_level} K")

  # Adjoint求解
  print("Step 3: Adjoint計算中...")
  lambda_all, adjoint_iters = multiple_time_step_solver_Adjoint(
    T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000
  )
  print(f"  Adjoint完了: 平均CG反復回数={np.mean(adjoint_iters):.1f} ± {np.std(adjoint_iters):.1f}")

  # 統計情報
  residual = T_cal[:, :, :, 0] - Y_obs
  print(f"\n結果統計:")
  print(f"  DHCP温度範囲: {np.min(T_cal):.2f} - {np.max(T_cal):.2f} K")
  print(f"  温度上昇: {np.max(T_cal) - T0:.2f} K")
  print(f"  残差範囲: {np.min(residual):.2e} - {np.max(residual):.2e} K")
  print(f"  残差RMS: {np.sqrt(np.mean(residual**2)):.2e} K")
  print(f"  随伴場範囲: {np.min(lambda_all):.2e} - {np.max(lambda_all):.2e}")

  # データ保存
  data = {
    'problem_type': 'adjoint_3D_temperature_dependent',
    'grid': {'ni': ni, 'nj': nj, 'nk': nk, 'dx': dx, 'dy': dy, 'dz': dz.tolist()},
    'z_coords': {'faces': z_faces.tolist(), 'centers': z_centers.tolist(), 'dz_b': dz_b.tolist(), 'dz_t': dz_t.tolist()},
    'properties': {'rho': rho, 'cp_coeffs': cp_coeffs.tolist(), 'k_coeffs': k_coeffs.tolist()},
    'boundary': {'q_surface_base': q_base, 'T0': T0, 'noise_level': noise_level},
    'time': {'nt': nt, 'dt': dt},
    'T_initial': T_initial.tolist(),
    'q_surface': q_surface.tolist(),
    'T_cal': T_cal.tolist(),
    'Y_obs': Y_obs.tolist(),
    'lambda_all': lambda_all.tolist(),
    'dhcp_iters': dhcp_iters.tolist(),
    'adjoint_iters': adjoint_iters.tolist(),
    'stats': {
      'T_min': float(np.min(T_cal)),
      'T_max': float(np.max(T_cal)),
      'dT': float(np.max(T_cal) - T0),
      'residual_min': float(np.min(residual)),
      'residual_max': float(np.max(residual)),
      'residual_rms': float(np.sqrt(np.mean(residual**2))),
      'lambda_min': float(np.min(lambda_all)),
      'lambda_max': float(np.max(lambda_all)),
      'dhcp_cg_mean': float(np.mean(dhcp_iters)),
      'dhcp_cg_std': float(np.std(dhcp_iters)),
      'adjoint_cg_mean': float(np.mean(adjoint_iters)),
      'adjoint_cg_std': float(np.std(adjoint_iters))
    }
  }

  return data


def main():
  """メイン処理"""
  print("="*60)
  print("Phase 3参照データ生成スクリプト（Adjoint随伴ソルバー）")
  print("="*60)

  # 出力ディレクトリ確認
  output_dir = os.path.dirname(os.path.abspath(__file__))
  print(f"出力ディレクトリ: {output_dir}")

  # Adjoint問題1（1D小規模）
  data1 = generate_simple_adjoint_problem()
  output_file1 = os.path.join(output_dir, 'phase3_reference_adjoint_1D.json')
  with open(output_file1, 'w') as f:
    json.dump(data1, f, indent=2)
  print(f"✓ 保存完了: {output_file1}")

  # Adjoint問題2（3D小規模）
  data2 = generate_3d_adjoint_problem()
  output_file2 = os.path.join(output_dir, 'phase3_reference_adjoint_3D.json')
  with open(output_file2, 'w') as f:
    json.dump(data2, f, indent=2)
  print(f"✓ 保存完了: {output_file2}")

  print("\n" + "="*60)
  print("参照データ生成完了")
  print("="*60)
  print("\nJuliaテストでの使用方法:")
  print("  julia> using JSON")
  print("  julia> data = JSON.parsefile(\"phase3_reference_adjoint_1D.json\")")
  print("  julia> T_cal = data[\"T_cal\"]")
  print("  julia> Y_obs = data[\"Y_obs\"]")
  print("  julia> lambda_all = data[\"lambda_all\"]")
  print("="*60)


if __name__ == '__main__':
  main()