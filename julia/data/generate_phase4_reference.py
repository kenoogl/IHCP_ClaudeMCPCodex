#!/usr/bin/env python3
"""
Phase 4参照データ生成スクリプト（CGM共役勾配法）

このスクリプトは、Julia実装の検証用に小規模なCGM逆問題の参照データを生成します：
1. 真の熱流束を設定（既知値）
2. DHCP直接解を計算（疑似観測データ生成）
3. 観測データにノイズを付加
4. 初期推定値からCGM最適化を実行
5. 収束履歴と逆解析結果を記録

生成データ：
- 真の熱流束（q_true）
- 疑似観測データ（Y_obs、ノイズ付加）
- CGM初期推定値（q_init）
- CGM反復履歴（J_hist, q_hist, gradient_hist, beta_hist）
- 最終逆解析結果（q_final）
- 各反復での中間場（T_cal, λ_field, dT）

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- global_CGM_time() (1371-1553行)
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


def coeffs_and_rhs_building_DHCP(T_initial, q_surface, rho, cp, k_field, dx, dy, dz, dz_b, dz_t, dt):
  """
  直接熱伝導問題（DHCP）の係数とRHS構築

  Pythonオリジナル: 1087-1129行
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
  b = np.zeros(N)

  def idx(i, j, kk):
    return i + j * ni + kk * ni * nj

  for kk in range(nk):
    for j in range(nj):
      for i in range(ni):
        p = idx(i, j, kk)

        V = dx * dy * dz[kk]
        rho_cp_V_dt = rho * cp[i, j, kk] * V / dt

        # 西（i-1）
        if i > 0:
          k_w = 0.5 * (k_field[i, j, kk] + k_field[i-1, j, kk])
          a_w[p] = k_w * dy * dz[kk] / dx

        # 東（i+1）
        if i < ni - 1:
          k_e = 0.5 * (k_field[i, j, kk] + k_field[i+1, j, kk])
          a_e[p] = k_e * dy * dz[kk] / dx

        # 南（j-1）
        if j > 0:
          k_s = 0.5 * (k_field[i, j, kk] + k_field[i, j-1, kk])
          a_s[p] = k_s * dx * dz[kk] / dy

        # 北（j+1）
        if j < nj - 1:
          k_n = 0.5 * (k_field[i, j, kk] + k_field[i, j+1, kk])
          a_n[p] = k_n * dx * dz[kk] / dy

        # 下（k-1、底面）
        if kk > 0:
          k_b = 0.5 * (k_field[i, j, kk] + k_field[i, j, kk-1])
          a_b[p] = k_b * dx * dy / dz_b[kk]

        # 上（k+1、表面）
        if kk < nk - 1:
          k_t = 0.5 * (k_field[i, j, kk] + k_field[i, j, kk+1])
          a_t[p] = k_t * dx * dy / dz_t[kk]

        # 対角成分
        a_p[p] = rho_cp_V_dt + a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p]

        # RHS（前ステップ温度）
        b[p] = rho_cp_V_dt * T_initial[i, j, kk]

        # 表面境界条件（k=nk-1、熱流束注入）
        if kk == nk - 1:
          b[p] += q_surface[i, j] * dx * dy

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def assemble_A_DHCP(a_w, a_e, a_s, a_n, a_b, a_t, a_p, ni, nj, nk):
  """
  DHCP係数行列をCSR形式で組み立て

  Pythonオリジナル: 1131-1153行
  """
  N = ni * nj * nk

  def idx(i, j, kk):
    return i + j * ni + kk * ni * nj

  row = []
  col = []
  data = []

  for kk in range(nk):
    for j in range(nj):
      for i in range(ni):
        p = idx(i, j, kk)

        # 対角
        row.append(p)
        col.append(p)
        data.append(a_p[p])

        # 西
        if i > 0:
          row.append(p)
          col.append(idx(i-1, j, kk))
          data.append(-a_w[p])

        # 東
        if i < ni - 1:
          row.append(p)
          col.append(idx(i+1, j, kk))
          data.append(-a_e[p])

        # 南
        if j > 0:
          row.append(p)
          col.append(idx(i, j-1, kk))
          data.append(-a_s[p])

        # 北
        if j < nj - 1:
          row.append(p)
          col.append(idx(i, j+1, kk))
          data.append(-a_n[p])

        # 下
        if kk > 0:
          row.append(p)
          col.append(idx(i, j, kk-1))
          data.append(-a_b[p])

        # 上
        if kk < nk - 1:
          row.append(p)
          col.append(idx(i, j, kk+1))
          data.append(-a_t[p])

  A_csr = sp.csr_matrix((data, (row, col)), shape=(N, N))
  return A_csr


def multiple_time_step_solver_DHCP(T_init, q, nt, rho, cp_coeffs, k_coeffs,
                                    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-6, maxiter=20000):
  """
  複数時間ステップDHCP求解（ホットスタート対応）

  Pythonオリジナル: 1156-1219行

  Args:
    T_init: 初期温度場 (ni, nj, nk)
    q: 表面熱流束 (nt-1, ni, nj)
    nt: 時間ステップ数
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
    T_all: 温度場の時間発展 (nt, ni, nj, nk)
  """
  ni, nj, nk = T_init.shape
  T_all = np.zeros((nt, ni, nj, nk))
  T_all[0] = T_init.copy()

  x0 = None  # ホットスタート用

  for t in range(1, nt):
    # 熱物性値計算
    cp, k = thermal_properties_calculator_simple(T_all[t-1], cp_coeffs, k_coeffs)

    # 係数とRHS構築
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
      T_all[t-1], q[t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 行列組み立て
    A_csr = assemble_A_DHCP(a_w, a_e, a_s, a_n, a_b, a_t, a_p, ni, nj, nk)

    # 対角前処理
    diag = A_csr.diagonal()
    inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)
    M = LinearOperator(A_csr.shape, matvec=lambda z: inv_diag * z)

    # CG求解
    x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)

    if info > 0:
      print(f"[警告][t={t}] DHCP CG未収束 (info={info})")
    elif info < 0:
      print(f"[エラー][t={t}] DHCP CG入力不正 (info={info})")

    # 温度場更新
    T_all[t] = x.reshape((ni, nj, nk), order='F')
    x0 = x  # ホットスタート

  return T_all


def coeffs_and_rhs_building_Adjoint(lambda_initial, T_cal_bottom, Y_obs,
                                      rho, cp, k_field, dx, dy, dz, dz_b, dz_t, dt):
  """
  随伴方程式の係数とRHS構築

  Pythonオリジナル: 1228-1270行
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
  b = np.zeros(N)

  def idx(i, j, kk):
    return i + j * ni + kk * ni * nj

  for kk in range(nk):
    for j in range(nj):
      for i in range(ni):
        p = idx(i, j, kk)

        V = dx * dy * dz[kk]
        rho_cp_V_dt = rho * cp[i, j, kk] * V / dt

        # DHCPと同じ係数（対称性）
        if i > 0:
          k_w = 0.5 * (k_field[i, j, kk] + k_field[i-1, j, kk])
          a_w[p] = k_w * dy * dz[kk] / dx

        if i < ni - 1:
          k_e = 0.5 * (k_field[i, j, kk] + k_field[i+1, j, kk])
          a_e[p] = k_e * dy * dz[kk] / dx

        if j > 0:
          k_s = 0.5 * (k_field[i, j, kk] + k_field[i, j-1, kk])
          a_s[p] = k_s * dx * dz[kk] / dy

        if j < nj - 1:
          k_n = 0.5 * (k_field[i, j, kk] + k_field[i, j+1, kk])
          a_n[p] = k_n * dx * dz[kk] / dy

        if kk > 0:
          k_b = 0.5 * (k_field[i, j, kk] + k_field[i, j, kk-1])
          a_b[p] = k_b * dx * dy / dz_b[kk]

        if kk < nk - 1:
          k_t = 0.5 * (k_field[i, j, kk] + k_field[i, j, kk+1])
          a_t[p] = k_t * dx * dy / dz_t[kk]

        a_p[p] = rho_cp_V_dt + a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p]

        # RHS（次ステップ随伴場）
        b[p] = rho_cp_V_dt * lambda_initial[i, j, kk]

        # 残差注入（底面 k=0）
        if kk == 0:
          residual = T_cal_bottom[i, j] - Y_obs[i, j]
          b[p] += 2.0 * residual * dx * dy

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def multiple_time_step_solver_Adjoint(T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
                                        dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=20000):
  """
  複数時間ステップAdjoint求解（後退時間積分）

  Pythonオリジナル: 1297-1364行

  Args:
    T_cal: DHCP計算温度場 (nt, ni, nj, nk)
    Y_obs: 観測温度（底面） (nt, ni, nj)
    nt: 時間ステップ数
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
    lambda_field: 随伴場の時間発展 (nt, ni, nj, nk)
  """
  ni, nj, nk = T_cal.shape[1:]
  lambda_field = np.zeros((nt, ni, nj, nk))

  # 最終時刻の随伴場はゼロ
  lambda_field[-1] = 0.0

  x0 = None  # ホットスタート用

  # 後退時間積分
  for t in range(nt - 2, -1, -1):
    # 熱物性値計算（t時刻の温度場を使用）
    cp, k = thermal_properties_calculator_simple(T_cal[t], cp_coeffs, k_coeffs)

    # 係数とRHS構築
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_Adjoint(
      lambda_field[t+1], T_cal[t+1, :, :, 0], Y_obs[t+1],
      rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    # 行列組み立て（DHCPと同じassemble関数を流用）
    A_csr = assemble_A_DHCP(a_w, a_e, a_s, a_n, a_b, a_t, a_p, ni, nj, nk)

    # 対角前処理
    diag = A_csr.diagonal()
    inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)
    M = LinearOperator(A_csr.shape, matvec=lambda z: inv_diag * z)

    # CG求解
    x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)

    if info > 0:
      print(f"[警告][t={t}] Adjoint CG未収束 (info={info})")
    elif info < 0:
      print(f"[エラー][t={t}] Adjoint CG入力不正 (info={info})")

    # 随伴場更新
    lambda_field[t] = x.reshape((ni, nj, nk), order='F')
    x0 = x  # ホットスタート

  return lambda_field


def global_CGM_time_reference(T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
                                rho, cp_coeffs, k_coeffs, CGM_iteration=20):
  """
  CGM最適化（参照データ生成用、簡略版）

  Pythonオリジナル: 1371-1553行

  この関数は、オリジナルのglobal_CGM_time()から以下を簡略化：
  - メモリ監視機能を削除
  - 数値異常検出を削除
  - プリント文を簡略化
  - 反復数を制限（CGM_iteration=20）

  Args:
    T_init: 初期温度場 (ni, nj, nk)
    Y_obs: 観測温度（底面） (nt, ni, nj)
    q_init: 初期熱流束推定 (nt-1, ni, nj)
    dx, dy: x, y方向格子幅 [m]
    dz: z方向格子幅配列 (nk,) [m]
    dz_b: 下側界面距離 (nk,) [m]
    dz_t: 上側界面距離 (nk,) [m]
    dt: 時間刻み [s]
    rho: 密度 [kg/m³]
    cp_coeffs: 比熱多項式係数
    k_coeffs: 熱伝導率多項式係数
    CGM_iteration: CGM最大反復数

  Returns:
    q: 最終逆解析熱流束 (nt-1, ni, nj)
    T_cal: 最終温度場 (nt, ni, nj, nk)
    J_hist: 目的関数履歴 (list)
    history: 詳細履歴（デバッグ用）
  """
  nt = Y_obs.shape[0]
  ni, nj, nk = T_init.shape
  q = q_init.copy()

  J_hist = []
  history = {
    'q_hist': [],
    'gradient_hist': [],
    'beta_hist': [],
    'gamma_hist': [],
    'search_direction_hist': []
  }

  M = ni * nj
  sigma = 1.8
  epsilon = M * (sigma ** 2) * (nt - 1)

  grad = np.zeros_like(q)
  grad_last = np.zeros_like(q)

  bottom_idx, top_idx = 0, -1
  eps = 1e-12
  dire_reset_every = 5

  p_n_last = np.zeros_like(q)

  # プラトー検出パラメータ
  P = 10
  eta = 1e-4
  min_iter = 10

  def _dot(a, b):
    return float(np.tensordot(a, b, axes=a.ndim))

  for it in range(CGM_iteration):
    print(f"\n=== CGM反復 {it} ===")

    # Step 1: 直接問題求解
    T_cal = multiple_time_step_solver_DHCP(
      T_init, q, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol=1e-6, maxiter=20000
    )

    # Step 2: 停止判定
    res_T = T_cal[1:, :, :, bottom_idx] - Y_obs[1:]
    delta_T = np.abs(res_T)
    J = float(np.tensordot(res_T, res_T, axes=res_T.ndim))
    J_hist.append(J)

    print(f"J = {J:.5e}")

    # Discrepancy停止判定
    if it >= min_iter and J < epsilon and delta_T.max() <= sigma:
      print(f"[停止] Discrepancy条件満足: J={J:.4e} < {epsilon:.4e}")
      break

    # プラトー停止判定
    rel_drop_avg = None
    if len(J_hist) >= P + 1:
      drops = []
      for i in range(-P, 0):
        prev_i, cur_i = J_hist[i-1], J_hist[i]
        drops.append(max(0.0, (prev_i - cur_i) / (abs(prev_i) + eps)))
      rel_drop_avg = sum(drops) / P

    if it >= min_iter and rel_drop_avg is not None and rel_drop_avg < eta:
      print(f"[停止] プラトー: rel_drop_avg={rel_drop_avg:.3e} < {eta:.1e}")
      break

    # Step 3: 随伴問題求解
    lambda_field = multiple_time_step_solver_Adjoint(
      T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol=1e-8, maxiter=20000
    )

    # Step 4: 勾配計算
    for n in range(nt - 1):
      grad[n] = lambda_field[n][:, :, top_idx]

    # Step 5: 共役勾配方向計算（ポラック・リビエール）
    if it == 0 or _dot(grad, p_n_last) <= 0 or it % dire_reset_every == 0:
      p_n = grad.copy()
      gamma = 0
    else:
      y = grad - grad_last
      denom = _dot(grad_last, grad_last) + eps
      gamma = max(0, _dot(grad, y) / denom)
      p_n_candidate = grad + gamma * p_n_last

      if _dot(grad, p_n_candidate) > 0:
        p_n = p_n_candidate
      else:
        p_n = grad.copy()
        gamma = 0

    p_n_last = p_n.copy()

    # Step 6: 感度問題求解（dT計算）
    dT_init = np.zeros_like(T_init)
    dT = multiple_time_step_solver_DHCP(
      dT_init, p_n, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol=1e-8, maxiter=20000
    )

    # Step 7: ステップサイズ計算
    Sp = dT[1:, :, :, bottom_idx]
    numerator = float(np.tensordot(res_T, Sp, axes=res_T.ndim))
    denominator = float(np.tensordot(Sp, Sp, axes=Sp.ndim))
    beta = numerator / (denominator + eps)

    # ステップサイズ制限
    beta_max = 1e8
    if it == 0 and abs(beta) > beta_max:
      print(f"[警告] beta制限: {beta:.2e} => {np.sign(beta)*beta_max:.2e}")
      beta = np.clip(beta, -beta_max, beta_max)

    print(f"beta = {beta:.4e}, gamma = {gamma:.4e}")

    # Step 8: 熱流束更新
    q = q - beta * p_n

    grad_last = grad.copy()

    # 履歴保存
    history['q_hist'].append(q.copy())
    history['gradient_hist'].append(grad.copy())
    history['beta_hist'].append(beta)
    history['gamma_hist'].append(gamma)
    history['search_direction_hist'].append(p_n.copy())

  return q, T_cal, J_hist, history


# ===== テストケース生成 =====

def generate_cgm_reference_1D():
  """
  1次元小規模CGM逆問題（3格子点 × 5時間ステップ）

  - 真の熱流束: 正弦波形（時間変化）
  - 観測ノイズ: σ=1.0K
  - CGM反復: 最大20回
  """
  print("\n" + "="*60)
  print("1次元小規模CGM逆問題（参照データ生成）")
  print("="*60)

  # 格子設定
  ni, nj, nk = 1, 1, 3  # 1×1×3（深さ方向のみ）
  nt = 6  # 時間ステップ数（0, 1, 2, 3, 4, 5）

  # 物理パラメータ
  dx = 0.001  # 1mm
  dy = 0.001  # 1mm
  dz = np.array([0.001, 0.001, 0.001])  # 各層1mm
  dz_b = np.array([0.0005, 0.001, 0.001])  # 下界面距離
  dz_t = np.array([0.001, 0.001, 0.0005])  # 上界面距離
  dt = 0.01  # 10ms

  rho = 7900.0  # SUS304密度 [kg/m³]

  # 簡略化熱物性値（定数）
  T_ref = 300.0  # 参照温度 [K]
  cp_coeffs = [0.0, 0.0, 0.0, 500.0]  # 定数比熱 500 J/(kg·K)
  k_coeffs = [0.0, 0.0, 0.0, 15.0]    # 定数熱伝導率 15 W/(m·K)

  # 初期温度場（一様）
  T_init = np.full((ni, nj, nk), T_ref)

  # 真の熱流束（正弦波）
  q_true = np.zeros((nt - 1, ni, nj))
  for t in range(nt - 1):
    q_true[t, 0, 0] = 5000.0 * np.sin(np.pi * t / (nt - 1))  # 0 ~ 5000 W/m²

  print("\n真の熱流束:")
  for t in range(nt - 1):
    print(f"  t={t}: q={q_true[t, 0, 0]:.2f} W/m²")

  # 疑似観測データ生成（ノイズ付加）
  T_cal_true = multiple_time_step_solver_DHCP(
    T_init, q_true, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt
  )

  # 観測データ（底面温度 + ノイズ）
  sigma_noise = 1.0  # ノイズ標準偏差 [K]
  np.random.seed(42)  # 再現性のため
  noise = np.random.normal(0, sigma_noise, (nt, ni, nj))
  Y_obs = T_cal_true[:, :, :, 0] + noise

  print("\n観測温度（底面、ノイズ付加）:")
  for t in range(nt):
    print(f"  t={t}: Y_obs={Y_obs[t, 0, 0]:.4f} K, T_true={T_cal_true[t, 0, 0, 0]:.4f} K")

  # 初期推定値（ゼロ）
  q_init = np.zeros((nt - 1, ni, nj))

  # CGM実行
  q_final, T_cal_final, J_hist, history = global_CGM_time_reference(
    T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
    rho, cp_coeffs, k_coeffs, CGM_iteration=20
  )

  print("\n最終逆解析結果:")
  for t in range(nt - 1):
    print(f"  t={t}: q_final={q_final[t, 0, 0]:.2f} W/m², q_true={q_true[t, 0, 0]:.2f} W/m²")

  print(f"\nCGM反復数: {len(J_hist)}")
  print(f"最終目的関数: J={J_hist[-1]:.4e}")

  # 参照データ保存
  reference_data = {
    # 問題設定
    'problem': {
      'ni': int(ni), 'nj': int(nj), 'nk': int(nk), 'nt': int(nt),
      'dx': float(dx), 'dy': float(dy),
      'dz': dz.tolist(), 'dz_b': dz_b.tolist(), 'dz_t': dz_t.tolist(), 'dt': float(dt),
      'rho': float(rho),
      'cp_coeffs': cp_coeffs, 'k_coeffs': k_coeffs,
      'sigma_noise': float(sigma_noise)
    },

    # 入力データ
    'input': {
      'T_init': T_init.tolist(),
      'Y_obs': Y_obs.tolist(),
      'q_init': q_init.tolist(),
      'q_true': q_true.tolist()
    },

    # CGM結果
    'output': {
      'q_final': q_final.tolist(),
      'T_cal_final': T_cal_final.tolist(),
      'J_hist': J_hist,
      'n_iterations': int(len(J_hist))
    },

    # 詳細履歴（デバッグ用）
    'history': {
      'beta_hist': history['beta_hist'],
      'gamma_hist': history['gamma_hist'],
      # q_histは大きすぎるので最初と最後のみ保存
      'q_hist_first': history['q_hist'][0].tolist() if history['q_hist'] else [],
      'q_hist_last': history['q_hist'][-1].tolist() if history['q_hist'] else [],
      'gradient_hist_first': history['gradient_hist'][0].tolist() if history['gradient_hist'] else [],
      'gradient_hist_last': history['gradient_hist'][-1].tolist() if history['gradient_hist'] else []
    }
  }

  output_file = "phase4_reference_cgm_1D.json"
  with open(output_file, 'w') as f:
    json.dump(reference_data, f, indent=2)

  print(f"\n参照データを保存: {output_file}")
  print(f"ファイルサイズ: {os.path.getsize(output_file) / 1024:.2f} KB")

  return reference_data


def generate_cgm_reference_3D_small():
  """
  3次元小規模CGM逆問題（3×3×4格子 × 5時間ステップ）

  - 真の熱流束: ガウス分布（空間変化） × ランプ関数（時間変化）
  - 観測ノイズ: σ=1.5K
  - CGM反復: 最大20回
  """
  print("\n" + "="*60)
  print("3次元小規模CGM逆問題（参照データ生成）")
  print("="*60)

  # 格子設定
  ni, nj, nk = 3, 3, 4  # 3×3×4
  nt = 6  # 時間ステップ数

  # 物理パラメータ
  dx = 0.001  # 1mm
  dy = 0.001  # 1mm
  dz = np.array([0.0005, 0.0005, 0.001, 0.001])  # 非均等格子（表面側細かく）
  dz_b = np.array([0.00025, 0.0005, 0.00075, 0.001])
  dz_t = np.array([0.0005, 0.00075, 0.001, 0.0005])
  dt = 0.01  # 10ms

  rho = 7900.0

  # 簡略化熱物性値（定数）
  T_ref = 300.0
  cp_coeffs = [0.0, 0.0, 0.0, 500.0]
  k_coeffs = [0.0, 0.0, 0.0, 15.0]

  # 初期温度場（一様）
  T_init = np.full((ni, nj, nk), T_ref)

  # 真の熱流束（ガウス分布 × ランプ関数）
  q_true = np.zeros((nt - 1, ni, nj))
  x_center, y_center = dx * (ni - 1) / 2, dy * (nj - 1) / 2

  for t in range(nt - 1):
    time_factor = (t + 1) / (nt - 1)  # 0.2, 0.4, 0.6, 0.8, 1.0
    for i in range(ni):
      for j in range(nj):
        x = i * dx
        y = j * dy
        r_sq = (x - x_center)**2 + (y - y_center)**2
        spatial_factor = np.exp(-r_sq / (0.0005**2))  # ガウス分布
        q_true[t, i, j] = 10000.0 * spatial_factor * time_factor

  print("\n真の熱流束（中心点）:")
  for t in range(nt - 1):
    print(f"  t={t}: q_center={q_true[t, 1, 1]:.2f} W/m²")

  # 疑似観測データ生成
  T_cal_true = multiple_time_step_solver_DHCP(
    T_init, q_true, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt
  )

  # 観測データ（ノイズ付加）
  sigma_noise = 1.5
  np.random.seed(123)
  noise = np.random.normal(0, sigma_noise, (nt, ni, nj))
  Y_obs = T_cal_true[:, :, :, 0] + noise

  print("\n観測温度（底面中心、ノイズ付加）:")
  for t in range(nt):
    print(f"  t={t}: Y_obs={Y_obs[t, 1, 1]:.4f} K, T_true={T_cal_true[t, 1, 1, 0]:.4f} K")

  # 初期推定値（一様小値）
  q_init = np.full((nt - 1, ni, nj), 1000.0)  # 1000 W/m²

  # CGM実行
  q_final, T_cal_final, J_hist, history = global_CGM_time_reference(
    T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
    rho, cp_coeffs, k_coeffs, CGM_iteration=20
  )

  print("\n最終逆解析結果（中心点）:")
  for t in range(nt - 1):
    print(f"  t={t}: q_final={q_final[t, 1, 1]:.2f} W/m², q_true={q_true[t, 1, 1]:.2f} W/m²")

  print(f"\nCGM反復数: {len(J_hist)}")
  print(f"最終目的関数: J={J_hist[-1]:.4e}")

  # 参照データ保存
  reference_data = {
    'problem': {
      'ni': int(ni), 'nj': int(nj), 'nk': int(nk), 'nt': int(nt),
      'dx': float(dx), 'dy': float(dy),
      'dz': dz.tolist(), 'dz_b': dz_b.tolist(), 'dz_t': dz_t.tolist(), 'dt': float(dt),
      'rho': float(rho),
      'cp_coeffs': cp_coeffs, 'k_coeffs': k_coeffs,
      'sigma_noise': float(sigma_noise)
    },

    'input': {
      'T_init': T_init.tolist(),
      'Y_obs': Y_obs.tolist(),
      'q_init': q_init.tolist(),
      'q_true': q_true.tolist()
    },

    'output': {
      'q_final': q_final.tolist(),
      'T_cal_final': T_cal_final.tolist(),
      'J_hist': J_hist,
      'n_iterations': int(len(J_hist))
    },

    'history': {
      'beta_hist': history['beta_hist'],
      'gamma_hist': history['gamma_hist'],
      'q_hist_first': history['q_hist'][0].tolist() if history['q_hist'] else [],
      'q_hist_last': history['q_hist'][-1].tolist() if history['q_hist'] else []
    }
  }

  output_file = "phase4_reference_cgm_3D.json"
  with open(output_file, 'w') as f:
    json.dump(reference_data, f, indent=2)

  print(f"\n参照データを保存: {output_file}")
  print(f"ファイルサイズ: {os.path.getsize(output_file) / 1024:.2f} KB")

  return reference_data


if __name__ == "__main__":
  print("="*60)
  print("Phase 4: CGM共役勾配法 参照データ生成")
  print("="*60)

  # 1次元小規模問題
  ref_1d = generate_cgm_reference_1D()

  # 3次元小規模問題
  ref_3d = generate_cgm_reference_3D_small()

  print("\n" + "="*60)
  print("参照データ生成完了")
  print("="*60)
  print("\n生成ファイル:")
  print("  - phase4_reference_cgm_1D.json")
  print("  - phase4_reference_cgm_3D.json")
  print("\nこれらのファイルをJuliaテストで使用してください。")
