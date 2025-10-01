#!/usr/bin/env python3
"""
Phase 5参照データ生成スクリプト（スライディングウィンドウCGM計算）

このスクリプトは、Julia実装の検証用にスライディングウィンドウ計算の参照データを生成します:
1. 真の熱流束を設定（既知値）
2. DHCP直接解を計算（疑似観測データ生成）
3. 観測データにノイズを付加
4. スライディングウィンドウでCGM最適化を実行
5. ウィンドウ間のオーバーラップ平均化処理
6. 全ウィンドウの結果を連結

生成データ:
- 真の熱流束（q_true）
- 疑似観測データ（Y_obs、ノイズ付加）
- 各ウィンドウのCGM結果
- オーバーラップ平均化後の最終結果
- ウィンドウ分割情報

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- sliding_window_CGM_q_saving() (1556-1626行)
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg, LinearOperator
import json
import os

# ===== Phase 4のCGM関数を再利用 =====
# （Phase 4生成スクリプトから必要な関数をコピー）

def polyval_numba_simple(coeffs, x):
  """多項式評価"""
  result = 0.0
  for i in range(len(coeffs)):
    result += coeffs[i] * x ** (len(coeffs) - i - 1)
  return result


def thermal_properties_calculator_simple(T, cp_coeffs, k_coeffs):
  """簡易版熱物性値計算"""
  ni, nj, nk = T.shape
  cp = np.zeros_like(T)
  k = np.zeros_like(T)

  for i in range(ni):
    for j in range(nj):
      for k_idx in range(nk):
        T_val = T[i, j, k_idx]
        cp[i, j, k_idx] = polyval_numba_simple(cp_coeffs, T_val)
        k[i, j, k_idx] = polyval_numba_simple(k_coeffs, T_val)

  return cp, k


def coeffs_and_rhs_building_DHCP(T_initial, q_surface, rho, cp, k_field, dx, dy, dz, dz_b, dz_t, dt):
  """直接熱伝導問題（DHCP）の係数とRHS構築"""
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
        b[p] = rho_cp_V_dt * T_initial[i, j, kk]

        if kk == nk - 1:
          b[p] += q_surface[i, j] * dx * dy

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def assemble_A_DHCP(a_w, a_e, a_s, a_n, a_b, a_t, a_p, ni, nj, nk):
  """DHCP係数行列をCSR形式で組み立て"""
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

        row.append(p)
        col.append(p)
        data.append(a_p[p])

        if i > 0:
          row.append(p)
          col.append(idx(i-1, j, kk))
          data.append(-a_w[p])

        if i < ni - 1:
          row.append(p)
          col.append(idx(i+1, j, kk))
          data.append(-a_e[p])

        if j > 0:
          row.append(p)
          col.append(idx(i, j-1, kk))
          data.append(-a_s[p])

        if j < nj - 1:
          row.append(p)
          col.append(idx(i, j+1, kk))
          data.append(-a_n[p])

        if kk > 0:
          row.append(p)
          col.append(idx(i, j, kk-1))
          data.append(-a_b[p])

        if kk < nk - 1:
          row.append(p)
          col.append(idx(i, j, kk+1))
          data.append(-a_t[p])

  A_csr = sp.csr_matrix((data, (row, col)), shape=(N, N))
  return A_csr


def multiple_time_step_solver_DHCP(T_init, q, nt, rho, cp_coeffs, k_coeffs,
                                    dx, dy, dz, dz_b, dz_t, dt, rtol=1e-6, maxiter=20000):
  """複数時間ステップDHCP求解"""
  ni, nj, nk = T_init.shape
  T_all = np.zeros((nt, ni, nj, nk))
  T_all[0] = T_init.copy()

  x0 = None

  for t in range(1, nt):
    cp, k = thermal_properties_calculator_simple(T_all[t-1], cp_coeffs, k_coeffs)

    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
      T_all[t-1], q[t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    A_csr = assemble_A_DHCP(a_w, a_e, a_s, a_n, a_b, a_t, a_p, ni, nj, nk)

    diag = A_csr.diagonal()
    inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)
    M = LinearOperator(A_csr.shape, matvec=lambda z: inv_diag * z)

    x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)

    if info != 0:
      print(f"[警告][t={t}] DHCP CG info={info}")

    T_all[t] = x.reshape((ni, nj, nk), order='F')
    x0 = x

  return T_all


def coeffs_and_rhs_building_Adjoint(lambda_initial, T_cal_bottom, Y_obs,
                                      rho, cp, k_field, dx, dy, dz, dz_b, dz_t, dt):
  """随伴方程式の係数とRHS構築"""
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
        b[p] = rho_cp_V_dt * lambda_initial[i, j, kk]

        if kk == 0:
          residual = T_cal_bottom[i, j] - Y_obs[i, j]
          b[p] += 2.0 * residual * dx * dy

  return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def multiple_time_step_solver_Adjoint(T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
                                        dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=20000):
  """複数時間ステップAdjoint求解（後退時間積分）"""
  ni, nj, nk = T_cal.shape[1:]
  lambda_field = np.zeros((nt, ni, nj, nk))

  lambda_field[-1] = 0.0

  x0 = None

  for t in range(nt - 2, -1, -1):
    cp, k = thermal_properties_calculator_simple(T_cal[t], cp_coeffs, k_coeffs)

    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_Adjoint(
      lambda_field[t+1], T_cal[t+1, :, :, 0], Y_obs[t+1],
      rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    A_csr = assemble_A_DHCP(a_w, a_e, a_s, a_n, a_b, a_t, a_p, ni, nj, nk)

    diag = A_csr.diagonal()
    inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)
    M = LinearOperator(A_csr.shape, matvec=lambda z: inv_diag * z)

    x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)

    if info != 0:
      print(f"[警告][t={t}] Adjoint CG info={info}")

    lambda_field[t] = x.reshape((ni, nj, nk), order='F')
    x0 = x

  return lambda_field


def global_CGM_time_reference(T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
                                rho, cp_coeffs, k_coeffs, CGM_iteration=20):
  """
  CGM最適化（参照データ生成用、簡略版）

  Phase 4と同じだが、収束判定を簡略化してテスト用に固定反復数で実行
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

    # Discrepancy停止判定
    if it >= min_iter and J < epsilon and delta_T.max() <= sigma:
      print(f"[停止] Discrepancy: J={J:.4e} < {epsilon:.4e}")
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
      beta = np.clip(beta, -beta_max, beta_max)

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


def sliding_window_CGM_q_saving_reference(
    Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, CGM_iteration=20
):
  """
  スライディングウィンドウCGM計算（参照データ生成用）

  Pythonオリジナル: 1556-1626行

  Args:
    Y_obs: 全時間の観測温度 (nt+1, ni, nj)
    T0: 初期温度場 (ni, nj, nk)
    dx, dy: x, y方向格子幅 [m]
    dz: z方向格子幅配列 (nk,) [m]
    dz_b: 下側界面距離 (nk,) [m]
    dz_t: 上側界面距離 (nk,) [m]
    dt: 時間刻み [s]
    rho: 密度 [kg/m³]
    cp_coeffs: 比熱多項式係数
    k_coeffs: 熱伝導率多項式係数
    window_size: ウィンドウサイズ（時間ステップ数）
    overlap: オーバーラップ領域（時間ステップ数）
    q_init_value: 初期熱流束の初期値
    CGM_iteration: CGM最大反復数

  Returns:
    q_global: 全時間の逆解析熱流束 (nt-1, ni, nj)
    windows_info: 各ウィンドウの詳細情報（デバッグ用）
  """
  nt = Y_obs.shape[0]
  T_init = T0.copy()
  ni, nj, nk = T_init.shape

  start_idx = 0
  q_total = []
  prev_q_win = None

  windows_info = []

  safety_counter = 0
  safety_limit = nt * 5

  print(f"\n=== スライディングウィンドウCGM計算開始 ===")
  print(f"全時間ステップ数: {nt}")
  print(f"ウィンドウサイズ: {window_size}")
  print(f"オーバーラップ: {overlap}")

  while start_idx < nt - 1:
    safety_counter += 1
    if safety_counter > safety_limit:
      print("[警告] 安全カウンタ超過")
      break

    # 当前可用窓長
    max_L = min(window_size, (nt - 1) - start_idx)
    end_idx = start_idx + max_L
    Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]

    print(f"\n--- ウィンドウ {len(windows_info)+1}: [{start_idx}, {end_idx}] (長さ={max_L}) ---")

    # 該窗的初始热流
    if prev_q_win is None:
      q_init_win = np.full((max_L, ni, nj), q_init_value, dtype=float)
      print(f"初期熱流束: 一定値 {q_init_value} W/m²")
    else:
      q_init_win = np.empty((max_L, ni, nj), dtype=float)
      L_overlap = min(overlap, max_L, prev_q_win.shape[0])
      if L_overlap > 0:
        q_init_win[:L_overlap] = prev_q_win[-L_overlap:]
        print(f"オーバーラップ継承: 前{L_overlap}ステップ")
      if L_overlap < max_L:
        edge = prev_q_win[-1]
        q_init_win[L_overlap:] = edge
        print(f"残り部分: 前ウィンドウ最終値で埋める")

    # CGM実行
    q_win, T_win_last, J_hist, history = global_CGM_time_reference(
      T_init, Y_obs_win, q_init_win, dx, dy, dz, dz_b, dz_t, dt,
      rho, cp_coeffs, k_coeffs, CGM_iteration=CGM_iteration
    )

    prev_q_win = q_win.copy()

    # 拼接 q（対重叠部分做平均）
    overlap_steps_actual = 0
    if len(q_total) == 0:
      q_total.append(q_win)
      print(f"第1ウィンドウ: そのまま追加")
    else:
      overlap_steps = min(overlap, q_win.shape[0], q_total[-1].shape[0])
      overlap_steps_actual = overlap_steps  # 記録用
      if overlap_steps > 0:
        # オーバーラップ平均化（Python原典1609行に合わせる: 0.5*old + new）
        q_total[-1][-overlap_steps:] = 0.5 * q_total[-1][-overlap_steps:] + q_win[:overlap_steps]
        remaining = q_win[overlap_steps:]
        if remaining.shape[0] > 0:
          q_total.append(remaining)
        print(f"オーバーラップ平均化: {overlap_steps}ステップ")
      else:
        q_total.append(q_win)
        print(f"オーバーラップなし: そのまま追加")

    # 温度場継承
    T_init = T_win_last.copy() if T_win_last.ndim == 3 else T_win_last[-1].copy()

    # ウィンドウ情報保存
    window_data = {
      'window_id': len(windows_info) + 1,
      'start_idx': int(start_idx),
      'end_idx': int(end_idx),
      'max_L': int(max_L),
      'n_iterations': int(len(J_hist)),
      'J_final': float(J_hist[-1]) if J_hist else 0.0,
      'q_win_shape': q_win.shape,
      'overlap_steps': int(overlap_steps_actual)
    }
    windows_info.append(window_data)

    print(f"CGM反復数: {len(J_hist)}, 最終J: {J_hist[-1]:.4e}")

    step = max(1, max_L - overlap)
    start_idx += step

  # 拼接为全局 q，并裁剪到 nt-1
  q_global = np.concatenate(q_total, axis=0)[:nt-1]

  print(f"\n=== スライディングウィンドウ計算完了 ===")
  print(f"総ウィンドウ数: {len(windows_info)}")
  print(f"最終q_global形状: {q_global.shape}")

  return q_global, windows_info


# ===== テストケース生成 =====

def generate_sliding_window_reference_1D():
  """
  1次元小規模スライディングウィンドウ（3格子点 × 15時間ステップ、2ウィンドウ）

  - 真の熱流束: ランプ関数（徐々に増加）
  - ウィンドウサイズ: 10
  - オーバーラップ: 3
  - CGM反復: 最大10回（高速化）
  """
  print("\n" + "="*60)
  print("1次元スライディングウィンドウ（参照データ生成）")
  print("="*60)

  # 格子設定
  ni, nj, nk = 1, 1, 3
  nt = 16  # 0~15（熱流束は0~14）

  # 物理パラメータ
  dx = 0.001
  dy = 0.001
  dz = np.array([0.001, 0.001, 0.001])
  dz_b = np.array([0.0005, 0.001, 0.001])
  dz_t = np.array([0.001, 0.001, 0.0005])
  dt = 0.01

  rho = 7900.0
  cp_coeffs = [0.0, 0.0, 0.0, 500.0]
  k_coeffs = [0.0, 0.0, 0.0, 15.0]

  T_ref = 300.0
  T_init = np.full((ni, nj, nk), T_ref)

  # 真の熱流束（ランプ関数）
  q_true = np.zeros((nt - 1, ni, nj))
  for t in range(nt - 1):
    q_true[t, 0, 0] = 1000.0 + 200.0 * t  # 1000 ~ 3800 W/m²

  print("\n真の熱流束:")
  for t in range(0, nt-1, 3):
    print(f"  t={t}: q={q_true[t, 0, 0]:.2f} W/m²")

  # 疑似観測データ生成
  T_cal_true = multiple_time_step_solver_DHCP(
    T_init, q_true, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt
  )

  sigma_noise = 1.0
  np.random.seed(42)
  noise = np.random.normal(0, sigma_noise, (nt, ni, nj))
  Y_obs = T_cal_true[:, :, :, 0] + noise

  # スライディングウィンドウ計算
  window_size = 10
  overlap = 3
  q_init_value = 0.0

  q_global, windows_info = sliding_window_CGM_q_saving_reference(
    Y_obs, T_init, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, CGM_iteration=10
  )

  print("\n最終逆解析結果（サンプリング）:")
  for t in range(0, nt-1, 3):
    print(f"  t={t}: q_final={q_global[t, 0, 0]:.2f} W/m², q_true={q_true[t, 0, 0]:.2f} W/m²")

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

    'sliding_window': {
      'window_size': int(window_size),
      'overlap': int(overlap),
      'q_init_value': float(q_init_value),
      'CGM_iteration': int(10)
    },

    'input': {
      'T_init': T_init.tolist(),
      'Y_obs': Y_obs.tolist(),
      'q_true': q_true.tolist()
    },

    'output': {
      'q_global': q_global.tolist(),
      'n_windows': int(len(windows_info)),
      'windows_info': windows_info
    }
  }

  output_file = "phase5_reference_sliding_window_1D.json"
  with open(output_file, 'w') as f:
    json.dump(reference_data, f, indent=2)

  print(f"\n参照データを保存: {output_file}")
  print(f"ファイルサイズ: {os.path.getsize(output_file) / 1024:.2f} KB")

  return reference_data


def generate_sliding_window_reference_2D_small():
  """
  2次元小規模スライディングウィンドウ（2×2×3格子 × 13時間ステップ、3ウィンドウ）

  - 真の熱流束: 空間変化（中心高、周辺低） × 時間変化（ステップ関数）
  - ウィンドウサイズ: 6
  - オーバーラップ: 2
  - CGM反復: 最大10回
  """
  print("\n" + "="*60)
  print("2次元スライディングウィンドウ（参照データ生成）")
  print("="*60)

  # 格子設定
  ni, nj, nk = 2, 2, 3
  nt = 14  # 0~13（熱流束は0~12）

  # 物理パラメータ
  dx = 0.001
  dy = 0.001
  dz = np.array([0.0008, 0.001, 0.001])
  dz_b = np.array([0.0004, 0.0009, 0.001])
  dz_t = np.array([0.0009, 0.001, 0.0005])
  dt = 0.01

  rho = 7900.0
  cp_coeffs = [0.0, 0.0, 0.0, 500.0]
  k_coeffs = [0.0, 0.0, 0.0, 15.0]

  T_ref = 300.0
  T_init = np.full((ni, nj, nk), T_ref)

  # 真の熱流束（空間変化 × 時間変化）
  q_true = np.zeros((nt - 1, ni, nj))
  for t in range(nt - 1):
    # 時間変化: 0~5は1000, 6~12は3000 W/m²
    base_q = 1000.0 if t < 6 else 3000.0
    for i in range(ni):
      for j in range(nj):
        # 空間変化: 中心 (0.5, 0.5) から距離に応じて減衰
        x = i * dx
        y = j * dy
        x_center = 0.0005
        y_center = 0.0005
        r = np.sqrt((x - x_center)**2 + (y - y_center)**2)
        spatial_factor = np.exp(-r / 0.0005)
        q_true[t, i, j] = base_q * spatial_factor

  print("\n真の熱流束（中心点）:")
  for t in range(0, nt-1, 3):
    print(f"  t={t}: q_center={q_true[t, 0, 0]:.2f} W/m²")

  # 疑似観測データ生成
  T_cal_true = multiple_time_step_solver_DHCP(
    T_init, q_true, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt
  )

  sigma_noise = 1.2
  np.random.seed(123)
  noise = np.random.normal(0, sigma_noise, (nt, ni, nj))
  Y_obs = T_cal_true[:, :, :, 0] + noise

  # スライディングウィンドウ計算
  window_size = 6
  overlap = 2
  q_init_value = 500.0

  q_global, windows_info = sliding_window_CGM_q_saving_reference(
    Y_obs, T_init, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, CGM_iteration=10
  )

  print("\n最終逆解析結果（中心点、サンプリング）:")
  for t in range(0, nt-1, 3):
    print(f"  t={t}: q_final={q_global[t, 0, 0]:.2f} W/m², q_true={q_true[t, 0, 0]:.2f} W/m²")

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

    'sliding_window': {
      'window_size': int(window_size),
      'overlap': int(overlap),
      'q_init_value': float(q_init_value),
      'CGM_iteration': int(10)
    },

    'input': {
      'T_init': T_init.tolist(),
      'Y_obs': Y_obs.tolist(),
      'q_true': q_true.tolist()
    },

    'output': {
      'q_global': q_global.tolist(),
      'n_windows': int(len(windows_info)),
      'windows_info': windows_info
    }
  }

  output_file = "phase5_reference_sliding_window_2D.json"
  with open(output_file, 'w') as f:
    json.dump(reference_data, f, indent=2)

  print(f"\n参照データを保存: {output_file}")
  print(f"ファイルサイズ: {os.path.getsize(output_file) / 1024:.2f} KB")

  return reference_data


if __name__ == "__main__":
  print("="*60)
  print("Phase 5: スライディングウィンドウCGM 参照データ生成")
  print("="*60)

  # 1次元小規模問題
  ref_1d = generate_sliding_window_reference_1D()

  # 2次元小規模問題
  ref_2d = generate_sliding_window_reference_2D_small()

  print("\n" + "="*60)
  print("参照データ生成完了")
  print("="*60)
  print("\n生成ファイル:")
  print("  - phase5_reference_sliding_window_1D.json")
  print("  - phase5_reference_sliding_window_2D.json")
  print("\nこれらのファイルをJuliaテストで使用してください。")
