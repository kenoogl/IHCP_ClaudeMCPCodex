#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
診断機能付きソルバー (Solvers with Diagnostics)

オリジナルのIHCP_CGM_Sliding_Window_Calculation_ver2.pyのソルバーを拡張し、
詳細な診断情報（反復回数、残差、実行時間）を返すバージョン。

Python-Julia比較検証用。
"""

import numpy as np
import time
from scipy.sparse.linalg import LinearOperator, cg
from scipy.sparse import diags


class SolverDiagnostics:
    """ソルバー診断情報を格納するクラス"""

    def __init__(self):
        self.iterations = []  # 各時間ステップの反復回数
        self.residuals = []   # 各時間ステップの最終残差
        self.initial_residuals = []  # 各時間ステップの初期残差
        self.elapsed_time = 0.0  # 総実行時間
        self.convergence_flags = []  # 収束フラグ（info値）

    def add_timestep(self, iters, init_resid, final_resid, info):
        """1時間ステップの診断情報を追加"""
        self.iterations.append(iters)
        self.initial_residuals.append(init_resid)
        self.residuals.append(final_resid)
        self.convergence_flags.append(info)

    def summary(self):
        """診断情報のサマリーを返す"""
        if not self.iterations:
            return {
                'total_iters': 0,
                'avg_iters': 0.0,
                'max_iters': 0,
                'min_iters': 0,
                'elapsed_time': self.elapsed_time,
                'convergence_rate': 0.0
            }

        total_iters = sum(self.iterations)
        avg_iters = total_iters / len(self.iterations)
        max_iters = max(self.iterations)
        min_iters = min(self.iterations)
        converged_count = sum(1 for info in self.convergence_flags if info == 0)
        convergence_rate = converged_count / len(self.convergence_flags) if self.convergence_flags else 0.0

        return {
            'total_iters': total_iters,
            'avg_iters': avg_iters,
            'max_iters': max_iters,
            'min_iters': min_iters,
            'elapsed_time': self.elapsed_time,
            'convergence_rate': convergence_rate,
            'timesteps': len(self.iterations)
        }


def cg_with_counter(A, b, M=None, x0=None, rtol=1e-5, maxiter=None):
    """
    反復回数と残差をカウントする共役勾配法ラッパー

    Returns:
        x: 解ベクトル
        info: 収束情報
        iters: 実際の反復回数
        init_resid: 初期残差ノルム
        final_resid: 最終残差ノルム
    """
    counter = {'iters': 0, 'residuals': []}

    def callback(xk):
        """各反復で呼ばれるコールバック"""
        counter['iters'] += 1
        # 残差計算: r = b - A*xk
        r = b - A @ xk
        resid_norm = np.linalg.norm(r)
        counter['residuals'].append(resid_norm)

    # 初期残差計算
    if x0 is not None:
        r0 = b - A @ x0
    else:
        r0 = b.copy()
    init_resid = np.linalg.norm(r0)
    counter['residuals'].append(init_resid)

    # scipy.cg実行
    x, info = cg(A, b, M=M, x0=x0, rtol=rtol, maxiter=maxiter, callback=callback)

    # 最終残差
    final_resid = counter['residuals'][-1] if counter['residuals'] else init_resid

    return x, info, counter['iters'], init_resid, final_resid


def multiple_time_step_solver_DHCP_with_diagnostics(
    T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol, maxiter,
    thermal_properties_calculator, coeffs_and_rhs_building_DHCP, assemble_A_DHCP,
    verbose=False
):
    """
    診断機能付きDHCPソルバー

    オリジナルのmultiple_time_step_solver_DHCPと同じ機能に加えて、
    各時間ステップの反復回数・残差情報を返す。

    Args:
        (オリジナルと同じ)
        thermal_properties_calculator: 熱物性値計算関数
        coeffs_and_rhs_building_DHCP: 係数行列構築関数
        assemble_A_DHCP: 行列組み立て関数

    Returns:
        T_all: 温度場 (nt, ni, nj, nk)
        diagnostics: SolverDiagnosticsオブジェクト
    """
    ni, nj, nk = T_initial.shape
    T_all = np.empty((nt, ni, nj, nk))
    T_all[0] = T_initial

    x0 = T_initial.ravel(order="F").copy()

    diagnostics = SolverDiagnostics()
    start_time = time.time()

    for t in range(1, nt):
        cp, k = thermal_properties_calculator(T_all[t-1], cp_coeffs, k_coeffs)

        a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
            T_all[t-1], q_surface[t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
        )

        A_csr = assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

        # Jacobi前処理
        diag = A_csr.diagonal()
        inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

        def Mv(z):
            return inv_diag * z

        M = LinearOperator(A_csr.shape, matvec=Mv)

        # 反復回数カウント付きCG
        x, info, iters, init_resid, final_resid = cg_with_counter(
            A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0
        )

        # 診断情報追加
        diagnostics.add_timestep(iters, init_resid, final_resid, info)

        # Julia版と同じ形式で詳細ログ出力
        if verbose:
            print(f"[t={t}/{nt-1}] converged: Iteration= {iters:2d} : Res_0= {init_resid:.15e} : time={time.time()-start_time:.4f}")

        if info != 0:
            if info > 0:
                print(f"[DHCP警告][t={t}] 収束せず (info={info}, iters={iters})")
            else:
                print(f"[DHCPエラー][t={t}] 入力不正 (info={info})")

        T_all[t] = x.reshape((ni, nj, nk), order="F")
        x0 = x.copy()

    diagnostics.elapsed_time = time.time() - start_time

    return T_all, diagnostics


def multiple_time_step_solver_Adjoint_with_diagnostics(
    T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt, rtol, maxiter,
    thermal_properties_calculator, coeffs_and_rhs_building_Adjoint, assemble_A_Adjoint
):
    """
    診断機能付きAdjointソルバー

    Args:
        (オリジナルと同じ)
        thermal_properties_calculator: 熱物性値計算関数
        coeffs_and_rhs_building_Adjoint: 係数行列構築関数（随伴）
        assemble_A_Adjoint: 行列組み立て関数（随伴）

    Returns:
        lambda_all: 随伴場 (nt, ni, nj, nk) ※最後はゼロ
        diagnostics: SolverDiagnosticsオブジェクト
    """
    ni, nj, nk = T_cal.shape[1:4]
    lambda_all = np.zeros((nt, ni, nj, nk))  # nt個確保、最後はゼロ初期化
    lambda_all[-1] = 0  # 最終ステップはゼロ

    diagnostics = SolverDiagnostics()
    start_time = time.time()

    x0 = None

    for n in range(nt - 2, -1, -1):
        cp, k = thermal_properties_calculator(T_cal[n], cp_coeffs, k_coeffs)

        lambda_initial = lambda_all[n + 1]

        a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_Adjoint(
            lambda_initial, T_cal[n][:,:,0], Y_obs[n], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
        )

        A_csr = assemble_A_Adjoint(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

        # Jacobi前処理
        diag = A_csr.diagonal()
        inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

        def Mv(z):
            return inv_diag * z

        M = LinearOperator(A_csr.shape, matvec=Mv)

        # 反復回数カウント付きCG
        x, info, iters, init_resid, final_resid = cg_with_counter(
            A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0
        )

        # 診断情報追加
        diagnostics.add_timestep(iters, init_resid, final_resid, info)

        if info != 0:
            if info > 0:
                print(f"[Adjoint警告][n={n}] 収束せず (info={info}, iters={iters})")
            else:
                print(f"[Adjointエラー][n={n}] 入力不正 (info={info})")

        lambda_all[n] = x.reshape((ni, nj, nk), order="F")
        x0 = x.copy()

    diagnostics.elapsed_time = time.time() - start_time

    return lambda_all, diagnostics


def global_CGM_time_with_diagnostics(
    T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
    rho, cp_coeffs, k_coeffs, CGM_iteration,
    thermal_properties_calculator,
    coeffs_and_rhs_building_DHCP, assemble_A_DHCP,
    coeffs_and_rhs_building_Adjoint, assemble_A_Adjoint,
    verbose=True
):
    """
    診断機能付きCGMソルバー

    オリジナルのglobal_CGM_timeと同じ機能に加えて、
    各CGM反復での詳細診断情報を返す。

    Returns:
        q: 最適化された熱流束 (nt-1, ni, nj)
        T_cal: 最終温度場 (nt, ni, nj, nk)
        J_hist: 目的関数履歴
        cgm_diagnostics: 各反復の診断情報リスト
    """
    nt = Y_obs.shape[0]
    ni, nj, nk = T_init.shape
    q = q_init.copy()

    J_hist = []
    cgm_diagnostics = []  # 各CGM反復の診断情報

    M = ni * nj
    sigma = 1.8
    epsilon = M * (sigma ** 2) * (nt - 1)

    grad = np.zeros_like(q)
    grad_last = np.zeros_like(q)

    bottom_idx, top_idx = 0, -1
    eps = 1e-12
    dire_reset_every = 5

    p_n_last = np.zeros_like(q)

    P = 10
    eta = 1e-4
    min_iter = 10

    def _dot(a, b):
        return float(np.tensordot(a, b, axes=a.ndim))

    for it in range(CGM_iteration):
        iter_diagnostics = {
            'iteration': it,
            'dhcp': None,
            'adjoint': None,
            'J': None,
            'beta': None,
            'gamma': None,
            'wall_time': 0.0
        }

        t0 = time.time()

        # DHCP実行（診断情報付き）
        T_cal, dhcp_diag = multiple_time_step_solver_DHCP_with_diagnostics(
            T_init, q, nt, rho, cp_coeffs, k_coeffs,
            dx, dy, dz, dz_b, dz_t, dt,
            rtol=1e-6, maxiter=20000,
            thermal_properties_calculator=thermal_properties_calculator,
            coeffs_and_rhs_building_DHCP=coeffs_and_rhs_building_DHCP,
            assemble_A_DHCP=assemble_A_DHCP,
            verbose=verbose
        )
        iter_diagnostics['dhcp'] = dhcp_diag.summary()

        # 目的関数計算
        res_T = T_cal[1:, :, :, bottom_idx] - Y_obs[1:]
        delta_T = np.abs(res_T)
        J = float(np.tensordot(res_T, res_T, axes=res_T.ndim))
        J_hist.append(J)
        iter_diagnostics['J'] = J

        # 停止判定
        if it >= min_iter and J < epsilon and delta_T.max() <= sigma:
            if verbose:
                print(f"[停止] Discrepancy条件満足: J={J:.4e} < {epsilon:.4e}")
            cgm_diagnostics.append(iter_diagnostics)
            break

        # プラトー検出
        rel_drop_avg = None
        if len(J_hist) >= P + 1:
            drops = []
            for i in range(-P, 0):
                prev_i, cur_i = J_hist[i-1], J_hist[i]
                drops.append(max(0.0, (prev_i - cur_i) / (abs(prev_i) + eps)))
            rel_drop_avg = sum(drops) / P

        if it >= min_iter and rel_drop_avg is not None and rel_drop_avg < eta:
            if verbose:
                print(f"[停止] プラトー: rel_drop_avg={rel_drop_avg:.3e} < {eta:.1e}")
            cgm_diagnostics.append(iter_diagnostics)
            break

        # Adjoint実行（診断情報付き）
        lambda_field, adjoint_diag = multiple_time_step_solver_Adjoint_with_diagnostics(
            T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
            dx, dy, dz, dz_b, dz_t, dt,
            rtol=1e-8, maxiter=20000,
            thermal_properties_calculator=thermal_properties_calculator,
            coeffs_and_rhs_building_Adjoint=coeffs_and_rhs_building_Adjoint,
            assemble_A_Adjoint=assemble_A_Adjoint
        )
        iter_diagnostics['adjoint'] = adjoint_diag.summary()

        # 勾配計算
        for n in range(nt - 1):
            grad[n] = lambda_field[n][:, :, top_idx]

        # 共役勾配方向
        if it == 0 or _dot(grad, p_n_last) <= 0 or it % dire_reset_every == 0:
            p_n = grad.copy()
            gamma = 0
        else:
            y = grad - grad_last
            denom = _dot(grad_last, grad_last) + eps
            gamma = max(0, _dot(grad, y) / denom)
            p_n = grad + gamma * p_n_last

        iter_diagnostics['gamma'] = gamma

        # 感度解析（簡易版 - オリジナルコードにはないが、完全性のため）
        # 注: オリジナルコードではdTは陽的に計算されていない
        dT = np.zeros((nt, ni, nj, nk))

        # ステップサイズ計算
        p_dot_p = _dot(p_n, p_n)
        dT_dot_dT = _dot(dT[1:, :, :, bottom_idx], dT[1:, :, :, bottom_idx])

        beta_max = 1.0e8
        if dT_dot_dT > eps:
            beta = p_dot_p / dT_dot_dT
            beta = min(beta, beta_max)
        else:
            beta = beta_max

        iter_diagnostics['beta'] = beta

        # 熱流束更新
        q = q - beta * p_n

        grad_last = grad.copy()
        p_n_last = p_n.copy()

        wall_s = time.time() - t0
        iter_diagnostics['wall_time'] = wall_s

        # 詳細出力
        if verbose:
            rel_drop = None
            if len(J_hist) >= 2:
                rel_drop = abs(J_hist[-1] - J_hist[-2]) / (J_hist[-2] + eps)

            print(f"@ ___ Iter {it:3d} ___ @ wall_s = {wall_s:.3f}s")
            print(f"J = {J:.5e}, beta = {beta:.4e}, rel_drop = {rel_drop if rel_drop else 'N/A'}")
            print(f"|T - Y|: max={delta_T.max():.3e}, min={delta_T.min():.3e}, mean={delta_T.mean():.3e}")
            print(f"grad: min={grad.min():.4e}, max={grad.max():.4e}, mean={grad.mean():.4e}")
            print(f"q:    min={q.min():.4e}, max={q.max():.4e}, mean={q.mean():.4e}")
            print(f"DHCP:    {dhcp_diag.summary()['total_iters']:4d} iters, {dhcp_diag.elapsed_time:.2f}s")
            print(f"Adjoint: {adjoint_diag.summary()['total_iters']:4d} iters, {adjoint_diag.elapsed_time:.2f}s")
            print()

        cgm_diagnostics.append(iter_diagnostics)

    return q, T_cal, J_hist, cgm_diagnostics
