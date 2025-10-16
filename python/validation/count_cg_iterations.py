#!/usr/bin/env python3
"""Python CG反復回数カウント版（10ステップベンチマーク）"""

import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import sys
import time
from pathlib import Path

import numpy as np

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]
ORIGINAL_DIR = THIS_FILE.parents[1] / "original"
if str(ORIGINAL_DIR) not in sys.path:
    sys.path.insert(0, str(ORIGINAL_DIR))

# SciPyのcg関数を独自実装でラップ
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import diags
import scipy.sparse.linalg as spla

# グローバルカウンター
cg_iteration_counts = []
current_cg_count = [0]  # リストで参照渡し

def cg_with_counter(A, b, M=None, rtol=1e-6, maxiter=20000, x0=None):
    """反復回数をカウントするCG法ラッパー"""
    counter = [0]  # ローカルカウンター

    def callback(xk):
        counter[0] += 1

    x, info = spla.cg(A, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0, callback=callback)

    return x, info, counter[0]

# オリジナルの実装をインポート
from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    thermal_properties_calculator,
    coeffs_and_rhs_building_DHCP,
    assemble_A_DHCP,
    coeffs_and_rhs_building_Adjoint,
    assemble_A_Adjoint,
    rho,
    cp_coeffs,
    k_coeffs,
    dx,
    dy,
    dz,
    dz_b,
    dz_t,
)

NT = 10
DT = 1.0e-3

# 反復回数記録用
dhcp_iter_counts = []
adjoint_iter_counts = []
sensitivity_iter_counts = []


def multiple_time_step_solver_DHCP_counted(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
                                           dx, dy, dz, dz_b, dz_t, dt, rtol, maxiter):
    """DHCP solver with iteration counting"""
    ni, nj, nk = T_initial.shape
    N = ni * nj * nk

    T_all = np.empty((nt, ni, nj, nk))
    T_all[0] = T_initial

    x0 = T_initial.ravel(order="F").copy()

    step_iters = []

    for t in range(1, nt):
        cp, k = thermal_properties_calculator(T_all[t-1], cp_coeffs, k_coeffs)

        a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
            T_all[t-1], q_surface[t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt)

        A_csr = assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

        diag = A_csr.diagonal()
        inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

        def Mv(z):
            return inv_diag * z

        M = LinearOperator(A_csr.shape, matvec=Mv)

        x, info, iter_count = cg_with_counter(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)
        step_iters.append(iter_count)

        if info > 0:
            print(f"[警告][DHCP t={t}] Krylov法が収束しませんでした (info={info}, iters={iter_count})")
        elif info < 0:
            print(f"[エラー][DHCP t={t}] Krylov法の入力が不正です (info={info})")

        T_all[t] = x.reshape((ni, nj, nk), order="F")
        x0 = x

    return T_all, step_iters


def multiple_time_step_solver_Adjoint_counted(T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
                                               dx, dy, dz, dz_b, dz_t, dt, rtol, maxiter):
    """Adjoint solver with iteration counting"""
    ni, nj, nk = T_cal[0].shape
    N = ni * nj * nk

    lambda_field = np.empty_like(T_cal)
    lambda_field[-1] = 0

    lambda_initial = lambda_field[-1]
    x0 = lambda_initial.ravel(order="F").copy()

    step_iters = []

    for t in range(nt-2, -1, -1):
        lambda_initial = lambda_field[t+1]

        cp, k = thermal_properties_calculator(T_cal[t], cp_coeffs, k_coeffs)

        a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_Adjoint(
            lambda_initial, T_cal[t][:,:,0], Y_obs[t], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt)

        A_csr = assemble_A_Adjoint(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

        diag = A_csr.diagonal()
        inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

        def Mv(z):
            return inv_diag * z

        M = LinearOperator(A_csr.shape, matvec=Mv)

        x, info, iter_count = cg_with_counter(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)
        step_iters.append(iter_count)

        if info > 0:
            print(f"[警告][Adjoint t={t}] Krylov法が収束しませんでした (info={info}, iters={iter_count})")
        elif info < 0:
            print(f"[エラー][Adjoint t={t}] Krylov法の入力が不正です (info={info})")

        lambda_field[t] = x.reshape((ni, nj, nk), order="F")
        x0 = x

    return lambda_field, step_iters


def main() -> int:
    print("="*80)
    print("Python 10-step CGM with CG iteration counting")
    print("="*80)
    print(f"NUMBA_NUM_THREADS={os.environ.get('NUMBA_NUM_THREADS')}")

    data_path = PROJECT_ROOT / "shared" / "data" / "T_measure_700um_1ms.npy"
    print(f"\nLoading measurement data: {data_path}")
    mmapped = np.load(data_path, mmap_mode="r")
    try:
        if mmapped.shape[0] < NT:
            raise ValueError(f"Need {NT} time steps")
        Y_obs = np.array(mmapped[:NT], dtype=np.float64, copy=True)
    finally:
        del mmapped

    # 初期温度場
    T_init = np.repeat(Y_obs[0][:, :, np.newaxis], dz.shape[0], axis=2).astype(np.float64, copy=False)
    ni, nj = Y_obs.shape[1:]
    q_init = np.zeros((NT - 1, ni, nj), dtype=np.float64)

    print(f"\nGrid size: {ni}x{nj}x{dz.shape[0]} = {ni*nj*dz.shape[0]} points")
    print(f"Time steps: {NT}, dt={DT}s")

    print("\n[1/3] Running CGM (DHCP + Adjoint + Sensitivity)...")
    cgm_start = time.perf_counter()

    # Step 1: DHCP
    print("  [DHCP] Running...")
    dhcp_start = time.perf_counter()
    T_cal, dhcp_iters = multiple_time_step_solver_DHCP_counted(
        T_init, q_init, NT, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, DT,
        rtol=1.0e-6, maxiter=20000
    )
    dhcp_elapsed = time.perf_counter() - dhcp_start
    print(f"  [DHCP] Complete: {dhcp_elapsed:.2f}s")
    print(f"  [DHCP] Iterations per step: {dhcp_iters}")
    print(f"  [DHCP] Total iterations: {sum(dhcp_iters)}")
    print(f"  [DHCP] Average: {sum(dhcp_iters)/len(dhcp_iters):.1f} iters/step")

    # Step 2: Adjoint
    print("  [Adjoint] Running...")
    adjoint_start = time.perf_counter()
    lambda_field, adjoint_iters = multiple_time_step_solver_Adjoint_counted(
        T_cal, Y_obs, NT, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, DT,
        rtol=1.0e-8, maxiter=20000
    )
    adjoint_elapsed = time.perf_counter() - adjoint_start
    print(f"  [Adjoint] Complete: {adjoint_elapsed:.2f}s")
    print(f"  [Adjoint] Iterations per step: {adjoint_iters}")
    print(f"  [Adjoint] Total iterations: {sum(adjoint_iters)}")
    print(f"  [Adjoint] Average: {sum(adjoint_iters)/len(adjoint_iters):.1f} iters/step")

    # Step 3: Sensitivity (dT)
    print("  [Sensitivity] Running...")
    sensitivity_start = time.perf_counter()
    dT_init = np.zeros_like(T_init)
    p_n = lambda_field[0][:, :, -1]  # gradient on top surface (simplified)
    p_n_full = np.zeros((NT-1, ni, nj))
    for n in range(NT-1):
        p_n_full[n] = lambda_field[n][:, :, -1]

    dT, sensitivity_iters = multiple_time_step_solver_DHCP_counted(
        dT_init, p_n_full, NT, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, DT,
        rtol=1.0e-8, maxiter=20000
    )
    sensitivity_elapsed = time.perf_counter() - sensitivity_start
    print(f"  [Sensitivity] Complete: {sensitivity_elapsed:.2f}s")
    print(f"  [Sensitivity] Iterations per step: {sensitivity_iters}")
    print(f"  [Sensitivity] Total iterations: {sum(sensitivity_iters)}")
    print(f"  [Sensitivity] Average: {sum(sensitivity_iters)/len(sensitivity_iters):.1f} iters/step")

    cgm_elapsed = time.perf_counter() - cgm_start

    print("\n" + "="*80)
    print("Summary")
    print("="*80)
    print(f"DHCP time:        {dhcp_elapsed:.2f}s")
    print(f"Adjoint time:     {adjoint_elapsed:.2f}s")
    print(f"Sensitivity time: {sensitivity_elapsed:.2f}s")
    print(f"Total CGM time:   {cgm_elapsed:.2f}s")
    print()
    print(f"DHCP total iters:        {sum(dhcp_iters)} ({sum(dhcp_iters)/len(dhcp_iters):.1f}/step)")
    print(f"Adjoint total iters:     {sum(adjoint_iters)} ({sum(adjoint_iters)/len(adjoint_iters):.1f}/step)")
    print(f"Sensitivity total iters: {sum(sensitivity_iters)} ({sum(sensitivity_iters)/len(sensitivity_iters):.1f}/step)")
    print(f"Grand total iters:       {sum(dhcp_iters) + sum(adjoint_iters) + sum(sensitivity_iters)}")
    print("="*80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
