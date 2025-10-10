#!/usr/bin/env python3
"""Python 10-step full-size CGM validation runner."""

import os
os.environ.setdefault("NUMBA_NUM_THREADS", "8")

import sys
import time
from pathlib import Path

import numpy as np

# Ensure the legacy implementation is importable before touching heavy modules.
THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]
ORIGINAL_DIR = THIS_FILE.parents[1] / "original"
if str(ORIGINAL_DIR) not in sys.path:
    sys.path.insert(0, str(ORIGINAL_DIR))

from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    global_CGM_time,
    multiple_time_step_solver_DHCP,
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
RESULT_PATH = PROJECT_ROOT / "shared" / "results" / "python_10steps_fullsize.npz"


def load_measurement(nt: int) -> np.ndarray:
    """Load the first ``nt`` measurement frames from the shared dataset."""
    data_path = PROJECT_ROOT / "shared" / "data" / "T_measure_700um_1ms.npy"
    if not data_path.exists():
        raise FileNotFoundError(f"Measurement file not found: {data_path}")

    print(f"Loading measurement data: {data_path}")
    mmapped = np.load(data_path, mmap_mode="r")
    try:
        if mmapped.shape[0] < nt:
            raise ValueError(
                f"Measurement file only has {mmapped.shape[0]} time steps; need {nt}."
            )
        subset = np.array(mmapped[:nt], dtype=np.float64, copy=True)
    finally:
        del mmapped
    print(f"Loaded measurement subset shape: {subset.shape}")
    return subset


def build_initial_temperature(first_frame: np.ndarray) -> np.ndarray:
    """Tile the first measurement frame through the depth to form ``T_init``."""
    nk = dz.shape[0]
    return np.repeat(first_frame[:, :, np.newaxis], nk, axis=2).astype(np.float64, copy=False)


def main() -> int:
    print("=" * 80)
    print("Python 10-step full-size CGM validation")
    print("=" * 80)
    print(f"NUMBA_NUM_THREADS={os.environ.get('NUMBA_NUM_THREADS')}")

    start_total = time.perf_counter()

    # ------------------------------------------------------------------
    # Prepare input data and initial conditions.
    # ------------------------------------------------------------------
    Y_obs = load_measurement(NT)
    T_init = build_initial_temperature(Y_obs[0])
    ni, nj = Y_obs.shape[1:]
    q_init = np.zeros((NT - 1, ni, nj), dtype=np.float64)

    # ------------------------------------------------------------------
    # Run CGM inversion for a single full window.
    # ------------------------------------------------------------------
    print("\n[1/3] Running global CGM (single window)...")
    cgm_start = time.perf_counter()
    q_result, _, J_history = global_CGM_time(
        T_init,
        Y_obs,
        q_init,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        DT,
        rho,
        cp_coeffs,
        k_coeffs,
        CGM_iteration=1,
    )
    cgm_elapsed = time.perf_counter() - cgm_start
    print(f"CGM complete. Iterations: {len(J_history)} | Elapsed: {cgm_elapsed:.2f} s")

    # ------------------------------------------------------------------
    # Forward verification: recompute temperature field with the estimated q.
    # ------------------------------------------------------------------
    print("\n[2/3] Verifying via forward solver...")
    verify_start = time.perf_counter()
    T_forward = multiple_time_step_solver_DHCP(
        T_init,
        q_result,
        NT,
        rho,
        cp_coeffs,
        k_coeffs,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        DT,
        rtol=1.0e-6,
        maxiter=20000,
    )
    verify_elapsed = time.perf_counter() - verify_start
    T_surface = T_forward[:, :, :, 0]
    residual = T_surface - Y_obs
    rms_error = float(np.sqrt(np.mean(residual**2)))
    max_error = float(np.max(np.abs(residual)))
    print(f"Forward verification complete. Elapsed: {verify_elapsed:.2f} s")
    print(f"Residual RMS: {rms_error:.4e} K | Residual max: {max_error:.4e} K")

    # ------------------------------------------------------------------
    # Persist results for downstream comparison and analysis.
    # ------------------------------------------------------------------
    print("\n[3/3] Saving results...")
    RESULT_PATH.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        RESULT_PATH,
        T=T_forward,
        q=q_result,
        Y_obs=Y_obs,
        T_init=T_init,
        dt=DT,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_b=dz_b,
        dz_t=dz_t,
        rho=rho,
        cp_coeffs=cp_coeffs,
        k_coeffs=k_coeffs,
        rms_error=rms_error,
        max_error=max_error,
        elapsed_cgm=cgm_elapsed,
        elapsed_verify=verify_elapsed,
        J_history=np.array(J_history, dtype=np.float64),
    )
    print(f"Saved bundled results to: {RESULT_PATH}")

    total_elapsed = time.perf_counter() - start_total
    print("\n" + "=" * 80)
    print("Execution summary")
    print("=" * 80)
    print(f"CGM time:      {cgm_elapsed:.2f} s")
    print(f"Forward time:  {verify_elapsed:.2f} s")
    print(f"Total runtime: {total_elapsed:.2f} s")
    print(f"Heat flux range: {q_result.min():.4e} ~ {q_result.max():.4e} W/m^2")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
