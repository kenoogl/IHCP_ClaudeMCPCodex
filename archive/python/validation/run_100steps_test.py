#!/usr/bin/env python3
"""Python 100-step CGM validation runner."""

import os
os.environ.setdefault("NUMBA_NUM_THREADS", "8")

import sys
import time
from pathlib import Path

import numpy as np

# Add original implementation to import path before importing heavy modules.
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


def load_measurement(nt: int) -> np.ndarray:
    """Load the first ``nt`` steps from the shared measurement file."""
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
        # Copy the window we actually need to keep resident in memory.
        subset = np.array(mmapped[:nt], dtype=np.float64, copy=True)
    finally:
        # Explicitly close the memmap file handle.
        del mmapped
    print(f"Loaded measurement subset shape: {subset.shape}")
    return subset


def build_initial_temperature(first_frame: np.ndarray) -> np.ndarray:
    """Construct the initial 3D temperature field from the first measurement frame."""
    nk = dz.shape[0]
    # Repeat along the z-axis to create the volumetric initial condition.
    return np.repeat(first_frame[:, :, np.newaxis], nk, axis=2).astype(np.float64, copy=False)


def main() -> int:
    print("=" * 80)
    print("Python 100-step CGM validation")
    print("=" * 80)
    print(f"NUMBA_NUM_THREADS={os.environ.get('NUMBA_NUM_THREADS')}")

    nt = 100
    dt = 1.0e-3  # 1 ms time step as documented.

    start_total = time.perf_counter()

    # ------------------------------------------------------------------
    # Prepare input data and initial conditions.
    # ------------------------------------------------------------------
    Y_obs = load_measurement(nt)
    T_init = build_initial_temperature(Y_obs[0])
    ni, nj = Y_obs.shape[1:]
    q_init = np.zeros((nt - 1, ni, nj), dtype=np.float64)

    # ------------------------------------------------------------------
    # Run CGM inversion to estimate surface heat flux.
    # ------------------------------------------------------------------
    print("\n[1/3] Running global CGM (1 iteration for comparison)...")
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
        dt,
        rho,
        cp_coeffs,
        k_coeffs,
        CGM_iteration=1,
    )
    cgm_elapsed = time.perf_counter() - cgm_start
    print(f"CGM complete. Iterations: {len(J_history)} | Elapsed: {cgm_elapsed:.2f} s")

    # ------------------------------------------------------------------
    # Forward verification: recompute temperature field with estimated q.
    # ------------------------------------------------------------------
    print("\n[2/3] Verifying via forward solver...")
    verify_start = time.perf_counter()
    T_forward = multiple_time_step_solver_DHCP(
        T_init,
        q_result,
        nt,
        rho,
        cp_coeffs,
        k_coeffs,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        dt,
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
    results_dir = PROJECT_ROOT / "shared" / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    q_path = results_dir / "python_100steps_cgm1_q.npy"
    T_path = results_dir / "python_100steps_cgm1_T.npy"
    npz_path = results_dir / "python_100steps_cgm1.npz"

    np.save(q_path, q_result)
    np.save(T_path, T_forward)
    np.savez_compressed(
        npz_path,
        Y_obs=Y_obs,
        T_init=T_init,
        dt=dt,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_t=dz_t,
        dz_b=dz_b,
        rho=rho,
        cp_coeffs=cp_coeffs,
        k_coeffs=k_coeffs,
        q_result=q_result,
        T_forward=T_forward,
        J_history=np.array(J_history, dtype=np.float64),
        elapsed_cgm=cgm_elapsed,
        elapsed_verify=verify_elapsed,
        rms_error=rms_error,
        max_error=max_error,
    )

    print(f"Saved heat flux to: {q_path}")
    print(f"Saved temperature field to: {T_path}")
    print(f"Saved bundled results to: {npz_path}")

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
