#!/usr/bin/env python3
"""
Standalone DHCP solver validation runner (Python version).

Replicates the behaviour of `julia/scripts/test_dhcp_solver.jl` by
solving the direct heat conduction problem (DHCP) against the standard
measurement dataset, using the matrix-form CG solver from the legacy
Python implementation.
"""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path
from typing import Tuple

import numpy as np

USE_NUMBA = os.environ.get("IHCP_USE_NUMBA", "0") == "1"

if USE_NUMBA:
    from numba import njit, prange  # type: ignore

    NUMBA_AVAILABLE = True
else:  # pragma: no cover - fallback when Numba is disabled or unavailable
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore
        if args and callable(args[0]) and not kwargs:
            return args[0]

        def wrapper(func):
            return func

        return wrapper

    def prange(*args):  # type: ignore
        if len(args) == 1:
            return range(args[0])
        if len(args) == 2:
            return range(args[0], args[1])
        return range(args[0], args[1], args[2])

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]


DEFAULT_DATA_PATH = PROJECT_ROOT / "shared" / "data" / "T_measure_700um_1ms.npy"
RESULTS_DIR = PROJECT_ROOT / "python" / "validation" / "shared" / "results"


# ---------------------------------------------------------------------------#
# Helpers
# ---------------------------------------------------------------------------#

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Solve the DHCP forward problem for benchmarking.",
    )
    parser.add_argument(
        "--nt",
        type=int,
        default=10,
        help="Number of time steps to simulate (default: 10).",
    )
    parser.add_argument(
        "--nk",
        type=int,
        default=20,
        help="Number of depth-wise grid cells (default: 20).",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=1.0e-3,
        help="Time step size in seconds (default: 1.0e-3).",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=1.0e-6,
        help="Relative residual tolerance for CG (default: 1e-6).",
    )
    parser.add_argument(
        "--maxiter",
        type=int,
        default=20_000,
        help="Maximum CG iterations per time step (default: 20000).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=RESULTS_DIR / "dhcp_validation.npz",
        help="Path to persist solver outputs (default: validation/shared/results/dhcp_validation.npz).",
    )
    return parser.parse_args()


# Material coefficients (aligned with Julia implementation)
CP_COEFFS = np.array(
    [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02], dtype=np.float64
)
K_COEFFS = np.array(
    [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00], dtype=np.float64
)
RHO = np.float64(7_823.493962874829)


# ---------------------------------------------------------------------------#
# Legacy core routines (re-implemented locally with Numba)
# ---------------------------------------------------------------------------#


@njit
def polyval_numba(coeffs: np.ndarray, x: float) -> float:
    result = 0.0
    degree = coeffs.shape[0]
    for idx in range(degree):
        result += coeffs[idx] * x ** (degree - idx - 1)
    return result


@njit(parallel=True)
def thermal_properties_calculator(
    Temperature: np.ndarray, cp_coeffs: np.ndarray, k_coeffs: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    ni, nj, nk = Temperature.shape
    cp = np.empty((ni, nj, nk), dtype=np.float64)
    k = np.empty((ni, nj, nk), dtype=np.float64)
    for i in prange(ni):
        for j in range(nj):
            for k_ijk in range(nk):
                T_current = Temperature[i, j, k_ijk]
                cp[i, j, k_ijk] = polyval_numba(cp_coeffs, T_current)
                k[i, j, k_ijk] = polyval_numba(k_coeffs, T_current)
    return cp, k


@njit(parallel=True, fastmath=False)
def coeffs_and_rhs_building_DHCP(
    T_initial: np.ndarray,
    q_surface: np.ndarray,
    rho: float,
    cp: np.ndarray,
    k: np.ndarray,
    dx: float,
    dy: float,
    dz: np.ndarray,
    dz_b: np.ndarray,
    dz_t: np.ndarray,
    dt: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ni, nj, nk = T_initial.shape
    N = ni * nj * nk

    a_w = np.zeros(N, dtype=np.float64)
    a_e = np.zeros(N, dtype=np.float64)
    a_s = np.zeros(N, dtype=np.float64)
    a_n = np.zeros(N, dtype=np.float64)
    a_b = np.zeros(N, dtype=np.float64)
    a_t = np.zeros(N, dtype=np.float64)
    a_p = np.zeros(N, dtype=np.float64)
    b = np.zeros(N, dtype=np.float64)

    for p in prange(N):
        i = p % ni
        j = (p // ni) % nj
        k_ijk = p // (ni * nj)

        dz_k = dz[k_ijk]
        dz_t_k = dz_t[k_ijk]
        dz_b_k = dz_b[k_ijk]

        k_p = k[i, j, k_ijk]

        a_p_0 = rho * cp[i, j, k_ijk] * dx * dy * dz_k / dt

        a_w[p] = 0.0 if j == 0 else (2 * k_p * k[i, j - 1, k_ijk] / (k_p + k[i, j - 1, k_ijk])) * dy * dz_k / dx
        a_e[p] = 0.0 if j == nj - 1 else (2 * k_p * k[i, j + 1, k_ijk] / (k_p + k[i, j + 1, k_ijk])) * dy * dz_k / dx
        a_s[p] = 0.0 if i == 0 else (2 * k_p * k[i - 1, j, k_ijk] / (k_p + k[i - 1, j, k_ijk])) * dx * dz_k / dy
        a_n[p] = 0.0 if i == ni - 1 else (2 * k_p * k[i + 1, j, k_ijk] / (k_p + k[i + 1, j, k_ijk])) * dx * dz_k / dy
        a_b[p] = 0.0 if k_ijk == 0 else (2 * k_p * k[i, j, k_ijk - 1] / (k_p + k[i, j, k_ijk - 1])) * dx * dy / dz_b_k
        a_t[p] = 0.0 if k_ijk == nk - 1 else (2 * k_p * k[i, j, k_ijk + 1] / (k_p + k[i, j, k_ijk + 1])) * dx * dy / dz_t_k

        a_p[p] = a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p] + a_p_0

        rhs = a_p_0 * T_initial[i, j, k_ijk]
        if k_ijk == nk - 1:
            rhs += q_surface[i, j] * dx * dy
        b[p] = rhs

    return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


# ---------------------------------------------------------------------------#
# Diagnostics container and CG solver utilities
# ---------------------------------------------------------------------------#


class SolverDiagnostics:
    """Collect per-timestep solver telemetry."""

    def __init__(self) -> None:
        self.iterations = []
        self.initial_residuals = []
        self.residuals = []
        self.convergence_flags = []
        self.elapsed_time = 0.0

    def add_timestep(self, iters: int, init_resid: float, final_resid: float, info: int) -> None:
        self.iterations.append(iters)
        self.initial_residuals.append(init_resid)
        self.residuals.append(final_resid)
        self.convergence_flags.append(info)

    def summary(self) -> dict:
        if not self.iterations:
            return {
                "total_iters": 0,
                "avg_iters": 0.0,
                "max_iters": 0,
                "min_iters": 0,
                "elapsed_time": self.elapsed_time,
                "convergence_rate": 0.0,
                "timesteps": 0,
            }
        total_iters = sum(self.iterations)
        max_iters = max(self.iterations)
        min_iters = min(self.iterations)
        avg_iters = total_iters / len(self.iterations)
        converged = sum(1 for info in self.convergence_flags if info == 0)
        return {
            "total_iters": total_iters,
            "avg_iters": avg_iters,
            "max_iters": max_iters,
            "min_iters": min_iters,
            "elapsed_time": self.elapsed_time,
            "convergence_rate": converged / len(self.iterations),
            "timesteps": len(self.iterations),
        }


@njit(parallel=True)
def matvec_stencil(
    x: np.ndarray,
    ni: int,
    nj: int,
    nk: int,
    a_w: np.ndarray,
    a_e: np.ndarray,
    a_s: np.ndarray,
    a_n: np.ndarray,
    a_b: np.ndarray,
    a_t: np.ndarray,
    a_p: np.ndarray,
) -> np.ndarray:
    N = x.size
    result = np.empty_like(x)
    stride_j = ni
    stride_k = ni * nj
    for p in prange(N):
        i = p % ni
        j = (p // ni) % nj
        k = p // (ni * nj)
        val = a_p[p] * x[p]
        if j > 0:
            val -= a_w[p] * x[p - stride_j]
        if j < nj - 1:
            val -= a_e[p] * x[p + stride_j]
        if i > 0:
            val -= a_s[p] * x[p - 1]
        if i < ni - 1:
            val -= a_n[p] * x[p + 1]
        if k > 0:
            val -= a_b[p] * x[p - stride_k]
        if k < nk - 1:
            val -= a_t[p] * x[p + stride_k]
        result[p] = val
    return result


@njit
def apply_jacobi(r: np.ndarray, a_p: np.ndarray) -> np.ndarray:
    z = np.empty_like(r)
    for idx in range(r.size):
        diag = a_p[idx]
        z[idx] = r[idx] / diag if diag != 0.0 else r[idx]
    return z


@njit
def cg_solve(
    ni: int,
    nj: int,
    nk: int,
    a_w: np.ndarray,
    a_e: np.ndarray,
    a_s: np.ndarray,
    a_n: np.ndarray,
    a_b: np.ndarray,
    a_t: np.ndarray,
    a_p: np.ndarray,
    b: np.ndarray,
    x0: np.ndarray,
    rtol: float,
    maxiter: int,
) -> Tuple[np.ndarray, int, int, float, float]:
    x = x0.copy()
    Ax = matvec_stencil(x, ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)
    r = b - Ax
    z = apply_jacobi(r, a_p)
    p = z.copy()

    rz_old = np.dot(r, z)
    res_norm = np.sqrt(np.dot(r, r))
    init_resid = res_norm
    norm_b = np.sqrt(np.dot(b, b))
    if norm_b == 0.0:
        norm_b = 1.0
    tol = rtol * norm_b

    if res_norm <= tol:
        return x, 0, 0, init_resid, res_norm

    for it in range(1, maxiter + 1):
        Ap = matvec_stencil(p, ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)
        denom = np.dot(p, Ap)
        if denom == 0.0:
            return x, -1, it - 1, init_resid, res_norm
        alpha = rz_old / denom
        x = x + alpha * p
        r = r - alpha * Ap
        res_norm = np.sqrt(np.dot(r, r))
        if res_norm <= tol:
            final_resid = res_norm
            return x, 0, it, init_resid, final_resid
        z = apply_jacobi(r, a_p)
        rz_new = np.dot(r, z)
        if rz_old == 0.0:
            return x, -2, it, init_resid, res_norm
        beta = rz_new / rz_old
        p = z + beta * p
        rz_old = rz_new

    final_resid = np.sqrt(np.dot(r, r))
    return x, 1, maxiter, init_resid, final_resid


def multiple_time_step_solver_dhcp(
    T_initial: np.ndarray,
    q_surface: np.ndarray,
    nt: int,
    rho: float,
    cp_coeffs: np.ndarray,
    k_coeffs: np.ndarray,
    dx: float,
    dy: float,
    dz: np.ndarray,
    dz_b: np.ndarray,
    dz_t: np.ndarray,
    dt: float,
    rtol: float,
    maxiter: int,
    verbose: bool = False,
) -> Tuple[np.ndarray, SolverDiagnostics]:
    ni, nj, nk = T_initial.shape
    T_all = np.empty((nt, ni, nj, nk), dtype=np.float64)
    T_all[0] = T_initial

    x0 = T_initial.ravel(order="F").copy()
    diagnostics = SolverDiagnostics()
    start_time = time.perf_counter()

    for t in range(1, nt):
        cp, k = thermal_properties_calculator(T_all[t - 1], cp_coeffs, k_coeffs)
        a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
            T_all[t - 1],
            q_surface[t - 1],
            rho,
            cp,
            k,
            dx,
            dy,
            dz,
            dz_b,
            dz_t,
            dt,
        )

        x, info, iters, init_resid, final_resid = cg_solve(
            ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p, b, x0, rtol, maxiter
        )
        diagnostics.add_timestep(iters, init_resid, final_resid, info)

        if verbose:
            print(
                f"  [DHCP][t={t}] iterations={iters:4d} info={info:2d} "
                f"res0={init_resid:.3e} res={final_resid:.3e}"
            )

        T_all[t] = x.reshape((ni, nj, nk), order="F")
        x0 = x

    diagnostics.elapsed_time = time.perf_counter() - start_time
    return T_all, diagnostics


def build_z_grid(nk: int, Lz: float, stretch_factor: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Construct stretched grid spacing arrays matching the legacy solver."""
    z_faces_normalized = np.linspace(1.0, 0.0, nk + 1, dtype=np.float64)
    exp_factor = np.exp(stretch_factor)
    z_faces = Lz - (Lz / (exp_factor - 1.0)) * (np.exp(stretch_factor * z_faces_normalized) - 1.0)

    dz = np.diff(z_faces).astype(np.float64)

    z_centers = np.empty(nk, dtype=np.float64)
    z_centers[0] = z_faces[0]
    z_centers[-1] = z_faces[-1]
    if nk > 2:
        z_centers[1:-1] = 0.5 * (z_faces[1:nk-1] + z_faces[2:nk])

    dz_t = np.empty(nk, dtype=np.float64)
    dz_t[:-1] = z_centers[1:] - z_centers[:-1]
    dz_t[-1] = np.inf

    dz_b = np.empty(nk, dtype=np.float64)
    dz_b[0] = np.inf
    dz_b[1:] = z_centers[1:] - z_centers[:-1]

    return dz, dz_b, dz_t


def load_measurement(nt: int) -> np.ndarray:
    """Load the first nt frames from the shared measurement dataset."""
    data_path = DEFAULT_DATA_PATH
    if not data_path.exists():
        raise FileNotFoundError(f"Measurement file not found: {data_path}")

    print(f"  loading measurement data from: {data_path}")
    load_start = time.perf_counter()
    mmapped = np.load(data_path, mmap_mode="r")
    try:
        if mmapped.shape[0] < nt:
            raise ValueError(
                f"Measurement file only has {mmapped.shape[0]} time steps; need {nt}."
            )
        subset = np.array(mmapped[:nt], dtype=np.float64, copy=True)
    finally:
        del mmapped
    elapsed = time.perf_counter() - load_start
    print(f"  measurement subset shape: {subset.shape} (t, i, j)")
    print(f"  load time: {elapsed:.2f} s")
    return subset


def build_initial_temperature(first_frame: np.ndarray, nk: int) -> np.ndarray:
    """Tile the first measurement frame through depth to initialise T."""
    return np.repeat(first_frame[:, :, np.newaxis], nk, axis=2).astype(np.float64, copy=False)


def analyze_residuals(T_result: np.ndarray, Y_obs: np.ndarray) -> Tuple[float, float]:
    """Compute RMS and max residuals between simulated surface and observation."""
    T_surface = T_result[:, :, :, 0]
    residual = T_surface - Y_obs
    rms_error = float(np.sqrt(np.mean(residual**2)))
    max_error = float(np.max(np.abs(residual)))
    print("  Mode: Measurement data")
    print(f"  RMS residual: {rms_error:.4e} K")
    print(f"  Max residual: {max_error:.4e} K")
    print(f"  Mean temperature: {float(np.mean(T_result)):.2f} K")
    print(f"  Temperature range: {float(np.min(T_result)):.2f} ~ {float(np.max(T_result)):.2f} K")
    return rms_error, max_error


def ensure_num_threads():
    """Warn if NUMBA_NUM_THREADS is unset, mirroring Julia helper."""
    if not NUMBA_AVAILABLE:
        print("  Numba disabled (set IHCP_USE_NUMBA=1 to enable acceleration).")
        return
    if "NUMBA_NUM_THREADS" not in os.environ:
        os.environ.setdefault("NUMBA_NUM_THREADS", "1")
        print("  NUMBA_NUM_THREADS not set; defaulting to 1. "
              "Set this environment variable to control parallelism.")


# ---------------------------------------------------------------------------#
# Main execution
# ---------------------------------------------------------------------------#

def main() -> int:
    args = parse_args()

    ensure_num_threads()
    print("=" * 80)
    print("Python DHCP solver validation (matrix CG)")
    print("=" * 80)
    print(f"Project root: {PROJECT_ROOT}")
    print(f"NUMBA_NUM_THREADS={os.environ.get('NUMBA_NUM_THREADS')}")

    if args.nt < 2:
        raise ValueError("--nt must be >= 2.")
    if args.nk < 1:
        raise ValueError("--nk must be positive.")

    # Problem constants (aligned with Julia test script)
    dx = 0.12e-3
    dy = 0.12e-3 * np.sin(np.deg2rad(80.0)) / np.sin(np.deg2rad(45.0))
    Lz = 0.5e-3
    stretch_factor = 3.0

    print("\n[Configuration]")
    print(f"  Time steps: {args.nt}")
    print(f"  Depth cells: {args.nk}")
    print(f"  Time step (dt): {args.dt:.3e} s")
    print(f"  Grid spacing: dx={dx:.3e} m, dy={dy:.3e} m")

    total_start = time.perf_counter()

    print("\n[1/5] Grid and material parameters")
    dz, dz_b, dz_t = build_z_grid(args.nk, Lz, stretch_factor)
    print(f"  dz range: {dz.min():.3e} ~ {dz.max():.3e} m")
    print(f"  rho (reference): {RHO:.6f} kg/m^3")

    print("\n[2/5] Loading measurement data")
    Y_obs = load_measurement(args.nt)
    ni, nj = Y_obs.shape[1:]
    print(f"  Grid size (ni, nj, nk): {ni} × {nj} × {args.nk} (N={ni * nj * args.nk})")
    T_init = build_initial_temperature(Y_obs[0], args.nk)

    # Zero heat flux for the standalone forward solve
    q_surface = np.zeros((args.nt - 1, ni, nj), dtype=np.float64)
    print("\n[3/5] Heat flux pattern")
    print(f"  q_surface shape: {q_surface.shape}")
    print(f"  q_surface range: {q_surface.min():.4e} ~ {q_surface.max():.4e} W/m^2")

    print("\n[4/5] Running DHCP solver")
    solve_start = time.perf_counter()
    T_result, diagnostics = multiple_time_step_solver_dhcp(
        T_init,
        q_surface,
        args.nt,
        RHO,
        CP_COEFFS,
        K_COEFFS,
        dx,
        dy,
        dz,
        dz_b,
        dz_t,
        args.dt,
        args.rtol,
        args.maxiter,
        verbose=True,
    )
    dhcp_elapsed = time.perf_counter() - solve_start
    print(f"  DHCP elapsed time: {dhcp_elapsed:.2f} s")
    if isinstance(diagnostics, SolverDiagnostics):
        summary = diagnostics.summary()
        print("  Solver diagnostics:")
        print(f"    Total iterations: {summary['total_iters']}")
        print(f"    Avg iterations:   {summary['avg_iters']:.2f}")
        print(f"    Convergence rate: {summary['convergence_rate'] * 100:.1f}%")

    print("\n[5/5] Residual analysis")
    rms_error, max_error = analyze_residuals(T_result, Y_obs)

    total_elapsed = time.perf_counter() - total_start
    print("\nSummary")
    print(f"  Total runtime: {total_elapsed:.2f} s")
    print(f"  DHCP share: {100.0 * dhcp_elapsed / total_elapsed:.1f}%")
    print("=" * 80)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output,
        T=T_result,
        Y_obs=Y_obs,
        q=q_surface,
        dt=args.dt,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_b=dz_b,
        dz_t=dz_t,
        rho=RHO,
        cp_coeffs=CP_COEFFS,
        k_coeffs=K_COEFFS,
        rms_error=rms_error,
        max_error=max_error,
        diagnostics=(
            diagnostics.summary() if isinstance(diagnostics, SolverDiagnostics) else None
        ),
    )
    print(f"Saved bundled results to: {args.output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
