#!/usr/bin/env python3
"""Compare 100-step outputs from the Python and Julia solvers."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import numpy as np

# Force a non-interactive backend for headless environments before importing pyplot.
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

EPS = 1e-10


@dataclass
class ErrorStats:
    """Container for absolute and relative error statistics."""

    abs_max: float
    abs_mean: float
    abs_rms: float
    rel_max: float
    rel_mean: float
    rel_rms: float


def load_results() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load Python and Julia results from numpy files."""
    results_dir = Path(__file__).resolve().parent.parent.parent / "shared" / "results"

    python_q_path = results_dir / "python_100steps_q.npy"
    python_T_path = results_dir / "python_100steps_T.npy"
    julia_q_path = results_dir / "julia_100steps_q.npy"
    julia_T_path = results_dir / "julia_100steps_T.npy"

    for path in (python_q_path, python_T_path, julia_q_path, julia_T_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing result file: {path}")

    q_python = np.load(python_q_path)
    T_python = np.load(python_T_path)
    q_julia = np.load(julia_q_path)
    T_julia = np.load(julia_T_path)

    # Align Julia arrays to the Python ordering.
    q_julia = np.transpose(q_julia, (2, 0, 1))
    T_julia = np.transpose(T_julia, (3, 0, 1, 2))

    if q_python.shape != q_julia.shape:
        raise ValueError(
            "Heat-flux arrays have incompatible shapes: "
            f"Python {q_python.shape} vs Julia {q_julia.shape}"
        )
    if T_python.shape != T_julia.shape:
        raise ValueError(
            "Temperature arrays have incompatible shapes: "
            f"Python {T_python.shape} vs Julia {T_julia.shape}"
        )

    return q_python, T_python, q_julia, T_julia


def compute_error_stats(reference: np.ndarray, comparison: np.ndarray) -> tuple[ErrorStats, np.ndarray, np.ndarray]:
    """Compute absolute and relative error statistics."""
    diff = reference - comparison
    abs_diff = np.abs(diff)
    rel_diff = abs_diff / (np.abs(reference) + EPS)

    stats = ErrorStats(
        abs_max=float(np.max(abs_diff)),
        abs_mean=float(np.mean(abs_diff)),
        abs_rms=float(np.sqrt(np.mean(abs_diff**2))),
        rel_max=float(np.max(rel_diff)),
        rel_mean=float(np.mean(rel_diff)),
        rel_rms=float(np.sqrt(np.mean(rel_diff**2))),
    )
    return stats, abs_diff, rel_diff


def print_stats(label: str, stats: ErrorStats, unit: str) -> None:
    """Pretty-print statistics for a given quantity."""
    print(f"{label} absolute error max   : {stats.abs_max:.6e} {unit}")
    print(f"{label} absolute error mean  : {stats.abs_mean:.6e} {unit}")
    print(f"{label} absolute error RMS   : {stats.abs_rms:.6e} {unit}")
    print(f"{label} relative error max   : {stats.rel_max:.6e}")
    print(f"{label} relative error mean  : {stats.rel_mean:.6e}")
    print(f"{label} relative error RMS   : {stats.rel_rms:.6e}")


def plot_heat_flux_timeseries(q_python: np.ndarray, q_julia: np.ndarray, output_path: Path) -> None:
    """Plot the spatially averaged heat flux over time."""
    q_python_mean = q_python.reshape(q_python.shape[0], -1).mean(axis=1)
    q_julia_mean = q_julia.reshape(q_julia.shape[0], -1).mean(axis=1)

    time_index = np.arange(q_python.shape[0])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(time_index, q_python_mean, label="Python", linewidth=2)
    ax.plot(time_index, q_julia_mean, label="Julia", linestyle="--", linewidth=2)
    ax.set_xlabel("Time step")
    ax.set_ylabel("q (W/m^2)")
    ax.set_title("Heat flux evolution (spatial average)")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_temperature_snapshot(T_python: np.ndarray, T_julia: np.ndarray, output_path: Path) -> None:
    """Plot a spatial slice of the temperature field at the final time step."""
    # Use the last time slice and the first depth index (assuming layout nt, ni, nj, nk).
    t_idx = -1
    k_idx = 0
    T_py_slice = T_python[t_idx, :, :, k_idx]
    T_ju_slice = T_julia[t_idx, :, :, k_idx]

    vmin = min(float(T_py_slice.min()), float(T_ju_slice.min()))
    vmax = max(float(T_py_slice.max()), float(T_ju_slice.max()))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    for ax, data, title in zip(axes, (T_py_slice, T_ju_slice), ("Python", "Julia")):
        im = ax.imshow(data, origin="lower", cmap="inferno", vmin=vmin, vmax=vmax)
        ax.set_title(f"Temperature snapshot ({title})")
        ax.set_xlabel("j-index")
        ax.set_ylabel("i-index")
        fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.046, pad=0.04)

    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_error_histograms(
    abs_q: np.ndarray,
    rel_q: np.ndarray,
    abs_T: np.ndarray,
    rel_T: np.ndarray,
    output_path: Path,
) -> None:
    """Plot histograms of absolute and relative errors for q and T."""
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    bins = 50

    axes[0, 0].hist(abs_q.ravel(), bins=bins, color="#1f77b4", alpha=0.8)
    axes[0, 0].set_title("q absolute error")
    axes[0, 0].set_xlabel("|Python - Julia| [W/m^2]")
    axes[0, 0].set_ylabel("Count")

    axes[0, 1].hist(rel_q.ravel(), bins=bins, color="#ff7f0e", alpha=0.8)
    axes[0, 1].set_title("q relative error")
    axes[0, 1].set_xlabel("Relative error")
    axes[0, 1].set_ylabel("Count")

    axes[1, 0].hist(abs_T.ravel(), bins=bins, color="#2ca02c", alpha=0.8)
    axes[1, 0].set_title("T absolute error")
    axes[1, 0].set_xlabel("|Python - Julia| [K]")
    axes[1, 0].set_ylabel("Count")

    axes[1, 1].hist(rel_T.ravel(), bins=bins, color="#d62728", alpha=0.8)
    axes[1, 1].set_title("T relative error")
    axes[1, 1].set_xlabel("Relative error")
    axes[1, 1].set_ylabel("Count")

    for ax in axes.flat:
        ax.grid(True, linestyle="--", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main() -> int:
    print("=" * 80)
    print("Python vs Julia comparison (100 steps)")
    print("=" * 80)

    print("Loading result files...")
    q_python, T_python, q_julia, T_julia = load_results()
    print(f"  q shape     : {q_python.shape}")
    print(f"  T shape     : {T_python.shape}")

    print("Computing error statistics...")
    q_stats, q_abs_diff, q_rel_diff = compute_error_stats(q_python, q_julia)
    T_stats, T_abs_diff, T_rel_diff = compute_error_stats(T_python, T_julia)

    print()
    print_stats("Heat flux (q)", q_stats, "W/m^2")
    print()
    print_stats("Temperature (T)", T_stats, "K")

    output_dir = Path(__file__).resolve().parent.parent.parent / "shared" / "results"
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Generating plots...")
    heat_flux_plot = output_dir / "comparison_100steps_heat_flux.png"
    temperature_plot = output_dir / "comparison_100steps_temperature.png"
    error_hist_plot = output_dir / "comparison_100steps_error_hist.png"

    plot_heat_flux_timeseries(
        q_python,
        q_julia,
        heat_flux_plot,
    )
    plot_temperature_snapshot(
        T_python,
        T_julia,
        temperature_plot,
    )
    plot_error_histograms(
        q_abs_diff,
        q_rel_diff,
        T_abs_diff,
        T_rel_diff,
        error_hist_plot,
    )
    print(f"  Saved plots to {heat_flux_plot}, {temperature_plot}, {error_hist_plot}")

    print("Done.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - defensive logging for CLI use
        print(f"Error: {exc}")
        sys.exit(1)
