#!/usr/bin/env python3
"""前回の100ステップ結果のパラメータを確認"""

import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parents[2]
results_dir = PROJECT_ROOT / "shared" / "results"

print("=" * 80)
print("前回の計算結果のパラメータ確認")
print("=" * 80)

# Python版の結果を確認
print("\n[Python版] python_100steps.npz")
python_npz = np.load(results_dir / "python_100steps.npz")
print("\n含まれるキー:")
for key in sorted(python_npz.keys()):
    print(f"  - {key}: {python_npz[key].shape if hasattr(python_npz[key], 'shape') else type(python_npz[key])}")

print("\nパラメータ詳細:")
print(f"  Y_obs shape: {python_npz['Y_obs'].shape}")
print(f"  T_init shape: {python_npz['T_init'].shape}")
print(f"  q_result shape: {python_npz['q_result'].shape}")
T_key = 'T_verify' if 'T_verify' in python_npz else 'T_forward'
print(f"  {T_key} shape: {python_npz[T_key].shape}")
print(f"  dt: {python_npz['dt']}")
print(f"  dx: {python_npz['dx']}")
print(f"  dy: {python_npz['dy']}")
print(f"  dz shape: {python_npz['dz'].shape}")
print(f"  dz: {python_npz['dz']}")
print(f"  rho: {python_npz['rho']}")
print(f"  cp_coeffs: {python_npz['cp_coeffs']}")
print(f"  k_coeffs: {python_npz['k_coeffs']}")
if 'cgm_iteration' in python_npz:
    print(f"  cgm_iteration: {python_npz['cgm_iteration']}")
if 'elapsed_time_cgm' in python_npz:
    print(f"  elapsed_time_cgm: {python_npz['elapsed_time_cgm']:.2f} s")

# Julia版の結果を確認
print("\n" + "=" * 80)
print("[Julia版] julia_100steps.npz")
julia_npz = np.load(results_dir / "julia_100steps.npz")
print("\n含まれるキー:")
for key in sorted(julia_npz.keys()):
    print(f"  - {key}: {julia_npz[key].shape if hasattr(julia_npz[key], 'shape') else type(julia_npz[key])}")

print("\nパラメータ詳細:")
print(f"  Y_obs shape: {julia_npz['Y_obs'].shape}")
print(f"  T_init shape: {julia_npz['T_init'].shape}")
print(f"  q_result shape: {julia_npz['q_result'].shape}")
julia_T_key = 'T_verify' if 'T_verify' in julia_npz else 'T_result_python'
print(f"  {julia_T_key} shape: {julia_npz[julia_T_key].shape}")
print(f"  dt: {julia_npz['dt']}")
print(f"  dx: {julia_npz['dx']}")
print(f"  dy: {julia_npz['dy']}")
print(f"  dz shape: {julia_npz['dz'].shape}")
print(f"  dz: {julia_npz['dz']}")
print(f"  rho: {julia_npz['rho']}")
if 'rho_coeffs' in julia_npz:
    print(f"  rho_coeffs: {julia_npz['rho_coeffs']}")
print(f"  cp_coeffs: {julia_npz['cp_coeffs']}")
print(f"  k_coeffs: {julia_npz['k_coeffs']}")
if 'cgm_iteration' in julia_npz:
    print(f"  cgm_iteration: {julia_npz['cgm_iteration']}")
if 'elapsed_time_cgm' in julia_npz:
    print(f"  elapsed_time_cgm: {julia_npz['elapsed_time_cgm']:.2f} s")

# グリッドサイズの比較
print("\n" + "=" * 80)
print("グリッドサイズ比較")
print("=" * 80)
python_shape = python_npz['Y_obs'].shape
julia_shape = julia_npz['Y_obs'].shape
print(f"Python Y_obs: {python_shape} → (nt={python_shape[0]}, ni={python_shape[1]}, nj={python_shape[2]})")
print(f"Julia Y_obs:  {julia_shape} → (ni={julia_shape[0]}, nj={julia_shape[1]}, nt={julia_shape[2]})")

python_T_key = 'T_verify' if 'T_verify' in python_npz else 'T_forward'
julia_T_key_cmp = 'T_verify' if 'T_verify' in julia_npz else 'T_result_python'
python_T_shape = python_npz[python_T_key].shape
julia_T_shape = julia_npz[julia_T_key_cmp].shape
print(f"Python T:     {python_T_shape} → (nt={python_T_shape[0]}, ni={python_T_shape[1]}, nj={python_T_shape[2]}, nk={python_T_shape[3]})")
print(f"Julia T:      {julia_T_shape} → (ni={julia_T_shape[0]}, nj={julia_T_shape[1]}, nk={julia_T_shape[2]}, nt={julia_T_shape[3]})")

print("\n" + "=" * 80)
