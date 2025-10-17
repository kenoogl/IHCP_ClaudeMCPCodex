#!/usr/bin/env python3
"""Python版の係数行列が対称かどうかを確認"""

import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")

import sys
from pathlib import Path

import numpy as np

THIS_FILE = Path(__file__).resolve()
PROJECT_ROOT = THIS_FILE.parents[2]
ORIGINAL_DIR = THIS_FILE.parents[1] / "original"
if str(ORIGINAL_DIR) not in sys.path:
    sys.path.insert(0, str(ORIGINAL_DIR))

from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    thermal_properties_calculator,
    coeffs_and_rhs_building_DHCP,
    assemble_A_DHCP,
    rho,
    cp_coeffs,
    k_coeffs,
    dx,
    dy,
    dz,
    dz_b,
    dz_t,
)

def check_matrix_symmetry():
    """係数行列の対称性をチェック"""

    # 小さなテスト問題
    ni, nj, nk = 10, 10, 5
    N = ni * nj * nk

    # 初期温度場
    T_init = np.ones((ni, nj, nk), dtype=np.float64) * 500.0
    q_surface = np.zeros((ni, nj), dtype=np.float64)
    dt = 1.0e-3

    # 熱物性値計算
    cp, k = thermal_properties_calculator(T_init, cp_coeffs, k_coeffs)

    # 係数構築
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
        T_init, q_surface, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt)

    # 疎行列組み立て
    A_csr = assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    print("="*80)
    print("Python版係数行列の対称性チェック")
    print("="*80)
    print(f"行列サイズ: {N}x{N} (ni={ni}, nj={nj}, nk={nk})")
    print(f"非ゼロ要素数: {A_csr.nnz}")
    print()

    # 対称性チェック: A - A^T のノルムを計算
    A_dense = A_csr.toarray()
    A_T = A_dense.T
    diff = A_dense - A_T

    # 対称性の指標
    frobenius_norm_diff = np.linalg.norm(diff, 'fro')
    frobenius_norm_A = np.linalg.norm(A_dense, 'fro')
    relative_asymmetry = frobenius_norm_diff / frobenius_norm_A

    max_diff = np.max(np.abs(diff))
    max_A = np.max(np.abs(A_dense))

    print(f"||A - A^T||_F: {frobenius_norm_diff:.6e}")
    print(f"||A||_F:       {frobenius_norm_A:.6e}")
    print(f"相対非対称度:   {relative_asymmetry:.6e}")
    print()
    print(f"max|A - A^T|:  {max_diff:.6e}")
    print(f"max|A|:        {max_A:.6e}")
    print()

    # 判定
    if relative_asymmetry < 1e-14:
        print("✅ 判定: 係数行列は**対称**です（数値誤差範囲内）")
        print("   → CG法の使用は適切")
    elif relative_asymmetry < 1e-10:
        print("⚠️  判定: 係数行列はほぼ対称ですが、わずかに非対称です")
        print("   → CG法は収束するが、BICGSTABの方が理論的に正確")
    else:
        print("❌ 判定: 係数行列は**非対称**です")
        print("   → CG法の使用は不適切、BICGSTABを使うべき")

    print()
    print("="*80)
    print("非対称性の詳細分析")
    print("="*80)

    # 非対称要素の分布を確認
    nonzero_diff = np.abs(diff) > 1e-14
    num_asymmetric = np.sum(nonzero_diff)

    print(f"非対称要素数: {num_asymmetric} / {N*N} ({num_asymmetric/(N*N)*100:.2f}%)")

    # どの方向の結合で非対称が生じているか
    if num_asymmetric > 0:
        print("\n非対称が生じている位置を解析...")

        # Z方向の結合を確認
        off_k = ni * nj
        z_coupling_asymmetry = 0
        for i in range(N - off_k):
            z_coupling_asymmetry += abs(diff[i, i+off_k] - diff[i+off_k, i])

        # X, Y方向の結合を確認
        off_i = 1
        off_j = ni
        x_coupling_asymmetry = 0
        y_coupling_asymmetry = 0

        for i in range(N - off_i):
            if (i+1) % ni != 0:  # 境界をまたがない
                x_coupling_asymmetry += abs(diff[i, i+off_i] - diff[i+off_i, i])

        for i in range(N - off_j):
            if (i+off_j) < N:
                y_coupling_asymmetry += abs(diff[i, i+off_j] - diff[i+off_j, i])

        print(f"  X方向結合の非対称度: {x_coupling_asymmetry:.6e}")
        print(f"  Y方向結合の非対称度: {y_coupling_asymmetry:.6e}")
        print(f"  Z方向結合の非対称度: {z_coupling_asymmetry:.6e}")

        if z_coupling_asymmetry > x_coupling_asymmetry and z_coupling_asymmetry > y_coupling_asymmetry:
            print("\n  → Z方向（非均一格子）が非対称性の主要因")
        else:
            print("\n  → 非対称性は均等に分布（予期しない）")

    print()
    print("="*80)

    return relative_asymmetry


if __name__ == "__main__":
    asymmetry = check_matrix_symmetry()

    if asymmetry > 1e-10:
        print("\n⚠️  警告: Python版はCG法を使用していますが、行列は非対称です")
        print("   CG法は対称行列専用のアルゴリズムです")
        print("   非対称行列に対してはBICGSTABなどを使用すべきです")
