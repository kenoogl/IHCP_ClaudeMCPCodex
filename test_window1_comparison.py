#!/usr/bin/env python3
"""
ウィンドウ1のPython-Julia比較テスト

ウィンドウ1（start_idx=0, end_idx=70）のみを実行し、
入力・出力を詳細に比較する。
"""

import numpy as np
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "python" / "original"))
sys.path.insert(0, str(PROJECT_ROOT / "python" / "validation"))

from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    thermal_properties_calculator,
    coeffs_and_rhs_building_DHCP,
    assemble_A_DHCP,
    coeffs_and_rhs_building_Adjoint,
    assemble_A_Adjoint,
    cp_coeffs,
    k_coeffs,
    rho,
    nz,
    dx,
    dy,
    dz,
    dz_b,
    dz_t
)

from solvers_with_diagnostics import global_CGM_time_with_diagnostics

def main():
    print("="*80)
    print("ウィンドウ1のPython計算テスト")
    print("="*80)

    # データ読み込み
    data_path = PROJECT_ROOT / "shared" / "data" / "T_measure_700um_1ms.npy"
    Y_obs_full = np.load(str(data_path))
    print(f"データ形状: {Y_obs_full.shape}")

    # ウィンドウ1の設定
    start_idx = 0
    end_idx = 70
    max_L = 70

    # Y_obsのスライス（71ステップ: t=0からt=70）
    Y_obs_win = Y_obs_full[start_idx:end_idx+1, :, :]
    print(f"Y_obs_win形状: {Y_obs_win.shape}")
    print(f"Y_obs_win範囲: [{np.min(Y_obs_win):.2f}, {np.max(Y_obs_win):.2f}] K")

    # 初期温度場
    ni, nj = 80, 100
    nk = 20
    T_measure_init_K = Y_obs_win[0, :, :]
    T_init = np.repeat(T_measure_init_K[:, :, np.newaxis], nk, axis=2).astype(np.float64)
    print(f"\nT_init形状: {T_init.shape}")
    print(f"T_init範囲: [{np.min(T_init):.2f}, {np.max(T_init):.2f}] K")

    # 初期熱流束（全てゼロ）
    q_init = np.zeros((max_L, ni, nj), dtype=np.float64)
    print(f"\nq_init形状: {q_init.shape}")
    print(f"q_init値: 全て {np.unique(q_init)} W/m²")

    # CGM計算（3反復）
    dt = 0.001
    CGM_iteration = 3

    print(f"\n{'='*80}")
    print("CGM計算開始（3反復）")
    print(f"{'='*80}\n")

    q_result, T_result, J_hist, cgm_diag = global_CGM_time_with_diagnostics(
        T_init, Y_obs_win, q_init, dx, dy, dz, dz_b, dz_t, dt,
        rho, cp_coeffs, k_coeffs, CGM_iteration,
        thermal_properties_calculator=thermal_properties_calculator,
        coeffs_and_rhs_building_DHCP=coeffs_and_rhs_building_DHCP,
        assemble_A_DHCP=assemble_A_DHCP,
        coeffs_and_rhs_building_Adjoint=coeffs_and_rhs_building_Adjoint,
        assemble_A_Adjoint=assemble_A_Adjoint,
        verbose=True
    )

    print(f"\n{'='*80}")
    print("CGM計算完了")
    print(f"{'='*80}")
    print(f"q_result形状: {q_result.shape}")
    print(f"q_result範囲: [{np.min(q_result):.6e}, {np.max(q_result):.6e}] W/m²")
    print(f"q_result平均: {np.mean(q_result):.6e} W/m²")
    print(f"q_result標準偏差: {np.std(q_result):.6e} W/m²")
    print(f"\nT_result形状: {T_result.shape}")
    print(f"T_result範囲: [{np.min(T_result):.2f}, {np.max(T_result):.2f}] K")
    print(f"\nJ履歴: {J_hist}")

    # 結果保存
    output_path = PROJECT_ROOT / "shared" / "results" / "window1_python_test.npz"
    np.savez(
        str(output_path),
        Y_obs_win=Y_obs_win,
        T_init=T_init,
        q_init=q_init,
        q_result=q_result,
        T_result=T_result,
        J_hist=J_hist
    )
    print(f"\n結果保存: {output_path}")

if __name__ == "__main__":
    main()
