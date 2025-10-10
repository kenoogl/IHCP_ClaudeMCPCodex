#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phase 5参照データ生成スクリプト

Python版のIHCP-CGMを使用してスライディングウィンドウ計算の参照データを生成し、
JSON形式で保存します。

生成されるファイル:
- julia/data/phase5_reference_sliding_window_1D.json
- julia/data/phase5_reference_sliding_window_2D.json
"""

import sys
import json
import numpy as np
from pathlib import Path

# 親ディレクトリのモジュールをインポート
sys.path.insert(0, str(Path(__file__).parent / "original"))
from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    multiple_time_step_solver_DHCP,
    global_CGM_time
)

def generate_1D_reference():
    """1D小規模スライディングウィンドウ参照データ生成"""
    print("="*70)
    print("1D参照データ生成開始")
    print("="*70)

    # 問題設定（scipy.sparseのdiags制約回避のため、最小サイズを3x3に変更）
    ni, nj, nk = 3, 3, 3
    nt = 16
    dx, dy, dt = 0.001, 0.001, 0.01

    dz = np.array([0.001, 0.001, 0.001])
    dz_b = np.array([0.0005, 0.001, 0.001])
    dz_t = np.array([0.001, 0.001, 0.0005])

    rho = 7900.0
    cp_coeffs = np.array([0.0, 0.0, 0.0, 500.0])
    k_coeffs = np.array([0.0, 0.0, 0.0, 15.0])

    # スライディングウィンドウパラメータ
    window_size = 10
    overlap = 3
    q_init_value = 0.0
    cgm_iteration = 10
    sigma_noise = 1.0

    # 初期条件
    T_init = np.full((ni, nj, nk), 300.0, dtype=np.float64)

    # 真の熱流束（線形増加）
    q_true = np.zeros((nt-1, ni, nj), dtype=np.float64)
    for t in range(nt-1):
        q_true[t, :, :] = 1000.0 + 200.0 * t

    # 観測データ生成（DHCP + ノイズ）
    print("観測データ生成中...")
    T_all = multiple_time_step_solver_DHCP(
        T_init, q_true, nt, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, dt,
        rtol=1e-8, maxiter=20000
    )

    # 底面温度 + ノイズ
    Y_obs = T_all[:, :, :, 0].copy()  # (nt, ni, nj)
    np.random.seed(42)
    noise = np.random.normal(0, sigma_noise, Y_obs.shape)
    Y_obs += noise

    print(f"観測データ形状: {Y_obs.shape}")
    print(f"温度範囲: {Y_obs.min():.2f} - {Y_obs.max():.2f} K")

    # スライディングウィンドウCGM実行
    print("\nスライディングウィンドウCGM実行中...")
    start_idx = 0
    q_total = []
    windows_info = []
    T_init_win = T_init.copy()
    prev_q_win = None
    window_id = 0

    safety_counter = 0
    safety_limit = nt * 5

    while start_idx < nt - 1:
        safety_counter += 1
        if safety_counter > safety_limit:
            print("Safety break")
            break

        window_id += 1
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        Y_obs_win = Y_obs[start_idx:end_idx+1, :, :]

        print(f"\nウィンドウ {window_id}: [{start_idx}, {end_idx}], 長さ={max_L}")

        # 初期熱流束
        if prev_q_win is None:
            q_init_win = np.full((max_L, ni, nj), q_init_value, dtype=np.float64)
        else:
            q_init_win = np.empty((max_L, ni, nj), dtype=np.float64)
            L_overlap = min(overlap, max_L, prev_q_win.shape[0])
            if L_overlap > 0:
                q_init_win[:L_overlap] = prev_q_win[-L_overlap:]
            if L_overlap < max_L:
                q_init_win[L_overlap:] = prev_q_win[-1]

        # CGM実行
        q_win, T_win_last, J_hist = global_CGM_time(
            T_init_win, Y_obs_win, q_init_win,
            dx, dy, dz, dz_b, dz_t, dt,
            rho, cp_coeffs, k_coeffs,
            CGM_iteration=cgm_iteration
        )

        print(f"  反復回数: {len(J_hist)}, 最終J: {J_hist[-1]:.6e}")

        prev_q_win = q_win.copy()

        # オーバーラップ処理
        overlap_steps = 0
        if len(q_total) == 0:
            q_total.append(q_win)
        else:
            overlap_steps = min(overlap, q_win.shape[0], q_total[-1].shape[0])
            if overlap_steps > 0:
                q_total[-1][-overlap_steps:] = 0.5 * q_total[-1][-overlap_steps:] + 0.5 * q_win[:overlap_steps]
                q_total.append(q_win[overlap_steps:])
            else:
                q_total.append(q_win)
            print(f"  オーバーラップ平均化: {overlap_steps}ステップ")

        # ウィンドウ情報記録
        windows_info.append({
            "window_id": window_id,
            "start_idx": int(start_idx),
            "end_idx": int(end_idx),
            "max_L": int(max_L),
            "n_iterations": len(J_hist),
            "J_final": float(J_hist[-1]),
            "q_win_shape": list(q_win.shape),
            "overlap_steps": int(overlap_steps)
        })

        T_init_win = T_win_last.copy()
        step = max(1, max_L - overlap)
        start_idx += step

    # グローバル熱流束
    q_global = np.concatenate(q_total, axis=0)[:nt-1]

    print(f"\n結果:")
    print(f"  総ウィンドウ数: {len(windows_info)}")
    print(f"  最終q_global形状: {q_global.shape}")

    # JSON形式で保存
    reference_data = {
        "problem": {
            "ni": int(ni),
            "nj": int(nj),
            "nk": int(nk),
            "nt": int(nt),
            "dx": float(dx),
            "dy": float(dy),
            "dz": dz.tolist(),
            "dz_b": dz_b.tolist(),
            "dz_t": dz_t.tolist(),
            "dt": float(dt),
            "rho": float(rho),
            "cp_coeffs": cp_coeffs.tolist(),
            "k_coeffs": k_coeffs.tolist(),
            "sigma_noise": float(sigma_noise)
        },
        "sliding_window": {
            "window_size": int(window_size),
            "overlap": int(overlap),
            "q_init_value": float(q_init_value),
            "CGM_iteration": int(cgm_iteration)
        },
        "input": {
            "T_init": T_init.tolist(),
            "Y_obs": Y_obs.tolist(),
            "q_true": q_true.tolist()
        },
        "output": {
            "q_global": q_global.tolist(),
            "n_windows": len(windows_info),
            "windows_info": windows_info
        }
    }

    output_path = Path(__file__).parent.parent / "julia" / "data" / "phase5_reference_sliding_window_1D.json"
    with open(output_path, 'w') as f:
        json.dump(reference_data, f, indent=2)

    print(f"\n✅ 1D参照データ保存完了: {output_path}")
    return reference_data


def generate_2D_reference():
    """2D小規模スライディングウィンドウ参照データ生成"""
    print("\n" + "="*70)
    print("2D参照データ生成開始")
    print("="*70)

    # 問題設定（scipy.sparseのdiags制約回避のため、最小サイズを3x3に変更）
    ni, nj, nk = 3, 3, 3
    nt = 14
    dx, dy, dt = 0.001, 0.001, 0.01

    dz = np.array([0.001, 0.001, 0.001])
    dz_b = np.array([0.0005, 0.001, 0.001])
    dz_t = np.array([0.001, 0.001, 0.0005])

    rho = 7900.0
    cp_coeffs = np.array([0.0, 0.0, 0.0, 500.0])
    k_coeffs = np.array([0.0, 0.0, 0.0, 15.0])

    # スライディングウィンドウパラメータ
    window_size = 10
    overlap = 3
    q_init_value = 500.0
    cgm_iteration = 10
    sigma_noise = 1.0

    # 初期条件
    T_init = np.full((ni, nj, nk), 300.0, dtype=np.float64)

    # 真の熱流束（正弦波）
    q_true = np.zeros((nt-1, ni, nj), dtype=np.float64)
    for t in range(nt-1):
        q_base = 500.0 * (1.0 + 0.5 * np.sin(2 * np.pi * t / (nt-1)))
        q_true[t, :, :] = q_base

    # 観測データ生成（DHCP + ノイズ）
    print("観測データ生成中...")
    T_all = multiple_time_step_solver_DHCP(
        T_init, q_true, nt, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, dt,
        rtol=1e-8, maxiter=20000
    )

    # 底面温度 + ノイズ
    Y_obs = T_all[:, :, :, 0].copy()
    np.random.seed(42)
    noise = np.random.normal(0, sigma_noise, Y_obs.shape)
    Y_obs += noise

    print(f"観測データ形状: {Y_obs.shape}")
    print(f"温度範囲: {Y_obs.min():.2f} - {Y_obs.max():.2f} K")

    # スライディングウィンドウCGM実行
    print("\nスライディングウィンドウCGM実行中...")
    start_idx = 0
    q_total = []
    windows_info = []
    T_init_win = T_init.copy()
    prev_q_win = None
    window_id = 0

    safety_counter = 0
    safety_limit = nt * 5

    while start_idx < nt - 1:
        safety_counter += 1
        if safety_counter > safety_limit:
            print("Safety break")
            break

        window_id += 1
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        Y_obs_win = Y_obs[start_idx:end_idx+1, :, :]

        print(f"\nウィンドウ {window_id}: [{start_idx}, {end_idx}], 長さ={max_L}")

        # 初期熱流束
        if prev_q_win is None:
            q_init_win = np.full((max_L, ni, nj), q_init_value, dtype=np.float64)
        else:
            q_init_win = np.empty((max_L, ni, nj), dtype=np.float64)
            L_overlap = min(overlap, max_L, prev_q_win.shape[0])
            if L_overlap > 0:
                q_init_win[:L_overlap] = prev_q_win[-L_overlap:]
            if L_overlap < max_L:
                q_init_win[L_overlap:] = prev_q_win[-1]

        # CGM実行
        q_win, T_win_last, J_hist = global_CGM_time(
            T_init_win, Y_obs_win, q_init_win,
            dx, dy, dz, dz_b, dz_t, dt,
            rho, cp_coeffs, k_coeffs,
            CGM_iteration=cgm_iteration
        )

        print(f"  反復回数: {len(J_hist)}, 最終J: {J_hist[-1]:.6e}")

        prev_q_win = q_win.copy()

        # オーバーラップ処理
        overlap_steps = 0
        if len(q_total) == 0:
            q_total.append(q_win)
        else:
            overlap_steps = min(overlap, q_win.shape[0], q_total[-1].shape[0])
            if overlap_steps > 0:
                q_total[-1][-overlap_steps:] = 0.5 * q_total[-1][-overlap_steps:] + 0.5 * q_win[:overlap_steps]
                q_total.append(q_win[overlap_steps:])
            else:
                q_total.append(q_win)
            print(f"  オーバーラップ平均化: {overlap_steps}ステップ")

        # ウィンドウ情報記録
        windows_info.append({
            "window_id": window_id,
            "start_idx": int(start_idx),
            "end_idx": int(end_idx),
            "max_L": int(max_L),
            "n_iterations": len(J_hist),
            "J_final": float(J_hist[-1]),
            "q_win_shape": list(q_win.shape),
            "overlap_steps": int(overlap_steps)
        })

        T_init_win = T_win_last.copy()
        step = max(1, max_L - overlap)
        start_idx += step

    # グローバル熱流束
    q_global = np.concatenate(q_total, axis=0)[:nt-1]

    print(f"\n結果:")
    print(f"  総ウィンドウ数: {len(windows_info)}")
    print(f"  最終q_global形状: {q_global.shape}")

    # JSON形式で保存
    reference_data = {
        "problem": {
            "ni": int(ni),
            "nj": int(nj),
            "nk": int(nk),
            "nt": int(nt),
            "dx": float(dx),
            "dy": float(dy),
            "dz": dz.tolist(),
            "dz_b": dz_b.tolist(),
            "dz_t": dz_t.tolist(),
            "dt": float(dt),
            "rho": float(rho),
            "cp_coeffs": cp_coeffs.tolist(),
            "k_coeffs": k_coeffs.tolist(),
            "sigma_noise": float(sigma_noise)
        },
        "sliding_window": {
            "window_size": int(window_size),
            "overlap": int(overlap),
            "q_init_value": float(q_init_value),
            "CGM_iteration": int(cgm_iteration)
        },
        "input": {
            "T_init": T_init.tolist(),
            "Y_obs": Y_obs.tolist(),
            "q_true": q_true.tolist()
        },
        "output": {
            "q_global": q_global.tolist(),
            "n_windows": len(windows_info),
            "windows_info": windows_info
        }
    }

    output_path = Path(__file__).parent.parent / "julia" / "data" / "phase5_reference_sliding_window_2D.json"
    with open(output_path, 'w') as f:
        json.dump(reference_data, f, indent=2)

    print(f"\n✅ 2D参照データ保存完了: {output_path}")
    return reference_data


def main():
    """メイン実行"""
    print("="*70)
    print("Phase 5参照データ生成スクリプト")
    print("="*70)
    print()

    # 1D参照データ生成
    ref_1d = generate_1D_reference()

    # 2D参照データ生成
    ref_2d = generate_2D_reference()

    print("\n" + "="*70)
    print("✅ Phase 5参照データ生成完了")
    print("="*70)
    print("\n生成されたファイル:")
    print("  - julia/data/phase5_reference_sliding_window_1D.json")
    print("  - julia/data/phase5_reference_sliding_window_2D.json")
    print("\n次のステップ:")
    print("  julia --project=julia julia/test/runtests.jl")


if __name__ == "__main__":
    main()
