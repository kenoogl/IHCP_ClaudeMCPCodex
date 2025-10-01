#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python版縮小テストデータ実行スクリプト（簡易版）

Julia版との比較検証用にPython版を実行し、結果を保存する
CGM反復数を最小限（1-2回）に抑えて、短時間で実行可能

実行方法:
  cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/python/benchmark
  python run_test_1min_simple.py

出力:
  - shared/results/python_test_1min.npz: 計算結果
"""

import numpy as np
import time
import sys
from pathlib import Path
import warnings
import os

# 警告を抑制
warnings.filterwarnings('ignore')

# 標準出力を抑制するコンテキストマネージャ
from contextlib import contextmanager

@contextmanager
def suppress_stdout():
    """標準出力を一時的に抑制"""
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

# プロジェクトルート設定
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(PROJECT_ROOT / "python" / "original"))

# オリジナルコードから必要な関数をインポート（出力抑制）
with suppress_stdout():
    from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
        global_CGM_time,
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


def sliding_window_CGM_simple(
    Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, CGM_iteration=1
):
    """
    スライディングウィンドウCGM計算（簡易版）

    - 詳細な出力を抑制
    - CGM反復数を最小限に設定
    """
    nt = Y_obs.shape[0]
    T_init = T0.copy()
    ni, nj, nk = T_init.shape

    start_idx = 0
    q_total = []
    prev_q_win = None

    windows_info = []

    safety_counter = 0
    safety_limit = nt * 5

    print(f"\nPython版計算開始: nt={nt}, window={window_size}, overlap={overlap}, CGM={CGM_iteration}")

    while start_idx < nt - 1:
        safety_counter += 1
        if safety_counter > safety_limit:
            print("[警告] 安全カウンタ超過")
            break

        # 現在のウィンドウ長
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]

        # 初期熱流束の設定
        if prev_q_win is None:
            q_init_win = np.full((max_L, ni, nj), q_init_value, dtype=float)
        else:
            q_init_win = np.empty((max_L, ni, nj), dtype=float)
            L_overlap = min(overlap, max_L, prev_q_win.shape[0])
            if L_overlap > 0:
                q_init_win[:L_overlap] = prev_q_win[-L_overlap:]
            if L_overlap < max_L:
                edge = prev_q_win[-1]
                q_init_win[L_overlap:] = edge

        # CGM計算実行（出力を抑制）
        start_time_win = time.time()
        with suppress_stdout():
            q_win, T_win_last, J_hist = global_CGM_time(
                T_init, Y_obs_win, q_init_win, dx, dy, dz, dz_b, dz_t, dt,
                rho, cp_coeffs, k_coeffs, CGM_iteration=CGM_iteration
            )
        elapsed_win = time.time() - start_time_win

        # ウィンドウ情報を保存
        windows_info.append({
            'window_id': len(windows_info) + 1,
            'start_idx': start_idx,
            'end_idx': end_idx,
            'max_L': max_L,
            'n_iterations': len(J_hist),
            'J_final': J_hist[-1] if len(J_hist) > 0 else 0.0,
            'elapsed_time': elapsed_win
        })

        print(f"  Window {len(windows_info)}: [{start_idx:3d},{end_idx:3d}] J={J_hist[-1]:.3e} t={elapsed_win:.1f}s")

        prev_q_win = q_win.copy()

        # 熱流束の連結（オーバーラップ部分は平均化）
        if len(q_total) == 0:
            q_total.append(q_win)
        else:
            overlap_steps = min(overlap, q_win.shape[0], q_total[-1].shape[0])
            if overlap_steps > 0:
                q_total[-1][-overlap_steps:] = 0.5 * q_total[-1][-overlap_steps:] + 0.5 * q_win[:overlap_steps]
                q_total.append(q_win[overlap_steps:])
            else:
                q_total.append(q_win)

        # 温度場の継承
        T_init = T_win_last.copy() if T_win_last.ndim == 3 else T_win_last[-1].copy()

        # 次のウィンドウへ
        step = max(1, max_L - overlap)
        start_idx += step

    # 全時間の熱流束を連結
    q_global = np.concatenate(q_total, axis=0)[:nt-1]

    print(f"計算完了: q_global.shape={q_global.shape}, windows={len(windows_info)}")

    return q_global, windows_info


def main():
    """メイン実行関数"""

    print("=" * 80)
    print("Python版縮小テストデータ実行（簡易版）")
    print("=" * 80)

    # テストデータ読み込み
    data_path = PROJECT_ROOT / "shared/data/T_measure_test_1min.npy"
    print(f"データファイル: {data_path}")

    if not data_path.exists():
        print(f"エラー: データファイルが見つかりません")
        return

    Y_obs = np.load(data_path)
    print(f"データ形状: {Y_obs.shape}, サイズ: {Y_obs.nbytes / 1e6:.2f} MB")

    # パラメータ設定
    nt = Y_obs.shape[0]
    window_size = 71
    overlap = 17
    cgm_iteration = 1  # 最小限の反復数
    q_init_value = 0.0
    dt = 0.001

    print(f"\nパラメータ: nt={nt}, window={window_size}, overlap={overlap}, CGM={cgm_iteration}")

    # 初期温度場の設定
    T_measure_init_K = Y_obs[0, :, :]
    T0 = np.repeat(T_measure_init_K[:, :, np.newaxis], nz, axis=2).astype(np.float64)

    print(f"初期温度場: {T0.shape}, T_range=[{np.min(T0):.2f}, {np.max(T0):.2f}]K")

    # 計算実行
    print("\n" + "=" * 80)
    start_time = time.time()

    q_result, windows_info = sliding_window_CGM_simple(
        Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt,
        rho, cp_coeffs, k_coeffs,
        window_size, overlap, q_init_value, cgm_iteration
    )

    elapsed = time.time() - start_time

    print("=" * 80)
    print(f"総実行時間: {elapsed:.2f}秒")
    print(f"1時間ステップあたり: {elapsed/nt*1000:.1f}ms")

    # 結果保存
    output_path = PROJECT_ROOT / "shared/results/python_test_1min.npz"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    windows_count = len(windows_info)
    J_final = windows_info[-1]['J_final'] if windows_info else 0.0

    np.savez(
        output_path,
        q_result=q_result,
        elapsed_time=elapsed,
        windows_count=windows_count,
        J_final=J_final,
        # パラメータ
        nt=nt,
        window_size=window_size,
        overlap=overlap,
        cgm_iteration=cgm_iteration,
        dt=dt,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_b=dz_b,
        dz_t=dz_t,
        rho=rho,
        cp_coeffs=cp_coeffs,
        k_coeffs=k_coeffs
    )

    print(f"\n結果保存: {output_path}")
    print(f"  q_result: {q_result.shape}")
    print(f"  elapsed: {elapsed:.2f}s")
    print(f"  windows: {windows_count}")
    print(f"  J_final: {J_final:.6e}")
    print(f"  q_range: [{np.min(q_result):.3e}, {np.max(q_result):.3e}] W/m²")
    print(f"  q_mean: {np.mean(q_result):.3e} W/m²")
    print(f"  q_std: {np.std(q_result):.3e} W/m²")

    print("\n" + "=" * 80)
    print("完了")
    print("=" * 80)


if __name__ == "__main__":
    main()
