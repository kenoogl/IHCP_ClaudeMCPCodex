#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python版縮小テストデータ実行スクリプト（1分想定）

Julia版との比較検証用にPython版を実行し、結果を保存する

実行方法:
  cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/python/benchmark
  python run_test_1min.py

出力:
  - shared/results/python_test_1min.npz: 計算結果
    - q_result: 熱流束結果 (nt-1, ni, nj)
    - elapsed_time: 実行時間 [秒]
    - J_final: 最終目的関数値
    - windows_count: ウィンドウ数
"""

import numpy as np
import time
import sys
from pathlib import Path
import warnings

# 警告を抑制
warnings.filterwarnings('ignore')

# プロジェクトルート設定
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(PROJECT_ROOT / "python" / "original"))

# オリジナルコードから必要な関数をインポート
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


def sliding_window_CGM_q_saving_for_benchmark(
    Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, CGM_iteration=10
):
    """
    スライディングウィンドウCGM計算（ベンチマーク用）

    オリジナルのsliding_window_CGM_q_saving()をベースに、
    結果保存機能を追加してベンチマークに適した形式に変更
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

    print("\n" + "=" * 80)
    print("Python版スライディングウィンドウCGM計算開始")
    print("=" * 80)
    print(f"全時間ステップ数: {nt}")
    print(f"ウィンドウサイズ: {window_size}")
    print(f"オーバーラップ: {overlap}")
    print(f"CGM反復数: {CGM_iteration}")
    print(f"空間格子サイズ: {ni} × {nj} × {nk}")

    while start_idx < nt - 1:
        safety_counter += 1
        if safety_counter > safety_limit:
            print("[警告] 安全カウンタ超過")
            break

        # 現在のウィンドウ長
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]

        print(f"\n--- ウィンドウ {len(windows_info)+1}: [{start_idx}, {end_idx}] (長さ={max_L}) ---")

        # 初期熱流束の設定
        if prev_q_win is None:
            q_init_win = np.full((max_L, ni, nj), q_init_value, dtype=float)
            print(f"初期熱流束: 一定値 {q_init_value} W/m²")
        else:
            q_init_win = np.empty((max_L, ni, nj), dtype=float)
            L_overlap = min(overlap, max_L, prev_q_win.shape[0])
            if L_overlap > 0:
                q_init_win[:L_overlap] = prev_q_win[-L_overlap:]
                print(f"オーバーラップ継承: 前{L_overlap}ステップ")
            if L_overlap < max_L:
                edge = prev_q_win[-1]
                q_init_win[L_overlap:] = edge

        # CGM計算実行
        start_time_win = time.time()
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

        print(f"計算完了: J = {J_hist[-1]:.6e}, 実行時間 = {elapsed_win:.2f}秒")

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

    print("\n" + "=" * 80)
    print("Python版計算完了")
    print("=" * 80)
    print(f"熱流束形状: {q_global.shape}")
    print(f"ウィンドウ数: {len(windows_info)}")

    return q_global, windows_info


def main():
    """メイン実行関数"""

    print("=" * 80)
    print("Python版縮小テストデータ実行（1分想定）")
    print("=" * 80)

    # テストデータ読み込み
    data_path = PROJECT_ROOT / "shared/data/T_measure_test_1min.npy"
    print(f"\nデータファイル: {data_path}")

    if not data_path.exists():
        print(f"エラー: データファイルが見つかりません")
        return

    Y_obs = np.load(data_path)
    print(f"データ形状: {Y_obs.shape}")
    print(f"データサイズ: {Y_obs.nbytes / 1e6:.2f} MB")

    # パラメータ設定（load_test_data.jlと同じ）
    nt = Y_obs.shape[0]
    window_size = 71
    overlap = 17
    cgm_iteration = 10  # 1分想定での適切な反復数
    q_init_value = 0.0  # 初期熱流束
    dt = 0.001  # 1ms

    print(f"\nパラメータ設定:")
    print(f"  時間ステップ数: nt = {nt}")
    print(f"  ウィンドウサイズ: {window_size}")
    print(f"  オーバーラップ: {overlap}")
    print(f"  CGM反復数: {cgm_iteration}")
    print(f"  初期熱流束: {q_init_value} W/m²")
    print(f"  時間刻み: {dt*1000} ms")
    print(f"  総時間: {nt*dt} s = {nt*dt*1000} ms")

    # 初期温度場の設定
    T_measure_init_K = Y_obs[0, :, :]
    T0 = np.repeat(T_measure_init_K[:, :, np.newaxis], nz, axis=2).astype(np.float64)

    print(f"\n初期温度場:")
    print(f"  形状: {T0.shape}")
    print(f"  温度範囲: {np.min(T0):.2f}K ~ {np.max(T0):.2f}K")

    # 計算実行
    print("\n" + "=" * 80)
    print("計算開始")
    print("=" * 80)

    start_time = time.time()

    q_result, windows_info = sliding_window_CGM_q_saving_for_benchmark(
        Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt,
        rho, cp_coeffs, k_coeffs,
        window_size, overlap, q_init_value, cgm_iteration
    )

    elapsed = time.time() - start_time

    print("\n" + "=" * 80)
    print("計算完了")
    print("=" * 80)
    print(f"総実行時間: {elapsed:.2f}秒")
    print(f"1時間ステップあたり: {elapsed/nt*1000:.1f}ms")

    # 結果保存
    output_path = PROJECT_ROOT / "shared/results/python_test_1min.npz"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # ウィンドウ情報を配列に変換
    windows_count = len(windows_info)
    J_final = windows_info[-1]['J_final'] if windows_info else 0.0

    np.savez(
        output_path,
        q_result=q_result,
        elapsed_time=elapsed,
        windows_count=windows_count,
        J_final=J_final,
        # パラメータも保存
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
    print(f"  熱流束形状: {q_result.shape}")
    print(f"  実行時間: {elapsed:.2f}秒")
    print(f"  ウィンドウ数: {windows_count}")
    print(f"  最終J値: {J_final:.6e}")

    # 統計情報
    print(f"\n熱流束統計:")
    print(f"  最小値: {np.min(q_result):.6e} W/m²")
    print(f"  最大値: {np.max(q_result):.6e} W/m²")
    print(f"  平均値: {np.mean(q_result):.6e} W/m²")
    print(f"  標準偏差: {np.std(q_result):.6e} W/m²")

    print("\n" + "=" * 80)
    print("完了")
    print("=" * 80)


if __name__ == "__main__":
    main()
