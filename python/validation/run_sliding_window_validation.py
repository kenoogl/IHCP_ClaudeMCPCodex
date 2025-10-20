#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pythonスライディングウィンドウ検証スクリプト

オリジナルのPythonコードでスライディングウィンドウ計算を実行し、
詳細な診断情報を含む結果を保存する。

実行方法:
  python python/validation/run_sliding_window_validation.py --cgm-iter 3
  python python/validation/run_sliding_window_validation.py --cgm-iter 20000

出力:
  shared/results/python_sliding_window_cgm{N}.npz
"""

import numpy as np
import time
import sys
import argparse
from pathlib import Path
import warnings

# 警告を抑制
warnings.filterwarnings('ignore')

# プロジェクトルート設定
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "python" / "original"))
sys.path.insert(0, str(PROJECT_ROOT / "python" / "validation"))

# オリジナルコードから必要な関数・変数をインポート
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

# 診断機能付きソルバーをインポート
from solvers_with_diagnostics import global_CGM_time_with_diagnostics


def sliding_window_CGM_with_diagnostics(
    Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
    window_size, overlap, q_init_value, CGM_iteration=3, verbose=True
):
    """
    診断機能付きスライディングウィンドウCGM計算

    Args:
        Y_obs: 観測温度 (nt, ni, nj)
        T0: 初期温度場 (ni, nj, nk)
        dx, dy, dz, dz_b, dz_t: 格子パラメータ
        dt: 時間刻み
        rho, cp_coeffs, k_coeffs: 熱物性値
        window_size: ウィンドウサイズ
        overlap: オーバーラップ
        q_init_value: 初期熱流束値
        CGM_iteration: CGM最大反復数
        verbose: 詳細出力フラグ

    Returns:
        q_global: 全時間の熱流束 (nt-1, ni, nj)
        T_final: 最終温度場 (ni, nj, nk)
        windows_info: 各ウィンドウの診断情報リスト
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

    if verbose:
        print("\n" + "=" * 80)
        print("Python スライディングウィンドウCGM計算")
        print("=" * 80)
        print(f"総時間ステップ: {nt}")
        print(f"ウィンドウサイズ: {window_size}")
        print(f"オーバーラップ: {overlap}")
        print(f"CGM最大反復数: {CGM_iteration}")
        print(f"格子: {ni}×{nj}×{nk}")
        print("=" * 80)

    while start_idx < nt - 1:
        safety_counter += 1
        if safety_counter > safety_limit:
            print("[警告] 安全カウンタ超過")
            break

        # 現在のウィンドウ長
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]

        window_id = len(windows_info) + 1

        if verbose:
            print(f"\n{'='*80}")
            print(f"ウィンドウ {window_id}: [{start_idx}, {end_idx}] (長さ={max_L})")
            print(f"{'='*80}")

        # 初期熱流束の設定
        if prev_q_win is None:
            q_init_win = np.full((max_L, ni, nj), q_init_value, dtype=float)
            if verbose:
                print(f"初期熱流束: 一定値 {q_init_value} W/m²")
        else:
            q_init_win = np.empty((max_L, ni, nj), dtype=float)
            L_overlap = min(overlap, max_L, prev_q_win.shape[0])
            if L_overlap > 0:
                q_init_win[:L_overlap] = prev_q_win[-L_overlap:]
            if L_overlap < max_L:
                edge = prev_q_win[-1]
                q_init_win[L_overlap:] = edge
            if verbose:
                print(f"初期熱流束: 前ウィンドウから継承 (オーバーラップ={L_overlap})")

        # CGM計算実行
        start_time_win = time.time()

        q_win, T_win_last, J_hist, cgm_diag = global_CGM_time_with_diagnostics(
            T_init, Y_obs_win, q_init_win, dx, dy, dz, dz_b, dz_t, dt,
            rho, cp_coeffs, k_coeffs, CGM_iteration,
            thermal_properties_calculator=thermal_properties_calculator,
            coeffs_and_rhs_building_DHCP=coeffs_and_rhs_building_DHCP,
            assemble_A_DHCP=assemble_A_DHCP,
            coeffs_and_rhs_building_Adjoint=coeffs_and_rhs_building_Adjoint,
            assemble_A_Adjoint=assemble_A_Adjoint,
            verbose=verbose
        )

        elapsed_win = time.time() - start_time_win

        # ウィンドウ診断情報を保存
        win_info = {
            'window_id': window_id,
            'start_idx': start_idx,
            'end_idx': end_idx,
            'max_L': max_L,
            'cgm_iterations': len(J_hist),
            'J_final': J_hist[-1] if len(J_hist) > 0 else 0.0,
            'J_history': J_hist,
            'elapsed_window': elapsed_win,
            'cgm_diagnostics': cgm_diag,
            'q_range': [float(np.min(q_win)), float(np.max(q_win))],
            'q_mean': float(np.mean(q_win)),
            'q_std': float(np.std(q_win))
        }
        windows_info.append(win_info)

        if verbose:
            print(f"\nウィンドウ {window_id} 完了:")
            print(f"  CGM反復数: {len(J_hist)}")
            print(f"  最終J: {J_hist[-1]:.6e}")
            print(f"  実行時間: {elapsed_win:.2f}秒")
            print(f"  q範囲: [{np.min(q_win):.3e}, {np.max(q_win):.3e}] W/m²")

        prev_q_win = q_win.copy()

        # 熱流束の連結（オーバーラップ部分は平均化）
        if len(q_total) == 0:
            q_total.append(q_win)
        else:
            overlap_steps = min(overlap, q_win.shape[0], q_total[-1].shape[0])
            if overlap_steps > 0:
                # オリジナルと同じ平均化
                q_total[-1][-overlap_steps:] = 0.5 * q_total[-1][-overlap_steps:] + 0.5 * q_win[:overlap_steps]
                q_total.append(q_win[overlap_steps:])
            else:
                q_total.append(q_win)

        # 温度場の継承
        if T_win_last.ndim == 3:
            T_init = T_win_last.copy()
        else:
            T_init = T_win_last[-1].copy()

        # 次のウィンドウへ
        step = max(1, max_L - overlap)
        start_idx += step

    # 全時間の熱流束を連結
    q_global = np.concatenate(q_total, axis=0)[:nt-1]

    if verbose:
        print(f"\n{'='*80}")
        print("スライディングウィンドウ計算完了")
        print(f"{'='*80}")
        print(f"総ウィンドウ数: {len(windows_info)}")
        print(f"q_global形状: {q_global.shape}")
        print(f"q範囲: [{np.min(q_global):.3e}, {np.max(q_global):.3e}] W/m²")

    return q_global, T_init, windows_info


def main():
    """メイン実行関数"""
    parser = argparse.ArgumentParser(description='Pythonスライディングウィンドウ検証')
    parser.add_argument('--cgm-iter', type=int, default=3,
                        help='CGM最大反復数 (デフォルト: 3)')
    parser.add_argument('--nt', type=int, default=300,
                        help='時間ステップ数 (デフォルト: 300)')
    parser.add_argument('--window', type=int, default=71,
                        help='ウィンドウサイズ (デフォルト: 71)')
    parser.add_argument('--overlap', type=int, default=17,
                        help='オーバーラップ (デフォルト: 17)')
    parser.add_argument('--quiet', action='store_true',
                        help='詳細出力を抑制')

    args = parser.parse_args()

    print("=" * 80)
    print("Pythonスライディングウィンドウ検証スクリプト")
    print("=" * 80)
    print(f"CGM反復数: {args.cgm_iter}")
    print(f"時間ステップ数: {args.nt}")
    print(f"ウィンドウサイズ: {args.window}")
    print(f"オーバーラップ: {args.overlap}")
    print("=" * 80)

    # データファイル読み込み
    data_path = PROJECT_ROOT / "shared/data/T_measure_700um_1ms.npy"
    if not data_path.exists():
        print(f"❌ エラー: データファイルが見つかりません: {data_path}")
        return 1

    print(f"\nデータ読み込み中: {data_path}")
    T_measure_full = np.load(data_path)
    print(f"データ形状: {T_measure_full.shape}")

    # 時間ステップを切り出し
    nt = min(args.nt, T_measure_full.shape[0])
    Y_obs = T_measure_full[:nt, :, :]
    print(f"使用データ: {Y_obs.shape} (最初{nt}ステップ)")

    # パラメータ設定
    dt = 0.001
    q_init_value = 0.0

    # 初期温度場の設定
    T_measure_init_K = Y_obs[0, :, :]
    T0 = np.repeat(T_measure_init_K[:, :, np.newaxis], nz, axis=2).astype(np.float64)

    print(f"\n初期温度場: {T0.shape}")
    print(f"温度範囲: [{np.min(T0):.2f}, {np.max(T0):.2f}] K")

    # 計算実行
    print("\n" + "=" * 80)
    print("計算開始")
    print("=" * 80)

    start_time = time.time()

    q_result, T_final, windows_info = sliding_window_CGM_with_diagnostics(
        Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt,
        rho, cp_coeffs, k_coeffs,
        args.window, args.overlap, q_init_value, args.cgm_iter,
        verbose=not args.quiet
    )

    elapsed_total = time.time() - start_time

    print("\n" + "=" * 80)
    print("計算完了")
    print("=" * 80)
    print(f"総実行時間: {elapsed_total:.2f}秒")
    print(f"1ウィンドウあたり: {elapsed_total/len(windows_info):.2f}秒")

    # 結果保存
    output_path = PROJECT_ROOT / f"shared/results/python_sliding_window_cgm{args.cgm_iter}.npz"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"\n結果保存中: {output_path}")

    # 診断情報を保存用に整理
    windows_summary = []
    for win_info in windows_info:
        # cgm_diagnosticsからDHCP/Adjoint情報を抽出
        dhcp_iters = []
        adjoint_iters = []
        dhcp_times = []
        adjoint_times = []

        for cgm_iter in win_info['cgm_diagnostics']:
            if cgm_iter['dhcp'] is not None:
                dhcp_iters.append(cgm_iter['dhcp']['total_iters'])
                dhcp_times.append(cgm_iter['dhcp']['elapsed_time'])
            if cgm_iter['adjoint'] is not None:
                adjoint_iters.append(cgm_iter['adjoint']['total_iters'])
                adjoint_times.append(cgm_iter['adjoint']['elapsed_time'])

        win_summary = {
            'window_id': win_info['window_id'],
            'start_idx': win_info['start_idx'],
            'end_idx': win_info['end_idx'],
            'max_L': win_info['max_L'],
            'cgm_iterations': win_info['cgm_iterations'],
            'J_final': win_info['J_final'],
            'J_history': win_info['J_history'],
            'elapsed_window': win_info['elapsed_window'],
            'dhcp_total_iters': dhcp_iters,
            'adjoint_total_iters': adjoint_iters,
            'dhcp_times': dhcp_times,
            'adjoint_times': adjoint_times,
            'q_range': win_info['q_range'],
            'q_mean': win_info['q_mean'],
            'q_std': win_info['q_std']
        }
        windows_summary.append(win_summary)

    np.savez(
        output_path,
        q_global=q_result,
        T_final=T_final,
        elapsed_total=elapsed_total,
        windows_info=windows_summary,
        # パラメータ
        nt=nt,
        window_size=args.window,
        overlap=args.overlap,
        cgm_max_iter=args.cgm_iter,
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

    print(f"✅ 保存完了")
    print(f"\n結果サマリー:")
    print(f"  q_global: {q_result.shape}")
    print(f"  q範囲: [{np.min(q_result):.3e}, {np.max(q_result):.3e}] W/m²")
    print(f"  q平均: {np.mean(q_result):.3e} W/m²")
    print(f"  q標準偏差: {np.std(q_result):.3e} W/m²")
    print(f"  総ウィンドウ数: {len(windows_info)}")
    print(f"  総CGM反復数: {sum(w['cgm_iterations'] for w in windows_info)}")

    # ウィンドウ別サマリー
    print(f"\nウィンドウ別サマリー:")
    print(f"{'ID':>3} {'範囲':>15} {'CGM':>4} {'DHCP合計':>10} {'Adj合計':>10} {'時間[s]':>8}")
    print("-" * 70)
    for win in windows_summary:
        dhcp_sum = sum(win['dhcp_total_iters']) if win['dhcp_total_iters'] else 0
        adj_sum = sum(win['adjoint_total_iters']) if win['adjoint_total_iters'] else 0
        print(f"{win['window_id']:3d} "
              f"[{win['start_idx']:3d},{win['end_idx']:3d}] "
              f"{win['cgm_iterations']:4d} "
              f"{dhcp_sum:10d} "
              f"{adj_sum:10d} "
              f"{win['elapsed_window']:8.2f}")

    print("\n" + "=" * 80)
    print("完了")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
