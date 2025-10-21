#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python版スライディングウィンドウCGM実行スクリプト
Julia版run_sliding_window.jlと同等のコマンドライン引数対応
"""

import sys
import os
import time
import argparse
import numpy as np
from pathlib import Path

# プロジェクトルートをパスに追加
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / "python" / "original"))

# オリジナルコードをインポート
from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    sliding_window_CGM_q_saving,
    dx, dy, dz, dz_b, dz_t,
    rho_coeffs, cp_coeffs, k_coeffs,
    nz
)

def main():
    parser = argparse.ArgumentParser(description='Python版スライディングウィンドウCGM')
    parser.add_argument('--nt', type=int, default=10, help='時間ステップ数 (default: 10)')
    parser.add_argument('--cgm-iter', type=int, default=3, help='CGM最大反復数 (default: 3)')
    parser.add_argument('--window', type=int, default=71, help='ウィンドウサイズ (default: 71)')
    parser.add_argument('--overlap', type=int, default=17, help='オーバーラップ (default: 17)')
    parser.add_argument('--output', type=str, default='python_sliding_window',
                        help='出力ファイル名プレフィックス (default: python_sliding_window)')
    parser.add_argument('--dry-run', action='store_true',
                        help='ウィンドウ分割のみ表示して終了（計算は実行しない）')

    args = parser.parse_args()

    print("=" * 80)
    print("Python版スライディングウィンドウCGM")
    print("=" * 80)
    print(f"\n[Configuration]")
    print(f"  CGM max iter: {args.cgm_iter}")
    print(f"  nt: {args.nt}")
    print(f"  window size: {args.window}")
    print(f"  overlap: {args.overlap}")
    print(f"  dt: 1.000e-03 s")

    # データ読み込み
    data_path = project_root / "shared" / "data" / "T_measure_700um_1ms.npy"
    print(f"\n[1/4] Loading measurement data")
    print(f"  loading: {data_path}")
    T_measure_K = np.load(data_path)
    print(f"  full data shape: {T_measure_K.shape}")

    # 指定ステップ数のみ使用
    Y_obs = T_measure_K[:args.nt, :, :]
    ni, nj = Y_obs.shape[1], Y_obs.shape[2]
    nk = nz

    print(f"  measurement subset: ({args.nt}, {ni}, {nj})")
    print(f"  grid: ni={ni}, nj={nj}, nk={nk}")

    # 初期温度場
    T_measure_init_K = T_measure_K[0, :, :]
    T0 = np.repeat(T_measure_init_K[:, :, np.newaxis], nk, axis=2).astype(np.float64)

    print(f"  initial temperature range: {T0.min():.2f} ~ {T0.max():.2f} K")

    # タイムステップ
    dt = 0.001

    # rho値計算
    T_avg = np.mean(T0)
    rho_val = np.polyval(rho_coeffs, T_avg)
    print(f"  rho: {rho_val:.6f} kg/m^3")

    # ウィンドウ分割の計算（Juliaと同じロジック）
    print(f"\n[2/4] Window configuration")
    total_flux_steps = args.nt - 1
    print(f"  total time steps: {args.nt}")
    print(f"  flux steps: {total_flux_steps}")
    print(f"  window size: {args.window}")
    print(f"  overlap: {args.overlap}")

    # ウィンドウ数の計算（Julia版と同じロジック）
    windows = []
    start_idx = 0
    window_id = 1
    while start_idx < total_flux_steps:
        max_L = min(args.window, total_flux_steps - start_idx)
        end_idx = start_idx + max_L
        windows.append({
            'id': window_id,
            'start': start_idx,
            'end': end_idx,
            'length': max_L
        })
        print(f"  [Window {window_id}] range=[{start_idx},{end_idx}] length={max_L}")
        window_id += 1

        # 次のウィンドウの開始位置（Julia版と同じ）
        step = max(1, max_L - args.overlap)
        start_idx += step

    print(f"\n  Total windows: {len(windows)}")

    # ドライラン時はここで終了
    if args.dry_run:
        print("\n[Dry-run mode] ウィンドウ分割の確認が完了しました。計算は実行しません。")
        print("=" * 80)
        return 0

    # 出力ファイル名
    output_dir = project_root / "shared" / "results"
    temp_filename = str(output_dir / args.output)

    # 実行
    print(f"\n[3/4] Running sliding-window CGM")
    print("=" * 80)

    start_time = time.time()

    try:
        q_result = sliding_window_CGM_q_saving(
            Y_obs=Y_obs,
            T0=T0,
            dx=dx,
            dy=dy,
            dz=dz,
            dz_b=dz_b,
            dz_t=dz_t,
            dt=dt,
            rho=rho_val,
            cp_coeffs=cp_coeffs,
            k_coeffs=k_coeffs,
            window_size=args.window,
            overlap=args.overlap,
            q_init_value=0.0,
            filename=temp_filename,
            CGM_iteration=args.cgm_iter,
            save_intermediate=False,
            output_dir=str(output_dir)
        )

        elapsed_time = time.time() - start_time

        print("=" * 80)
        print("\n[4/4] Saving outputs")
        print(f"  Total runtime: {elapsed_time:.2f} s")
        print(f"  Heat-flux shape: {q_result.shape}")
        print(f"  Heat-flux range: {q_result.min():.6e} ~ {q_result.max():.6e} W/m^2")

        # 結果保存
        output_path = output_dir / f"{args.output}.npz"
        np.savez(
            output_path,
            q_result=q_result,
            T0=T0,
            Y_obs=Y_obs,
            elapsed_time=elapsed_time,
            cgm_max_iter=args.cgm_iter,
            window_size=args.window,
            overlap=args.overlap,
            nt=args.nt,
            ni=ni,
            nj=nj,
            nk=nk,
            dt=dt,
            dx=dx,
            dy=dy
        )
        print(f"  saved (.npz) → {output_path}")

        # メタデータ保存
        metadata_path = output_dir / f"{args.output}_metadata.txt"
        with open(metadata_path, 'w') as f:
            f.write(f"cgm_max_iter: {args.cgm_iter}\n")
            f.write(f"nt: {args.nt}\n")
            f.write(f"window_size: {args.window}\n")
            f.write(f"overlap: {args.overlap}\n")
            f.write(f"dt: {dt:e}\n")
            f.write(f"q_min: {q_result.min():e}\n")
            f.write(f"q_max: {q_result.max():e}\n")
            f.write(f"total_elapsed: {elapsed_time:.2f}\n")
        print(f"  saved metadata → {metadata_path}")

        print("=" * 80)
        return 0

    except Exception as e:
        print(f"\n❌ エラー発生: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
