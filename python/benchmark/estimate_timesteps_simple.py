#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
時間ステップ数推定スクリプト（簡易版）

Python版（Numba並列化）で約1分の実行時間となる適切な時間ステップ数を推定する
警告メッセージを抑制してクリーンな出力を提供
"""

import numpy as np
import time
from pathlib import Path
import sys
import warnings
from math import sin, radians

# 警告を抑制
warnings.filterwarnings('ignore')

# プロジェクトルート設定
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(PROJECT_ROOT / "python" / "original"))

# 標準出力を一時的に抑制するためのコンテキストマネージャ
import os
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

# オリジナルコードから必要な関数をインポート（警告を抑制）
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


def benchmark_timesteps(nt_list=[10, 20, 50, 100, 200]):
    """異なる時間ステップ数での実行時間を計測"""

    print("=" * 80)
    print("時間ステップ数推定ベンチマーク")
    print("=" * 80)

    # 実データ読み込み
    data_path = PROJECT_ROOT / "shared/data/T_measure_700um_1ms.npy"
    print(f"\nデータファイル: {data_path}")

    if not data_path.exists():
        print(f"エラー: データファイルが見つかりません")
        return None

    T_measure_full = np.load(data_path)

    print(f"実データ形状: {T_measure_full.shape}")
    print(f"データサイズ: {T_measure_full.nbytes / 1e9:.3f} GB")
    print(f"時間ステップ総数: {T_measure_full.shape[0]}")
    print(f"空間格子サイズ (ni x nj): {T_measure_full.shape[1]} x {T_measure_full.shape[2]}")

    # 時間刻み
    dt = 0.001  # 1ms

    # 初期温度場の設定
    T_measure_init_K = T_measure_full[0, :, :]
    T0 = np.repeat(T_measure_init_K[:, :, np.newaxis], nz, axis=2).astype(np.float64)

    results = []

    print("\n" + "=" * 80)
    print("ベンチマーク開始")
    print("=" * 80)

    for nt in nt_list:
        # データの時間方向を縮小
        Y_obs = T_measure_full[:nt, :, :]
        ni, nj = Y_obs.shape[1], Y_obs.shape[2]
        nk = T0.shape[2]

        # 初期熱流束（全てゼロ）
        q_init = np.zeros((nt - 1, ni, nj))

        print(f"\n時間ステップ数: nt = {nt:4d}  ", end='', flush=True)

        # CGM計算実行（反復回数は1回、警告抑制）
        with suppress_stdout():
            start_time = time.time()
            _, _, _ = global_CGM_time(
                T0, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, dt,
                rho, cp_coeffs, k_coeffs, CGM_iteration=1
            )
            elapsed = time.time() - start_time

        print(f"実行時間: {elapsed:6.2f}秒  ", end='')
        print(f"1時間ステップあたり: {elapsed/nt*1000:6.1f}ms", flush=True)

        # 結果を保存
        results.append({
            'nt': nt,
            'time': elapsed,
            'time_per_step': elapsed / nt
        })

        # 1分を大幅に超えたら終了
        if elapsed > 90:
            print(f"\n実行時間が90秒を超えたため、ベンチマークを終了します")
            break

    # 結果分析
    print(f"\n{'='*80}")
    print("ベンチマーク結果サマリー")
    print("=" * 80)
    print(f"{'nt':>6} {'実行時間':>12} {'1ステップ':>12} {'60秒比':>10}")
    print("-" * 80)

    for r in results:
        ratio = r['time'] / 60.0
        print(f"{r['nt']:6d} {r['time']:10.2f}秒 "
              f"{r['time_per_step']*1000:9.1f}ms {ratio:9.1%}")

    # 1分となる時間ステップ数を推定
    if len(results) >= 2:
        # 線形回帰で推定
        nt_array = np.array([r['nt'] for r in results])
        time_array = np.array([r['time'] for r in results])

        # 最小二乗法
        coeffs = np.polyfit(nt_array, time_array, 1)
        slope, intercept = coeffs[0], coeffs[1]

        # 60秒となるntを計算
        estimated_nt = int((60.0 - intercept) / slope)

        print(f"\n{'='*80}")
        print("推定結果")
        print("=" * 80)
        print(f"線形近似: 実行時間 = {slope:.6f} * nt + {intercept:.4f}")
        print(f"推定: 約1分（60秒）で実行可能な時間ステップ数 = {estimated_nt}")

        # 安全マージンを考慮した推奨値
        recommended_nt = int(estimated_nt * 0.9)  # 10%のマージン
        print(f"推奨値（10%マージン）: nt = {recommended_nt}")

        # ウィンドウサイズとオーバーラップの推奨値
        print(f"\n{'='*80}")
        print("スライディングウィンドウパラメータ推奨値")
        print("=" * 80)
        print(f"時間ステップ数: nt = {recommended_nt}")
        print(f"ウィンドウサイズ: window_size = {min(100, recommended_nt // 5)}")
        print(f"オーバーラップ: overlap = {min(20, recommended_nt // 20)}")

        return {
            'estimated_nt': estimated_nt,
            'recommended_nt': recommended_nt,
            'window_size': min(100, recommended_nt // 5),
            'overlap': min(20, recommended_nt // 20),
            'slope': slope,
            'intercept': intercept,
            'results': results
        }

    return None


def generate_test_data(recommended_nt):
    """推奨時間ステップ数に基づいてテストデータを作成"""

    print(f"\n{'='*80}")
    print("テストデータ作成")
    print("=" * 80)

    data_path = PROJECT_ROOT / "shared/data/T_measure_700um_1ms.npy"
    T_measure_full = np.load(data_path)

    # 推奨時間ステップ数でデータを切り出し
    T_measure_test = T_measure_full[:recommended_nt, :, :]

    # テストデータを保存
    output_path = PROJECT_ROOT / "shared/data/T_measure_test_1min.npy"
    np.save(output_path, T_measure_test)

    print(f"テストデータ作成完了")
    print(f"保存先: {output_path}")
    print(f"形状: {T_measure_test.shape}")
    print(f"サイズ: {T_measure_test.nbytes / 1e6:.2f} MB")

    return output_path


if __name__ == "__main__":
    # ベンチマーク実行
    result = benchmark_timesteps()

    if result:
        # テストデータ作成
        test_data_path = generate_test_data(result['recommended_nt'])

        print(f"\n{'='*80}")
        print("完了")
        print("=" * 80)
        print(f"推奨時間ステップ数: {result['recommended_nt']}")
        print(f"テストデータ: {test_data_path}")
        print("=" * 80)
    else:
        print("\nベンチマーク失敗")
