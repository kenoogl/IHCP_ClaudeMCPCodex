#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python版結果分析スクリプト

Python版の計算結果を詳細に分析し、可視化とレポート出力を行う
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def analyze_python_results():
    """Python版結果の詳細分析"""

    print("=" * 80)
    print("Python版結果分析")
    print("=" * 80)

    # 結果読み込み
    result_path = PROJECT_ROOT / "shared/results/python_test_1min.npz"
    if not result_path.exists():
        print(f"エラー: 結果ファイルが見つかりません: {result_path}")
        return

    data = np.load(result_path)

    # データ抽出
    q_result = data['q_result']
    elapsed_time = float(data['elapsed_time'])
    windows_count = int(data['windows_count'])
    J_final = float(data['J_final'])

    # パラメータ抽出
    nt = int(data['nt'])
    window_size = int(data['window_size'])
    overlap = int(data['overlap'])
    cgm_iteration = int(data['cgm_iteration'])
    dt = float(data['dt'])

    print(f"\n基本情報:")
    print(f"  データファイル: {result_path}")
    print(f"  熱流束形状: {q_result.shape}")
    print(f"  データサイズ: {q_result.nbytes / 1e6:.2f} MB")

    print(f"\n計算パラメータ:")
    print(f"  時間ステップ数: {nt}")
    print(f"  ウィンドウサイズ: {window_size}")
    print(f"  オーバーラップ: {overlap}")
    print(f"  CGM反復数: {cgm_iteration}")
    print(f"  時間刻み: {dt*1000} ms")
    print(f"  総時間: {nt*dt} s")

    print(f"\n計算結果:")
    print(f"  実行時間: {elapsed_time:.2f} 秒")
    print(f"  1時間ステップあたり: {elapsed_time/nt*1000:.1f} ms")
    print(f"  ウィンドウ数: {windows_count}")
    print(f"  最終目的関数値: {J_final:.6e}")

    # 統計情報
    print(f"\n熱流束統計:")
    print(f"  最小値: {np.min(q_result):.6e} W/m²")
    print(f"  最大値: {np.max(q_result):.6e} W/m²")
    print(f"  平均値: {np.mean(q_result):.6e} W/m²")
    print(f"  中央値: {np.median(q_result):.6e} W/m²")
    print(f"  標準偏差: {np.std(q_result):.6e} W/m²")
    print(f"  分散: {np.var(q_result):.6e} (W/m²)²")

    # 各時間ステップでの統計
    q_mean_t = np.mean(q_result, axis=(1, 2))
    q_std_t = np.std(q_result, axis=(1, 2))
    q_min_t = np.min(q_result, axis=(1, 2))
    q_max_t = np.max(q_result, axis=(1, 2))

    print(f"\n時系列統計:")
    print(f"  平均値の範囲: [{np.min(q_mean_t):.3e}, {np.max(q_mean_t):.3e}] W/m²")
    print(f"  標準偏差の範囲: [{np.min(q_std_t):.3e}, {np.max(q_std_t):.3e}] W/m²")

    # 空間分布の統計
    q_mean_xy = np.mean(q_result, axis=0)
    print(f"\n空間平均分布:")
    print(f"  最小値: {np.min(q_mean_xy):.6e} W/m²")
    print(f"  最大値: {np.max(q_mean_xy):.6e} W/m²")

    # JSON形式でサマリー保存
    summary = {
        'computation': {
            'elapsed_time': elapsed_time,
            'time_per_step_ms': elapsed_time / nt * 1000,
            'windows_count': windows_count,
            'J_final': J_final
        },
        'parameters': {
            'nt': nt,
            'window_size': window_size,
            'overlap': overlap,
            'cgm_iteration': cgm_iteration,
            'dt': dt
        },
        'statistics': {
            'q_min': float(np.min(q_result)),
            'q_max': float(np.max(q_result)),
            'q_mean': float(np.mean(q_result)),
            'q_median': float(np.median(q_result)),
            'q_std': float(np.std(q_result)),
            'q_var': float(np.var(q_result))
        }
    }

    summary_path = PROJECT_ROOT / "shared/results/python_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nサマリー保存: {summary_path}")

    # 簡易可視化用データ作成
    viz_data = {
        'time_series': {
            't': (np.arange(nt-1) * dt).tolist(),
            'q_mean': q_mean_t.tolist(),
            'q_std': q_std_t.tolist(),
            'q_min': q_min_t.tolist(),
            'q_max': q_max_t.tolist()
        }
    }

    viz_path = PROJECT_ROOT / "shared/results/python_viz_data.json"
    with open(viz_path, 'w') as f:
        json.dump(viz_data, f, indent=2)

    print(f"可視化データ保存: {viz_path}")

    # サンプル点での熱流束値
    print(f"\nサンプル点（中心付近）での熱流束:")
    ni, nj = q_result.shape[1], q_result.shape[2]
    sample_i, sample_j = ni // 2, nj // 2
    q_sample = q_result[:, sample_i, sample_j]

    print(f"  座標: ({sample_i}, {sample_j})")
    print(f"  初期値: {q_sample[0]:.6e} W/m²")
    print(f"  中間値: {q_sample[len(q_sample)//2]:.6e} W/m²")
    print(f"  最終値: {q_sample[-1]:.6e} W/m²")

    print("\n" + "=" * 80)
    print("分析完了")
    print("=" * 80)

    return summary


if __name__ == "__main__":
    analyze_python_results()
