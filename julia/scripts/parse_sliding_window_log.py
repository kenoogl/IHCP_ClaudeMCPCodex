#!/usr/bin/env python3
"""
スライディングウィンドウログファイルから統計情報を抽出

使用方法:
  python julia/scripts/parse_sliding_window_log.py shared/results/julia_sliding_window_10steps_test.log

出力:
  同じディレクトリに {logfile}_stats.npz を生成
"""

import re
import sys
import numpy as np
from pathlib import Path

def parse_solver_stats(log_content):
    """
    ログからソルバー統計情報を抽出

    Returns:
        dict: {
            'windows': list of window IDs,
            'cgm_iterations': list of CGM iterations per window,
            'dhcp_times': list of DHCP solve times per CGM iteration,
            'dhcp_iters': list of DHCP total iterations per CGM iteration,
            'gradient_times': list of Gradient solve times per CGM iteration,
            'gradient_iters': list of Gradient total iterations per CGM iteration,
            'sensitivity_times': list of Sensitivity solve times per CGM iteration,
            'sensitivity_iters': list of Sensitivity total iterations per CGM iteration
        }
    """

    # パターン定義
    window_pattern = r'\[Window (\d+)\] range=\[(\d+),(\d+)\] length=(\d+)'
    dhcp_pattern = r'DHCP solve time: ([\d.]+) s \(nt=(\d+) steps, total_iters=(\d+), avg=([\d.]+)\)'
    gradient_pattern = r'Gradient solve time: ([\d.]+) s \(nt=(\d+) steps, total_iters=(\d+), avg=([\d.]+)\)'
    sensitivity_pattern = r'Sensitivity solve time: ([\d.]+) s \(nt=(\d+) steps, total_iters=(\d+), avg=([\d.]+)\)'

    windows = []
    window_ranges = []
    window_lengths = []

    dhcp_times = []
    dhcp_iters = []
    dhcp_avg_iters = []

    gradient_times = []
    gradient_iters = []
    gradient_avg_iters = []

    sensitivity_times = []
    sensitivity_iters = []
    sensitivity_avg_iters = []

    current_window = None
    current_window_dhcp = []
    current_window_gradient = []
    current_window_sensitivity = []

    for line in log_content.split('\n'):
        # ウィンドウ開始検出
        window_match = re.search(window_pattern, line)
        if window_match:
            # 前のウィンドウのデータを保存
            if current_window is not None:
                dhcp_times.append(current_window_dhcp.copy())
                gradient_times.append(current_window_gradient.copy())
                sensitivity_times.append(current_window_sensitivity.copy())

            # 新しいウィンドウ
            current_window = int(window_match.group(1))
            windows.append(current_window)
            window_ranges.append((int(window_match.group(2)), int(window_match.group(3))))
            window_lengths.append(int(window_match.group(4)))

            current_window_dhcp = []
            current_window_gradient = []
            current_window_sensitivity = []
            continue

        # DHCP統計
        dhcp_match = re.search(dhcp_pattern, line)
        if dhcp_match and current_window is not None:
            time_val = float(dhcp_match.group(1))
            nt = int(dhcp_match.group(2))
            total_iters = int(dhcp_match.group(3))
            avg = float(dhcp_match.group(4))

            current_window_dhcp.append({
                'time': time_val,
                'nt': nt,
                'total_iters': total_iters,
                'avg': avg
            })
            continue

        # Gradient統計
        gradient_match = re.search(gradient_pattern, line)
        if gradient_match and current_window is not None:
            time_val = float(gradient_match.group(1))
            nt = int(gradient_match.group(2))
            total_iters = int(gradient_match.group(3))
            avg = float(gradient_match.group(4))

            current_window_gradient.append({
                'time': time_val,
                'nt': nt,
                'total_iters': total_iters,
                'avg': avg
            })
            continue

        # Sensitivity統計
        sensitivity_match = re.search(sensitivity_pattern, line)
        if sensitivity_match and current_window is not None:
            time_val = float(sensitivity_match.group(1))
            nt = int(sensitivity_match.group(2))
            total_iters = int(sensitivity_match.group(3))
            avg = float(sensitivity_match.group(4))

            current_window_sensitivity.append({
                'time': time_val,
                'nt': nt,
                'total_iters': total_iters,
                'avg': avg
            })
            continue

    # 最後のウィンドウのデータを保存
    if current_window is not None:
        dhcp_times.append(current_window_dhcp)
        gradient_times.append(current_window_gradient)
        sensitivity_times.append(current_window_sensitivity)

    return {
        'windows': np.array(windows),
        'window_ranges': np.array(window_ranges),
        'window_lengths': np.array(window_lengths),
        'dhcp_stats': dhcp_times,
        'gradient_stats': gradient_times,
        'sensitivity_stats': sensitivity_times
    }


def print_summary(stats):
    """統計情報のサマリーを表示"""
    print(f"\n=== 統計情報サマリー ===")
    print(f"ウィンドウ数: {len(stats['windows'])}")

    for i, window_id in enumerate(stats['windows']):
        print(f"\n[ウィンドウ {window_id}]")
        print(f"  範囲: {stats['window_ranges'][i]}, 長さ: {stats['window_lengths'][i]}")
        print(f"  CGM反復数: {len(stats['dhcp_stats'][i])}")

        if len(stats['dhcp_stats'][i]) > 0:
            dhcp_total_time = sum(s['time'] for s in stats['dhcp_stats'][i])
            dhcp_total_iters = sum(s['total_iters'] for s in stats['dhcp_stats'][i])
            print(f"    DHCP: {dhcp_total_time:.2f}s, {dhcp_total_iters}反復")

        if len(stats['gradient_stats'][i]) > 0:
            grad_total_time = sum(s['time'] for s in stats['gradient_stats'][i])
            grad_total_iters = sum(s['total_iters'] for s in stats['gradient_stats'][i])
            print(f"    Gradient: {grad_total_time:.2f}s, {grad_total_iters}反復")

        if len(stats['sensitivity_stats'][i]) > 0:
            sens_total_time = sum(s['time'] for s in stats['sensitivity_stats'][i])
            sens_total_iters = sum(s['total_iters'] for s in stats['sensitivity_stats'][i])
            print(f"    Sensitivity: {sens_total_time:.2f}s, {sens_total_iters}反復")


def save_stats_npz(stats, output_path):
    """統計情報をNPZファイルに保存"""
    # リスト形式のデータを保存可能な形式に変換
    save_dict = {
        'windows': stats['windows'],
        'window_ranges': stats['window_ranges'],
        'window_lengths': stats['window_lengths']
    }

    # 各ウィンドウの各ソルバー統計を配列に変換
    for i, window_id in enumerate(stats['windows']):
        # DHCP
        if len(stats['dhcp_stats'][i]) > 0:
            save_dict[f'window_{window_id}_dhcp_times'] = np.array([s['time'] for s in stats['dhcp_stats'][i]])
            save_dict[f'window_{window_id}_dhcp_iters'] = np.array([s['total_iters'] for s in stats['dhcp_stats'][i]])
            save_dict[f'window_{window_id}_dhcp_avg'] = np.array([s['avg'] for s in stats['dhcp_stats'][i]])

        # Gradient
        if len(stats['gradient_stats'][i]) > 0:
            save_dict[f'window_{window_id}_gradient_times'] = np.array([s['time'] for s in stats['gradient_stats'][i]])
            save_dict[f'window_{window_id}_gradient_iters'] = np.array([s['total_iters'] for s in stats['gradient_stats'][i]])
            save_dict[f'window_{window_id}_gradient_avg'] = np.array([s['avg'] for s in stats['gradient_stats'][i]])

        # Sensitivity
        if len(stats['sensitivity_stats'][i]) > 0:
            save_dict[f'window_{window_id}_sensitivity_times'] = np.array([s['time'] for s in stats['sensitivity_stats'][i]])
            save_dict[f'window_{window_id}_sensitivity_iters'] = np.array([s['total_iters'] for s in stats['sensitivity_stats'][i]])
            save_dict[f'window_{window_id}_sensitivity_avg'] = np.array([s['avg'] for s in stats['sensitivity_stats'][i]])

    np.savez_compressed(output_path, **save_dict)
    print(f"\n統計情報を保存しました: {output_path}")


def main():
    if len(sys.argv) < 2:
        print("使用方法: python parse_sliding_window_log.py <logfile>")
        sys.exit(1)

    log_file = Path(sys.argv[1])
    if not log_file.exists():
        print(f"エラー: ログファイルが見つかりません: {log_file}")
        sys.exit(1)

    print(f"ログファイルを解析中: {log_file}")

    # ログファイル読み込み
    with open(log_file, 'r') as f:
        log_content = f.read()

    # 統計情報抽出
    stats = parse_solver_stats(log_content)

    # サマリー表示
    print_summary(stats)

    # NPZファイルに保存
    output_path = log_file.parent / f"{log_file.stem}_stats.npz"
    save_stats_npz(stats, output_path)


if __name__ == '__main__':
    main()
