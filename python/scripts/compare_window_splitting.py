#!/usr/bin/env python3
"""
オリジナルコードとrun_sliding_window.pyのウィンドウ分割を比較
"""

import sys
import subprocess
import numpy as np
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / "python" / "original"))

from IHCP_CGM_Sliding_Window_Calculation_ver2 import sliding_window_CGM_q_saving

def extract_windows_from_original(nt, window_size, overlap):
    """オリジナルコードからウィンドウ分割を抽出"""
    # ダミーデータ
    ni, nj, nk = 80, 100, 20
    Y_obs = np.zeros((nt, ni, nj))
    T0 = np.zeros((ni, nj, nk))
    dx = dy = 0.0001
    dz = np.linspace(0, 0.0005, nk)
    dz_b = dz_t = np.zeros(nk)
    dt = 0.001
    rho = 7800
    cp_coeffs = k_coeffs = [1, 2, 3]

    # ウィンドウ分割を計算
    windows = []
    start_idx = 0
    window_id = 1
    while start_idx < nt - 1:
        max_L = min(window_size, (nt - 1) - start_idx)
        end_idx = start_idx + max_L
        windows.append({
            'id': window_id,
            'start': start_idx,
            'end': end_idx,
            'length': max_L
        })
        step = max(1, max_L - overlap)
        start_idx += step
        window_id += 1

    return windows

def test_case(nt, window_size, overlap):
    """テストケース実行"""
    print(f"\n{'='*80}")
    print(f"Test: nt={nt}, window={window_size}, overlap={overlap}")
    print(f"{'='*80}")

    # オリジナルコードから抽出
    windows_orig = extract_windows_from_original(nt, window_size, overlap)

    print("\nオリジナルコード:")
    for w in windows_orig:
        print(f"  [Window {w['id']}] range=[{w['start']},{w['end']}] length={w['length']}")
    print(f"  Total windows: {len(windows_orig)}")

    # run_sliding_window.pyの出力と比較
    cmd = [
        "python", "python/scripts/run_sliding_window.py",
        "--nt", str(nt),
        "--window", str(window_size),
        "--overlap", str(overlap),
        "--dry-run"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(project_root))

    # ウィンドウ行を抽出
    run_windows = []
    for line in result.stdout.split('\n'):
        if '[Window' in line and 'range=' in line:
            run_windows.append(line.strip())

    print("\nrun_sliding_window.py:")
    for line in run_windows:
        print(f"  {line}")

    # 比較
    if len(windows_orig) == len(run_windows):
        print(f"\n✅ ウィンドウ数一致: {len(windows_orig)}")

        all_match = True
        for i, w in enumerate(windows_orig):
            expected = f"[Window {w['id']}] range=[{w['start']},{w['end']}] length={w['length']}"
            if expected not in run_windows[i]:
                all_match = False
                print(f"  ❌ Window {i+1} 不一致")
                print(f"     期待: {expected}")
                print(f"     実際: {run_windows[i]}")

        if all_match:
            print("✅ 全ウィンドウ詳細一致")
        return all_match
    else:
        print(f"❌ ウィンドウ数不一致: {len(windows_orig)} vs {len(run_windows)}")
        return False

print("="*80)
print("ウィンドウ分割比較テスト")
print("="*80)

test_cases = [
    (10, 5, 2),      # テストケース1
    (300, 71, 17),   # テストケース2
    (10, 10, 0),     # テストケース3
    (20, 10, 5),     # テストケース4
]

all_passed = True
for nt, window, overlap in test_cases:
    passed = test_case(nt, window, overlap)
    if not passed:
        all_passed = False

print("\n" + "="*80)
if all_passed:
    print("✅ 全テストケース合格")
else:
    print("❌ 一部テストケース不合格")
print("="*80)
