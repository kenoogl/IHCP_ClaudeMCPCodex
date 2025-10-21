#!/usr/bin/env python3
"""
オリジナルコードのドライラン機能をテスト
"""

import sys
import numpy as np
from pathlib import Path

# プロジェクトルートをパスに追加
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / "python" / "original"))

from IHCP_CGM_Sliding_Window_Calculation_ver2 import sliding_window_CGM_q_saving

print("=" * 80)
print("オリジナルコード ドライランテスト")
print("=" * 80)

# テストケース1: nt=10, window=5, overlap=2
print("\n=== Test Case 1: nt=10, window=5, overlap=2 ===")

# ダミーデータ作成
nt = 10
ni, nj, nk = 80, 100, 20
Y_obs = np.random.rand(nt, ni, nj) * 50 + 550  # 温度データ
T0 = np.random.rand(ni, nj, nk) * 50 + 550

# その他のパラメータ（ダミー）
dx = 0.00012
dy = 0.0001671
dz = np.linspace(0, 0.0005, nk+1)
dz_diff = np.diff(dz)
dz_b = np.append(np.inf, (dz[:-1] + dz[1:]) / 2)
dz_t = np.append((dz[:-1] + dz[1:]) / 2, np.inf)
dt = 0.001
rho = 7800
cp_coeffs = [1, 2, 3]
k_coeffs = [1, 2, 3]

# ドライラン実行
result = sliding_window_CGM_q_saving(
    Y_obs=Y_obs,
    T0=T0,
    dx=dx,
    dy=dy,
    dz=dz_diff,
    dz_b=dz_b,
    dz_t=dz_t,
    dt=dt,
    rho=rho,
    cp_coeffs=cp_coeffs,
    k_coeffs=k_coeffs,
    window_size=5,
    overlap=2,
    q_init_value=0.0,
    filename="dummy.npy",
    CGM_iteration=3,
    dry_run=True
)

print("\n" + "=" * 80)
print("テスト完了")
print("=" * 80)
