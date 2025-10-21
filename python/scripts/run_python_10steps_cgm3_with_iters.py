#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python版 10ステップ CGM3回 線形ソルバー反復数カウント版

オリジナルのコードにcgコールバックを追加して反復数を記録
"""

import sys
import os
import time
import numpy as np
from pathlib import Path

# プロジェクトルートをパスに追加
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(project_root / "python" / "original"))

print("="*80)
print("Python版 10ステップ CGM3回 (線形ソルバー反復数カウント版)")
print("="*80)
print("\n⚠️ 注意: このスクリプトはオリジナルコードを直接実行します")
print("反復数をカウントするには、オリジナルコードにコールバック追加が必要です\n")

# データ読み込みとパラメータ設定
data_path = project_root / "shared" / "data" / "T_measure_700um_1ms.npy"
print(f"データ読み込み: {data_path}")
T_measure_K = np.load(data_path)

# 10ステップのみ
nt = 10
Y_obs = T_measure_K[:nt, :, :]
ni, nj = Y_obs.shape[1], Y_obs.shape[2]

print(f"使用データ: nt={nt}, ni={ni}, nj={nj}")
print("\n" + "="*80)
print("オリジナルコードの修正方法:")
print("="*80)
print("""
IHCP_CGM_Sliding_Window_Calculation_ver2.py の cg() 呼び出しを修正:

【変更前】
x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)

【変更後】
iter_count = [0]  # 反復数カウンター

def callback(xk):
    iter_count[0] += 1

x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0, callback=callback)

# 反復数を記録
print(f"  [t={t}] CG iterations: {iter_count[0]}, info={info}")
""")

print("\n" + "="*80)
print("推奨: Julia側の詳細ログと比較するため、以下の修正を実施してください:")
print("="*80)
print("1. coeffs_and_rhs_building_DHCP にコールバック追加")
print("2. coeffs_and_rhs_building_Adjoint にコールバック追加")  
print("3. coeffs_and_rhs_building_Sensitivity にコールバック追加")
print("4. 各時刻ステップの反復数を合計してログ出力")
print("\n現在の実装では反復数が取得できないため、ここで終了します。")

sys.exit(0)
