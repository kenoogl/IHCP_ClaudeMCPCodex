#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python版 10ステップ CGM3回 ベンチマーク実行スクリプト

オリジナルのIHCP_CGM_Sliding_Window_Calculation_ver2.pyを使用して
10ステップ、CGM反復3回の実行を行う。
"""

import sys
import os
import time
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
    print("=" * 80)
    print("Python版 10ステップ CGM3回 ベンチマーク")
    print("=" * 80)

    # データ読み込み
    data_path = project_root / "shared" / "data" / "T_measure_700um_1ms.npy"
    print(f"データ読み込み: {data_path}")
    T_measure_K = np.load(data_path)
    print(f"データ形状: {T_measure_K.shape}")

    # 10ステップのみ使用
    nt = 10
    Y_obs = T_measure_K[:nt, :, :]
    ni, nj = Y_obs.shape[1], Y_obs.shape[2]
    nk = nz

    print(f"使用データ: nt={nt}, ni={ni}, nj={nj}, nk={nk}")

    # 初期温度場
    T_measure_init_K = T_measure_K[0, :, :]
    T0 = np.repeat(T_measure_init_K[:, :, np.newaxis], nk, axis=2).astype(np.float64)

    print(f"初期温度範囲: {T0.min():.2f} ~ {T0.max():.2f} K")

    # 初期熱流束（ゼロ）
    q_init = np.zeros((nt - 1, ni, nj))

    # タイムステップ
    dt = 0.001

    # CGMパラメータ
    cgm_max_iter = 3
    window_size = 71
    overlap = 17

    print(f"\nCGM設定:")
    print(f"  最大反復数: {cgm_max_iter}")
    print(f"  ウィンドウサイズ: {window_size}")
    print(f"  オーバーラップ: {overlap}")
    print()

    # 実行
    print("=" * 80)
    print("計算開始")
    print("=" * 80)

    start_time = time.time()

    try:
        # オリジナル関数のシグネチャに合わせる
        # def sliding_window_CGM_q_saving(
        #     Y_obs, T0, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
        #     window_size, overlap, q_init_value, filename, CGM_iteration=20000
        # ):

        # rho は関数内で計算されるためrho_coeffsではなくrho値が必要
        # ここでは平均値を使用
        T_avg = np.mean(T0)
        rho_val = np.polyval(rho_coeffs, T_avg)

        # ファイル名は一時的なものを使用（結果は保存しない）
        temp_filename = str(project_root / "shared" / "results" / "temp_python_10steps")

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
            window_size=window_size,
            overlap=overlap,
            q_init_value=0.0,
            filename=temp_filename,
            CGM_iteration=cgm_max_iter
        )

        elapsed_time = time.time() - start_time

        print("=" * 80)
        print("計算完了")
        print("=" * 80)
        print(f"実行時間: {elapsed_time:.2f} 秒")
        print(f"熱流束形状: {q_result.shape}")
        print(f"熱流束範囲: {q_result.min():.6e} ~ {q_result.max():.6e} W/m²")

        # 結果保存
        output_path = project_root / "shared" / "results" / "python_10steps_cgm3.npz"
        np.savez(
            output_path,
            q_result=q_result,
            T0=T0,
            Y_obs=Y_obs,
            elapsed_time=elapsed_time,
            cgm_max_iter=cgm_max_iter,
            nt=nt,
            ni=ni,
            nj=nj,
            nk=nk
        )
        print(f"\n結果保存: {output_path}")

        return 0

    except Exception as e:
        print(f"\n❌ エラー発生: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
