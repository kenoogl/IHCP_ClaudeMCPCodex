#!/usr/bin/env python3
"""
Python版 100ステップ一致検証テストスクリプト

このスクリプトは、Julia版との計算結果一致確認用に、T_measure_700um_1ms.npyの
最初の100ステップを使用してスライディングウィンドウCGM計算を実行します。

改訂版の特徴:
- ステップ数: 300 → 100に削減（計算時間短縮）
- CGM反復数: 2 → 1に削減（既存検証と同じ条件）

実行方法:
  python python/validation/test_python_julia_match_100steps.py
"""

import os
# Numba並列化設定（インポート前に設定する必要がある）
os.environ['NUMBA_NUM_THREADS'] = '8'

import sys
from pathlib import Path
import numpy as np
import time

# オリジナルコードのディレクトリをパスに追加
ORIGINAL_DIR = Path(__file__).parent.parent / "original"
sys.path.insert(0, str(ORIGINAL_DIR))

# オリジナルコードから必要な関数をインポート
from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    thermal_properties_calculator,
    multiple_time_step_solver_DHCP,
    multiple_time_step_solver_Adjoint,
    global_CGM_time,
    sliding_window_CGM_q_saving,
    rho,
    cp_coeffs,
    k_coeffs,
    dx,
    dy,
    dz,
    dz_b,
    dz_t
)

def main():
    """メイン実行関数"""
    print("="*80)
    print("Python版 100ステップ一致検証テスト")
    print("="*80)
    print(f"Numbaスレッド数: {os.environ.get('NUMBA_NUM_THREADS', '不明')}")

    start_time_total = time.time()

    # ========================================
    # 1. パラメータ設定
    # ========================================
    print("\n[1/6] パラメータ設定...")

    # 時間パラメータ
    nt = 100  # 100ステップに削減
    dt = 0.001  # 1 ms

    # CGMパラメータ
    window_size = 71
    overlap = 17
    cgm_iteration = 1  # 1反復に削減（既存検証と同じ条件）
    q_init_value = 0.0

    # 数値ソルバーパラメータ（デフォルト）
    rtol_dhcp = 1.0e-6
    maxiter_dhcp = 20000
    rtol_adjoint = 1.0e-8
    maxiter_adjoint = 20000

    print(f"  時間ステップ数: nt={nt}")
    print(f"  時間刻み: dt={dt} s")
    print(f"  空間刻み: dx={dx:.6e} m, dy={dy:.6e} m")
    print(f"  材料密度: rho={rho:.6f} kg/m³")
    print(f"  CGMパラメータ: window_size={window_size}, overlap={overlap}, iteration={cgm_iteration}")

    # ========================================
    # 2. データ読み込み
    # ========================================
    print("\n[2/6] 測定データ読み込み...")
    data_path = Path(__file__).parent.parent.parent / "shared/data/T_measure_700um_1ms.npy"

    if not data_path.exists():
        raise FileNotFoundError(f"測定データが見つかりません: {data_path}")

    # 全データ読み込み
    T_measure_full = np.load(data_path)
    print(f"  全データ形状: {T_measure_full.shape}")

    # 100ステップに制限
    Y_obs = T_measure_full[:nt, :, :]
    print(f"  使用データ形状: {Y_obs.shape}")
    print(f"  データ型: {Y_obs.dtype}")
    print(f"  温度範囲: {Y_obs.min():.2f} ~ {Y_obs.max():.2f} K")

    ni, nj = Y_obs.shape[1], Y_obs.shape[2]
    nk = dz.shape[0]

    # ========================================
    # 3. 初期温度場の作成
    # ========================================
    print("\n[3/6] 初期温度場作成...")

    # 最初のフレームから初期温度を設定
    T_measure_init = Y_obs[0, :, :]
    T_init = np.repeat(T_measure_init[:, :, np.newaxis], nk, axis=2).astype(np.float64)

    print(f"  T_init形状: {T_init.shape}")
    print(f"  T_init型: {T_init.dtype}")
    print(f"  初期温度範囲: {T_init.min():.2f} ~ {T_init.max():.2f} K")

    # ========================================
    # 4. スライディングウィンドウCGM計算
    # ========================================
    print("\n[4/6] スライディングウィンドウCGM計算開始...")
    expected_windows = (nt - 1) // (window_size - overlap) + 1
    print(f"  予想ウィンドウ数: 約{expected_windows}個")

    # 一時ファイル名
    temp_filename = "temp_python_100steps_q.npy"

    start_time_cgm = time.time()

    q_result = sliding_window_CGM_q_saving(
        Y_obs=Y_obs,
        T0=T_init,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_b=dz_b,
        dz_t=dz_t,
        dt=dt,
        rho=rho,
        cp_coeffs=cp_coeffs,
        k_coeffs=k_coeffs,
        window_size=window_size,
        overlap=overlap,
        q_init_value=q_init_value,
        filename=temp_filename,
        CGM_iteration=cgm_iteration
    )

    elapsed_cgm = time.time() - start_time_cgm

    print(f"\n  CGM計算完了！")
    print(f"  計算時間: {elapsed_cgm:.2f}秒")
    print(f"  結果形状: {q_result.shape}")
    print(f"  熱流束範囲: {q_result.min():.2e} ~ {q_result.max():.2e} W/m²")

    # ========================================
    # 5. 結果検証（DHCP forward計算）
    # ========================================
    print("\n[5/6] 結果検証（DHCP順解析）...")
    start_time_verify = time.time()

    T_verify = multiple_time_step_solver_DHCP(
        T_initial=T_init,
        q_surface=q_result,
        nt=nt,
        rho=rho,
        cp_coeffs=cp_coeffs,
        k_coeffs=k_coeffs,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_b=dz_b,
        dz_t=dz_t,
        dt=dt,
        rtol=rtol_dhcp,
        maxiter=maxiter_dhcp
    )

    elapsed_verify = time.time() - start_time_verify

    # 検証誤差計算（底面温度）
    # Python版: T_verify.shape = (nt, ni, nj, nk)
    # 底面はk=0
    T_calc_surface = T_verify[:, :, :, 0]
    residual = T_calc_surface - Y_obs
    rms_error = np.sqrt(np.mean(residual**2))
    max_error = np.abs(residual).max()

    print(f"  検証計算時間: {elapsed_verify:.2f}秒")
    print(f"  温度誤差（RMS）: {rms_error:.4e} K")
    print(f"  温度誤差（最大）: {max_error:.4e} K")

    # ========================================
    # 6. 結果保存
    # ========================================
    print("\n[6/6] 結果保存...")
    output_dir = Path(__file__).parent.parent.parent / "shared/results"
    output_dir.mkdir(parents=True, exist_ok=True)

    # 個別に保存（.npy形式、非圧縮）
    q_path = output_dir / "python_100steps_q.npy"
    T_path = output_dir / "python_100steps_T.npy"

    np.save(q_path, q_result)
    np.save(T_path, T_verify)

    # メタデータ付き圧縮版も保存
    npz_path = output_dir / "python_100steps.npz"
    np.savez_compressed(
        npz_path,
        # 入力データ
        Y_obs=Y_obs,
        T_init=T_init,
        # パラメータ
        dt=dt,
        dx=dx,
        dy=dy,
        dz=dz,
        dz_t=dz_t,
        dz_b=dz_b,
        rho=rho,
        cp_coeffs=cp_coeffs,
        k_coeffs=k_coeffs,
        # CGMパラメータ
        window_size=window_size,
        overlap=overlap,
        cgm_iteration=cgm_iteration,
        # 結果
        q_result=q_result,
        T_verify=T_verify,
        # 統計情報
        elapsed_time_cgm=elapsed_cgm,
        elapsed_time_verify=elapsed_verify,
        rms_error=rms_error,
        max_error=max_error
    )

    elapsed_total = time.time() - start_time_total

    print(f"  保存完了:")
    print(f"    - {q_path} ({q_path.stat().st_size / 1024 / 1024:.2f} MB)")
    print(f"    - {T_path} ({T_path.stat().st_size / 1024 / 1024:.2f} MB)")
    print(f"    - {npz_path} ({npz_path.stat().st_size / 1024 / 1024:.2f} MB)")

    # ========================================
    # サマリー表示
    # ========================================
    print("\n" + "="*80)
    print("実行完了サマリー")
    print("="*80)
    print(f"総実行時間: {elapsed_total:.2f}秒")
    print(f"  CGM計算:  {elapsed_cgm:.2f}秒 ({elapsed_cgm/elapsed_total*100:.1f}%)")
    print(f"  検証計算: {elapsed_verify:.2f}秒 ({elapsed_verify/elapsed_total*100:.1f}%)")
    print(f"\n結果統計:")
    print(f"  熱流束範囲: {q_result.min():.4e} ~ {q_result.max():.4e} W/m²")
    print(f"  温度誤差（RMS）: {rms_error:.4e} K")
    print(f"  温度誤差（最大）: {max_error:.4e} K")
    print(f"\n次のステップ:")
    print(f"  Julia版を実行: JULIA_NUM_THREADS=1 julia julia/examples/test_python_julia_match_100steps.jl")
    print("="*80)

    # 一時ファイルの削除
    if Path(temp_filename).exists():
        Path(temp_filename).unlink()

    return q_result, T_verify


if __name__ == "__main__":
    try:
        q_result, T_verify = main()
    except Exception as e:
        print(f"\n❌ エラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
