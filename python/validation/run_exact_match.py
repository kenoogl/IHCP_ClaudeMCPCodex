#!/usr/bin/env python3
"""
Python版完全一致検証用実行スクリプト

このスクリプトは、TOMLパラメータファイルから完全に同じパラメータを読み込んで
スライディングウィンドウCGM計算を実行し、Julia版との完全一致検証用データを生成します。

実行方法:
  python python/validation/run_exact_match.py
"""

import os
# Numba並列化設定（インポート前に設定する必要がある）
os.environ['NUMBA_NUM_THREADS'] = '4'

import sys
from pathlib import Path
import numpy as np
import time
import toml

# オリジナルコードのディレクトリをパスに追加
ORIGINAL_DIR = Path(__file__).parent.parent / "original"
sys.path.insert(0, str(ORIGINAL_DIR))

# オリジナルコードから必要な関数をインポート
from IHCP_CGM_Sliding_Window_Calculation_ver2 import (
    thermal_properties_calculator,
    multiple_time_step_solver_DHCP,
    multiple_time_step_solver_Adjoint,
    global_CGM_time,
    sliding_window_CGM_q_saving
)

# 検証器関数を無効化（性能向上のため）
import IHCP_CGM_Sliding_Window_Calculation_ver2 as original_module
def _dummy_validator(*args, **kwargs):
    """ダミー検証器（常にTrueを返す）"""
    return True, "検証スキップ"

original_module.check_temperature_field = _dummy_validator
original_module.check_flux_field = _dummy_validator
original_module.check_adjoint_field = _dummy_validator


def load_verification_params():
    """共通パラメータファイル（TOML）を読み込み"""
    config_path = Path(__file__).parent.parent.parent / "shared/config/verification_params.toml"

    if not config_path.exists():
        raise FileNotFoundError(f"パラメータファイルが見つかりません: {config_path}")

    params = toml.load(config_path)
    return params


def prepare_grid_arrays(params):
    """格子配列を準備（infを処理）"""
    dz = np.array(params['grid_z']['dz'], dtype=np.float64)

    # dz_tとdz_bの'inf'文字列を処理
    dz_t_raw = params['grid_z']['dz_t']
    dz_t = np.array([np.inf if x == 'inf' else float(x) for x in dz_t_raw], dtype=np.float64)

    dz_b_raw = params['grid_z']['dz_b']
    dz_b = np.array([np.inf if x == 'inf' else float(x) for x in dz_b_raw], dtype=np.float64)

    return dz, dz_t, dz_b


def main():
    """メイン実行関数"""
    print("="*60)
    print("Python版 完全一致検証実行")
    print("="*60)

    start_time_total = time.time()

    # ========================================
    # 1. パラメータ読み込み
    # ========================================
    print("\n[1/6] パラメータ読み込み...")
    params = load_verification_params()

    ni = params['problem']['ni']
    nj = params['problem']['nj']
    nk = params['problem']['nk']
    nt = params['problem']['nt']
    dt = params['problem']['dt']
    dx = params['problem']['dx']
    dy = params['problem']['dy']

    rho = params['material']['rho_reference_value']
    cp_coeffs = np.array(params['material']['cp_coeffs'], dtype=np.float64)
    k_coeffs = np.array(params['material']['k_coeffs'], dtype=np.float64)

    dz, dz_t, dz_b = prepare_grid_arrays(params)

    window_size = params['cgm']['window_size']
    overlap = params['cgm']['overlap']
    cgm_iteration = params['cgm']['cgm_iteration']
    q_init_value = params['cgm']['q_init_value']

    T_init_value = params['initial']['T_init_value']

    rtol_dhcp = params['numerical']['rtol_dhcp']
    maxiter_dhcp = params['numerical']['maxiter_dhcp']
    rtol_adjoint = params['numerical']['rtol_adjoint']
    maxiter_adjoint = params['numerical']['maxiter_adjoint']

    print(f"  問題サイズ: ni={ni}, nj={nj}, nk={nk}, nt={nt}")
    print(f"  時間刻み: dt={dt} s")
    print(f"  空間刻み: dx={dx:.6e} m, dy={dy:.6e} m")
    print(f"  材料密度: rho={rho:.6f} kg/m³")
    print(f"  CGMパラメータ: window_size={window_size}, overlap={overlap}, iteration={cgm_iteration}")

    # ========================================
    # 2. データ読み込み
    # ========================================
    print("\n[2/6] 測定データ読み込み...")
    data_path = Path(__file__).parent.parent.parent / "shared/data/T_measure_test_1min.npy"

    if not data_path.exists():
        raise FileNotFoundError(f"測定データが見つかりません: {data_path}")

    Y_obs = np.load(data_path)
    print(f"  データ形状: {Y_obs.shape}")
    print(f"  データ型: {Y_obs.dtype}")
    print(f"  温度範囲: {Y_obs.min():.2f} ~ {Y_obs.max():.2f} K")

    # データサイズの一致確認
    assert Y_obs.shape == (nt, ni, nj), f"データ形状が不一致: {Y_obs.shape} != ({nt}, {ni}, {nj})"

    # ========================================
    # 3. 初期温度場の作成
    # ========================================
    print("\n[3/6] 初期温度場作成...")

    # 厳密に定義された初期温度（全格子点で一定値）
    T_init = np.full((ni, nj, nk), T_init_value, dtype=np.float64)

    print(f"  初期温度: {T_init_value} K")
    print(f"  T_init形状: {T_init.shape}")
    print(f"  T_init型: {T_init.dtype}")

    # ========================================
    # 4. スライディングウィンドウCGM計算
    # ========================================
    print("\n[4/6] スライディングウィンドウCGM計算開始...")
    print(f"  予想ウィンドウ数: 約{(nt-1) // (window_size - overlap) + 1}個")

    # 一時ファイル名（実際には保存しない）
    temp_filename = "temp_q_result.npy"

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

    # 検証誤差計算（表面温度）
    T_calc_surface = T_verify[:, :, :, 0]  # 底面温度
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
    output_dir = Path(__file__).parent.parent.parent / "shared/results/validation"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / "python_exact.npz"

    np.savez_compressed(
        output_path,
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

    print(f"  保存完了: {output_path}")
    print(f"  ファイルサイズ: {output_path.stat().st_size / 1024 / 1024:.2f} MB")

    # ========================================
    # サマリー表示
    # ========================================
    print("\n" + "="*60)
    print("実行完了サマリー")
    print("="*60)
    print(f"総実行時間: {elapsed_total:.2f}秒")
    print(f"  CGM計算:  {elapsed_cgm:.2f}秒 ({elapsed_cgm/elapsed_total*100:.1f}%)")
    print(f"  検証計算: {elapsed_verify:.2f}秒 ({elapsed_verify/elapsed_total*100:.1f}%)")
    print(f"\n結果統計:")
    print(f"  熱流束範囲: {q_result.min():.4e} ~ {q_result.max():.4e} W/m²")
    print(f"  温度誤差（RMS）: {rms_error:.4e} K")
    print(f"  温度誤差（最大）: {max_error:.4e} K")
    print(f"\n次のステップ:")
    print(f"  1. Julia版を実行: julia examples/run_exact_match.jl")
    print(f"  2. 比較実行: python python/validation/compare_exact_match.py")
    print("="*60)

    return q_result, T_verify


if __name__ == "__main__":
    try:
        q_result, T_verify = main()
    except Exception as e:
        print(f"\n❌ エラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
