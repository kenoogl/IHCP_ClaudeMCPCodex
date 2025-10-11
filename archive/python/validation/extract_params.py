#!/usr/bin/env python3
"""
Pythonオリジナルコードからパラメータを抽出してTOML形式で保存

このスクリプトは完全一致検証のため、Pythonオリジナルコードで使用されている
すべてのパラメータを抽出し、PythonとJuliaで同じ値を使用できるようにします。
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from math import sin, radians

# オリジナルコードのディレクトリをパスに追加
ORIGINAL_DIR = Path(__file__).parent.parent / "original"
sys.path.insert(0, str(ORIGINAL_DIR))

def extract_parameters():
    """Pythonオリジナルコードからパラメータを抽出"""

    print("="*60)
    print("パラメータ抽出開始")
    print("="*60)

    # ========================================
    # 熱物性値係数の取得
    # ========================================
    thermal_props_path = ORIGINAL_DIR / "metal_thermal_properties.csv"
    sus304_data = pd.read_csv(thermal_props_path)

    sus304_temp = sus304_data['Temperature/K'].values
    sus304_rho = sus304_data['Density'].values
    sus304_cp = sus304_data['Specific_Heat'].values
    sus304_k = sus304_data['Thermal_Conductivity'].values

    # 3次多項式フィッティング
    rho_coeffs = np.polyfit(sus304_temp, sus304_rho, 3)
    cp_coeffs = np.polyfit(sus304_temp, sus304_cp, 3)
    k_coeffs = np.polyfit(sus304_temp, sus304_k, 3)

    # 参照温度での密度（オリジナルコード80行目）
    def polyval_numba(coeffs, x):
        result = 0.0
        for i in range(len(coeffs)):
            result += coeffs[i] * x ** (len(coeffs) - i - 1)
        return result

    rho_value = polyval_numba(rho_coeffs, 225 + 273.15)

    print(f"\n熱物性値多項式係数:")
    print(f"  密度 (rho @ 498.15K): {rho_value:.6f} kg/m³")
    print(f"  rho_coeffs: {rho_coeffs.tolist()}")
    print(f"  cp_coeffs:  {cp_coeffs.tolist()}")
    print(f"  k_coeffs:   {k_coeffs.tolist()}")

    # ========================================
    # 格子パラメータ（オリジナルコード1043-1077行）
    # ========================================
    nz = 20
    dx = 0.12e-3
    dy = 0.12e-3 * sin(radians(80)) / sin(radians(45))
    Lz = 0.5e-3

    stretch_factor = 3

    # z方向の非均等格子生成
    z_faces = np.linspace(1, 0, nz + 1)
    z_faces = Lz - (Lz / (np.exp(stretch_factor) - 1)) * (np.exp(stretch_factor * z_faces) - 1)

    z_centers = np.zeros(nz)
    z_centers[0] = z_faces[0]
    z_centers[-1] = z_faces[-1]
    z_centers[1:-1] = (z_faces[1:-2] + z_faces[2:-1]) / 2

    dz = np.diff(z_faces)

    dz_t = np.zeros(nz)
    dz_t[-1] = np.inf
    dz_t[:-1] = z_centers[1:] - z_centers[:-1]

    dz_b = np.zeros(nz)
    dz_b[0] = np.inf
    dz_b[1:] = z_centers[1:] - z_centers[:-1]

    print(f"\n格子パラメータ:")
    print(f"  nz: {nz}")
    print(f"  dx: {dx:.10e} m")
    print(f"  dy: {dy:.10e} m")
    print(f"  Lz: {Lz:.10e} m")
    print(f"  stretch_factor: {stretch_factor}")
    print(f"  dz: {dz.tolist()}")
    print(f"  dz_t: {[float('inf') if np.isinf(x) else x for x in dz_t.tolist()]}")
    print(f"  dz_b: {[float('inf') if np.isinf(x) else x for x in dz_b.tolist()]}")

    # ========================================
    # テストデータ設定
    # ========================================
    # T_measure_test_1min.npy のサイズを確認
    test_data_path = Path(__file__).parent.parent.parent / "shared/data/T_measure_test_1min.npy"
    if test_data_path.exists():
        Y_obs = np.load(test_data_path)
        nt, ni, nj = Y_obs.shape
        print(f"\nテストデータ:")
        print(f"  データパス: {test_data_path}")
        print(f"  形状: ({nt}, {ni}, {nj})")
    else:
        # デフォルト値（357ステップ = 1分相当）
        nt = 357
        ni = 80
        nj = 100
        print(f"\nテストデータ（デフォルト値）:")
        print(f"  nt: {nt}")
        print(f"  ni: {ni}")
        print(f"  nj: {nj}")

    # 時間刻み（オリジナルコード1663行）
    dt = 0.001  # [s]

    print(f"  dt: {dt} s")

    # ========================================
    # CGMパラメータ
    # ========================================
    # オリジナルコードでは明示的に定義されていないため、
    # 既存のphase5参照データ生成スクリプトから推定
    window_size = 71
    overlap = 17
    cgm_iteration = 20  # テスト用（実際は20000だが検証では短縮）
    q_init_value = 0.0

    print(f"\nCGMパラメータ:")
    print(f"  window_size: {window_size}")
    print(f"  overlap: {overlap}")
    print(f"  cgm_iteration: {cgm_iteration}")
    print(f"  q_init_value: {q_init_value} W/m²")

    # ========================================
    # 初期条件
    # ========================================
    T_init_value = 300.0  # [K]（室温相当）

    print(f"\n初期条件:")
    print(f"  T_init_value: {T_init_value} K")

    # ========================================
    # TOML形式で保存
    # ========================================
    import toml

    params = {
        'problem': {
            'ni': int(ni),
            'nj': int(nj),
            'nk': int(nz),
            'nt': int(nt),
            'dt': float(dt),
            'dx': float(dx),
            'dy': float(dy),
            'Lz': float(Lz),
            'stretch_factor': int(stretch_factor),
            'description': '1分間テストデータ用パラメータ（完全一致検証）'
        },

        'grid_z': {
            'dz': [float(x) for x in dz],
            'dz_t': [float(x) if not np.isinf(x) else 'inf' for x in dz_t],
            'dz_b': [float(x) if not np.isinf(x) else 'inf' for x in dz_b],
            'z_faces': [float(x) for x in z_faces],
            'z_centers': [float(x) for x in z_centers]
        },

        'material': {
            'rho_reference_value': float(rho_value),
            'rho_coeffs': [float(x) for x in rho_coeffs],
            'cp_coeffs': [float(x) for x in cp_coeffs],
            'k_coeffs': [float(x) for x in k_coeffs],
            'description': 'SUS304熱物性値（3次多項式係数）'
        },

        'cgm': {
            'window_size': int(window_size),
            'overlap': int(overlap),
            'cgm_iteration': int(cgm_iteration),
            'q_init_value': float(q_init_value),
            'description': '共役勾配法パラメータ'
        },

        'initial': {
            'T_init_value': float(T_init_value),
            'description': '初期温度（全格子点で一定）'
        },

        'numerical': {
            'rtol_dhcp': 1e-6,
            'maxiter_dhcp': 20000,
            'rtol_adjoint': 1e-8,
            'maxiter_adjoint': 20000,
            'description': '数値ソルバーのパラメータ'
        }
    }

    # TOMLファイルに保存
    output_path = Path(__file__).parent.parent.parent / "shared/config/verification_params.toml"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        toml.dump(params, f)

    print(f"\n{'='*60}")
    print(f"パラメータを保存: {output_path}")
    print(f"{'='*60}")

    return params


if __name__ == "__main__":
    params = extract_parameters()

    print("\n検証用パラメータファイルが生成されました。")
    print("次のステップ:")
    print("  1. Python版実行スクリプトで計算実行")
    print("  2. Julia版実行スクリプトで計算実行")
    print("  3. 比較スクリプトで一致性確認")
