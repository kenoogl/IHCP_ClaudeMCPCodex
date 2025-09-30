"""
Phase 1参照データ生成スクリプト
TDD用にPythonの計算結果をJSON形式で保存する
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path

# Numba実装（オリジナルと同じ）
from numba import njit, prange

@njit
def polyval_numba(coeffs, x):
  result = 0.0
  for i in range(len(coeffs)):
    result += coeffs[i] * x ** (len(coeffs) - i - 1)
  return result

@njit(parallel=True)
def thermal_properties_calculator(Temperature, cp_coeffs, k_coeffs):
  ni, nj, nk = Temperature.shape
  cp = np.empty((ni, nj, nk))
  k = np.empty((ni, nj, nk))

  for i in prange(ni):
    for j in range(nj):
      for k_ijk in range(nk):
        T_current = Temperature[i, j, k_ijk]
        cp[i, j, k_ijk] = polyval_numba(cp_coeffs, T_current)
        k[i, j, k_ijk] = polyval_numba(k_coeffs, T_current)

  return cp, k


# データ読み込み
BASE_DIR = Path(__file__).resolve().parent.parent.parent
Thermal_properties_file_path = BASE_DIR / "shared/data/metal_thermal_properties.csv"
sus304_data = pd.read_csv(Thermal_properties_file_path)

sus304_temp = sus304_data['Temperature/K'].values
sus304_rho = sus304_data['Density'].values
sus304_cp = sus304_data['Specific_Heat'].values
sus304_k = sus304_data['Thermal_Conductivity'].values

# 多項式フィッティング（3次）
rho_coeffs = np.polyfit(sus304_temp, sus304_rho, 3)
cp_coeffs = np.polyfit(sus304_temp, sus304_cp, 3)
k_coeffs = np.polyfit(sus304_temp, sus304_k, 3)

# ρ計算（225℃ + 273.15K = 498.15K）
rho_value = polyval_numba(rho_coeffs, 225 + 273.15)

# テスト用温度配列生成（小規模: 3x3x3）
test_temp_small = np.array([
  [[300.0, 400.0, 500.0],
   [600.0, 700.0, 800.0],
   [900.0, 1000.0, 1100.0]],

  [[350.0, 450.0, 550.0],
   [650.0, 750.0, 850.0],
   [950.0, 1050.0, 1150.0]],

  [[375.0, 475.0, 575.0],
   [675.0, 775.0, 875.0],
   [975.0, 1075.0, 1175.0]]
])

# 熱物性値計算
cp_small, k_small = thermal_properties_calculator(test_temp_small, cp_coeffs, k_coeffs)

# 参照データ辞書作成
reference_data = {
  # 元データ
  "sus304_temp": sus304_temp.tolist(),
  "sus304_rho": sus304_rho.tolist(),
  "sus304_cp": sus304_cp.tolist(),
  "sus304_k": sus304_k.tolist(),

  # 多項式係数
  "rho_coeffs": rho_coeffs.tolist(),
  "cp_coeffs": cp_coeffs.tolist(),
  "k_coeffs": k_coeffs.tolist(),

  # polyval_numbaテスト用
  "polyval_test_temp": 498.15,
  "polyval_test_rho": float(rho_value),

  # thermal_properties_calculatorテスト用（3x3x3）
  "test_temp_small": test_temp_small.tolist(),
  "cp_small": cp_small.tolist(),
  "k_small": k_small.tolist(),

  # 追加の個別テストケース
  "single_temp_tests": {
    "300K": {
      "temp": 300.0,
      "cp": float(polyval_numba(cp_coeffs, 300.0)),
      "k": float(polyval_numba(k_coeffs, 300.0))
    },
    "800K": {
      "temp": 800.0,
      "cp": float(polyval_numba(cp_coeffs, 800.0)),
      "k": float(polyval_numba(k_coeffs, 800.0))
    },
    "1600K": {
      "temp": 1600.0,
      "cp": float(polyval_numba(cp_coeffs, 1600.0)),
      "k": float(polyval_numba(k_coeffs, 1600.0))
    }
  }
}

# JSON保存
output_path = BASE_DIR / "julia/data/phase1_reference_data.json"
output_path.parent.mkdir(parents=True, exist_ok=True)

with open(output_path, 'w') as f:
  json.dump(reference_data, f, indent=2)

print(f"参照データ保存完了: {output_path}")
print(f"\n主要な値:")
print(f"  ρ(498.15K) = {rho_value:.10f} kg/m³")
print(f"  cp係数 = {cp_coeffs}")
print(f"  k係数 = {k_coeffs}")
print(f"\nテスト配列形状: {test_temp_small.shape}")
print(f"  cp範囲: {cp_small.min():.2f} - {cp_small.max():.2f}")
print(f"  k範囲: {k_small.min():.2f} - {k_small.max():.2f}")