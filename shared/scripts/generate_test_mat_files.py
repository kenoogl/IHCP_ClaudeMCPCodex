#!/usr/bin/env python3
"""
Phase 6 C-1: テスト用MATファイル生成スクリプト

TDDテスト用の小規模IRカメラデータ（MATLABファイル）を生成
Juliaテストの参照データとして使用
"""

import numpy as np
from scipy.io import savemat
import os
from pathlib import Path

def generate_test_mat_files(output_dir: str):
  """
  テスト用MATファイルを生成

  Args:
    output_dir: 出力ディレクトリパス
  """
  # ディレクトリ作成
  os.makedirs(output_dir, exist_ok=True)

  print(f"テスト用MATファイル生成: {output_dir}")

  # テストデータ: 10x8の温度配列（ケルビン単位）
  # 各ファイルで異なる温度パターンを設定

  # SUS1.MAT: 基準温度300K
  mat1_data = np.full((10, 8), 300.0, dtype=np.float64)
  file1 = os.path.join(output_dir, "SUS1.MAT")
  savemat(file1, {"SUS1": mat1_data})
  print(f"  生成: SUS1.MAT (一定温度 300K)")

  # SUS2.MAT: 温度勾配あり
  mat2_data = np.array([[300.0 + i + j for j in range(1, 9)] for i in range(1, 11)],
                       dtype=np.float64)
  file2 = os.path.join(output_dir, "SUS2.MAT")
  savemat(file2, {"SUS2": mat2_data})
  print(f"  生成: SUS2.MAT (温度勾配パターン)")

  # SUS10.MAT: 自然順ソートテスト用（10 > 2）
  mat10_data = np.full((10, 8), 350.0, dtype=np.float64)
  file10 = os.path.join(output_dir, "SUS10.MAT")
  savemat(file10, {"SUS10": mat10_data})
  print(f"  生成: SUS10.MAT (一定温度 350K)")

  # SUS20.MAT: さらに大きい番号
  mat20_data = np.full((10, 8), 400.0, dtype=np.float64)
  file20 = os.path.join(output_dir, "SUS20.MAT")
  savemat(file20, {"SUS20": mat20_data})
  print(f"  生成: SUS20.MAT (一定温度 400K)")

  # ノイズファイル（マッチしないファイル）
  readme_file = os.path.join(output_dir, "README.txt")
  with open(readme_file, "w") as f:
    f.write("Test files for Phase 6 C-1 data loading tests\n")
  print(f"  生成: README.txt (ノイズファイル)")

  # 小文字拡張子（マッチしない）
  mat_lower_data = np.full((10, 8), 999.0, dtype=np.float64)
  file_lower = os.path.join(output_dir, "SUS99.mat")
  savemat(file_lower, {"SUS99": mat_lower_data})
  print(f"  生成: SUS99.mat (小文字拡張子、マッチしない)")

  print("\n生成完了！")

  # 検証用データ出力
  print("\n=== 検証用データサマリー ===")
  print(f"SUS1.MAT: shape={mat1_data.shape}, 値範囲=[{mat1_data.min():.1f}, {mat1_data.max():.1f}]K")
  print(f"SUS2.MAT: shape={mat2_data.shape}, 値範囲=[{mat2_data.min():.1f}, {mat2_data.max():.1f}]K")
  print(f"  例: SUS2[0,0]={mat2_data[0,0]:.1f}K, SUS2[2,1]={mat2_data[2,1]:.1f}K")
  print(f"SUS10.MAT: shape={mat10_data.shape}, 値範囲=[{mat10_data.min():.1f}, {mat10_data.max():.1f}]K")
  print(f"SUS20.MAT: shape={mat20_data.shape}, 値範囲=[{mat20_data.min():.1f}, {mat20_data.max():.1f}]K")


def generate_custom_test_files(output_dir: str):
  """
  カスタムテスト用MATファイルを生成（prefix/extension変更テスト用）

  Args:
    output_dir: 出力ディレクトリパス
  """
  temp_dir = os.path.join(output_dir, "custom_prefix")
  os.makedirs(temp_dir, exist_ok=True)

  print(f"\nカスタムテストファイル生成: {temp_dir}")

  # TEMP接頭辞のファイル
  temp1_data = np.full((5, 5), 300.0, dtype=np.float64)
  file_temp1 = os.path.join(temp_dir, "TEMP1.MAT")
  savemat(file_temp1, {"TEMP1": temp1_data})
  print(f"  生成: TEMP1.MAT")

  temp5_data = np.full((5, 5), 350.0, dtype=np.float64)
  file_temp5 = os.path.join(temp_dir, "TEMP5.MAT")
  savemat(file_temp5, {"TEMP5": temp5_data})
  print(f"  生成: TEMP5.MAT")


def generate_boundary_test_files(output_dir: str):
  """
  境界値テスト用MATファイルを生成（ゼロパディング番号）

  Args:
    output_dir: 出力ディレクトリパス
  """
  pad_dir = os.path.join(output_dir, "zero_padding")
  os.makedirs(pad_dir, exist_ok=True)

  print(f"\n境界値テストファイル生成: {pad_dir}")

  # ゼロパディング版
  sus001_data = np.full((5, 5), 300.0, dtype=np.float64)
  file001 = os.path.join(pad_dir, "SUS001.MAT")
  savemat(file001, {"SUS001": sus001_data})
  print(f"  生成: SUS001.MAT")

  sus010_data = np.full((5, 5), 350.0, dtype=np.float64)
  file010 = os.path.join(pad_dir, "SUS010.MAT")
  savemat(file010, {"SUS010": sus010_data})
  print(f"  生成: SUS010.MAT")

  sus100_data = np.full((5, 5), 400.0, dtype=np.float64)
  file100 = os.path.join(pad_dir, "SUS100.MAT")
  savemat(file100, {"SUS100": sus100_data})
  print(f"  生成: SUS100.MAT")


def generate_nan_test_file(output_dir: str):
  """
  NaN/Inf検証用MATファイルを生成

  Args:
    output_dir: 出力ディレクトリパス
  """
  nan_dir = os.path.join(output_dir, "nan_test")
  os.makedirs(nan_dir, exist_ok=True)

  print(f"\nNaN検証テストファイル生成: {nan_dir}")

  # NaNを含むデータ
  nan_data = np.full((10, 8), 300.0, dtype=np.float64)
  nan_data[4, 4] = np.nan  # (5, 5)位置にNaN（1-indexed）
  file_nan = os.path.join(nan_dir, "SUS1.MAT")
  savemat(file_nan, {"SUS1": nan_data})
  print(f"  生成: SUS1.MAT (NaNを含む)")
  print(f"  NaN位置: [4, 4] (0-indexed)")


if __name__ == "__main__":
  # 基準出力ディレクトリ
  script_dir = Path(__file__).parent
  base_output_dir = script_dir.parent.parent / "julia" / "test" / "data" / "mat_files"

  print("=" * 60)
  print("Phase 6 C-1: テスト用MATファイル生成")
  print("=" * 60)

  # メインテストファイル生成
  generate_test_mat_files(str(base_output_dir))

  # カスタムテストファイル生成
  generate_custom_test_files(str(base_output_dir))

  # 境界値テストファイル生成
  generate_boundary_test_files(str(base_output_dir))

  # NaN検証テストファイル生成
  generate_nan_test_file(str(base_output_dir))

  print("\n" + "=" * 60)
  print("全ファイル生成完了")
  print(f"出力先: {base_output_dir}")
  print("=" * 60)
