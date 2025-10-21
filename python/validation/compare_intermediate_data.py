#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
シナリオ3: 中間データ数値比較スクリプト

Python版とJulia版のCGM中間データ（T_cal, λ, q, grad, dT）を
要素ごとに比較し、配列形状の対応確認と数値誤差を定量評価する。

Usage:
  python python/validation/compare_intermediate_data.py [--cgm-iter N]

Options:
  --cgm-iter N    比較するCGM反復数（デフォルト: 0-3）
"""

import sys
import argparse
import numpy as np
from pathlib import Path
from typing import Dict, Tuple
import pandas as pd

# プロジェクトルート
PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent


def load_data(cgm_iter: int) -> Tuple[Dict, Dict]:
    """
    指定されたCGM反復のPython/Juliaデータを読み込み

    Args:
        cgm_iter: CGM反復番号

    Returns:
        python_data: Python版データ辞書
        julia_data: Julia版データ辞書
    """
    results_dir = PROJECT_ROOT / "shared" / "results"

    py_path = results_dir / f"python_cgm{cgm_iter}_data.npz"
    jl_path = results_dir / f"julia_cgm{cgm_iter}_data.npz"

    if not py_path.exists():
        raise FileNotFoundError(f"Python data not found: {py_path}")
    if not jl_path.exists():
        raise FileNotFoundError(f"Julia data not found: {jl_path}")

    print(f"\n読み込み: CGM反復{cgm_iter}")
    print(f"  Python: {py_path.name}")
    print(f"  Julia:  {jl_path.name}")

    py_data = dict(np.load(py_path))
    jl_data = dict(np.load(jl_path))

    return py_data, jl_data


def check_array_shapes(py_data: Dict, jl_data: Dict, cgm_iter: int) -> bool:
    """
    配列形状の一致確認

    Args:
        py_data: Python版データ
        jl_data: Julia版データ
        cgm_iter: CGM反復番号

    Returns:
        すべての配列形状が一致すればTrue
    """
    print(f"\n[CGM反復{cgm_iter}] 配列形状確認")
    print("-" * 80)

    arrays = ['T_cal', 'lambda_field', 'q', 'grad', 'dT']
    all_match = True

    for name in arrays:
        if name not in py_data or name not in jl_data:
            print(f"  ❌ {name}: データが見つかりません")
            all_match = False
            continue

        py_shape = py_data[name].shape
        jl_shape = jl_data[name].shape

        match = py_shape == jl_shape
        symbol = "✅" if match else "❌"

        print(f"  {symbol} {name:15s}: Python {py_shape} vs Julia {jl_shape}")

        if not match:
            all_match = False

    return all_match


def compute_errors(py_array: np.ndarray, jl_array: np.ndarray,
                   name: str, cgm_iter: int) -> Dict[str, float]:
    """
    要素ごとの誤差統計を計算

    Args:
        py_array: Python版配列
        jl_array: Julia版配列
        name: 配列名
        cgm_iter: CGM反復番号

    Returns:
        誤差統計辞書
    """
    if py_array.shape != jl_array.shape:
        raise ValueError(f"Shape mismatch for {name}: {py_array.shape} vs {jl_array.shape}")

    # 絶対誤差
    abs_diff = np.abs(py_array - jl_array)

    # 相対誤差（ゼロ除算回避）
    denominator = np.abs(py_array) + 1e-10
    rel_diff = abs_diff / denominator

    stats = {
        'name': name,
        'cgm_iter': cgm_iter,
        'max_abs_err': abs_diff.max(),
        'mean_abs_err': abs_diff.mean(),
        'median_abs_err': np.median(abs_diff),
        'max_rel_err': rel_diff.max(),
        'mean_rel_err': rel_diff.mean(),
        'median_rel_err': np.median(rel_diff),
        'py_min': py_array.min(),
        'py_max': py_array.max(),
        'py_mean': py_array.mean(),
        'jl_min': jl_array.min(),
        'jl_max': jl_array.max(),
        'jl_mean': jl_array.mean(),
    }

    return stats


def print_error_stats(stats: Dict):
    """誤差統計を表形式で出力"""
    print(f"\n{stats['name']} (CGM反復{stats['cgm_iter']})")
    print("-" * 80)

    print(f"  Python値範囲:  {stats['py_min']:+.6e} ~ {stats['py_max']:+.6e} (平均: {stats['py_mean']:+.6e})")
    print(f"  Julia値範囲:   {stats['jl_min']:+.6e} ~ {stats['jl_max']:+.6e} (平均: {stats['jl_mean']:+.6e})")
    print()
    print(f"  絶対誤差:")
    print(f"    最大: {stats['max_abs_err']:.6e}")
    print(f"    平均: {stats['mean_abs_err']:.6e}")
    print(f"    中央: {stats['median_abs_err']:.6e}")
    print()
    print(f"  相対誤差:")
    print(f"    最大: {stats['max_rel_err']:.6e} ({stats['max_rel_err']*100:.4f}%)")
    print(f"    平均: {stats['mean_rel_err']:.6e} ({stats['mean_rel_err']*100:.6f}%)")
    print(f"    中央: {stats['median_rel_err']:.6e} ({stats['median_rel_err']*100:.6f}%)")


def compare_cgm_iteration(cgm_iter: int) -> Dict:
    """
    指定されたCGM反復の中間データを比較

    Args:
        cgm_iter: CGM反復番号

    Returns:
        すべての配列の誤差統計リスト
    """
    # データ読み込み
    py_data, jl_data = load_data(cgm_iter)

    # 形状確認
    shapes_ok = check_array_shapes(py_data, jl_data, cgm_iter)
    if not shapes_ok:
        print("\n❌ 配列形状が一致しません。配列の対応を確認してください。")
        return None

    # 各配列の誤差計算
    arrays = ['T_cal', 'lambda_field', 'q', 'grad', 'dT']
    all_stats = []

    print(f"\n[CGM反復{cgm_iter}] 数値誤差比較")
    print("=" * 80)

    for name in arrays:
        if name not in py_data or name not in jl_data:
            continue

        stats = compute_errors(py_data[name], jl_data[name], name, cgm_iter)
        print_error_stats(stats)
        all_stats.append(stats)

    # 目的関数値の比較
    if 'J' in py_data and 'J' in jl_data:
        J_py = float(py_data['J'])
        J_jl = float(jl_data['J'])
        J_diff = abs(J_py - J_jl)
        J_rel = J_diff / (abs(J_py) + 1e-10)

        print(f"\n目的関数 J")
        print("-" * 80)
        print(f"  Python: {J_py:.6e}")
        print(f"  Julia:  {J_jl:.6e}")
        print(f"  差分:   {J_diff:.6e} (相対: {J_rel:.6e}, {J_rel*100:.4f}%)")

    return all_stats


def create_summary_report(all_iterations_stats: Dict[int, list], output_path: Path):
    """
    全反復の誤差統計サマリーレポートを作成

    Args:
        all_iterations_stats: {cgm_iter: [stats, ...]}
        output_path: 出力ファイルパス
    """
    # DataFrameに変換
    rows = []
    for cgm_iter, stats_list in all_iterations_stats.items():
        for stats in stats_list:
            rows.append(stats)

    df = pd.DataFrame(rows)

    # Markdown形式でレポート作成
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# シナリオ3: 中間データ数値比較レポート\n\n")
        f.write(f"**生成日時**: {pd.Timestamp.now()}\n\n")
        f.write("---\n\n")

        f.write("## サマリー\n\n")

        # CGM反復ごとに表作成
        for cgm_iter in sorted(all_iterations_stats.keys()):
            f.write(f"### CGM反復{cgm_iter}\n\n")

            df_iter = df[df['cgm_iter'] == cgm_iter]

            f.write("| 配列 | 最大絶対誤差 | 平均絶対誤差 | 最大相対誤差 | 平均相対誤差 |\n")
            f.write("|------|--------------|--------------|--------------|---------------|\n")

            for _, row in df_iter.iterrows():
                f.write(f"| {row['name']:15s} | "
                       f"{row['max_abs_err']:.3e} | "
                       f"{row['mean_abs_err']:.3e} | "
                       f"{row['max_rel_err']*100:.4f}% | "
                       f"{row['mean_rel_err']*100:.6f}% |\n")

            f.write("\n")

        # 詳細統計
        f.write("---\n\n")
        f.write("## 詳細統計\n\n")

        for cgm_iter in sorted(all_iterations_stats.keys()):
            f.write(f"### CGM反復{cgm_iter}\n\n")

            df_iter = df[df['cgm_iter'] == cgm_iter]

            for _, row in df_iter.iterrows():
                f.write(f"#### {row['name']}\n\n")
                f.write(f"- **Python値範囲**: {row['py_min']:+.6e} ~ {row['py_max']:+.6e} (平均: {row['py_mean']:+.6e})\n")
                f.write(f"- **Julia値範囲**: {row['jl_min']:+.6e} ~ {row['jl_max']:+.6e} (平均: {row['jl_mean']:+.6e})\n")
                f.write(f"- **絶対誤差**: 最大 {row['max_abs_err']:.6e}, 平均 {row['mean_abs_err']:.6e}, 中央 {row['median_abs_err']:.6e}\n")
                f.write(f"- **相対誤差**: 最大 {row['max_rel_err']*100:.4f}%, 平均 {row['mean_rel_err']*100:.6f}%, 中央 {row['median_rel_err']*100:.6f}%\n\n")

    print(f"\n✅ サマリーレポート保存: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Python-Julia中間データ数値比較（シナリオ3）"
    )
    parser.add_argument(
        "--cgm-iter",
        type=int,
        nargs='+',
        default=[0, 1, 2, 3],
        help="比較するCGM反復番号（デフォルト: 0 1 2 3）"
    )
    args = parser.parse_args()

    print("=" * 80)
    print("シナリオ3: Python-Julia中間データ数値比較")
    print("=" * 80)

    all_iterations_stats = {}

    for cgm_iter in args.cgm_iter:
        try:
            stats = compare_cgm_iteration(cgm_iter)
            if stats:
                all_iterations_stats[cgm_iter] = stats
        except FileNotFoundError as e:
            print(f"\n⚠️  CGM反復{cgm_iter}: データファイルが見つかりません")
            print(f"   {e}")
            continue
        except Exception as e:
            print(f"\n❌ CGM反復{cgm_iter}: エラー発生")
            print(f"   {e}")
            import traceback
            traceback.print_exc()
            continue

    # サマリーレポート作成
    if all_iterations_stats:
        output_path = PROJECT_ROOT / "shared" / "results" / "scenario3_comparison_report.md"
        create_summary_report(all_iterations_stats, output_path)

    print("\n" + "=" * 80)
    print("比較完了")
    print("=" * 80)

    return 0


if __name__ == "__main__":
    sys.exit(main())
