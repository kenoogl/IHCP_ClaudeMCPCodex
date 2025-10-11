#!/usr/bin/env python3
"""
PythonとJuliaの完全一致検証比較スクリプト

このスクリプトは、Python版とJulia版の計算結果を詳細に比較し、
数値的な一致性を検証します。

実行方法:
  python python/validation/compare_exact_match.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # GUIなし環境対応
import matplotlib.pyplot as plt
from pathlib import Path
import json


def load_results():
    """Python版とJulia版の結果を読み込み"""
    results_dir = Path(__file__).parent.parent.parent / "shared/results/validation"

    python_path = results_dir / "python_exact.npz"
    julia_path = results_dir / "julia_exact.npz"

    if not python_path.exists():
        raise FileNotFoundError(f"Python結果が見つかりません: {python_path}")

    if not julia_path.exists():
        raise FileNotFoundError(f"Julia結果が見つかりません: {julia_path}")

    python_data = np.load(python_path)
    julia_data = np.load(julia_path)

    return python_data, julia_data


def compare_arrays(name, arr_python, arr_julia, verbose=True):
    """配列の詳細比較"""
    if arr_python.shape != arr_julia.shape:
        print(f"⚠️  {name}: 形状が異なります！ Python={arr_python.shape}, Julia={arr_julia.shape}")
        return None

    diff = arr_julia - arr_python
    abs_diff = np.abs(diff)
    rel_diff = abs_diff / (np.abs(arr_python) + 1e-15)

    stats = {
        'name': name,
        'shape': arr_python.shape,
        'max_abs_error': np.max(abs_diff),
        'mean_abs_error': np.mean(abs_diff),
        'rms_error': np.sqrt(np.mean(diff**2)),
        'max_rel_error': np.max(rel_diff),
        'mean_rel_error': np.mean(rel_diff),
        'python_range': (np.min(arr_python), np.max(arr_python)),
        'julia_range': (np.min(arr_julia), np.max(arr_julia))
    }

    if verbose:
        print(f"\n  {name}:")
        print(f"    形状: {stats['shape']}")
        print(f"    最大絶対誤差: {stats['max_abs_error']:.6e}")
        print(f"    平均絶対誤差: {stats['mean_abs_error']:.6e}")
        print(f"    RMS誤差:      {stats['rms_error']:.6e}")
        print(f"    最大相対誤差: {stats['max_rel_error']:.6e}")
        print(f"    平均相対誤差: {stats['mean_rel_error']:.6e}")
        print(f"    Python範囲: {stats['python_range'][0]:.4e} ~ {stats['python_range'][1]:.4e}")
        print(f"    Julia範囲:  {stats['julia_range'][0]:.4e} ~ {stats['julia_range'][1]:.4e}")

    return stats


def analyze_temporal_error(q_python, q_julia, nt):
    """時間ステップ毎の誤差分析"""
    print("\n時間ステップ毎の誤差:")

    temporal_stats = []

    # サンプリング間隔
    sample_indices = [0, nt//4, nt//2, 3*nt//4, nt-1]

    for t in sample_indices:
        if t >= q_python.shape[0]:
            continue

        diff_t = np.abs(q_julia[t] - q_python[t])
        max_err = np.max(diff_t)
        mean_err = np.mean(diff_t)

        temporal_stats.append({
            't': t,
            'max_error': max_err,
            'mean_error': mean_err
        })

        print(f"  t={t:3d}: 最大誤差={max_err:.6e}, 平均誤差={mean_err:.6e}")

    return temporal_stats


def analyze_spatial_error(q_python, q_julia):
    """空間分布の誤差分析"""
    print("\n空間分布の誤差（最終ステップ）:")

    diff_final = np.abs(q_julia[-1] - q_python[-1])

    percentiles = [0, 25, 50, 75, 95, 100]
    values = np.percentile(diff_final, percentiles)

    print(f"  最大誤差: {values[-1]:.6e}")
    for p, v in zip(percentiles[:-1], values[:-1]):
        print(f"  {p}パーセンタイル: {v:.6e}")

    spatial_stats = {
        'percentiles': percentiles,
        'values': values.tolist()
    }

    return spatial_stats


def create_comparison_plots(q_python, q_julia, output_dir):
    """比較可視化"""
    print("\n可視化作成中...")

    nt = q_python.shape[0]
    ni, nj = q_python.shape[1:3]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle('Python vs Julia 完全一致検証', fontsize=16)

    # 1. 時間変化（中心点）
    ax = axes[0, 0]
    i_c, j_c = ni//2, nj//2
    ax.plot(q_python[:, i_c, j_c], 'b-', label='Python', linewidth=1)
    ax.plot(q_julia[:, i_c, j_c], 'r--', label='Julia', linewidth=1)
    ax.set_xlabel('Time step')
    ax.set_ylabel('Heat flux [W/m²]')
    ax.set_title(f'Time series (center: i={i_c}, j={j_c})')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. 誤差の時間変化
    ax = axes[0, 1]
    diff_time = np.abs(q_julia - q_python)
    max_err_time = np.max(diff_time.reshape(nt, -1), axis=1)
    mean_err_time = np.mean(diff_time.reshape(nt, -1), axis=1)
    ax.semilogy(max_err_time, 'r-', label='Max error', linewidth=1)
    ax.semilogy(mean_err_time, 'b-', label='Mean error', linewidth=1)
    ax.set_xlabel('Time step')
    ax.set_ylabel('Absolute error [W/m²]')
    ax.set_title('Temporal error evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 3. 誤差ヒストグラム
    ax = axes[0, 2]
    diff_all = diff_time.flatten()
    ax.hist(diff_all, bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Absolute error [W/m²]')
    ax.set_ylabel('Frequency')
    ax.set_title('Error distribution')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3)

    # 4. 空間分布（Python、最終ステップ）
    ax = axes[1, 0]
    im = ax.imshow(q_python[-1], cmap='hot', aspect='auto')
    ax.set_title('Python (final step)')
    ax.set_xlabel('j')
    ax.set_ylabel('i')
    plt.colorbar(im, ax=ax, label='Heat flux [W/m²]')

    # 5. 空間分布（Julia、最終ステップ）
    ax = axes[1, 1]
    im = ax.imshow(q_julia[-1], cmap='hot', aspect='auto')
    ax.set_title('Julia (final step)')
    ax.set_xlabel('j')
    ax.set_ylabel('i')
    plt.colorbar(im, ax=ax, label='Heat flux [W/m²]')

    # 6. 誤差マップ（最終ステップ）
    ax = axes[1, 2]
    diff_final = np.abs(q_julia[-1] - q_python[-1])
    im = ax.imshow(diff_final, cmap='viridis', aspect='auto')
    ax.set_title('Absolute error (final step)')
    ax.set_xlabel('j')
    ax.set_ylabel('i')
    plt.colorbar(im, ax=ax, label='Absolute error [W/m²]')

    plt.tight_layout()

    output_path = output_dir / "comparison_plots.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  保存完了: {output_path}")


def generate_report(all_stats, output_dir):
    """詳細レポート生成"""
    print("\n詳細レポート作成中...")

    report_path = output_dir / "exact_match_comparison_report.md"

    with open(report_path, 'w') as f:
        f.write("# Python vs Julia 完全一致検証レポート\n\n")
        f.write("## 実行日時\n")
        from datetime import datetime
        f.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        f.write("## 概要\n")
        f.write("PythonオリジナルコードとJulia移植版で完全に同じパラメータを使用して計算を実行し、\n")
        f.write("結果の数値的一致性を検証しました。\n\n")

        f.write("## 比較結果サマリー\n\n")

        # 配列比較結果
        f.write("### 配列比較\n\n")
        f.write("| 項目 | 形状 | 最大絶対誤差 | RMS誤差 | 最大相対誤差 |\n")
        f.write("|------|------|--------------|---------|-------------|\n")

        for stats in all_stats['arrays']:
            if stats is None:
                continue
            f.write(f"| {stats['name']} | {stats['shape']} | ")
            f.write(f"{stats['max_abs_error']:.4e} | ")
            f.write(f"{stats['rms_error']:.4e} | ")
            f.write(f"{stats['max_rel_error']:.4e} |\n")

        f.write("\n")

        # 実行時間比較
        f.write("### 実行時間比較\n\n")
        py_time = all_stats['timing']['python_cgm']
        jl_time = all_stats['timing']['julia_cgm']
        speedup = py_time / jl_time if jl_time > 0 else float('inf')

        f.write(f"- Python CGM計算時間: {py_time:.2f}秒\n")
        f.write(f"- Julia CGM計算時間: {jl_time:.2f}秒\n")
        f.write(f"- 速度比（Python/Julia）: {speedup:.2f}x\n\n")

        # 判定基準
        f.write("## 一致性判定\n\n")

        q_stats = next(s for s in all_stats['arrays'] if s and s['name'] == 'q_result')

        criteria = {
            '初期条件': all_stats['arrays'][0]['max_abs_error'] < 1e-14 if all_stats['arrays'][0] else False,
            '熱流束（最大絶対誤差）': q_stats['max_abs_error'] < 1e-6,
            '熱流束（最大相対誤差）': q_stats['max_rel_error'] < 1e-8,
        }

        all_passed = all(criteria.values())

        f.write("| 判定項目 | 基準 | 結果 |\n")
        f.write("|---------|------|------|\n")

        for name, passed in criteria.items():
            status = "✅ 合格" if passed else "❌ 不合格"
            f.write(f"| {name} | - | {status} |\n")

        f.write("\n")

        f.write(f"### 総合判定: {'✅ 完全一致確認' if all_passed else '❌ 要調査'}\n\n")

        # 詳細統計
        f.write("## 詳細統計\n\n")
        f.write("### 時間ステップ毎の誤差\n\n")
        f.write("| Time step | 最大誤差 | 平均誤差 |\n")
        f.write("|-----------|----------|----------|\n")

        for t_stat in all_stats['temporal']:
            f.write(f"| {t_stat['t']:3d} | {t_stat['max_error']:.6e} | {t_stat['mean_error']:.6e} |\n")

        f.write("\n")

        # 可視化
        f.write("## 可視化\n\n")
        f.write("![比較プロット](comparison_plots.png)\n\n")

        # 結論
        f.write("## 結論\n\n")
        if all_passed:
            f.write("PythonオリジナルコードとJulia移植版は、数値計算上完全に一致することが確認されました。\n")
            f.write("Julia移植は正確に行われています。\n")
        else:
            f.write("一部の判定基準で不一致が検出されました。以下の点を確認してください:\n\n")
            for name, passed in criteria.items():
                if not passed:
                    f.write(f"- {name}\n")

    print(f"  保存完了: {report_path}")

    # JSON形式でも保存
    json_path = output_dir / "exact_match_comparison_stats.json"
    with open(json_path, 'w') as f:
        # NumPy配列をリストに変換
        json_stats = {}
        for key, value in all_stats.items():
            if key == 'arrays':
                json_stats[key] = [
                    {k: (list(v) if isinstance(v, (np.ndarray, tuple)) else v)
                     for k, v in item.items()} if item else None
                    for item in value
                ]
            else:
                json_stats[key] = value

        json.dump(json_stats, f, indent=2)

    print(f"  JSON保存: {json_path}")


def main():
    """メイン比較関数"""
    print("="*60)
    print("Python vs Julia 完全一致検証比較")
    print("="*60)

    # ========================================
    # 1. データ読み込み
    # ========================================
    print("\n[1/5] 結果データ読み込み...")
    python_data, julia_data = load_results()

    print(f"  Python結果: {len(python_data.files)} 項目")
    print(f"  Julia結果: {len(julia_data.files)} 項目")

    # ========================================
    # 2. 初期条件の一致確認
    # ========================================
    print("\n[2/5] 初期条件の一致確認...")

    T_init_python = python_data['T_init']
    T_init_julia = julia_data['T_init']

    T_init_stats = compare_arrays('T_init', T_init_python, T_init_julia)

    if T_init_stats and T_init_stats['max_abs_error'] == 0.0:
        print("  ✅ 初期温度完全一致")
    else:
        print("  ⚠️  初期温度に差異があります")

    # ========================================
    # 3. 主要結果の比較
    # ========================================
    print("\n[3/5] 主要結果の比較...")

    q_python = python_data['q_result']
    q_julia = julia_data['q_result']

    q_stats = compare_arrays('q_result', q_python, q_julia)

    T_verify_python = python_data['T_verify']
    T_verify_julia = julia_data['T_verify']

    T_verify_stats = compare_arrays('T_verify', T_verify_python, T_verify_julia)

    # ========================================
    # 4. 時間・空間誤差分析
    # ========================================
    print("\n[4/5] 時間・空間誤差分析...")

    nt = q_python.shape[0]
    temporal_stats = analyze_temporal_error(q_python, q_julia, nt)
    spatial_stats = analyze_spatial_error(q_python, q_julia)

    # 実行時間比較
    print("\n実行時間比較:")
    py_cgm_time = float(python_data['elapsed_time_cgm'])
    jl_cgm_time = float(julia_data['elapsed_time_cgm'])
    speedup = py_cgm_time / jl_cgm_time if jl_cgm_time > 0 else float('inf')

    print(f"  Python CGM計算: {py_cgm_time:.2f}秒")
    print(f"  Julia CGM計算:  {jl_cgm_time:.2f}秒")
    print(f"  速度比（Python/Julia）: {speedup:.2f}x")

    # ========================================
    # 5. 可視化とレポート生成
    # ========================================
    print("\n[5/5] 可視化とレポート生成...")

    output_dir = Path(__file__).parent.parent.parent / "shared/results/validation"
    create_comparison_plots(q_python, q_julia, output_dir)

    # 統計情報まとめ
    all_stats = {
        'arrays': [T_init_stats, q_stats, T_verify_stats],
        'temporal': temporal_stats,
        'spatial': spatial_stats,
        'timing': {
            'python_cgm': py_cgm_time,
            'julia_cgm': jl_cgm_time,
            'speedup': speedup
        }
    }

    generate_report(all_stats, output_dir)

    # ========================================
    # サマリー表示
    # ========================================
    print("\n" + "="*60)
    print("比較完了サマリー")
    print("="*60)

    if q_stats:
        print(f"\n熱流束（q_result）の一致性:")
        print(f"  最大絶対誤差: {q_stats['max_abs_error']:.6e}")
        print(f"  RMS誤差:      {q_stats['rms_error']:.6e}")
        print(f"  最大相対誤差: {q_stats['max_rel_error']:.6e}")

        # 判定
        if q_stats['max_abs_error'] < 1e-6 and q_stats['max_rel_error'] < 1e-8:
            print("\n✅ 判定: 完全一致確認")
        elif q_stats['max_abs_error'] < 1e-3:
            print("\n⚠️  判定: 高精度一致（実用上問題なし）")
        else:
            print("\n❌ 判定: 要調査")

    print("\n生成ファイル:")
    print(f"  - {output_dir}/exact_match_comparison_report.md")
    print(f"  - {output_dir}/exact_match_comparison_stats.json")
    print(f"  - {output_dir}/comparison_plots.png")
    print("="*60)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ エラーが発生しました: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
