#!/usr/bin/env python3
"""
Python版とJulia版の10ステップCGM3反復結果を比較

目的:
- 熱流束(q_result)の数値的一致性確認
- 温度場の比較（もしあれば）
- 詳細な統計情報の出力
"""

import numpy as np
import sys

def load_python_results():
    """Python版の結果を読み込み"""
    try:
        data = np.load('shared/results/python_10steps_cgm3.npz')
        return data
    except FileNotFoundError:
        print("❌ Python版の結果ファイルが見つかりません: shared/results/python_10steps_cgm3.npz")
        return None

def load_julia_results():
    """Julia版の結果を読み込み（複数候補を試行）"""
    candidates = [
        'shared/results/julia_10steps_cgm3.npz',
        'shared/results/julia_sliding_window_10steps_test_stats.npz',
        'shared/results/julia_10steps_fullsize.npz',
    ]

    for path in candidates:
        try:
            data = np.load(path)
            print(f"✅ Julia版結果読み込み成功: {path}")
            return data, path
        except FileNotFoundError:
            continue

    print("❌ Julia版の結果ファイルが見つかりません")
    print("   候補:")
    for path in candidates:
        print(f"   - {path}")
    return None, None

def compare_arrays(py_array, jl_array, name, rtol=1e-5, atol=1e-8):
    """2つの配列を比較"""
    print(f"\n{'='*80}")
    print(f"{name}の比較")
    print(f"{'='*80}")

    # 形状確認
    print(f"Python版形状: {py_array.shape}")
    print(f"Julia版形状:  {jl_array.shape}")

    if py_array.shape != jl_array.shape:
        print(f"⚠️  形状が一致しません")
        return False

    # 統計情報
    print(f"\nPython版:")
    print(f"  範囲: [{py_array.min():.6e}, {py_array.max():.6e}]")
    print(f"  平均: {py_array.mean():.6e}")
    print(f"  標準偏差: {py_array.std():.6e}")

    print(f"\nJulia版:")
    print(f"  範囲: [{jl_array.min():.6e}, {jl_array.max():.6e}]")
    print(f"  平均: {jl_array.mean():.6e}")
    print(f"  標準偏差: {jl_array.std():.6e}")

    # 差分解析
    diff = py_array - jl_array
    abs_diff = np.abs(diff)
    rel_diff = abs_diff / (np.abs(py_array) + atol)

    print(f"\n差分統計:")
    print(f"  絶対誤差 最大: {abs_diff.max():.6e}")
    print(f"  絶対誤差 平均: {abs_diff.mean():.6e}")
    print(f"  絶対誤差 RMS:  {np.sqrt((diff**2).mean()):.6e}")
    print(f"  相対誤差 最大: {rel_diff.max():.6e}")
    print(f"  相対誤差 平均: {rel_diff.mean():.6e}")

    # 一致性判定
    is_close = np.allclose(py_array, jl_array, rtol=rtol, atol=atol)

    if is_close:
        print(f"\n✅ 一致: rtol={rtol}, atol={atol}で一致しています")
    else:
        print(f"\n❌ 不一致: rtol={rtol}, atol={atol}で不一致です")
        # 不一致点の数
        mismatch = ~np.isclose(py_array, jl_array, rtol=rtol, atol=atol)
        print(f"   不一致要素数: {mismatch.sum()} / {py_array.size} ({100*mismatch.sum()/py_array.size:.2f}%)")

    return is_close

def main():
    print("="*80)
    print("Python-Julia 10ステップ CGM3反復 計算結果比較")
    print("="*80)

    # Python版読み込み
    py_data = load_python_results()
    if py_data is None:
        return 1

    print("\n=== Python版データ ===")
    print(f"キー: {list(py_data.keys())}")

    # Julia版読み込み
    jl_data, jl_path = load_julia_results()
    if jl_data is None:
        print("\n⚠️  Julia版のnpzファイルが見つかりません")
        print("   Julia版を実行して結果を保存してください")
        print("\n📊 Python版の結果のみ表示:")
        print(f"\n熱流束 (q_result):")
        q_py = py_data['q_result']
        print(f"  形状: {q_py.shape}")
        print(f"  範囲: [{q_py.min():.6e}, {q_py.max():.6e}]")
        print(f"  平均: {q_py.mean():.6e}")
        return 0

    print(f"\n=== Julia版データ ({jl_path}) ===")
    print(f"キー: {list(jl_data.keys())}")

    # 熱流束比較
    if 'q_result' in py_data and 'q_result' in jl_data:
        q_py = py_data['q_result']
        q_jl = jl_data['q_result']

        # メモリレイアウト確認（Python: (nt-1, ni, nj), Julia: (ni, nj, nt-1) の可能性）
        if q_py.shape != q_jl.shape:
            print(f"\n⚠️  形状不一致を検出。Juliaデータの転置を試みます...")
            if q_py.shape == q_jl.T.shape or q_py.shape == tuple(reversed(q_jl.shape)):
                q_jl_orig = q_jl
                # Julia版が(ni, nj, nt-1)なら(nt-1, ni, nj)に変換
                if len(q_jl.shape) == 3:
                    q_jl = np.transpose(q_jl, (2, 0, 1))
                    print(f"   転置後のJulia形状: {q_jl.shape}")

        compare_arrays(q_py, q_jl, "熱流束 (q_result)", rtol=1e-5, atol=1e-8)
    else:
        print(f"\n⚠️  熱流束データが見つかりません")
        if 'q_result' not in py_data:
            print(f"   Python版にq_resultが存在しません")
        if 'q_result' not in jl_data:
            print(f"   Julia版にq_resultが存在しません")

    # 温度場比較（あれば）
    if 'T0' in py_data and 'T0' in jl_data:
        T_py = py_data['T0']
        T_jl = jl_data['T0']
        compare_arrays(T_py, T_jl, "初期温度場 (T0)", rtol=1e-6, atol=1e-8)

    print("\n" + "="*80)
    print("比較完了")
    print("="*80)

    return 0

if __name__ == "__main__":
    sys.exit(main())
