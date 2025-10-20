#!/usr/bin/env python3
"""
ウィンドウ1のPython-Julia結果比較
"""

import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent

def main():
    print("="*80)
    print("ウィンドウ1結果比較")
    print("="*80)

    # データ読み込み
    python_path = PROJECT_ROOT / "shared" / "results" / "window1_python_test.npz"
    julia_path = PROJECT_ROOT / "shared" / "results" / "window1_julia_test.npz"

    if not python_path.exists():
        print(f"❌ Python結果が見つかりません: {python_path}")
        return

    if not julia_path.exists():
        print(f"❌ Julia結果が見つかりません: {julia_path}")
        return

    python_data = np.load(str(python_path))
    julia_data = np.load(str(julia_path))

    print("\n" + "="*80)
    print("入力データの確認")
    print("="*80)

    # Y_obsの比較
    Y_obs_py = python_data['Y_obs_win']
    Y_obs_jl = julia_data['Y_obs_win']
    print(f"\nY_obs形状:")
    print(f"  Python: {Y_obs_py.shape}")
    print(f"  Julia:  {Y_obs_jl.shape}")
    print(f"Y_obs一致: {np.allclose(Y_obs_py, Y_obs_jl)}")
    if not np.allclose(Y_obs_py, Y_obs_jl):
        diff = np.abs(Y_obs_py - Y_obs_jl)
        print(f"  最大差: {np.max(diff):.6e} K")

    # T_initの比較
    T_init_py = python_data['T_init']
    T_init_jl = julia_data['T_init']
    print(f"\nT_init形状:")
    print(f"  Python: {T_init_py.shape}")
    print(f"  Julia:  {T_init_jl.shape}")
    print(f"T_init一致: {np.allclose(T_init_py, T_init_jl)}")
    if not np.allclose(T_init_py, T_init_jl):
        diff = np.abs(T_init_py - T_init_jl)
        print(f"  最大差: {np.max(diff):.6e} K")

    # q_initの比較
    q_init_py = python_data['q_init']
    q_init_jl = julia_data['q_init']
    print(f"\nq_init形状:")
    print(f"  Python: {q_init_py.shape}")
    print(f"  Julia:  {q_init_jl.shape}")
    print(f"q_init一致: {np.allclose(q_init_py, q_init_jl)}")

    print("\n" + "="*80)
    print("出力結果の比較")
    print("="*80)

    # J履歴の比較
    J_py = python_data['J_hist']
    J_jl = julia_data['J_hist']
    print(f"\nJ履歴:")
    print(f"  Python: {J_py}")
    print(f"  Julia:  {J_jl}")
    if len(J_py) == len(J_jl):
        J_diff = np.abs(np.array(J_py) - np.array(J_jl))
        J_rel_diff = J_diff / np.abs(J_py)
        print(f"  絶対差: {J_diff}")
        print(f"  相対差: {J_rel_diff * 100} %")
    else:
        print(f"  ⚠️ 反復回数が異なります")

    # q_resultの比較
    q_py = python_data['q_result']
    q_jl = julia_data['q_result']
    print(f"\nq_result形状:")
    print(f"  Python: {q_py.shape}")
    print(f"  Julia:  {q_jl.shape}")

    if q_py.shape == q_jl.shape:
        q_diff = np.abs(q_py - q_jl)
        q_rel_diff = q_diff / (np.abs(q_py) + 1e-10)

        print(f"\nq_result統計:")
        print(f"  Python範囲: [{np.min(q_py):.6e}, {np.max(q_py):.6e}] W/m²")
        print(f"  Julia範囲:  [{np.min(q_jl):.6e}, {np.max(q_jl):.6e}] W/m²")
        print(f"\n  絶対差:")
        print(f"    最大: {np.max(q_diff):.6e} W/m²")
        print(f"    平均: {np.mean(q_diff):.6e} W/m²")
        print(f"    RMS:  {np.sqrt(np.mean(q_diff**2)):.6e} W/m²")
        print(f"\n  相対差:")
        print(f"    最大: {np.max(q_rel_diff) * 100:.6f} %")
        print(f"    平均: {np.mean(q_rel_diff) * 100:.6f} %")
        print(f"    RMS:  {np.sqrt(np.mean(q_rel_diff**2)) * 100:.6f} %")

        # 許容範囲チェック
        if np.allclose(q_py, q_jl, rtol=1e-4, atol=1e-3):
            print(f"\n✅ q_resultは許容範囲内で一致 (rtol=1e-4, atol=1e-3)")
        else:
            print(f"\n❌ q_resultに許容範囲を超える差があります")

            # 差が大きい点を特定
            large_diff_mask = q_rel_diff > 0.01  # 1%以上の相対差
            if np.any(large_diff_mask):
                n_large = np.sum(large_diff_mask)
                total = q_py.size
                print(f"\n  相対差1%以上の点: {n_large}/{total} ({n_large/total*100:.2f}%)")

    else:
        print(f"  ⚠️ 形状が異なります")

    # T_resultの比較
    T_py = python_data['T_result']
    T_jl = julia_data['T_result']
    print(f"\nT_result形状:")
    print(f"  Python: {T_py.shape}")
    print(f"  Julia:  {T_jl.shape}")

    if T_py.shape == T_jl.shape:
        T_diff = np.abs(T_py - T_jl)

        print(f"\nT_result統計:")
        print(f"  Python範囲: [{np.min(T_py):.2f}, {np.max(T_py):.2f}] K")
        print(f"  Julia範囲:  [{np.min(T_jl):.2f}, {np.max(T_jl):.2f}] K")
        print(f"\n  絶対差:")
        print(f"    最大: {np.max(T_diff):.6e} K")
        print(f"    平均: {np.mean(T_diff):.6e} K")
        print(f"    RMS:  {np.sqrt(np.mean(T_diff**2)):.6e} K")

        if np.allclose(T_py, T_jl, rtol=1e-6, atol=1e-4):
            print(f"\n✅ T_resultは許容範囲内で一致 (rtol=1e-6, atol=1e-4)")
        else:
            print(f"\n❌ T_resultに許容範囲を超える差があります")
    else:
        print(f"  ⚠️ 形状が異なります")

    print("\n" + "="*80)
    print("比較完了")
    print("="*80)

if __name__ == "__main__":
    main()
