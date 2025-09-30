#!/usr/bin/env python3
"""
Phase 2参照データ生成スクリプト（DHCP直接ソルバー）

このスクリプトは、Julia実装の検証用に小規模な参照データを生成します：
1. 製作解問題（解析解が既知）
2. 温度依存熱物性値問題（実問題に近い条件）

生成データ：
- 係数行列（疎行列、IJV形式）
- RHSベクトル
- 温度場の時間発展
- CG収束履歴
"""

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import cg, LinearOperator
from scipy.sparse import diags
import json
import os

# ===== Pythonオリジナルコードからの関数抽出 =====

def thermal_properties_calculator_simple(T, cp_coeffs, k_coeffs):
    """
    簡易版熱物性値計算（Numba依存を除去）

    Args:
        T: 温度場 (ni, nj, nk)
        cp_coeffs: 比熱多項式係数 [c0, c1, c2, c3]
        k_coeffs: 熱伝導率多項式係数 [k0, k1, k2, k3]

    Returns:
        cp: 比熱場 (ni, nj, nk)
        k: 熱伝導率場 (ni, nj, nk)
    """
    ni, nj, nk = T.shape
    cp = np.zeros_like(T)
    k = np.zeros_like(T)

    for i in range(ni):
        for j in range(nj):
            for k_idx in range(nk):
                T_val = T[i, j, k_idx]
                # Hornerの方法で多項式評価
                cp[i, j, k_idx] = cp_coeffs[0] + T_val * (cp_coeffs[1] + T_val * (cp_coeffs[2] + T_val * cp_coeffs[3]))
                k[i, j, k_idx] = k_coeffs[0] + T_val * (k_coeffs[1] + T_val * (k_coeffs[2] + T_val * k_coeffs[3]))

    return cp, k


def coeffs_and_rhs_building_DHCP(T_initial, q_surface, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt):
    """
    直接熱伝導問題（DHCP）の係数とRHS構築

    Pythonオリジナルコード（1087-1129行）のNumba最適化を除去したバージョン

    陰解法（後退差分）:
        (I - α∇²)T^(n+1) = T^n + q_boundary
        α = (k·dt) / (ρ·cp·dx²)

    Args:
        T_initial: 初期温度場 (ni, nj, nk)
        q_surface: 表面熱流束 (ni, nj)
        rho: 密度 [kg/m³]
        cp: 比熱場 (ni, nj, nk)
        k: 熱伝導率場 (ni, nj, nk)
        dx, dy: x, y方向格子幅 [m]
        dz: z方向格子幅配列 (nk,) [m]
        dz_b: 下側界面距離 (nk,) [m]
        dz_t: 上側界面距離 (nk,) [m]
        dt: 時間刻み [s]

    Returns:
        a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列 (N,)
        b: RHSベクトル (N,)
    """
    ni, nj, nk = T_initial.shape
    N = ni * nj * nk

    a_w = np.zeros(N)
    a_e = np.zeros(N)
    a_s = np.zeros(N)
    a_n = np.zeros(N)
    a_b = np.zeros(N)
    a_t = np.zeros(N)
    a_p = np.zeros(N)
    b   = np.zeros(N)

    for p in range(N):
        # Fortran順序（列優先）のインデックス変換
        i = p % ni
        j = (p // ni) % nj
        k_idx = p // (ni * nj)

        dz_k = dz[k_idx]
        dz_t_k = dz_t[k_idx]
        dz_b_k = dz_b[k_idx]

        k_p = k[i, j, k_idx]

        # 時間項（蓄熱項）
        a_p_0 = rho * cp[i, j, k_idx] * dx * dy * dz_k / dt

        # 6方向の熱伝導係数（調和平均）
        a_w[p] = 0.0 if j == 0      else (2 * k_p * k[i, j-1, k_idx] / (k_p + k[i, j-1, k_idx])) * dy * dz_k / dx
        a_e[p] = 0.0 if j == nj-1   else (2 * k_p * k[i, j+1, k_idx] / (k_p + k[i, j+1, k_idx])) * dy * dz_k / dx
        a_s[p] = 0.0 if i == 0      else (2 * k_p * k[i-1, j, k_idx] / (k_p + k[i-1, j, k_idx])) * dx * dz_k / dy
        a_n[p] = 0.0 if i == ni-1   else (2 * k_p * k[i+1, j, k_idx] / (k_p + k[i+1, j, k_idx])) * dx * dz_k / dy
        a_b[p] = 0.0 if k_idx == 0  else (2 * k_p * k[i, j, k_idx-1] / (k_p + k[i, j, k_idx-1])) * dx * dy / dz_b_k
        a_t[p] = 0.0 if k_idx == nk-1 else (2 * k_p * k[i, j, k_idx+1] / (k_p + k[i, j, k_idx+1])) * dx * dy / dz_t_k

        # 対角項
        a_p[p] = a_w[p] + a_e[p] + a_s[p] + a_n[p] + a_b[p] + a_t[p] + a_p_0

        # RHS（初期温度項 + 境界条件）
        rhs = a_p_0 * T_initial[i, j, k_idx]

        # 表面（上端）での熱流束境界条件
        if k_idx == nk - 1:
            rhs += q_surface[i, j] * dx * dy

        b[p] = rhs

    return a_w, a_e, a_s, a_n, a_b, a_t, a_p, b


def assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p):
    """
    DHCPの疎行列を組み立て

    Pythonオリジナルコード（1131-1153行）と同一

    7点ステンシル（対角+6方向）のCSR疎行列を構築
    COO形式で構築してCSRに変換（小規模格子での重複オフセット問題を回避）

    Args:
        ni, nj, nk: 格子点数
        a_w, a_e, a_s, a_n, a_b, a_t, a_p: 係数配列

    Returns:
        A_csr: CSR形式疎行列 (N, N)
    """
    N = ni * nj * nk

    # オフセット計算
    off_i = 1           # i方向（南北）
    off_j = ni          # j方向（西東）
    off_k = ni * nj     # k方向（上下）

    # COO形式で係数を収集
    I_list = []
    J_list = []
    V_list = []

    for p in range(N):
        # 対角成分
        I_list.append(p)
        J_list.append(p)
        V_list.append(a_p[p])

        # i-1 (南)
        if p >= off_i:
            I_list.append(p)
            J_list.append(p - off_i)
            V_list.append(-a_s[p])

        # i+1 (北)
        if p < N - off_i:
            I_list.append(p)
            J_list.append(p + off_i)
            V_list.append(-a_n[p])

        # j-1 (西)
        if p >= off_j:
            I_list.append(p)
            J_list.append(p - off_j)
            V_list.append(-a_w[p])

        # j+1 (東)
        if p < N - off_j:
            I_list.append(p)
            J_list.append(p + off_j)
            V_list.append(-a_e[p])

        # k-1 (下)
        if p >= off_k:
            I_list.append(p)
            J_list.append(p - off_k)
            V_list.append(-a_b[p])

        # k+1 (上)
        if p < N - off_k:
            I_list.append(p)
            J_list.append(p + off_k)
            V_list.append(-a_t[p])

    # COO→CSR変換
    A_csr = sp.coo_matrix((V_list, (I_list, J_list)), shape=(N, N)).tocsr()

    return A_csr


def multiple_time_step_solver_DHCP(T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
                                   dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000):
    """
    複数時間ステップDHCPソルバー

    Pythonオリジナルコード（1156-1219行）の簡易版

    Args:
        T_initial: 初期温度場 (ni, nj, nk)
        q_surface: 表面熱流束時系列 (nt-1, ni, nj)
        nt: 時間ステップ数
        rho: 密度 [kg/m³]
        cp_coeffs, k_coeffs: 熱物性値多項式係数
        dx, dy, dz, dz_b, dz_t: 格子情報
        dt: 時間刻み [s]
        rtol: CG相対許容誤差
        maxiter: CG最大反復回数

    Returns:
        T_all: 温度場時系列 (nt, ni, nj, nk)
        cg_iters: CG反復回数履歴 (nt-1,)
    """
    ni, nj, nk = T_initial.shape
    N = ni * nj * nk

    T_all = np.zeros((nt, ni, nj, nk))
    T_all[0] = T_initial.copy()

    cg_iters = []

    # ホットスタート用初期推定値
    x0 = T_initial.ravel(order='F').copy()

    for t in range(1, nt):
        # 前ステップ温度から熱物性値計算
        cp, k = thermal_properties_calculator_simple(T_all[t-1], cp_coeffs, k_coeffs)

        # 係数とRHS構築
        a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = coeffs_and_rhs_building_DHCP(
            T_all[t-1], q_surface[t-1], rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
        )

        # 疎行列組み立て
        A_csr = assemble_A_DHCP(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

        # 対角前処理器
        diag = A_csr.diagonal()
        inv_diag = np.where(diag != 0.0, 1.0/diag, 0.0)

        def Mv(z):
            return inv_diag * z

        M = LinearOperator(A_csr.shape, matvec=Mv)

        # CG法求解
        iter_count = [0]  # クロージャでカウント
        def callback(xk):
            iter_count[0] += 1

        x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0, callback=callback)

        if info > 0:
            print(f"[警告][t={t}] CG法が収束しませんでした (info={info}, iter={iter_count[0]})")
        elif info < 0:
            print(f"[エラー][t={t}] CG法の入力が不正です (info={info})")

        # 結果保存とホットスタート更新
        T_all[t] = x.reshape((ni, nj, nk), order='F')
        x0 = x.copy()
        cg_iters.append(iter_count[0])

    return T_all, np.array(cg_iters)


# ===== 参照データ生成関数 =====

def generate_manufactured_solution():
    """
    製作解問題1: 1D定常熱伝導（解析解が既知）

    問題設定:
        - 1×1×10の1D問題
        - 定数熱物性値（ρ=8000, cp=500, k=15）
        - 初期温度: T(z) = 300 + 10*z/Lz
        - 表面熱流束: q = 1000 W/m²
        - 10ステップ実行

    解析解（定常状態）:
        T(z) = T0 + (q/k) * z
    """
    print("\n" + "="*60)
    print("製作解問題1: 1D定常熱伝導")
    print("="*60)

    # 格子設定
    ni, nj, nk = 1, 1, 10
    dx, dy = 1e-3, 1e-3
    Lz = 5e-3

    # 等間隔格子
    z_faces = np.linspace(0, Lz, nk+1)
    z_centers = (z_faces[:-1] + z_faces[1:]) / 2
    dz = np.diff(z_faces)

    # 界面距離（等間隔なので単純）
    dz_t = np.full(nk, dz[0])
    dz_t[-1] = np.inf
    dz_b = np.full(nk, dz[0])
    dz_b[0] = np.inf

    # 物性値（定数）
    rho = 8000.0  # kg/m³
    cp_const = 500.0  # J/(kg·K)
    k_const = 15.0    # W/(m·K)

    # 多項式係数（定数近似）
    cp_coeffs = np.array([cp_const, 0.0, 0.0, 0.0])
    k_coeffs = np.array([k_const, 0.0, 0.0, 0.0])

    # 初期条件（線形温度分布）
    T0 = 300.0
    T_initial = np.zeros((ni, nj, nk))
    for k_idx in range(nk):
        T_initial[0, 0, k_idx] = T0 + 10 * z_centers[k_idx] / Lz

    # 境界条件（一定熱流束）
    q = 1000.0  # W/m²
    nt = 11
    dt = 0.1  # s
    q_surface = np.full((nt-1, ni, nj), q)

    # 求解
    T_all, cg_iters = multiple_time_step_solver_DHCP(
        T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, dt, rtol=1e-10, maxiter=1000
    )

    # 解析解（定常状態）
    T_analytical = T0 + (q / k_const) * z_centers
    T_final = T_all[-1, 0, 0, :]

    error = np.abs(T_final - T_analytical)
    rel_error = error / np.abs(T_analytical)

    print(f"格子: {ni}×{nj}×{nk}, 時間ステップ: {nt}")
    print(f"CG平均反復回数: {np.mean(cg_iters):.1f}")
    print(f"最終温度場誤差: max={np.max(error):.2e}, mean={np.mean(error):.2e}")
    print(f"相対誤差: max={np.max(rel_error):.2e}, mean={np.mean(rel_error):.2e}")

    # データ保存
    data = {
        'problem_type': 'manufactured_1D_steady',
        'grid': {'ni': ni, 'nj': nj, 'nk': nk, 'dx': dx, 'dy': dy, 'dz': dz.tolist()},
        'z_coords': {'faces': z_faces.tolist(), 'centers': z_centers.tolist(), 'dz_b': dz_b.tolist(), 'dz_t': dz_t.tolist()},
        'properties': {'rho': rho, 'cp_coeffs': cp_coeffs.tolist(), 'k_coeffs': k_coeffs.tolist()},
        'boundary': {'q_surface': q, 'T0': T0},
        'time': {'nt': nt, 'dt': dt},
        'T_initial': T_initial.tolist(),
        'T_all': T_all.tolist(),
        'T_analytical': T_analytical.tolist(),
        'cg_iters': cg_iters.tolist(),
        'error': {'max': float(np.max(error)), 'mean': float(np.mean(error)), 'max_rel': float(np.max(rel_error))}
    }

    return data


def generate_small_3d_problem():
    """
    製作解問題2: 小規模3D問題（温度依存熱物性値）

    問題設定:
        - 5×5×10の3D問題
        - 温度依存熱物性値（SUS304近似）
        - 初期温度: T = 300K
        - 表面熱流束: q = 5000 * sin(π*x/Lx) * sin(π*y/Ly) W/m²
        - 20ステップ実行
    """
    print("\n" + "="*60)
    print("製作解問題2: 小規模3D問題（温度依存熱物性値）")
    print("="*60)

    # 格子設定
    ni, nj, nk = 5, 5, 10
    dx, dy = 0.12e-3, 0.12e-3
    Lz = 0.5e-3

    # 非等間隔格子（表面集中）
    stretch_factor = 3
    z_faces_norm = np.linspace(1, 0, nk+1)
    z_faces = Lz - (Lz / (np.exp(stretch_factor) - 1)) * (np.exp(stretch_factor * z_faces_norm) - 1)

    z_centers = np.zeros(nk)
    z_centers[0] = z_faces[0]
    z_centers[-1] = z_faces[-1]
    z_centers[1:-1] = (z_faces[1:-2] + z_faces[2:-1]) / 2

    dz = np.diff(z_faces)

    dz_t = np.zeros(nk)
    dz_t[-1] = np.inf
    dz_t[:-1] = z_centers[1:] - z_centers[:-1]

    dz_b = np.zeros(nk)
    dz_b[0] = np.inf
    dz_b[1:] = z_centers[1:] - z_centers[:-1]

    # 物性値（SUS304近似、温度依存）
    rho = 7900.0  # kg/m³
    # cp(T) = 462 + 0.134*T [J/(kg·K)]
    # k(T) = 14.6 + 0.0127*T [W/(m·K)]
    cp_coeffs = np.array([462.0, 0.134, 0.0, 0.0])
    k_coeffs = np.array([14.6, 0.0127, 0.0, 0.0])

    # 初期条件（一定温度）
    T0 = 300.0
    T_initial = np.full((ni, nj, nk), T0)

    # 境界条件（空間変動熱流束）
    nt = 21
    dt = 1e-3  # 1 ms
    q_surface = np.zeros((nt-1, ni, nj))

    Lx, Ly = dx * ni, dy * nj
    x = np.linspace(0, Lx, ni)
    y = np.linspace(0, Ly, nj)
    X, Y = np.meshgrid(x, y, indexing='ij')

    q_base = 5000.0  # W/m²
    for t in range(nt-1):
        q_surface[t] = q_base * np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)

    # 求解
    T_all, cg_iters = multiple_time_step_solver_DHCP(
        T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
        dx, dy, dz, dz_b, dz_t, dt, rtol=1e-8, maxiter=1000
    )

    print(f"格子: {ni}×{nj}×{nk}, 時間ステップ: {nt}")
    print(f"CG平均反復回数: {np.mean(cg_iters):.1f} ± {np.std(cg_iters):.1f}")
    print(f"CG最大反復回数: {np.max(cg_iters)}")
    print(f"温度範囲: {np.min(T_all):.2f} - {np.max(T_all):.2f} K")
    print(f"温度上昇: {np.max(T_all) - T0:.2f} K")

    # データ保存
    data = {
        'problem_type': 'small_3D_temperature_dependent',
        'grid': {'ni': ni, 'nj': nj, 'nk': nk, 'dx': dx, 'dy': dy, 'dz': dz.tolist()},
        'z_coords': {'faces': z_faces.tolist(), 'centers': z_centers.tolist(), 'dz_b': dz_b.tolist(), 'dz_t': dz_t.tolist()},
        'properties': {'rho': rho, 'cp_coeffs': cp_coeffs.tolist(), 'k_coeffs': k_coeffs.tolist()},
        'boundary': {'q_surface_base': q_base, 'T0': T0},
        'time': {'nt': nt, 'dt': dt},
        'T_initial': T_initial.tolist(),
        'q_surface': q_surface.tolist(),
        'T_all': T_all.tolist(),
        'cg_iters': cg_iters.tolist(),
        'stats': {
            'T_min': float(np.min(T_all)),
            'T_max': float(np.max(T_all)),
            'dT': float(np.max(T_all) - T0),
            'cg_mean': float(np.mean(cg_iters)),
            'cg_std': float(np.std(cg_iters)),
            'cg_max': int(np.max(cg_iters))
        }
    }

    return data


def main():
    """メイン処理"""
    print("="*60)
    print("Phase 2参照データ生成スクリプト（DHCP直接ソルバー）")
    print("="*60)

    # 出力ディレクトリ確認
    output_dir = os.path.dirname(os.path.abspath(__file__))
    print(f"出力ディレクトリ: {output_dir}")

    # 製作解問題1
    data1 = generate_manufactured_solution()
    output_file1 = os.path.join(output_dir, 'phase2_reference_1D_steady.json')
    with open(output_file1, 'w') as f:
        json.dump(data1, f, indent=2)
    print(f"✓ 保存完了: {output_file1}")

    # 製作解問題2
    data2 = generate_small_3d_problem()
    output_file2 = os.path.join(output_dir, 'phase2_reference_3D_small.json')
    with open(output_file2, 'w') as f:
        json.dump(data2, f, indent=2)
    print(f"✓ 保存完了: {output_file2}")

    print("\n" + "="*60)
    print("参照データ生成完了")
    print("="*60)
    print("\nJuliaテストでの使用方法:")
    print("  julia> using JSON")
    print("  julia> data = JSON.parsefile(\"phase2_reference_1D_steady.json\")")
    print("  julia> T_all = data[\"T_all\"]")
    print("="*60)


if __name__ == '__main__':
    main()