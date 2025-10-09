"""
SlidingWindowSolver.jl

Phase 5: スライディングウィンドウCGM計算

長時間の逆熱伝導問題を複数のウィンドウに分割し、各ウィンドウでCGM最適化を実行:
1. ウィンドウ分割（時間方向の切り出し）
2. オーバーラップ領域の初期値継承
3. 各ウィンドウでのCGM実行（Phase 4ソルバー呼び出し）
4. ウィンドウ間のオーバーラップ平均化
5. 温度場の継承（次ウィンドウの初期値として使用）
6. 全ウィンドウ結果の連結

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
- sliding_window_CGM_q_saving() (1556-1626行)

主要関数:
- solve_sliding_window_cgm: スライディングウィンドウCGMメインループ
"""

module SlidingWindowSolver

using Printf

# 共通モジュール
using ..Commons: WorkBuffers

# Phase 4モジュールは親モジュールで既にinclude済み
using ..CGMSolver

export solve_sliding_window_cgm, WindowInfo


"""
ウィンドウ情報構造体

各ウィンドウの計算情報を保持
"""
struct WindowInfo
  window_id::Int           # ウィンドウID（1始まり）
  start_idx::Int           # 開始インデックス（0始まり）
  end_idx::Int             # 終了インデックス（0始まり）
  max_L::Int               # ウィンドウ長（時間ステップ数）
  n_iterations::Int        # CGM反復数
  J_final::Float64         # 最終目的関数値
  overlap_steps::Int       # オーバーラップステップ数
end


"""
    solve_sliding_window_cgm(
      Y_obs, T0, dx, dy, dt, rho, cp_coeffs, k_coeffs,
      window_size, overlap, q_init_value, cgm_iteration
    ) -> (q_global, windows_info)

スライディングウィンドウCGM計算

Pythonオリジナル: sliding_window_CGM_q_saving() (1556-1626行)

アルゴリズム:
1. ウィンドウループ:
   - 開始インデックスから観測データを切り出し
   - 前ウィンドウの結果を初期値として継承
   - CGM最適化を実行
   - 結果を保存

2. オーバーラップ平均化:
   - 重複部分を0.5*old + 0.5*newで平均化
   - 連続性を確保

3. 温度場の継承:
   - 各ウィンドウの最終温度場を次ウィンドウの初期値とする

Args:
  Y_obs: 観測温度（底面） (ni, nj, nt)
  T0: 初期温度場 (ni, nj, nk)
  work: Heat3d用配列群
  dx, dy: x, y方向格子幅 [m]
  Z: z方向格子配列 (nk,) [m]
  ΔZ: CV幅 (nk,) [m]
  dt: 時間刻み [s]
  rho: 密度 [kg/m³]
  cp_coeffs: 比熱多項式係数 (4,)
  k_coeffs: 熱伝導率多項式係数 (4,)
  window_size: ウィンドウサイズ（時間ステップ数）
  overlap: オーバーラップ領域（時間ステップ数）
  q_init_value: 初期熱流束の初期値 [W/m²]
  cgm_iteration: CGM最大反復数

Returns:
  q_global: 全時間の逆解析熱流束 (ni, nj, nt-1)
  windows_info: 各ウィンドウの詳細情報
"""
function solve_sliding_window_cgm(
  Y_obs::Array{Float64,3}, T0::Array{Float64,3}, work::WorkBuffers,
  dx::Float64, dy::Float64, Z::Vector{Float64}, ΔZ::Vector{Float64},
  dt::Float64,
  rho::Float64, cp_coeffs::Vector{Float64}, k_coeffs::Vector{Float64},
  window_size::Int, overlap::Int, q_init_value::Float64, cgm_iteration::Int;
  rtol_dhcp::Float64=1e-6,
  maxiter_dhcp::Int=20000,
  rtol_adjoint::Float64=1e-8,
  maxiter_adjoint::Int=20000,
  verbose::Bool=false
)

  ni, nj, nt = size(Y_obs)  # Y_obs形状: (ni, nj, nt)
  T_init = copy(T0)
  nk = size(T_init, 3)

  start_idx = 0  # 0始まり（Pythonと同じ）
  q_total = []   # Vector{Array{Float64,3}}型 CHECK: メモリ確保必要？？
  prev_q_win = nothing

  windows_info = WindowInfo[]

  safety_counter = 0
  safety_limit = nt * 5

  if verbose
    println("\n=== Starting Sliding Window CGM calculation ===")
    println("Total time steps: $nt")
    println("Window size: $window_size")
    println("Overlap: $overlap")
  end

  while start_idx < nt - 1
    safety_counter += 1
    if safety_counter > safety_limit
      if verbose
        println("[Warning] Safety counter exceeded") # CHECK 何のため？？
      end
      break
    end

    # 現在のウィンドウ長（Pythonオリジナル: 1577-1580行）
    max_L = min(window_size, (nt - 1) - start_idx)
    end_idx = start_idx + max_L
    Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]  # Y_obs形状: (ni, nj, nt)

    if verbose
      println("\n--- Window $(length(windows_info)+1): [$start_idx, $end_idx] (length=$max_L) ---")
    end

    # 初期熱流束の設定（Pythonオリジナル: 1583-1592行、メモリレイアウト最適化: Phase 2.2）
    if isnothing(prev_q_win)
      # 第1ウィンドウ: 一定値
      q_init_win = fill(q_init_value, (ni, nj, max_L))
      if verbose
        println("Initial heat flux: constant value $q_init_value W/m²")
      end
    else
      # 第2ウィンドウ以降: 前ウィンドウから継承
      q_init_win = zeros(Float64, ni, nj, max_L) # CHECK : 毎回確保している
      L_overlap = min(overlap, max_L, size(prev_q_win, 3))

      if L_overlap > 0
        # オーバーラップ部分: 前ウィンドウの最後L_overlapステップを継承
        q_init_win[:, :, 1:L_overlap] = prev_q_win[:, :, end-L_overlap+1:end]
        if verbose
          println("Overlap inheritance: previous $(L_overlap) steps")
        end
      end

      if L_overlap < max_L
        # 残り部分: 前ウィンドウの最終値で埋める
        edge = prev_q_win[:, :, end]
        for t in L_overlap+1:max_L
          q_init_win[:, :, t] = edge
        end
        if verbose
          println("Remaining part: filled with previous window's final value")
        end
      end
    end

    # CGM実行（Pythonオリジナル: 1595-1598行）
    # paramsでCGM反復数を指定
    cgm_params = (
      max_iter=cgm_iteration,
      verbose=false,
      rtol_dhcp=rtol_dhcp,
      maxiter_dhcp=maxiter_dhcp,
      rtol_adjoint=rtol_adjoint,
      maxiter_adjoint=maxiter_adjoint
    )
    q_win, T_cal_win, J_hist = solve_cgm!(
      T_init, Y_obs_win, q_init_win, work,
      dx, dy, Z, ΔZ,
      dt, rho, cp_coeffs, k_coeffs; params=cgm_params
    )

    prev_q_win = copy(q_win)

    # 結果の拼接（オーバーラップ平均化）（Pythonオリジナル: 1603-1612行、メモリレイアウト最適化: Phase 2.2）
    overlap_steps_actual = 0

    if length(q_total) == 0
      # 第1ウィンドウ: そのまま追加
      push!(q_total, q_win)
      if verbose
        println("First window: added as-is")
      end
    else
      # 第2ウィンドウ以降: オーバーラップ平均化
      overlap_steps_actual = min(overlap, size(q_win, 3), size(q_total[end], 3))

      if overlap_steps_actual > 0
        # オーバーラップ部分を平均化（Python原典に合わせた重み付け）
        # 注: Python原典は 0.5*old + new（数学的には1.5倍だが、意図的な重み付けと判断）
        # q_total[end]の最後overlap_steps_actualステップと
        # q_winの最初overlap_steps_actualステップを組み合わせ
        q_total[end][:, :, end-overlap_steps_actual+1:end] =
          0.5 * q_total[end][:, :, end-overlap_steps_actual+1:end] +
          q_win[:, :, 1:overlap_steps_actual]

        # q_winの残り部分を追加
        if overlap_steps_actual < size(q_win, 3)
          push!(q_total, q_win[:, :, overlap_steps_actual+1:end])
        end

        if verbose
          println("Overlap averaging: $(overlap_steps_actual) steps")
        end
      else
        # オーバーラップなし
        push!(q_total, q_win)
        if verbose
          println("No overlap: added as-is")
        end
      end
    end

    # 温度場の継承（Pythonオリジナル: 1614行、メモリレイアウト最適化: Phase 2.2）
    # solve_cgm!は(ni, nj, nk, nt)を返すので、最終時刻を取得
    T_init = copy(T_cal_win[:, :, :, end])

    # ウィンドウ情報の保存
    win_info = WindowInfo(
      length(windows_info) + 1,
      start_idx,
      end_idx,
      max_L,
      length(J_hist),
      J_hist[end],
      overlap_steps_actual
    )
    push!(windows_info, win_info)

    if verbose
      println("CGM iterations: $(length(J_hist)), final J: $(J_hist[end])")
    end

    # インデックス進行（Pythonオリジナル: 1619-1620行）
    step = max(1, max_L - overlap)
    start_idx += step
  end

  # 全ウィンドウ結果の連結（Pythonオリジナル: 1622-1623行、メモリレイアウト最適化: Phase 2.2）
  q_global = cat(q_total..., dims=3) # CHECK: in-placeにできないか

  # nt-1に切り詰め
  if size(q_global, 3) > nt - 1
    q_global = q_global[:, :, 1:nt-1]
  end

  if verbose
    println("\n=== Sliding Window calculation completed ===")
    println("Total windows: $(length(windows_info))")
    println("Final q_global shape: $(size(q_global))")
  end

  return q_global, windows_info
end

end  # module SlidingWindowSolver
