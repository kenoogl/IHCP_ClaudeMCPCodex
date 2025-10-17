#!/usr/bin/env julia
"""
DHCP対流境界条件テスト

Phase 6: 実際のDHCPソルバーで対流境界条件の動作確認

テストシナリオ:
- 小規模格子: 10×10×5
- Z-plus面（上面）: 対流境界条件 (h=10 W/m²·K, T_∞=300K)
- 他の面: 断熱条件
- 初期温度: 550K（均一）
- 時間ステップ: 5ステップ

期待される動作:
1. PBICGSTAB!が収束すること
2. Z-plus面から熱が逃げる（対流により冷却）
3. 温度場が物理的に妥当であること（300K〜550Kの範囲）
4. NaN/Infが発生しないこと

実行方法:
  julia --project=julia julia/scripts/test_dhcp_convection.jl
"""

using Printf
using Statistics

# プロジェクトのソースコードをロード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM
using .IHCP_CGM.BoundaryConditions
using .IHCP_CGM.DHCPSolver: solve_dhcp!, set_dhcp_bc_parameters_convection
using .IHCP_CGM.Commons

"""
z方向格子生成（非均等格子）
"""
function build_z_grid_simple(nk::Int, Lz::Float64, stretch_factor::Float64)
  # Step 1: z_faces生成
  z_faces_normalized = range(1.0, stop = 0.0, length = nk + 1)
  z_faces = Lz .- (Lz / (exp(stretch_factor) - 1.0)) .* (exp.(stretch_factor .* z_faces_normalized) .- 1.0)

  # Step 2: dzの計算
  dz = diff(z_faces)

  # Step 3: ΔZ配列（ガイドセル込み、nk+2個）
  ΔZ = zeros(nk + 2)
  ΔZ[2:nk+1] = dz[1:nk]
  ΔZ[1] = ΔZ[2]
  ΔZ[nk+2] = ΔZ[nk+1]

  # Step 4: z_centers（ガイドセル込み、nk+2個）
  z_centers = zeros(nk + 2)
  # ガイドセル（k=1, k=nk+2）
  z_centers[1] = z_faces[1]
  z_centers[nk+2] = z_faces[nk+1]

  # 実セル（k=2からk=nk+1）のセンター座標
  for k in 2:nk+1
    z_centers[k] = (z_faces[k-1] + z_faces[k]) * 0.5
  end

  return z_centers, ΔZ, z_faces, dz
end

"""
対流境界条件ありのDHCPテスト
"""
function test_dhcp_with_convection(h_conv::Float64 = 10.0, nt_steps::Int = 20)
  println("="^60)
  println("DHCP対流境界条件テスト（Phase 6）")
  println("="^60)

  # パラメータ設定
  ni, nj, nk = 10, 10, 5
  nt = nt_steps

  dx = 0.12e-3
  dy = 0.12e-3 * sin(deg2rad(80.0)) / sin(deg2rad(45.0))
  dt = 1e-3
  rho = 7823.493962874829

  # z方向格子（非均等）
  Lz = 0.5e-3
  stretch_factor = 3.0
  z_centers, ΔZ, z_faces, dz = build_z_grid_simple(nk, Lz, stretch_factor)

  # 熱物性値係数（SUS304、多項式）
  cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]

  # 初期温度場（550K均一）
  T_init = fill(550.0, ni, nj, nk)

  # Z上面の熱流束（ゼロ、対流境界条件のみ）
  q_surface = zeros(Float64, ni, nj, nt - 1)

  # WorkBuffers
  work = WorkBuffers(ni + 2, nj + 2, nk + 2)

  println("\n【テスト条件】")
  println("格子: $(ni)×$(nj)×$(nk)")
  println("時間ステップ: $(nt), dt=$(dt)s")
  println("初期温度: 550K（均一）")
  println("dx=$(dx*1e3) mm, dy=$(dy*1e3) mm, Lz=$(Lz*1e3) mm")
  println("密度: $(rho) kg/m³")

  println("\n【境界条件設定】")
  println("X-minus面: 断熱")
  println("X-plus面:  断熱")
  println("Y-minus面: 断熱")
  println("Y-plus面:  断熱")
  println("Z-minus面: 断熱")
  println("Z-plus面:  対流 (h=$(h_conv) W/m²·K, T_∞=300K)")

  println("\n【期待される動作】")
  println("1. PBICGSTAB!が収束")
  println("2. Z-plus面から熱が逃げる（対流冷却）")
  println("3. 温度が300K〜550Kの範囲内")
  println("4. NaN/Infが発生しない")

  println("\n" * "="^60)
  println("DHCP実行開始")
  println("="^60)

  # 対流境界条件を設定
  bc_set = set_dhcp_bc_parameters_convection(h_conv, 300.0)

  # DHCPソルバー実行
  start_time = time()

  T_all, iter_counts, residuals = solve_dhcp!(
    T_init,
    q_surface,
    work,
    nt,
    rho,
    cp_coeffs,
    k_coeffs,
    dx,
    dy,
    z_centers,
    ΔZ,  # ガイドセル込み配列を使用（nk+2要素）
    dt;
    rtol=1.0e-6,
    maxiter=20000,
    verbose=true,
    par="sequential",
    use_previous_solution=true,  # 前ステップを初期推定値として使用
    extrapolation=:linear,  # 線形外挿
    smoother=:gs,
    solver=:pbicgstab,
    bc_set=bc_set
  )

  elapsed_time = time() - start_time

  println("\n" * "="^60)
  println("DHCP実行完了")
  println("="^60)
  @printf("実行時間: %.6f 秒\n", elapsed_time)

  # 結果検証
  println("\n" * "="^60)
  println("結果検証")
  println("="^60)

  # 1. 収束性チェック
  println("\n【1. 収束性チェック】")
  println("各時間ステップの反復回数と残差:")
  for t in 1:nt
    if t == 1
      println("  t=$(t): 初期条件（反復なし）")
    else
      @printf("  t=%d: %d回, 残差=%.3e\n", t, iter_counts[t], residuals[t])
    end
  end

  max_iter = maximum(iter_counts[2:end])
  max_residual = maximum(residuals[2:end])
  avg_iter = sum(iter_counts[2:end]) / (nt - 1)
  @printf("\n統計情報:\n")
  @printf("  最大反復回数: %d回\n", max_iter)
  @printf("  平均反復回数: %.2f回\n", avg_iter)
  @printf("  最大残差: %.3e\n", max_residual)

  if max_iter < 20000
    println("✓ 全時間ステップで収束")
  else
    println("✗ 収束失敗")
  end

  # 2. 温度範囲チェック
  println("\n【2. 温度場の妥当性チェック】")
  for t in 1:nt
    T_min = minimum(T_all[:, :, :, t])
    T_max = maximum(T_all[:, :, :, t])
    T_mean = mean(T_all[:, :, :, t])
    @printf("  t=%d: 温度範囲 %.2f - %.2f K (平均 %.2f K)\n", t, T_min, T_max, T_mean)
  end

  T_min_all = minimum(T_all)
  T_max_all = maximum(T_all)

  println("\n全体の温度範囲: $(T_min_all) - $(T_max_all) K")

  if T_min_all >= 300.0 && T_max_all <= 550.0
    println("✓ 温度が物理的に妥当な範囲（300K〜550K）")
  else
    println("✗ 温度が範囲外")
  end

  # 3. NaN/Infチェック
  println("\n【3. 数値異常チェック】")
  has_nan = any(isnan.(T_all))
  has_inf = any(isinf.(T_all))

  if !has_nan && !has_inf
    println("✓ NaN/Inf検出なし")
  else
    println("✗ 数値異常検出")
    println("  NaN: $(has_nan)")
    println("  Inf: $(has_inf)")
  end

  # 4. 対流冷却効果チェック
  println("\n【4. 対流冷却効果チェック】")
  println("Z-plus面（k=$(nk)）の温度変化:")

  # 中央点(i=5, j=5, k=nk)の温度履歴
  i_center, j_center = div(ni, 2), div(nj, 2)
  T_surface = [T_all[i_center, j_center, nk, t] for t in 1:nt]

  println("  中央点(i=$(i_center), j=$(j_center), k=$(nk))の温度履歴:")
  for t in 1:nt
    @printf("    t=%d: %.2f K\n", t, T_surface[t])
  end

  T_drop = T_surface[1] - T_surface[end]
  @printf("\n温度降下: %.2f K (%.2f%%)\n", T_drop, T_drop / T_surface[1] * 100)

  if T_drop > 0
    println("✓ 対流による冷却効果を確認")
  else
    println("✗ 冷却効果が見られない")
  end

  # 5. Z方向温度分布（最終ステップ）
  println("\n【5. Z方向温度分布（t=$(nt)）】")
  println("  中央点(i=$(i_center), j=$(j_center))のZ方向温度:")
  for k in nk:-1:1
    @printf("    k=%d (z=%.4f mm): %.2f K\n", k, z_centers[k+1]*1e3, T_all[i_center, j_center, k, nt])
  end

  # 総合判定
  println("\n" * "="^60)
  println("総合判定")
  println("="^60)

  all_pass = (max_iter < 20000) &&
             (T_min_all >= 300.0) &&
             (T_max_all <= 550.0) &&
             (!has_nan) &&
             (!has_inf) &&
             (T_drop > 0)

  if all_pass
    println("✅ Phase 6テスト成功")
    println("   対流境界条件が正しく動作していることを確認")
  else
    println("❌ Phase 6テスト失敗")
    println("   対流境界条件の実装に問題がある可能性")
  end

  println("="^60)

  return T_all, iter_counts, all_pass
end

# メイン実行
if abspath(PROGRAM_FILE) == @__FILE__
  try
    # コマンドライン引数から h値を取得（デフォルト: 10.0）
    h_value = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 10.0
    nt_steps = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20

    println("\n実行パラメータ:")
    println("  h = $(h_value) W/m²·K")
    println("  時間ステップ数 = $(nt_steps)")
    println()

    T_all, iter_counts, success = test_dhcp_with_convection(h_value, nt_steps)

    if success
      println("\n✓ テスト完了: 成功")
      exit(0)
    else
      println("\n✗ テスト完了: 失敗")
      exit(1)
    end
  catch e
    println("\n❌ エラーが発生しました: ", e)
    Base.showerror(stdout, e, catch_backtrace())
    println()
    exit(1)
  end
end
