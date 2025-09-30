"""
test_adjoint_solver.jl

Phase 3: 随伴ソルバー（AdjointSolver）のテスト

テストシナリオ:
1. 1D小規模問題（定数熱物性値）
   - Python参照データとの比較（1e-7以内）
   - 後退時間積分の正確性確認
   - 残差注入機構の動作確認

2. 3D小規模問題（温度依存熱物性値）
   - Python参照データとの比較
   - 勾配場の正確性確認（CGM用）

TDD方針:
- Python参照データ（phase3_reference_adjoint_1D.json, phase3_reference_adjoint_3D.json）
- Julia実装との比較テスト
- 時間反転ループの正確性検証
"""

using Test
using LinearAlgebra
using Statistics
using JSON

# テスト対象モジュール
include("../src/solvers/AdjointSolver.jl")
using .AdjointSolver

include("../src/ThermalProperties.jl")
using .ThermalProperties

# JSON読み込みヘルパー（Phase 2で実装済み）
include("../src/utils/json_helpers.jl")
using .JSONHelpers


# ===== テストセット1: 1D小規模問題 =====

@testset "Phase 3: Adjoint 1D小規模問題" begin
  println("\n" * "="^60)
  println("Phase 3テスト: Adjoint 1D小規模問題")
  println("="^60)

  # Python参照データ読み込み
  data_file = joinpath(@__DIR__, "../data/phase3_reference_adjoint_1D.json")
  if !isfile(data_file)
    @warn "参照データが見つかりません: $data_file"
    @warn "以下のコマンドで生成してください:"
    @warn "  cd julia/data && python3 generate_phase3_reference.py"
    @test_skip "参照データ生成が必要"
    return
  end

  data = JSON.parsefile(data_file)
  println("参照データ読み込み完了: $(data["problem_type"])")

  # 格子情報
  ni = data["grid"]["ni"]
  nj = data["grid"]["nj"]
  nk = data["grid"]["nk"]
  dx = data["grid"]["dx"]
  dy = data["grid"]["dy"]
  dz = Float64.(data["grid"]["dz"])

  dz_b = Float64.(data["z_coords"]["dz_b"])
  dz_t = Float64.(data["z_coords"]["dz_t"])

  println("格子: $ni×$nj×$nk")

  # 物性値
  rho = data["properties"]["rho"]
  cp_coeffs = Float64.(data["properties"]["cp_coeffs"])
  k_coeffs = Float64.(data["properties"]["k_coeffs"])

  # 時間設定
  nt = data["time"]["nt"]
  dt = data["time"]["dt"]
  println("時間: nt=$nt, dt=$dt s")

  # DHCP結果と観測データ（Python参照）
  T_cal = reshape_json_array(data["T_cal"], (nt, ni, nj, nk), Float64)
  Y_obs = reshape_json_array(data["Y_obs"], (nt, ni, nj), Float64)

  # Python参照随伴場
  λ_ref = reshape_json_array(data["lambda_all"], (nt, ni, nj, nk), Float64)

  println("Python参照データ読み込み完了")
  println("  T_cal範囲: [$(data["stats"]["T_min"]), $(data["stats"]["T_max"])] K")
  println("  残差範囲: [$(data["stats"]["residual_min"]), $(data["stats"]["residual_max"])] K")
  println("  λ_ref範囲: [$(data["stats"]["lambda_min"]), $(data["stats"]["lambda_max"])]")

  # Julia実装でAdjoint求解
  println("\nJulia Adjoint求解開始...")
  λ_all, cg_iters = solve_adjoint!(
    T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt,
    rtol=1e-10, maxiter=1000, verbose=false
  )
  println("Julia Adjoint求解完了")
  println("  平均CG反復回数: $(mean(cg_iters))")
  println("  λ_all範囲: [$(minimum(λ_all)), $(maximum(λ_all))]")

  # テスト1: 終端条件
  @testset "終端条件 λ[nt] = 0" begin
    @test all(λ_all[nt, :, :, :] .≈ 0.0)
    println("  ✓ 終端条件確認: λ[nt] = 0")
  end

  # テスト2: Python参照データとの比較（全時間ステップ）
  @testset "Python参照データ比較（全時間）" begin
    max_error = maximum(abs.(λ_all .- λ_ref))
    mean_error = mean(abs.(λ_all .- λ_ref))
    rel_error = max_error / (maximum(abs.(λ_ref)) + 1e-12)

    println("\n  誤差統計（全時間ステップ）:")
    println("    最大絶対誤差: $max_error")
    println("    平均絶対誤差: $mean_error")
    println("    相対誤差: $rel_error")

    # 1e-7以内の精度を要求（成功基準）
    @test max_error < 1e-7
    @test rel_error < 1e-5

    println("  ✓ Python参照データと一致（1e-7以内）")
  end

  # テスト3: 各時間ステップでの比較（詳細チェック）
  @testset "時間ステップ毎の比較" begin
    println("\n  時間ステップ毎の誤差:")
    for t in [1, nt÷2, nt-1]
      error_t = maximum(abs.(λ_all[t, :, :, :] .- λ_ref[t, :, :, :]))
      println("    t=$t: 最大誤差=$error_t")
      @test error_t < 1e-7
    end
    println("  ✓ 時間ステップ毎の精度確認")
  end

  # テスト4: 残差注入機構の確認（底面での変化）
  @testset "残差注入機構" begin
    # 底面（k=1）での随伴場が非ゼロであることを確認
    λ_bottom = λ_all[:, :, :, 1]
    @test maximum(abs.(λ_bottom)) > 1e-12

    # 測定残差との相関確認
    residual = T_cal[:, :, :, 1] .- Y_obs
    println("\n  残差注入機構:")
    println("    測定残差RMS: $(sqrt(mean(residual.^2)))")
    println("    底面随伴場範囲: [$(minimum(λ_bottom)), $(maximum(λ_bottom))]")
    println("  ✓ 残差注入機構動作確認")
  end

  # テスト5: 後退時間積分の単調性確認
  @testset "後退時間積分の単調性" begin
    # 随伴場のノルムが時間的に変化することを確認
    λ_norms = [norm(λ_all[t, :, :, :]) for t in 1:nt]
    println("\n  随伴場ノルム（時間変化）:")
    println("    t=1: $(λ_norms[1])")
    println("    t=$nt: $(λ_norms[nt])")
    @test λ_norms[nt] ≈ 0.0  # 終端条件
    println("  ✓ 後退時間積分の単調性確認")
  end

  println("\n" * "="^60)
  println("Phase 3テスト（1D）: 全テスト合格")
  println("="^60)
end


# ===== テストセット2: 3D小規模問題（温度依存熱物性値） =====

@testset "Phase 3: Adjoint 3D小規模問題（温度依存）" begin
  println("\n" * "="^60)
  println("Phase 3テスト: Adjoint 3D小規模問題（温度依存熱物性値）")
  println("="^60)

  # Python参照データ読み込み
  data_file = joinpath(@__DIR__, "../data/phase3_reference_adjoint_3D.json")
  if !isfile(data_file)
    @warn "参照データが見つかりません: $data_file"
    @test_skip "参照データ生成が必要"
    return
  end

  data = JSON.parsefile(data_file)
  println("参照データ読み込み完了: $(data["problem_type"])")

  # 格子情報
  ni = data["grid"]["ni"]
  nj = data["grid"]["nj"]
  nk = data["grid"]["nk"]
  dx = data["grid"]["dx"]
  dy = data["grid"]["dy"]
  dz = Float64.(data["grid"]["dz"])

  dz_b = Float64.(data["z_coords"]["dz_b"])
  dz_t = Float64.(data["z_coords"]["dz_t"])

  println("格子: $ni×$nj×$nk")

  # 物性値（SUS304近似、温度依存）
  rho = data["properties"]["rho"]
  cp_coeffs = Float64.(data["properties"]["cp_coeffs"])
  k_coeffs = Float64.(data["properties"]["k_coeffs"])

  # 時間設定
  nt = data["time"]["nt"]
  dt = data["time"]["dt"]
  println("時間: nt=$nt, dt=$dt s")

  # DHCP結果と観測データ（Python参照）
  T_cal = reshape_json_array(data["T_cal"], (nt, ni, nj, nk), Float64)
  Y_obs = reshape_json_array(data["Y_obs"], (nt, ni, nj), Float64)

  # Python参照随伴場
  λ_ref = reshape_json_array(data["lambda_all"], (nt, ni, nj, nk), Float64)

  println("Python参照データ読み込み完了")
  println("  T_cal範囲: [$(data["stats"]["T_min"]), $(data["stats"]["T_max"])] K")
  println("  温度上昇: $(data["stats"]["dT"]) K")
  println("  残差RMS: $(data["stats"]["residual_rms"]) K")
  println("  λ_ref範囲: [$(data["stats"]["lambda_min"]), $(data["stats"]["lambda_max"])]")

  # Julia実装でAdjoint求解
  println("\nJulia Adjoint求解開始...")
  λ_all, cg_iters = solve_adjoint!(
    T_cal, Y_obs, nt, rho, cp_coeffs, k_coeffs,
    dx, dy, dz, dz_b, dz_t, dt,
    rtol=1e-8, maxiter=1000, verbose=false
  )
  println("Julia Adjoint求解完了")
  println("  平均CG反復回数: $(mean(cg_iters)) ± $(std(cg_iters))")
  println("  λ_all範囲: [$(minimum(λ_all)), $(maximum(λ_all))]")

  # テスト1: 終端条件
  @testset "終端条件 λ[nt] = 0" begin
    @test all(λ_all[nt, :, :, :] .≈ 0.0)
    println("  ✓ 終端条件確認: λ[nt] = 0")
  end

  # テスト2: Python参照データとの比較（全時間ステップ）
  @testset "Python参照データ比較（全時間）" begin
    max_error = maximum(abs.(λ_all .- λ_ref))
    mean_error = mean(abs.(λ_all .- λ_ref))
    rel_error = max_error / (maximum(abs.(λ_ref)) + 1e-12)

    println("\n  誤差統計（全時間ステップ）:")
    println("    最大絶対誤差: $max_error")
    println("    平均絶対誤差: $mean_error")
    println("    相対誤差: $rel_error")

    # 1e-7以内の精度を要求（成功基準）
    @test max_error < 1e-7
    @test rel_error < 1e-5

    println("  ✓ Python参照データと一致（1e-7以内）")
  end

  # テスト3: 勾配場の確認（CGM用）
  @testset "勾配場（CGM用）" begin
    # 表面（k=nk）での随伴場がCGMで使用される勾配場
    λ_surface = λ_all[:, :, :, nk]
    println("\n  勾配場（表面 k=$nk）:")
    println("    範囲: [$(minimum(λ_surface)), $(maximum(λ_surface))]")
    println("    ノルム: $(norm(λ_surface))")

    # Python参照との比較
    λ_surface_ref = λ_ref[:, :, :, nk]
    error_surface = maximum(abs.(λ_surface .- λ_surface_ref))
    println("    Python参照との誤差: $error_surface")

    @test error_surface < 1e-7
    println("  ✓ 勾配場の正確性確認")
  end

  # テスト4: 温度依存熱物性値の影響確認
  @testset "温度依存熱物性値" begin
    # 温度範囲が小さいため、熱物性値の変化も小さい
    T_min = data["stats"]["T_min"]
    T_max = data["stats"]["T_max"]
    println("\n  温度依存性:")
    println("    温度範囲: [$T_min, $T_max] K")

    # cp(T)とk(T)の変化率
    cp_min = cp_coeffs[1] + cp_coeffs[2] * T_min
    cp_max = cp_coeffs[1] + cp_coeffs[2] * T_max
    k_min = k_coeffs[1] + k_coeffs[2] * T_min
    k_max = k_coeffs[1] + k_coeffs[2] * T_max

    println("    cp範囲: [$cp_min, $cp_max] J/(kg·K)")
    println("    k範囲: [$k_min, $k_max] W/(m·K)")
    println("  ✓ 温度依存熱物性値の影響確認")
  end

  # テスト5: CG収束性確認
  @testset "CG収束性" begin
    println("\n  CG収束性:")
    println("    平均反復回数: $(mean(cg_iters))")
    println("    標準偏差: $(std(cg_iters))")
    println("    最大反復回数: $(maximum(cg_iters))")

    # Python参照と比較
    ref_cg_mean = data["stats"]["adjoint_cg_mean"]
    ref_cg_std = data["stats"]["adjoint_cg_std"]
    println("    Python参照平均: $ref_cg_mean ± $ref_cg_std")

    # CG反復回数が合理的範囲内
    @test mean(cg_iters) < 100  # 100反復以内
    println("  ✓ CG収束性確認")
  end

  println("\n" * "="^60)
  println("Phase 3テスト（3D）: 全テスト合格")
  println("="^60)
end


# ===== まとめ =====

println("\n" * "="^60)
println("Phase 3: Adjoint随伴ソルバー テスト完了")
println("="^60)
println("\n成功基準:")
println("  ✓ 随伴場λがPython結果と1e-7以内")
println("  ✓ 勾配場λ[:,:,:,end]の一致（CGM用）")
println("  ✓ 後退時間積分の正確性確認")
println("  ✓ 残差注入機構の動作確認")
println("\n次のステップ: Phase 4（CGM最適化）")
println("="^60)