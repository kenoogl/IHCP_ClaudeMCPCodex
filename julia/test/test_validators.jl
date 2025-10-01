"""
test_validators.jl

検証器関数群のテストスイート
"""

using Test
using IHCP_CGM: Validators

@testset "Validators Module Tests" begin

  # ============================================================================
  # 1. check_field_finite: NaN/Inf検出
  # ============================================================================
  @testset "check_field_finite" begin
    # Test 1: 正常値（すべて有限）
    @test Validators.check_field_finite([1.0, 2.0, 3.0])
    @test Validators.check_field_finite(zeros(10, 10))
    @test Validators.check_field_finite(ones(5, 5, 5))

    # Test 2: NaN検出
    @test !Validators.check_field_finite([1.0, NaN, 3.0])
    @test !Validators.check_field_finite([NaN, 2.0, 3.0])
    @test !Validators.check_field_finite([1.0, 2.0, NaN])

    # Test 3: Inf検出
    @test !Validators.check_field_finite([1.0, Inf, 3.0])
    @test !Validators.check_field_finite([1.0, -Inf, 3.0])
    @test !Validators.check_field_finite([Inf, 2.0, 3.0])

    # Test 4: 多次元配列でのNaN/Inf検出
    field_with_nan = zeros(3, 3)
    field_with_nan[2, 2] = NaN
    @test !Validators.check_field_finite(field_with_nan)

    field_with_inf = ones(5, 5, 5)
    field_with_inf[3, 3, 3] = Inf
    @test !Validators.check_field_finite(field_with_inf)
  end

  # ============================================================================
  # 2. check_temperature_range: 温度範囲チェック
  # ============================================================================
  @testset "check_temperature_range" begin
    # Test 1: 正常範囲（150K ~ 3000K）
    is_valid, min_val, max_val = Validators.check_temperature_range([300.0, 500.0, 800.0], 150.0, 3000.0)
    @test is_valid
    @test min_val ≈ 300.0
    @test max_val ≈ 800.0

    # Test 2: 正常範囲（デフォルト引数）
    is_valid, min_val, max_val = Validators.check_temperature_range([300.0, 500.0, 800.0])
    @test is_valid
    @test min_val ≈ 300.0
    @test max_val ≈ 800.0

    # Test 3: 異常低温検出（100K < 150K）
    is_valid, min_val, max_val = Validators.check_temperature_range([100.0, 500.0, 800.0], 150.0, 3000.0)
    @test !is_valid
    @test min_val ≈ 100.0
    @test max_val ≈ 800.0

    # Test 4: 異常高温検出（3500K > 3000K）
    is_valid, min_val, max_val = Validators.check_temperature_range([300.0, 500.0, 3500.0], 150.0, 3000.0)
    @test !is_valid
    @test min_val ≈ 300.0
    @test max_val ≈ 3500.0

    # Test 5: 両端異常（低温 & 高温）
    is_valid, min_val, max_val = Validators.check_temperature_range([100.0, 500.0, 3500.0], 150.0, 3000.0)
    @test !is_valid
    @test min_val ≈ 100.0
    @test max_val ≈ 3500.0

    # Test 6: 3次元配列での範囲チェック
    T_3d = fill(500.0, (3, 3, 3))
    T_3d[1, 1, 1] = 200.0
    T_3d[3, 3, 3] = 1000.0
    is_valid, min_val, max_val = Validators.check_temperature_range(T_3d, 150.0, 3000.0)
    @test is_valid
    @test min_val ≈ 200.0
    @test max_val ≈ 1000.0

    # Test 7: カスタム範囲（200K ~ 1500K）
    is_valid, min_val, max_val = Validators.check_temperature_range([300.0, 500.0, 800.0], 200.0, 1500.0)
    @test is_valid

    is_valid, min_val, max_val = Validators.check_temperature_range([300.0, 500.0, 1800.0], 200.0, 1500.0)
    @test !is_valid
  end

  # ============================================================================
  # 3. check_flux_range: 熱流束範囲チェック
  # ============================================================================
  @testset "check_flux_range" begin
    # Test 1: 正常範囲（絶対値 ≤ 1e7 W/m²）
    is_valid, min_val, max_val = Validators.check_flux_range([1e5, 5e5, 1e6], 1e7)
    @test is_valid
    @test min_val ≈ 1e5
    @test max_val ≈ 1e6

    # Test 2: 正常範囲（負の値を含む）
    is_valid, min_val, max_val = Validators.check_flux_range([-1e5, 0.0, 1e5], 1e7)
    @test is_valid
    @test min_val ≈ -1e5
    @test max_val ≈ 1e5

    # Test 3: 異常な正の熱流束（1.5e7 > 1e7）
    is_valid, min_val, max_val = Validators.check_flux_range([1e5, 5e5, 1.5e7], 1e7)
    @test !is_valid
    @test min_val ≈ 1e5
    @test max_val ≈ 1.5e7

    # Test 4: 異常な負の熱流束（-1.5e7 < -1e7）
    is_valid, min_val, max_val = Validators.check_flux_range([-1.5e7, 0.0, 1e5], 1e7)
    @test !is_valid
    @test min_val ≈ -1.5e7
    @test max_val ≈ 1e5

    # Test 5: 2次元配列（表面熱流束）
    q_2d = fill(5e5, (10, 10))
    is_valid, min_val, max_val = Validators.check_flux_range(q_2d, 1e7)
    @test is_valid
    @test min_val ≈ 5e5
    @test max_val ≈ 5e5

    # Test 6: カスタム閾値（1e6 W/m²）
    is_valid, min_val, max_val = Validators.check_flux_range([5e5, 8e5, 1.2e6], 1e6)
    @test !is_valid
    @test max_val ≈ 1.2e6

    # Test 7: ゼロ熱流束
    is_valid, min_val, max_val = Validators.check_flux_range(zeros(5, 5), 1e7)
    @test is_valid
    @test min_val ≈ 0.0
    @test max_val ≈ 0.0
  end

  # ============================================================================
  # 4. check_gradient_magnitude: 勾配急変検出
  # ============================================================================
  @testset "check_gradient_magnitude" begin
    dx = 0.12e-3  # 0.12mm
    dy = 0.12e-3 * 0.866  # y方向（調整済み）
    dz = fill(1.0e-4, 5)  # 均等格子（0.1mm）

    # Test 1: 2次元フィールド（勾配なし）
    q_2d_flat = fill(5e5, (10, 10))
    is_valid, max_grad = Validators.check_gradient_magnitude(q_2d_flat, dx, dy, dz, 5000.0, 1e8)
    @test is_valid
    @test max_grad ≈ 0.0

    # Test 2: 2次元フィールド（正常な勾配）
    q_2d_normal = zeros(10, 10)
    for i in 1:10
      for j in 1:10
        q_2d_normal[i, j] = 1e5 + (i-1) * 1e4  # x方向に線形増加
      end
    end
    is_valid, max_grad = Validators.check_gradient_magnitude(q_2d_normal, dx, dy, dz, 5000.0, 1e8)
    # 勾配: 1e4 / dy ≈ 1e4 / (0.12e-3 * 0.866) ≈ 96.2e6 < 1e8
    @test is_valid
    @test max_grad < 1e8

    # Test 3: 2次元フィールド（急勾配）
    q_2d_steep = zeros(10, 10)
    q_2d_steep[1, :] .= 0.0
    q_2d_steep[2, :] .= 1e7  # 1セル間で1e7の変化
    is_valid, max_grad = Validators.check_gradient_magnitude(q_2d_steep, dx, dy, dz, 5000.0, 1e8)
    # 勾配: 1e7 / dy ≈ 1e7 / (0.12e-3 * 0.866) ≈ 9.62e10 > 1e8
    @test !is_valid
    @test max_grad > 1e8

    # Test 4: 3次元フィールド（温度場、勾配なし）
    T_3d_flat = fill(500.0, (5, 5, 5))
    is_valid, max_grad = Validators.check_gradient_magnitude(T_3d_flat, dx, dy, dz, 5000.0, 1e8)
    @test is_valid
    @test max_grad ≈ 0.0

    # Test 5: 3次元フィールド（z方向正常勾配）
    T_3d_z = zeros(5, 5, 5)
    for k in 1:5
      T_3d_z[:, :, k] .= 300.0 + (k-1) * 50.0  # z方向に50K/層
    end
    is_valid, max_grad = Validators.check_gradient_magnitude(T_3d_z, dx, dy, dz, 5000.0, 1e8)
    # 勾配: 50 / 1e-4 = 5e5 K/m > 5000 K/m
    @test !is_valid
    @test max_grad > 5000.0

    # Test 6: 3次元フィールド（緩やかな勾配）
    T_3d_gentle = zeros(5, 5, 5)
    for k in 1:5
      T_3d_gentle[:, :, k] .= 300.0 + (k-1) * 0.1  # z方向に0.1K/層
    end
    is_valid, max_grad = Validators.check_gradient_magnitude(T_3d_gentle, dx, dy, dz, 5000.0, 1e8)
    # 勾配: 0.1 / 1e-4 = 1000 K/m < 5000 K/m
    @test is_valid
    @test max_grad < 5000.0

    # Test 7: 非均等格子でのz方向勾配
    dz_nonuniform = [0.5e-4, 0.5e-4, 1.0e-4, 1.5e-4, 2.0e-4]  # 表面側に細かい格子
    T_3d_nonuniform = zeros(3, 3, 5)
    T_3d_nonuniform[:, :, 1] .= 300.0
    T_3d_nonuniform[:, :, 2] .= 310.0  # 10K変化、dz=0.5e-4m → 2e5 K/m
    T_3d_nonuniform[:, :, 3] .= 320.0  # 10K変化、dz=0.5e-4m → 2e5 K/m
    T_3d_nonuniform[:, :, 4] .= 330.0  # 10K変化、dz=1.0e-4m → 1e5 K/m
    T_3d_nonuniform[:, :, 5] .= 340.0  # 10K変化、dz=1.5e-4m → 6.67e4 K/m
    is_valid, max_grad = Validators.check_gradient_magnitude(T_3d_nonuniform, dx, dy, dz_nonuniform, 5000.0, 1e8)
    # 最大勾配: 10 / 0.5e-4 = 2e5 K/m > 5000 K/m
    @test !is_valid
    @test max_grad > 5000.0
  end

  # ============================================================================
  # 5. detect_numerical_anomalies: 包括的異常検出
  # ============================================================================
  @testset "detect_numerical_anomalies" begin
    dx = 0.12e-3
    dy = 0.12e-3 * 0.866
    dz = fill(1.0e-4, 5)

    # Test 1: 正常な温度場
    T_normal = fill(500.0, (5, 5, 5))
    has_anomaly, anomalies = Validators.detect_numerical_anomalies(
      T_normal, "温度場",
      iteration=1, timestep=10,
      temperature_range=(150.0, 3000.0),
      flux_range=(-1e7, 1e7),
      dx=dx, dy=dy, dz=dz
    )
    @test !has_anomaly
    @test length(anomalies) == 0

    # Test 2: NaN検出
    T_with_nan = fill(500.0, (5, 5, 5))
    T_with_nan[3, 3, 3] = NaN
    has_anomaly, anomalies = Validators.detect_numerical_anomalies(
      T_with_nan, "温度場",
      temperature_range=(150.0, 3000.0),
      dx=dx, dy=dy, dz=dz
    )
    @test has_anomaly
    @test any(occursin("NaN/Inf", a) for a in anomalies)

    # Test 3: 異常低温検出
    T_low = fill(500.0, (5, 5, 5))
    T_low[1, 1, 1] = 100.0  # 100K < 150K
    has_anomaly, anomalies = Validators.detect_numerical_anomalies(
      T_low, "温度場",
      temperature_range=(150.0, 3000.0),
      dx=dx, dy=dy, dz=dz
    )
    @test has_anomaly
    @test any(occursin("異常低温", a) for a in anomalies)

    # Test 4: 異常高温検出
    T_high = fill(500.0, (5, 5, 5))
    T_high[5, 5, 5] = 3500.0  # 3500K > 3000K
    has_anomaly, anomalies = Validators.detect_numerical_anomalies(
      T_high, "Temperature",
      temperature_range=(150.0, 3000.0),
      dx=dx, dy=dy, dz=dz
    )
    @test has_anomaly
    @test any(occursin("異常高温", a) for a in anomalies)

    # Test 5: 異常な熱流束
    q_high = fill(5e5, (10, 10))
    q_high[5, 5] = 1.5e7  # 1.5e7 > 1e7
    has_anomaly, anomalies = Validators.detect_numerical_anomalies(
      q_high, "熱流束",
      flux_range=(-1e7, 1e7),
      dx=dx, dy=dy
    )
    @test has_anomaly
    @test any(occursin("異常な熱流束", a) for a in anomalies)

    # Test 6: 複数の異常検出
    T_multi = fill(500.0, (5, 5, 5))
    T_multi[1, 1, 1] = NaN
    T_multi[2, 2, 2] = 100.0  # 異常低温
    T_multi[5, 5, 5] = 3500.0  # 異常高温
    has_anomaly, anomalies = Validators.detect_numerical_anomalies(
      T_multi, "温度場",
      temperature_range=(150.0, 3000.0),
      dx=dx, dy=dy, dz=dz
    )
    @test has_anomaly
    @test length(anomalies) >= 2  # NaN + 低温 or 高温
  end

  # ============================================================================
  # 6. check_temperature_field: 温度場ラッパー関数
  # ============================================================================
  @testset "check_temperature_field" begin
    dx = 0.12e-3
    dy = 0.12e-3 * 0.866
    dz = fill(1.0e-4, 5)

    # Test 1: 正常な温度場
    T_normal = fill(500.0, (5, 5, 5))
    is_valid, status_msg = Validators.check_temperature_field(
      T_normal, timestep=10, dx=dx, dy=dy, dz=dz
    )
    @test is_valid
    @test status_msg == "正常"

    # Test 2: 異常な温度場
    T_abnormal = fill(500.0, (5, 5, 5))
    T_abnormal[3, 3, 3] = 3500.0
    is_valid, status_msg = Validators.check_temperature_field(
      T_abnormal, timestep=10, dx=dx, dy=dy, dz=dz
    )
    @test !is_valid
    @test occursin("異常検出", status_msg)
  end

  # ============================================================================
  # 7. check_flux_field: 熱流束場ラッパー関数
  # ============================================================================
  @testset "check_flux_field" begin
    dx = 0.12e-3
    dy = 0.12e-3 * 0.866

    # Test 1: 正常な熱流束場
    q_normal = fill(5e5, (10, 10))
    is_valid, status_msg = Validators.check_flux_field(
      q_normal, iteration=5, timestep=10, dx=dx, dy=dy
    )
    @test is_valid
    @test status_msg == "正常"

    # Test 2: 異常な熱流束場
    q_abnormal = fill(5e5, (10, 10))
    q_abnormal[5, 5] = 1.5e7
    is_valid, status_msg = Validators.check_flux_field(
      q_abnormal, iteration=5, timestep=10, dx=dx, dy=dy
    )
    @test !is_valid
    @test occursin("異常検出", status_msg)
  end

  # ============================================================================
  # 8. check_adjoint_field: 随伴場ラッパー関数
  # ============================================================================
  @testset "check_adjoint_field" begin
    dx = 0.12e-3
    dy = 0.12e-3 * 0.866
    dz = fill(1.0e-4, 5)

    # Test 1: 正常な随伴場
    lambda_normal = fill(1e5, (5, 5, 5))
    is_valid, status_msg = Validators.check_adjoint_field(
      lambda_normal, iteration=5, timestep=10, dx=dx, dy=dy, dz=dz
    )
    @test is_valid
    @test status_msg == "正常"

    # Test 2: 異常な随伴場（オーバーフロー危険）
    lambda_abnormal = fill(1e5, (5, 5, 5))
    lambda_abnormal[3, 3, 3] = 1e20  # 極端に大きい値
    is_valid, status_msg = Validators.check_adjoint_field(
      lambda_abnormal, iteration=5, timestep=10, dx=dx, dy=dy, dz=dz
    )
    @test !is_valid
    @test occursin("異常検出", status_msg)
  end

end
