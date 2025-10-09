"""
test_rhs_core.jl

RHSCore.jl の単体テスト

calRHS_core! 関数のコア機能（初期化 + 6面一様境界条件）を検証
"""

using Test
using LinearAlgebra

# プロジェクトルートからの相対パスでモジュールをインクルード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

import .Commons
using .Commons: WorkBuffers

import .RHSCore
using .RHSCore: calRHS_core!

@testset "RHSCore.jl 単体テスト" begin

  @testset "基本動作: ゼロクリアと戻り値" begin
    # 小規模グリッド
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    # wk.bに初期値設定（非ゼロ）
    fill!(wk.b, 99.0)

    # パラメータ設定
    HF = zeros(6)  # すべての面で熱流束ゼロ
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]  # ガイドセル座標
    par = "sequential"

    # 実行
    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # 検証1: 戻り値の正確性
    @test SZ == (ni+2, nj+2, nk+2)
    @test dx1 ≈ 10.0
    @test dy1 ≈ 10.0
    @test z_st == 2
    @test z_ed == nk+1
    @test ddt ≈ 1000.0

    # 検証2: wk.bがゼロクリアされている
    @test all(wk.b .== 0.0)
  end

  @testset "X_minus面 境界条件（HF[1]）" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[1] = 100.0  # X_minus面のみ熱流束
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # X_minus面（i=2）のみに値が設定される
    expected_value = HF[1] / dx  # 100.0 / 0.1 = 1000.0

    # i=2, j=2:SZ[2]-1, k=z_st:z_ed の範囲
    for k in z_st:z_ed
      for j in 2:SZ[2]-1
        @test wk.b[2, j, k] ≈ expected_value
      end
    end

    # それ以外はゼロ
    @test all(wk.b[1, :, :] .== 0.0)
    @test all(wk.b[3:end, :, :] .== 0.0)
  end

  @testset "X_plus面 境界条件（HF[2]）" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[2] = 200.0  # X_plus面のみ熱流束
    dx = 0.2
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # X_plus面（i=SZ[1]-1=4）に負の値が設定される
    expected_value = -HF[2] / dx  # -200.0 / 0.2 = -1000.0

    i_plus = SZ[1] - 1
    for k in z_st:z_ed
      for j in 2:SZ[2]-1
        @test wk.b[i_plus, j, k] ≈ expected_value
      end
    end
  end

  @testset "Y_minus面 境界条件（HF[3]）" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[3] = 150.0  # Y_minus面のみ熱流束
    dx = 0.1
    dy = 0.15
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # Y_minus面（j=2）のみに値が設定される
    expected_value = HF[3] / dy  # 150.0 / 0.15 = 1000.0

    for k in z_st:z_ed
      for i in 2:SZ[1]-1
        @test wk.b[i, 2, k] ≈ expected_value
      end
    end
  end

  @testset "Y_plus面 境界条件（HF[4]）" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[4] = 180.0  # Y_plus面のみ熱流束
    dx = 0.1
    dy = 0.18
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # Y_plus面（j=SZ[2]-1=4）に負の値が設定される
    expected_value = -HF[4] / dy  # -180.0 / 0.18 = -1000.0

    j_plus = SZ[2] - 1
    for k in z_st:z_ed
      for i in 2:SZ[1]-1
        @test wk.b[i, j_plus, k] ≈ expected_value
      end
    end
  end

  @testset "Z_minus面 境界条件（HF[5]）" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[5] = 250.0  # Z_minus面のみ熱流束
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # Z_minus面（k=z_st=2）のみに値が設定される
    k_minus = z_st
    expected_value = HF[5] / ΔZ[k_minus]  # 250.0 / 0.05 = 5000.0

    for j in 2:SZ[2]-1
      for i in 2:SZ[1]-1
        @test wk.b[i, j, k_minus] ≈ expected_value
      end
    end

    # 他のk層はゼロ
    for k in (k_minus+1):(z_ed)
      @test all(wk.b[2:SZ[1]-1, 2:SZ[2]-1, k] .== 0.0)
    end
  end

  @testset "Z_plus面 境界条件（HF[6]）" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[6] = 300.0  # Z_plus面のみ熱流束
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # Z_plus面（k=z_ed=6）に負の値が設定される
    k_plus = z_ed
    expected_value = -HF[6] / ΔZ[k_plus]  # -300.0 / 0.05 = -6000.0

    for j in 2:SZ[2]-1
      for i in 2:SZ[1]-1
        @test wk.b[i, j, k_plus] ≈ expected_value
      end
    end
  end

  @testset "複数面同時設定" begin
    ni, nj, nk = 4, 4, 6
    wk = WorkBuffers(ni, nj, nk)

    # 6面すべてに異なる熱流束を設定
    HF = [100.0, 200.0, 150.0, 180.0, 250.0, 300.0]
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # 各面の値を検証
    # X_minus (i=2)
    @test wk.b[2, 3, 3] ≈ HF[1] / dx

    # X_plus (i=SZ[1]-1)
    @test wk.b[SZ[1]-1, 3, 3] ≈ -HF[2] / dx

    # Y_minus (j=2)
    @test wk.b[3, 2, 3] ≈ HF[3] / dy

    # Y_plus (j=SZ[2]-1)
    @test wk.b[3, SZ[2]-1, 3] ≈ -HF[4] / dy

    # Z_minus (k=z_st)
    @test wk.b[3, 3, z_st] ≈ HF[5] / ΔZ[z_st]

    # Z_plus (k=z_ed)
    @test wk.b[3, 3, z_ed] ≈ -HF[6] / ΔZ[z_ed]
  end

  @testset "ゼロ熱流束スキップ動作" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    # すべての面でゼロ熱流束
    HF = zeros(6)
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # wk.bがすべてゼロのまま
    @test all(wk.b .== 0.0)
  end

  @testset "非均等z格子" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = zeros(6)
    HF[5] = 100.0  # Z_minus面
    HF[6] = 200.0  # Z_plus面
    dx = 0.1
    dy = 0.1
    Δt = 0.001

    # 非均等格子
    ΔZ = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06]
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # Z_minus (k=2): ΔZ[2] = 0.02
    expected_z_minus = HF[5] / ΔZ[2]  # 100.0 / 0.02 = 5000.0
    @test wk.b[2, 2, 2] ≈ expected_z_minus

    # Z_plus (k=6): ΔZ[6] = 0.06
    expected_z_plus = -HF[6] / ΔZ[6]  # -200.0 / 0.06 = -3333.333...
    @test wk.b[2, 2, 6] ≈ expected_z_plus
  end

  @testset "z_range範囲外は処理されない" begin
    ni, nj, nk = 3, 3, 8
    wk = WorkBuffers(ni, nj, nk)

    HF = fill(100.0, 6)  # 全面に熱流束
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)

    # z方向の一部のみ計算対象
    z_range = [3, 6]  # k=3~6のみ
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # k=2（z_range外）はゼロ
    @test all(wk.b[:, :, 2] .== 0.0)

    # k=7（z_range外）はゼロ
    @test all(wk.b[:, :, 7] .== 0.0)

    # k=3~6（z_range内）は非ゼロ（X/Y面の影響）
    # X_minus (i=2, k=3~6)は非ゼロであるべき
    @test wk.b[2, 2, 3] != 0.0
  end

  @testset "ガイドセル境界は処理されない" begin
    ni, nj, nk = 3, 3, 5
    wk = WorkBuffers(ni, nj, nk)

    HF = fill(100.0, 6)  # 全面に熱流束
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]
    par = "sequential"

    SZ, dx1, dy1, z_st, z_ed, ddt = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par)

    # ガイドセル（i=1, i=SZ[1], j=1, j=SZ[2], k=1, k=SZ[3]）はゼロ
    @test all(wk.b[1, :, :] .== 0.0)  # i=1
    @test all(wk.b[SZ[1], :, :] .== 0.0)  # i=SZ[1]
    @test all(wk.b[:, 1, :] .== 0.0)  # j=1
    @test all(wk.b[:, SZ[2], :] .== 0.0)  # j=SZ[2]
    @test all(wk.b[:, :, 1] .== 0.0)  # k=1
    @test all(wk.b[:, :, SZ[3]] .== 0.0)  # k=SZ[3]
  end

  @testset "並列化backend=\"thread\"動作確認" begin
    ni, nj, nk = 10, 10, 10
    wk = WorkBuffers(ni, nj, nk)

    HF = fill(100.0, 6)
    dx = 0.1
    dy = 0.1
    Δt = 0.001
    ΔZ = fill(0.05, nk+1)
    z_range = [2, nk+1]

    # sequential版
    par_seq = "sequential"
    SZ1, dx1_1, dy1_1, z_st_1, z_ed_1, ddt_1 = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par_seq)
    b_seq = copy(wk.b)

    # thread版
    par_thread = "thread"
    SZ2, dx1_2, dy1_2, z_st_2, z_ed_2, ddt_2 = calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, z_range, par_thread)
    b_thread = copy(wk.b)

    # 結果が一致すること
    @test b_seq ≈ b_thread
    @test SZ1 == SZ2
    @test dx1_1 ≈ dx1_2
    @test dy1_1 ≈ dy1_2
    @test z_st_1 == z_st_2
    @test z_ed_1 == z_ed_2
    @test ddt_1 ≈ ddt_2
  end

end
