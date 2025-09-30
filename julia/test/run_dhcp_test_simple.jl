"""
簡易DHCP機能テスト（Phase 2）

JSONパース問題を回避して基本機能のみテスト
"""

using Test
using LinearAlgebra
using SparseArrays

# Phase 1モジュールのインクルード
include("../src/ThermalProperties.jl")
using .ThermalProperties

# Phase 2モジュールのインクルード
include("../src/solvers/DHCPSolver.jl")
using .DHCPSolver


println("="^60)
println("Phase 2簡易機能テスト")
println("="^60)

@testset "DHCP基本機能" begin

  @testset "1. 係数構築" begin
    println("\n[テスト1] 係数構築")

    # 小規模問題設定
    ni, nj, nk = 3, 3, 5
    N = ni * nj * nk

    dx, dy = 1e-3, 1e-3
    dz = fill(0.1e-3, nk)
    dz_b = fill(0.1e-3, nk); dz_b[1] = Inf
    dz_t = fill(0.1e-3, nk); dz_t[end] = Inf
    dt = 0.01

    rho = 8000.0
    cp = fill(500.0, ni, nj, nk)
    k = fill(15.0, ni, nj, nk)

    T_initial = fill(300.0, ni, nj, nk)
    q_surface = fill(1000.0, ni, nj)

    # 係数構築実行
    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_dhcp_system!(
      T_initial, q_surface, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    @test length(a_w) == N
    @test length(a_p) == N
    @test all(a_p .> 0)
    @test all(b .!= 0)

    println("  ✓ 係数配列が正しく構築された")
  end


  @testset "2. 疎行列組み立て" begin
    println("\n[テスト2] 疎行列組み立て")

    ni, nj, nk = 5, 5, 10
    N = ni * nj * nk

    dx, dy = 0.12e-3, 0.12e-3
    dz = fill(0.05e-3, nk)
    dz_b = fill(0.05e-3, nk); dz_b[1] = Inf
    dz_t = fill(0.05e-3, nk); dz_t[end] = Inf
    dt = 1e-3

    rho = 7900.0
    cp = fill(500.0, ni, nj, nk)
    k = fill(15.0, ni, nj, nk)

    T_initial = fill(300.0, ni, nj, nk)
    q_surface = fill(5000.0, ni, nj)

    a_w, a_e, a_s, a_n, a_b, a_t, a_p, b = build_dhcp_system!(
      T_initial, q_surface, rho, cp, k, dx, dy, dz, dz_b, dz_t, dt
    )

    A = assemble_dhcp_matrix(ni, nj, nk, a_w, a_e, a_s, a_n, a_b, a_t, a_p)

    @test size(A) == (N, N)
    @test issparse(A)
    @test all(diag(A) .> 0)

    println("  ✓ 疎行列が正しく組み立てられた")
    println("  　サイズ: $(size(A)), 非ゼロ要素: $(nnz(A))")
  end


  @testset "3. 時間積分ソルバー" begin
    println("\n[テスト3] 時間積分ソルバー")

    ni, nj, nk = 3, 3, 5
    dx, dy = 1e-3, 1e-3
    dz = fill(0.1e-3, nk)
    dz_b = fill(0.1e-3, nk); dz_b[1] = Inf
    dz_t = fill(0.1e-3, nk); dz_t[end] = Inf
    dt = 0.01

    rho = 8000.0
    cp_coeffs = [500.0, 0.0, 0.0, 0.0]
    k_coeffs = [15.0, 0.0, 0.0, 0.0]

    T_initial = fill(300.0, ni, nj, nk)

    nt = 6
    q_surface = fill(1000.0, nt-1, ni, nj)

    # 求解実行
    T_all = solve_dhcp!(
      T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol=1e-8, maxiter=1000, verbose=true
    )

    @test size(T_all) == (nt, ni, nj, nk)
    @test !any(isnan.(T_all))
    @test !any(isinf.(T_all))
    @test all(T_all .>= 0)

    # 温度上昇確認
    @test maximum(T_all[end, :, :, :]) > T_initial[1, 1, 1]

    println("  ✓ 時間積分が正常に完了")
    println("  　温度範囲: $(minimum(T_all)) - $(maximum(T_all)) K")
  end


  @testset "4. 温度依存熱物性値" begin
    println("\n[テスト4] 温度依存熱物性値")

    ni, nj, nk = 5, 5, 10
    dx, dy = 0.12e-3, 0.12e-3
    dz = fill(0.05e-3, nk)
    dz_b = fill(0.05e-3, nk); dz_b[1] = Inf
    dz_t = fill(0.05e-3, nk); dz_t[end] = Inf
    dt = 1e-3

    rho = 7900.0
    # SUS304近似（温度依存）
    cp_coeffs = [462.0, 0.134, 0.0, 0.0]
    k_coeffs = [14.6, 0.0127, 0.0, 0.0]

    T_initial = fill(300.0, ni, nj, nk)

    nt = 11
    q_surface = fill(5000.0, nt-1, ni, nj)

    T_all = solve_dhcp!(
      T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
      dx, dy, dz, dz_b, dz_t, dt,
      rtol=1e-8, maxiter=1000, verbose=false
    )

    @test size(T_all) == (nt, ni, nj, nk)
    @test !any(isnan.(T_all))

    println("  ✓ 温度依存熱物性値問題が正常に完了")
  end

end

println("\n" * "="^60)
println("Phase 2簡易機能テスト完了")
println("="^60)