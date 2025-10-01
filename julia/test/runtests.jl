"""
runtests.jl

IHCP_CGMプロジェクトのメインテストランナー

実行方法:
  julia test/runtests.jl

またはJuliaパッケージマネージャーから:
  ]test
"""

using Test

println("\n" * "="^70)
println("IHCP_CGM.jl テストスイート")
println("="^70)

# Phase 1: 熱物性値計算のテスト
println("\n[Phase 1] 熱物性値計算モジュールのテスト")
println("-"^70)
include("test_thermal_properties.jl")

# Phase 2: 直接解法（DHCP）のテスト
println("\n[Phase 2] 直接解法（DHCP）のテスト")
println("-"^70)
include("test_dhcp_solver.jl")

# Phase 3: 随伴解法（Adjoint）のテスト
println("\n[Phase 3] 随伴解法（Adjoint）のテスト")
println("-"^70)
include("test_adjoint_solver.jl")

# Phase 4: 共役勾配法（CGM）のテスト
println("\n[Phase 4] 共役勾配法（CGM）のテスト")
println("-"^70)
include("test_cgm_solver.jl")

# Phase 5: スライディングウィンドウ計算のテスト
println("\n[Phase 5] スライディングウィンドウ計算のテスト")
println("-"^70)
include("test_sliding_window.jl")

# Phase 6: 基本機能完成のテスト
println("\n[Phase 6 A-1] 検証器関数群のテスト")
println("-"^70)
include("test_validators.jl")

println("\n[Phase 6 C-1] 実データ読込機能のテスト")
println("-"^70)
include("test_data_loaders.jl")

println("\n" * "="^70)
println("全テスト完了")
println("="^70)