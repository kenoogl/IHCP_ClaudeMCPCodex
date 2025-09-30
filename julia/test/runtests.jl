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

# 今後のPhase用プレースホルダー
# println("\n[Phase 2] 直接解法（DHCP）のテスト")
# include("test_dhcp.jl")

# println("\n[Phase 3] 随伴解法（Adjoint）のテスト")
# include("test_adjoint.jl")

# println("\n[Phase 4] 共役勾配法（CGM）のテスト")
# include("test_cgm.jl")

# println("\n[Phase 5] スライディングウィンドウ計算のテスト")
# include("test_sliding_window.jl")

println("\n" * "="^70)
println("全テスト完了")
println("="^70)