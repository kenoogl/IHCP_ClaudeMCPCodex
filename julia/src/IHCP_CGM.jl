"""
IHCP_CGM.jl

逆熱伝導問題（IHCP）を共役勾配法（CGM）で解くメインモジュール

Phase 1（基盤構築）:
- 熱物性値計算
- データ読み込み

Phase 2（直接解法）:
- DHCP（直接熱伝導問題）ソルバー

Phase 3（随伴解法）:
- Adjoint（随伴方程式）ソルバー

Phase 4（共役勾配法）:
- CGM（共役勾配法）最適化

Phase 5（スライディングウィンドウ計算）:
- 長時間計算の分割処理

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
"""

module IHCP_CGM

# Phase 1モジュールのインクルード
include("ThermalProperties.jl")
include("DataLoaders.jl")

# Phase 2モジュールのインクルード
include("solvers/DHCPSolver.jl")

# Phase 3モジュールのインクルード
include("solvers/AdjointSolver.jl")

# Phase 4モジュールのインクルード
include("solvers/StoppingCriteria.jl")
include("solvers/CGMSolver.jl")

# Phase 5モジュールのインクルード
include("solvers/SlidingWindowSolver.jl")

# Phase 6モジュールのインクルード
include("utils/validators.jl")

# 再エクスポート
using .ThermalProperties
using .DataLoaders
using .DHCPSolver
using .AdjointSolver
using .StoppingCriteria
using .CGMSolver
using .SlidingWindowSolver
using .Validators

export polyval_numba, thermal_properties_calculator
export load_sus304_thermal_properties, polyfit, fit_sus304_coefficients
export build_dhcp_system!, assemble_dhcp_matrix, solve_dhcp!, dhcp_index
export build_adjoint_system!, assemble_adjoint_matrix, solve_adjoint!, adjoint_index
export check_discrepancy, check_plateau, check_stopping_criteria, StoppingStatus
export solve_cgm!, compute_gradient!, compute_sensitivity!, compute_step_size
export solve_sliding_window_cgm, WindowInfo
export check_field_finite, check_temperature_range, check_flux_range
export check_gradient_magnitude, detect_numerical_anomalies
export check_temperature_field, check_flux_field, check_adjoint_field

# バージョン情報
const VERSION = v"0.5.0"
const PHASE = "Phase 5: スライディングウィンドウ計算"

"""
  version_info()

モジュールのバージョン情報を表示
"""
function version_info()
  println("="^60)
  println("IHCP_CGM.jl - 逆熱伝導問題 共役勾配法ソルバー")
  println("="^60)
  println("Version: $(VERSION)")
  println("Current Phase: $(PHASE)")
  println("")
  println("実装済み機能:")
  println("  ✓ 熱物性値計算 (ThermalProperties.jl)")
  println("  ✓ データ読み込み (DataLoaders.jl)")
  println("  ✓ 直接解法 DHCP (DHCPSolver.jl)")
  println("  ✓ 随伴解法 Adjoint (AdjointSolver.jl)")
  println("  ✓ 停止判定 (StoppingCriteria.jl)")
  println("  ✓ 共役勾配法 CGM (CGMSolver.jl)")
  println("  ✓ スライディングウィンドウ計算 (SlidingWindowSolver.jl)")
  println("")
  println("全Phase完了!")
  println("="^60)
end

end # module IHCP_CGM