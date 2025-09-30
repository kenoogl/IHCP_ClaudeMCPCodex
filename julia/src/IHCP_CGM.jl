"""
IHCP_CGM.jl

逆熱伝導問題（IHCP）を共役勾配法（CGM）で解くメインモジュール

Phase 1（基盤構築）:
- 熱物性値計算
- データ読み込み

今後のPhaseで追加予定:
- Phase 2: 直接解法（DHCP）
- Phase 3: 随伴解法（Adjoint）
- Phase 4: 共役勾配法（CGM）
- Phase 5: スライディングウィンドウ計算

対応Pythonコード:
/python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
"""

module IHCP_CGM

# Phase 1モジュールのインクルード
include("ThermalProperties.jl")
include("DataLoaders.jl")

# 再エクスポート
using .ThermalProperties
using .DataLoaders

export polyval_numba, thermal_properties_calculator
export load_sus304_thermal_properties, polyfit, fit_sus304_coefficients

# バージョン情報
const VERSION = v"0.1.0"
const PHASE = "Phase 1: 基盤構築"

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
  println("")
  println("今後の実装予定:")
  println("  ⏳ Phase 2: 直接解法（DHCP）")
  println("  ⏳ Phase 3: 随伴解法（Adjoint）")
  println("  ⏳ Phase 4: 共役勾配法（CGM）")
  println("  ⏳ Phase 5: スライディングウィンドウ計算")
  println("="^60)
end

end # module IHCP_CGM