"""
コミット4c87fcf用の参照結果生成スクリプト

旧API（workパラメータなし）で実行し、結果を保存する。
"""

using LinearAlgebra
using Printf

# テスト用にIHCP_CGMモジュールをインクルード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

println("="^70)
println("コミット4c87fcf 中規模問題の参照結果生成")
println("="^70)

# 問題設定
ni, nj, nk = 5, 5, 10
dx, dy = 0.12e-3, 0.12e-3
dz = fill(0.05e-3, nk)
dz_b = fill(0.05e-3, nk)
dz_b[1] = Inf
dz_t = fill(0.05e-3, nk)
dz_t[end] = Inf
dt = 1e-3

rho = 7900.0
cp_coeffs = [462.0, 0.134, 0.0, 0.0]
k_coeffs = [14.6, 0.0127, 0.0, 0.0]

T_initial = fill(300.0, ni, nj, nk)

nt = 11
q_surface = fill(5000.0, ni, nj, nt-1)

println("\n問題設定:")
println("  格子: $(ni)×$(nj)×$(nk)")
println("  時間ステップ: $(nt), dt=$(dt)s")

# 旧API（workなし）で実行
println("\n実行中（旧API）...")

T_all = solve_dhcp!(
  T_initial, q_surface, nt, rho, cp_coeffs, k_coeffs,
  dx, dy, dz, dz_b, dz_t, dt;
  rtol=1e-8, maxiter=1000, verbose=true
)

println("\n実行完了")
println("  温度範囲: $(minimum(T_all)) - $(maximum(T_all)) K")

# 結果を保存
result_file = joinpath(@__DIR__, "medium_problem_4c87fcf.txt")

open(result_file, "w") do io
  println(io, "# コミット4c87fcf 中規模問題計算結果")
  println(io, "# 格子: $(ni)×$(nj)×$(nk), nt=$(nt)")
  println(io, "# 温度場 T_all[i,j,k,t] の全データ")
  for t in 1:nt
    for k_idx in 1:nk
      for j in 1:nj
        for i in 1:ni
          @printf(io, "%.16e\n", T_all[i,j,k_idx,t])
        end
      end
    end
  end
end

println("\n参照結果を保存しました: $(result_file)")
println("\n" * "="^70)
println("参照結果生成完了")
println("="^70)
