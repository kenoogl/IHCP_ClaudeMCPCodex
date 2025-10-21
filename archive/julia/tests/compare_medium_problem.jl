"""
中規模問題の比較スクリプト

コミット4c87fcfと現在のコードで中規模問題を実行し、
計算結果が一致していることを確認する。

問題設定:
- 格子: 5×5×10
- 時間ステップ: 11 (dt=1ms)
- 温度依存熱物性値（SUS304）
- 一定表面熱流束: 5000 W/m²
"""

using LinearAlgebra
using Printf

# テスト用にIHCP_CGMモジュールをインクルード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

println("="^70)
println("中規模問題の比較テスト")
println("="^70)

# 問題設定（テスト5と同じ）
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
println("  初期温度: $(T_initial[1,1,1]) K")
println("  表面熱流束: $(q_surface[1,1,1]) W/m²")
println("  密度: $(rho) kg/m³")
println("  cp係数: $(cp_coeffs)")
println("  k係数: $(k_coeffs)")

# 現在のコードで実行
println("\n" * "="^70)
println("現在のコードで実行中...")
println("="^70)

work = WorkBuffers(ni+2, nj+2, nk+2)
Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

T_all_current = solve_dhcp!(
  T_initial, q_surface, work,
  nt, rho, cp_coeffs, k_coeffs,
  dx, dy, Z, ΔZ, dz, dz_b, dz_t, dt;
  rtol=1e-8, maxiter=1000, verbose=true
)

println("\n現在のコード実行完了")
println("  温度範囲: $(minimum(T_all_current)) - $(maximum(T_all_current)) K")

# コミット4c87fcfの結果を保存しておく（初回実行時に生成）
result_file = joinpath(@__DIR__, "medium_problem_4c87fcf.txt")

if !isfile(result_file)
  println("\n" * "="^70)
  println("参照結果ファイルが存在しません")
  println("現在の結果を参照結果として保存します")
  println("="^70)

  # 現在の結果を保存
  open(result_file, "w") do io
    println(io, "# コミット4c87fcf相当の中規模問題計算結果")
    println(io, "# 格子: $(ni)×$(nj)×$(nk), nt=$(nt)")
    println(io, "# 温度場 T_all[i,j,k,t] の全データ")
    for t in 1:nt
      for k_idx in 1:nk
        for j in 1:nj
          for i in 1:ni
            @printf(io, "%.16e\n", T_all_current[i,j,k_idx,t])
          end
        end
      end
    end
  end

  println("参照結果を保存しました: $(result_file)")
  println("\nコミット4c87fcfでこのスクリプトを実行して参照結果を生成してください:")
  println("  git stash")
  println("  git checkout 4c87fcf")
  println("  julia --project=. test/compare_medium_problem.jl")
  println("  git checkout tuning2")
  println("  git stash pop")

else
  # 参照結果を読み込み
  println("\n" * "="^70)
  println("参照結果を読み込み中...")
  println("="^70)

  T_all_ref = zeros(Float64, ni, nj, nk, nt)
  open(result_file, "r") do io
    # 全てのデータを順番に読み込み
    for t in 1:nt
      for k_idx in 1:nk
        for j in 1:nj
          for i in 1:ni
            # 空行とコメント行をスキップ
            line = ""
            while isempty(line) || startswith(line, "#")
              line = strip(readline(io))
            end
            T_all_ref[i,j,k_idx,t] = parse(Float64, line)
          end
        end
      end
    end
  end

  println("参照結果読み込み完了")
  println("  温度範囲: $(minimum(T_all_ref)) - $(maximum(T_all_ref)) K")

  # 結果の比較
  println("\n" * "="^70)
  println("結果の比較")
  println("="^70)

  # 絶対誤差
  abs_error = abs.(T_all_current .- T_all_ref)
  max_abs_error = maximum(abs_error)
  mean_abs_error = sum(abs_error) / length(abs_error)

  # 相対誤差
  rel_error = abs.(T_all_current .- T_all_ref) ./ (abs.(T_all_ref) .+ 1e-10)
  max_rel_error = maximum(rel_error)
  mean_rel_error = sum(rel_error) / length(rel_error)

  println("絶対誤差:")
  @printf("  最大: %.6e K\n", max_abs_error)
  @printf("  平均: %.6e K\n", mean_abs_error)

  println("\n相対誤差:")
  @printf("  最大: %.6e (%.6f%%)\n", max_rel_error, max_rel_error*100)
  @printf("  平均: %.6e (%.6f%%)\n", mean_rel_error, mean_rel_error*100)

  # 最大誤差の位置を特定
  max_error_idx = argmax(abs_error)
  i_max, j_max, k_max, t_max = Tuple(max_error_idx)

  println("\n最大誤差の位置:")
  println("  (i,j,k,t) = ($(i_max), $(j_max), $(k_max), $(t_max))")
  @printf("  参照値: %.16e K\n", T_all_ref[i_max, j_max, k_max, t_max])
  @printf("  現在値: %.16e K\n", T_all_current[i_max, j_max, k_max, t_max])
  @printf("  差分:   %.16e K\n", abs_error[i_max, j_max, k_max, t_max])

  # 判定基準
  tolerance_abs = 1e-10  # 絶対誤差許容値 [K]
  tolerance_rel = 1e-8   # 相対誤差許容値

  println("\n" * "="^70)
  println("判定結果")
  println("="^70)

  abs_pass = max_abs_error < tolerance_abs
  rel_pass = max_rel_error < tolerance_rel

  println("絶対誤差判定 (< $(tolerance_abs) K): $(abs_pass ? "✅ 合格" : "❌ 不合格")")
  println("相対誤差判定 (< $(tolerance_rel)): $(rel_pass ? "✅ 合格" : "❌ 不合格")")

  if abs_pass || rel_pass
    println("\n✅ 計算結果が一致しています")
    exit(0)
  else
    println("\n❌ 計算結果に差異があります")
    exit(1)
  end
end

println("\n" * "="^70)
println("比較テスト完了")
println("="^70)
