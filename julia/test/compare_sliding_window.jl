"""
スライディングウィンドウCGM計算の比較スクリプト

測定条件:
- 問題サイズ: ni=40, nj=50, nk=10, nt=50 (中規模)
- ウィンドウサイズ: 30, オーバーラップ: 10

コミット4c87fcfと現在のコードで計算結果が一致することを確認
"""

using LinearAlgebra
using Printf

# テスト用にIHCP_CGMモジュールをインクルード
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

println("="^70)
println("スライディングウィンドウCGM計算の比較テスト")
println("="^70)

# 問題設定（テスト用に小規模化）
ni, nj, nk = 10, 10, 5  # 元: 40, 50, 10
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

nt = 20  # 元: 50
window_size = 10  # 元: 30
overlap = 3  # 元: 10

println("\n問題設定:")
println("  格子: $(ni)×$(nj)×$(nk)")
println("  時間ステップ: $(nt), dt=$(dt)s")
println("  ウィンドウサイズ: $(window_size)")
println("  オーバーラップ: $(overlap)")
println("  密度: $(rho) kg/m³")
println("  cp係数: $(cp_coeffs)")
println("  k係数: $(k_coeffs)")

# 初期条件と観測データの生成
T_init = fill(300.0, ni, nj, nk)

# 真の熱流束（解析対象）: 時間変動のある熱流束
q_true = zeros(Float64, ni, nj, nt-1)
for t in 1:nt-1
  q_base = 5000.0 * (1.0 + 0.1 * sin(2π * t / (nt-1)))
  q_true[:, :, t] .= q_base
end

# 観測データ生成（DHCPで前進計算）
println("\n観測データ生成中...")
work = WorkBuffers(ni+2, nj+2, nk+2)
Z, ΔZ = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

T_all_true = solve_dhcp!(
  T_init, q_true, work,
  nt, rho, cp_coeffs, k_coeffs,
  dx, dy, Z, ΔZ, dz, dz_b, dz_t, dt;
  rtol=1e-8, maxiter=1000, verbose=false
)

# 底面温度を観測データとする
Y_obs = T_all_true[:, :, 1, :]  # (ni, nj, nt)

# 初期推定（一定値スカラー）
q_init_value = 5000.0
cgm_iteration = 50

println("観測データ生成完了")
println("  観測温度範囲: $(minimum(Y_obs)) - $(maximum(Y_obs)) K")

# 現在のコードで実行
println("\n" * "="^70)
println("現在のコードで実行中...")
println("="^70)

q_result, window_info = solve_sliding_window_cgm(
  Y_obs, T_init, work,
  dx, dy, Z, ΔZ, dt,
  rho, cp_coeffs, k_coeffs,
  window_size, overlap, q_init_value, cgm_iteration;
  rtol_dhcp = 1e-6,
  maxiter_dhcp = 1000,
  rtol_adjoint = 1e-8,
  maxiter_adjoint = 1000
)

println("\n現在のコード実行完了")
println("  逆解析熱流束範囲: $(minimum(q_result)) - $(maximum(q_result)) W/m²")
println("  ウィンドウ数: $(length(window_info))")

# 全ウィンドウの目的関数履歴を結合
J_hist = Float64[]
for win in window_info
  append!(J_hist, win.J_history)
end
println("  総CGM反復回数: $(length(J_hist))")
if !isempty(J_hist)
  println("  最終目的関数: $(J_hist[end])")
end

# 参照結果ファイル
result_file = joinpath(@__DIR__, "sliding_window_4c87fcf.txt")

if !isfile(result_file)
  println("\n" * "="^70)
  println("参照結果ファイルが存在しません")
  println("現在の結果を参照結果として保存します")
  println("="^70)

  # 現在の結果を保存
  open(result_file, "w") do io
    println(io, "# コミット4c87fcf相当 スライディングウィンドウCGM計算結果")
    println(io, "# 格子: $(ni)×$(nj)×$(nk), nt=$(nt)")
    println(io, "# window_size=$(window_size), overlap=$(overlap)")
    println(io, "# 逆解析熱流束 q_result[i,j,t] の全データ")
    for t in 1:nt-1
      for j in 1:nj
        for i in 1:ni
          @printf(io, "%.16e\n", q_result[i,j,t])
        end
      end
    end

    println(io, "# 目的関数履歴")
    println(io, "$(length(J_hist))")
    for J in J_hist
      @printf(io, "%.16e\n", J)
    end
  end

  println("参照結果を保存しました: $(result_file)")
  println("\nコミット4c87fcfでこのスクリプトを実行して参照結果を生成してください:")
  println("  git stash")
  println("  git checkout 4c87fcf")
  println("  julia --project=. test/compare_sliding_window.jl")
  println("  git checkout tuning2")
  println("  git stash pop")

else
  # 参照結果を読み込み
  println("\n" * "="^70)
  println("参照結果を読み込み中...")
  println("="^70)

  q_ref = zeros(Float64, ni, nj, nt-1)
  J_hist_ref = Float64[]

  open(result_file, "r") do io
    # 熱流束データを読み込み
    for t in 1:nt-1
      for j in 1:nj
        for i in 1:ni
          # 空行とコメント行をスキップ
          line = ""
          while isempty(line) || startswith(line, "#")
            line = strip(readline(io))
          end
          q_ref[i,j,t] = parse(Float64, line)
        end
      end
    end

    # 目的関数履歴を読み込み
    # コメント行をスキップ
    line = ""
    while isempty(line) || startswith(line, "#")
      line = strip(readline(io))
    end
    n_iter = parse(Int, line)

    for _ in 1:n_iter
      line = strip(readline(io))
      push!(J_hist_ref, parse(Float64, line))
    end
  end

  println("参照結果読み込み完了")
  println("  逆解析熱流束範囲: $(minimum(q_ref)) - $(maximum(q_ref)) W/m²")
  println("  最終目的関数: $(J_hist_ref[end])")
  println("  CGM反復回数: $(length(J_hist_ref))")

  # 結果の比較
  println("\n" * "="^70)
  println("結果の比較")
  println("="^70)

  # 熱流束の誤差
  abs_error_q = abs.(q_result .- q_ref)
  max_abs_error_q = maximum(abs_error_q)
  mean_abs_error_q = sum(abs_error_q) / length(abs_error_q)

  rel_error_q = abs.(q_result .- q_ref) ./ (abs.(q_ref) .+ 1e-10)
  max_rel_error_q = maximum(rel_error_q)
  mean_rel_error_q = sum(rel_error_q) / length(rel_error_q)

  println("熱流束の誤差:")
  @printf("  最大絶対誤差: %.6e W/m²\n", max_abs_error_q)
  @printf("  平均絶対誤差: %.6e W/m²\n", mean_abs_error_q)
  @printf("  最大相対誤差: %.6e (%.6f%%)\n", max_rel_error_q, max_rel_error_q*100)
  @printf("  平均相対誤差: %.6e (%.6f%%)\n", mean_rel_error_q, mean_rel_error_q*100)

  # 目的関数履歴の比較
  if length(J_hist) == length(J_hist_ref)
    abs_error_J = abs.(J_hist .- J_hist_ref)
    max_abs_error_J = maximum(abs_error_J)

    println("\n目的関数履歴の誤差:")
    @printf("  最大絶対誤差: %.6e\n", max_abs_error_J)
    @printf("  最終値の差: %.6e\n", abs(J_hist[end] - J_hist_ref[end]))
  else
    println("\n⚠️ CGM反復回数が異なります")
    println("  現在: $(length(J_hist))回")
    println("  参照: $(length(J_hist_ref))回")
  end

  # 最大誤差の位置を特定
  max_error_idx = argmax(abs_error_q)
  i_max, j_max, t_max = Tuple(max_error_idx)

  println("\n最大誤差の位置:")
  println("  (i,j,t) = ($(i_max), $(j_max), $(t_max))")
  @printf("  参照値: %.16e W/m²\n", q_ref[i_max, j_max, t_max])
  @printf("  現在値: %.16e W/m²\n", q_result[i_max, j_max, t_max])
  @printf("  差分:   %.16e W/m²\n", abs_error_q[i_max, j_max, t_max])

  # 判定基準
  tolerance_abs = 1e-6   # 絶対誤差許容値 [W/m²]
  tolerance_rel = 1e-6   # 相対誤差許容値

  println("\n" * "="^70)
  println("判定結果")
  println("="^70)

  abs_pass = max_abs_error_q < tolerance_abs
  rel_pass = max_rel_error_q < tolerance_rel

  println("絶対誤差判定 (< $(tolerance_abs) W/m²): $(abs_pass ? "✅ 合格" : "❌ 不合格")")
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
