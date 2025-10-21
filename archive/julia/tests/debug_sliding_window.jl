# Phase 5スライディングウィンドウテストのデバッグスクリプト
# Python参照データとの詳細な数値比較

using JSON

# 参照データ読み込み関数
function load_reference_data(filename::String)
  filepath = joinpath(@__DIR__, "..", "data", filename)
  return JSON.parsefile(filepath)
end

println("=== 1D問題の詳細診断 ===\n")

# 1D参照データ読み込み
ref_data_1d = load_reference_data("phase5_reference_sliding_window_1D.json")
q_ref_1d = ref_data_1d["output"]["q_global"]

println("Python参照データ（1D）:")
println("形状: $(length(q_ref_1d)) × $(length(q_ref_1d[1])) × $(length(q_ref_1d[1][1]))")
println("\n全15ステップの熱流束値:")
for t in 1:15
  val = q_ref_1d[t][1][1]
  println("  t=$(t-1): q = $val")
end

println("\n真値 q_true:")
q_true_1d = ref_data_1d["input"]["q_true"]
for t in 1:15
  val = q_true_1d[t][1][1]
  println("  t=$(t-1): q_true = $val")
end

println("\n誤差分析:")
for t in 1:15
  ref = q_ref_1d[t][1][1]
  true_val = q_true_1d[t][1][1]
  error = abs(ref - true_val)
  rel_error = error / abs(true_val)
  println("  t=$(t-1): |ref - true| = $error (相対誤差: $rel_error)")
end

println("\n\n=== 2D問題の詳細診断 ===\n")

# 2D参照データ読み込み
ref_data_2d = load_reference_data("phase5_reference_sliding_window_2D.json")
q_ref_2d = ref_data_2d["output"]["q_global"]

println("Python参照データ（2D）:")
println("形状: $(length(q_ref_2d)) × $(length(q_ref_2d[1])) × $(length(q_ref_2d[1][1]))")
println("\n全13ステップの熱流束値（中心点 [1,1]）:")
for t in 1:13
  val = q_ref_2d[t][1][1]
  println("  t=$(t-1): q = $val")
end

println("\n真値 q_true（中心点 [1,1]）:")
q_true_2d = ref_data_2d["input"]["q_true"]
for t in 1:13
  val = q_true_2d[t][1][1]
  println("  t=$(t-1): q_true = $val")
end

println("\n誤差分析:")
for t in 1:13
  ref = q_ref_2d[t][1][1]
  true_val = q_true_2d[t][1][1]
  error = abs(ref - true_val)
  rel_error = error / abs(true_val)
  println("  t=$(t-1): |ref - true| = $error (相対誤差: $rel_error)")
end

println("\n\n=== 診断結果まとめ ===")
println("\n1D問題:")
println("  - Python参照データは極めて小さい値（1e-7オーダー）")
println("  - 真値は1000-3800 W/m²")
println("  - Python参照データは真値から大きくずれている")
println("  - これはノイズ付き観測データを用いた逆解析の結果")

println("\n2D問題:")
println("  - Python参照データは500.0付近（q_init_value=500.0）")
println("  - 真値は243.1 または 729.3 W/m²")
println("  - Python参照データは初期値から大きく変化していない")
println("  - 逆解析が収束していない可能性")

println("\n重要な発見:")
println("  - テストは「Python参照データとの一致」を要求")
println("  - しかし参照データ自体が真値から大きくずれている")
println("  - 相対誤差の計算方法に問題がある可能性")
println("    * 1D: rel_error = max_abs_diff / max(|q_ref|)")
println("    * q_refが1e-7オーダーなので、わずかな差でも相対誤差が巨大になる")
