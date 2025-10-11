#!/usr/bin/env julia
"""
fill! vs myfill! 性能比較ベンチマーク

シングルスレッド環境での性能差を測定
"""

using FLoops
using Statistics

# Commons.get_backendの簡易版
function get_backend(par::String)
  if par == "thread"
    return ThreadedEx()
  else
    return SequentialEx()
  end
end

# myfill!の実装（CommonSolver.jlからコピー）
function myfill!(a::Array{Float64,3}, val::Float64, par::String)
  backend = get_backend(par)
  SZ = size(a)

  @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
    a[i,j,k] = val
  end
end

# 簡易ベンチマーク関数（複数回実行して平均を取る）
function benchmark_func(f, args...; n_runs=100)
  # ウォームアップ
  for _ in 1:10
    f(args...)
  end

  times = Float64[]
  for _ in 1:n_runs
    t = @elapsed f(args...)
    push!(times, t)
  end
  return minimum(times)  # 最小値を返す（キャッシュ効果などを排除）
end

# テストケース
function run_benchmark()
  println("="^80)
  println("fill! vs myfill! 性能比較ベンチマーク")
  println("="^80)
  println()

  # 小規模配列（10×10×10）
  println("[1/3] 小規模配列 (10×10×10)")
  a_small = zeros(Float64, 10, 10, 10)

  println("  fill!:")
  t_fill_small = benchmark_func(fill!, a_small, 0.0; n_runs=1000)
  println("    時間: $(t_fill_small * 1e9) ns")

  println("  myfill! (sequential):")
  t_myfill_small = benchmark_func(myfill!, a_small, 0.0, "sequential"; n_runs=1000)
  println("    時間: $(t_myfill_small * 1e9) ns")
  println("    オーバーヘッド: $(round((t_myfill_small / t_fill_small - 1) * 100, digits=1))%")
  println()

  # 中規模配列（80×100×20）
  println("[2/3] 中規模配列 (80×100×20)")
  a_medium = zeros(Float64, 80, 100, 20)

  println("  fill!:")
  t_fill_medium = benchmark_func(fill!, a_medium, 0.0; n_runs=500)
  println("    時間: $(t_fill_medium * 1e6) μs")

  println("  myfill! (sequential):")
  t_myfill_medium = benchmark_func(myfill!, a_medium, 0.0, "sequential"; n_runs=500)
  println("    時間: $(t_myfill_medium * 1e6) μs")
  println("    オーバーヘッド: $(round((t_myfill_medium / t_fill_medium - 1) * 100, digits=1))%")
  println()

  # 大規模配列（200×200×50）
  println("[3/3] 大規模配列 (200×200×50)")
  a_large = zeros(Float64, 200, 200, 50)

  println("  fill!:")
  t_fill_large = benchmark_func(fill!, a_large, 0.0; n_runs=100)
  println("    時間: $(t_fill_large * 1e3) ms")

  println("  myfill! (sequential):")
  t_myfill_large = benchmark_func(myfill!, a_large, 0.0, "sequential"; n_runs=100)
  println("    時間: $(t_myfill_large * 1e3) ms")
  println("    オーバーヘッド: $(round((t_myfill_large / t_fill_large - 1) * 100, digits=1))%")
  println()

  println("="^80)
  println("結論:")
  if t_myfill_medium > t_fill_medium * 1.1
    println("  ⚠️  myfill!はfill!より約$(round((t_myfill_medium / t_fill_medium - 1) * 100, digits=0))%遅い")
    println("  推奨: シングルスレッド環境ではfill!を使用すべき")
  else
    println("  ✅ myfill!とfill!の性能差は小さい（10%以内）")
  end
  println("="^80)

  return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
  try
    exit(run_benchmark())
  catch err
    println("\n❌ Error: ", err)
    Base.showerror(stdout, err, catch_backtrace())
    println()
    exit(1)
  end
end
