#!/usr/bin/env julia
"""
basesize効果測定テスト

異なるbasesizeの値での性能を比較し、
basesizeパラメータの効果を実証する。
"""

using FLoops
using BenchmarkTools

# Ref型グローバル変数方式
const GLOBAL_BASESIZE = Ref{Int}(1)

function set_config_ref(; basesize::Int=1)
    GLOBAL_BASESIZE[] = basesize
    return nothing
end

function get_backend_ref(par::String)
    if par == "thread"
        return ThreadedEx(basesize=GLOBAL_BASESIZE[])
    else
        return SequentialEx()
    end
end

# テスト用の計算関数
function compute_sum!(a::Array{Float64,3}, b::Array{Float64,3}, c::Array{Float64,3}, par::String)
    if par == "sequential" || Threads.nthreads() == 1
        @. a = b + c
    else
        backend = get_backend_ref(par)
        SZ = size(a)
        @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
            a[i,j,k] = b[i,j,k] + c[i,j,k]
        end
    end
    return nothing
end

# ベンチマーク実行
function run_benchmark(basesize::Int, size_tuple::Tuple{Int,Int,Int})
    set_config_ref(basesize=basesize)

    # テスト配列準備
    a = zeros(Float64, size_tuple...)
    b = rand(Float64, size_tuple...)
    c = rand(Float64, size_tuple...)

    # ウォームアップ
    compute_sum!(a, b, c, "thread")

    # ベンチマーク
    result = @benchmark compute_sum!($a, $b, $c, "thread") samples=10 seconds=5

    return result
end

# メイン処理
function main()
    println("="^80)
    println("basesize効果測定テスト")
    println("="^80)
    println("Threads: $(Threads.nthreads())")
    println()

    # テストサイズ
    test_sizes = [
        (20, 20, 20),
        (50, 50, 20),
        (80, 100, 20),
    ]

    # basesizeの値
    basesize_values = [1, 100, 1000, 10000]

    for test_size in test_sizes
        println("="^80)
        println("配列サイズ: $(test_size)")
        println("="^80)
        println()

        results = Dict()

        for bs in basesize_values
            println("basesize = $bs")
            println("-"^80)

            result = run_benchmark(bs, test_size)
            results[bs] = result

            med_time = median(result.times) / 1e6  # ns -> ms
            min_time = minimum(result.times) / 1e6
            println("  中央値: $(round(med_time, digits=3)) ms")
            println("  最小値: $(round(min_time, digits=3)) ms")
            println()
        end

        # 相対性能比較
        println("相対性能（basesize=1を基準）")
        println("-"^80)
        baseline = median(results[1].times)
        for bs in basesize_values
            speedup = baseline / median(results[bs].times)
            println("  basesize=$bs: $(round(speedup, digits=3))x")
        end
        println()
    end

    println("="^80)
    println("測定完了")
    println("="^80)
end

main()
