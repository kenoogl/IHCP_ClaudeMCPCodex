#!/usr/bin/env julia
"""
FLoops + Task-local storage テスト

Task-local storageで設定したbasesizeが
@floopマクロで正しく動作するかを確認する。
"""

using FLoops

# Task-local storage方式
function set_config(; basesize::Int=1)
    task_local_storage(:backend_basesize, basesize)
    return nothing
end

function get_config()
    basesize = get(task_local_storage(), :backend_basesize, 1)
    return (basesize=basesize,)
end

function get_backend_tls(par::String)
    if par == "thread"
        config = get_config()
        return ThreadedEx(basesize=config.basesize)
    else
        return SequentialEx()
    end
end

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

# テスト用のmyfill!関数（Task-local storage版）
function myfill_tls!(a::Array{Float64,3}, val::Float64, par::String)
    if par == "sequential" || Threads.nthreads() == 1
        fill!(a, val)
    else
        backend = get_backend_tls(par)
        SZ = size(a)
        @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
            a[i,j,k] = val
        end
    end
    return nothing
end

# テスト用のmyfill!関数（Ref型グローバル変数版）
function myfill_ref!(a::Array{Float64,3}, val::Float64, par::String)
    if par == "sequential" || Threads.nthreads() == 1
        fill!(a, val)
    else
        backend = get_backend_ref(par)
        SZ = size(a)
        @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
            a[i,j,k] = val
        end
    end
    return nothing
end

# テスト実行
function main()
    println("="^80)
    println("FLoops + Backend設定 テスト")
    println("="^80)
    println("Threads: $(Threads.nthreads())")
    println()

    # テスト配列
    a_tls = zeros(Float64, 20, 20, 20)
    a_ref = zeros(Float64, 20, 20, 20)

    # Task-local storage方式のテスト
    println("[Test 1] Task-local storage方式")
    println("-"^80)
    try
        set_config(basesize=1000)
        config = get_config()
        println("設定確認: basesize=$(config.basesize)")

        println("myfill_tls!実行中...")
        myfill_tls!(a_tls, 42.0, "thread")

        println("✅ 成功: sum=$(sum(a_tls)), expected=$(42.0 * 8000)")
    catch e
        println("❌ エラー発生:")
        println(e)
        println()
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
    println()

    # Ref型グローバル変数方式のテスト
    println("[Test 2] Ref型グローバル変数方式")
   println("-"^80)
    try
        set_config_ref(basesize=1000)
        println("設定確認: basesize=$(GLOBAL_BASESIZE[])")

        println("myfill_ref!実行中...")
        myfill_ref!(a_ref, 42.0, "thread")

        println("✅ 成功: sum=$(sum(a_ref)), expected=$(42.0 * 8000)")
    catch e
        println("❌ エラー発生:")
        println(e)
        println()
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
    end
    println()

    println("="^80)
    println("テスト完了")
    println("="^80)
end

main()
