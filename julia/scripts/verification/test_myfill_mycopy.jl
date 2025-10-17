#!/usr/bin/env julia

"""
myfill!とmycopy!の動作確認テスト
Commons.jlに移動した関数が正しく動作することを確認
"""

using IHCP_CGM

const Commons = IHCP_CGM.Commons

function main()
  println("=== myfill! と mycopy! の動作確認 ===\n")

  # テストデータ準備
  ni, nj, nk = 5, 5, 5
  a = zeros(Float64, ni, nj, nk)
  b = rand(Float64, ni, nj, nk)

  # テスト1: myfill!（sequential）
  println("テスト1: myfill! (sequential)")
  Commons.myfill!(a, 3.14, "sequential")
  all_filled = all(x -> x == 3.14, a)
  println("  全要素が3.14か: ", all_filled)
  println("  a[3,3,3] = ", a[3,3,3])
  @assert all_filled "myfill! failed (sequential)"
  println("  ✓ 成功\n")

  # テスト2: mycopy!（sequential）
  println("テスト2: mycopy! (sequential)")
  fill!(a, 0.0)
  Commons.mycopy!(a, b, "sequential")
  all_copied = all(a .== b)
  println("  全要素がコピーされたか: ", all_copied)
  println("  a[3,3,3] = ", a[3,3,3])
  println("  b[3,3,3] = ", b[3,3,3])
  @assert all_copied "mycopy! failed (sequential)"
  println("  ✓ 成功\n")

  # テスト3: Float32でmyfill!
  println("テスト3: myfill! (Float32)")
  a32 = zeros(Float32, ni, nj, nk)
  Commons.myfill!(a32, Float32(2.718), "sequential")
  all_filled32 = all(x -> x == Float32(2.718), a32)
  println("  全要素が2.718か: ", all_filled32)
  println("  a32[3,3,3] = ", a32[3,3,3])
  @assert all_filled32 "myfill! failed (Float32)"
  println("  ✓ 成功\n")

  # テスト4: Float32でmycopy!
  println("テスト4: mycopy! (Float32)")
  b32 = rand(Float32, ni, nj, nk)
  fill!(a32, Float32(0.0))
  Commons.mycopy!(a32, b32, "sequential")
  all_copied32 = all(a32 .== b32)
  println("  全要素がコピーされたか: ", all_copied32)
  println("  a32[3,3,3] = ", a32[3,3,3])
  println("  b32[3,3,3] = ", b32[3,3,3])
  @assert all_copied32 "mycopy! failed (Float32)"
  println("  ✓ 成功\n")

  # マルチスレッドテスト（スレッド数が1より大きい場合のみ）
  if Threads.nthreads() > 1
    println("テスト5: myfill! (thread)")
    fill!(a, 0.0)
    Commons.myfill!(a, 1.414, "thread")
    all_filled_thread = all(x -> x == 1.414, a)
    println("  全要素が1.414か: ", all_filled_thread)
    @assert all_filled_thread "myfill! failed (thread)"
    println("  ✓ 成功\n")

    println("テスト6: mycopy! (thread)")
    fill!(a, 0.0)
    Commons.mycopy!(a, b, "thread")
    all_copied_thread = all(a .== b)
    println("  全要素がコピーされたか: ", all_copied_thread)
    @assert all_copied_thread "mycopy! failed (thread)"
    println("  ✓ 成功\n")
  else
    println("スキップ: スレッドテスト（スレッド数=1）\n")
  end

  println("=== 全テスト完了 ✓ ===")
  println("myfill!とmycopy!は正常にCommons.jlに移動されました")
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
