#!/bin/bash

# basesize最適化測定スクリプト
# 問題サイズ: 80×50×20, 80×100×40
# basesizeパラメータ: 400, 800, 1000, 2000
# スレッド数: 4

echo "=== basesize最適化測定開始 ==="
echo "日時: $(date)"
echo ""

# 結果ディレクトリ作成
mkdir -p shared/results/basesize_optimization

# 問題サイズ1: 80×50×20 (N=80,000)
echo "=== 問題サイズ1: 80×50×20 (N=80,000) ==="
for basesize in 400 800 1000 2000; do
  echo ""
  echo "--- basesize=$basesize ---"
  JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
    --ni 80 --nj 50 --nk 20 --nt 10 --synthetic --basesize $basesize \
    2>&1 | tee shared/results/basesize_optimization/size_80x50x20_basesize_${basesize}.log

  echo "完了: basesize=$basesize"
  sleep 2
done

# 問題サイズ2: 80×100×40 (N=320,000)
echo ""
echo "=== 問題サイズ2: 80×100×40 (N=320,000) ==="
for basesize in 400 800 1000 2000; do
  echo ""
  echo "--- basesize=$basesize ---"
  JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
    --ni 80 --nj 100 --nk 40 --nt 10 --synthetic --basesize $basesize \
    2>&1 | tee shared/results/basesize_optimization/size_80x100x40_basesize_${basesize}.log

  echo "完了: basesize=$basesize"
  sleep 2
done

echo ""
echo "=== 全測定完了 ==="
echo "日時: $(date)"
echo ""
echo "結果ファイル一覧:"
ls -lh shared/results/basesize_optimization/
