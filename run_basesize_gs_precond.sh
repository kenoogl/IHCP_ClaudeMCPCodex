#!/bin/bash

# basesize最適化測定スクリプト（GS前処理）
# 前回のレポート（phase5_2_step1）と同じ条件で測定
# 問題サイズ: 80×100×20 (N=160,000)
# 前処理: Gauss-Seidel (GS)
# データ: 実測データモード
# basesizeパラメータ: 400, 800, 1000, 2000
# スレッド数: 4

echo "=== basesize最適化測定（GS前処理、実測データ） ==="
echo "日時: $(date)"
echo ""

# 結果ディレクトリ作成
mkdir -p shared/results/basesize_gs_precond

# 問題サイズ: 80×100×20 (N=160,000)
echo "=== 問題サイズ: 80×100×20 (N=160,000) ==="
echo "前処理: Gauss-Seidel (GS)"
echo "データ: 実測データモード"
echo ""

for basesize in 400 800 1000 2000; do
  echo ""
  echo "--- basesize=$basesize ---"
  JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
    --solver pbicgstab --precond gs --nt 10 --basesize $basesize \
    2>&1 | tee shared/results/basesize_gs_precond/basesize_${basesize}.log

  echo "完了: basesize=$basesize"
  sleep 2
done

echo ""
echo "=== 全測定完了 ==="
echo "日時: $(date)"
echo ""
echo "結果ファイル一覧:"
ls -lh shared/results/basesize_gs_precond/
