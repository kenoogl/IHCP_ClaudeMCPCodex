#!/bin/bash
# PBICGSTAB + GS smoother テスト実行スクリプト

echo "=================================================="
echo "PBICGSTAB + GS Smoother ベンチマーク"
echo "=================================================="
echo "開始時刻: $(date)"
echo ""

cd "$(dirname "$0")"

# ブランチ確認
BRANCH=$(git branch --show-current)
if [ "$BRANCH" != "tuning7" ]; then
    echo "警告: 現在のブランチは $BRANCH です（tuning7 を期待）"
    exit 1
fi

# コミット確認
COMMIT=$(git rev-parse --short HEAD)
echo "コミット: $COMMIT"
echo ""

# Julia実行
echo "Juliaテスト開始..."
time julia --project=julia julia/scripts/run_10steps_fullsize_test.jl

# 結果確認
if [ -f shared/results/julia_10steps_fullsize.npz ]; then
    echo ""
    echo "✓ 結果ファイル生成成功: shared/results/julia_10steps_fullsize.npz"
    ls -lh shared/results/julia_10steps_fullsize.npz
else
    echo ""
    echo "✗ エラー: 結果ファイルが生成されませんでした"
    exit 1
fi

echo ""
echo "終了時刻: $(date)"
echo "=================================================="
