#!/bin/bash
# 全ソルバー・前処理組み合わせベンチマーク実行スクリプト
#
# test_dhcp_solver.jlを使って全ての組み合わせをテストします
#
# 実行方法:
#   bash julia/benchmarks/run_all_solver_tests.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTPUT_DIR="$PROJECT_ROOT/shared/results/solver_comparison"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

mkdir -p "$OUTPUT_DIR"

echo "========================================================================"
echo "全ソルバー・前処理組み合わせベンチマーク"
echo "========================================================================"
echo "実行日時: $(date)"
echo "出力先: $OUTPUT_DIR"
echo ""

# テスト条件
SOLVERS=("pbicgstab" "pcg")
PRECONDITIONERS=("none" "diagonal" "gs")
NT=10

# 結果ファイル
SUMMARY_FILE="$OUTPUT_DIR/benchmark_results_$TIMESTAMP.md"

# ヘッダー作成
cat > "$SUMMARY_FILE" <<EOF
# ソルバー・前処理組み合わせベンチマーク結果

**実行日時**: $(date)
**格子サイズ**: 80×100×20
**時間ステップ**: $NT
**テストスクリプト**: test_dhcp_solver.jl

---

## 結果サマリー

| ソルバー | 前処理 | DHCP時間 | 平均反復 | 最大反復 | 最終ステップ初期残差 | RMS残差 |
|---------|--------|---------|---------|---------|---------------------|---------|
EOF

# 各組み合わせでテスト実行
for solver in "${SOLVERS[@]}"; do
  for precond in "${PRECONDITIONERS[@]}"; do
    echo ""
    echo "========================================================================"
    echo "テスト: $solver + $precond"
    echo "========================================================================"

    LOG_FILE="$OUTPUT_DIR/${solver}_${precond}_$TIMESTAMP.log"

    # test_dhcp_solver.jl実行
    julia --project=julia "$PROJECT_ROOT/julia/scripts/test_dhcp_solver.jl" \
      --solver "$solver" \
      --precond "$precond" \
      --nt "$NT" \
      2>&1 | tee "$LOG_FILE"

    # 結果を抽出してサマリーに追加
    if grep -q "DHCP elapsed time" "$LOG_FILE"; then
      DHCP_TIME=$(grep "DHCP elapsed time" "$LOG_FILE" | awk '{print $4}')
      RMS_ERROR=$(grep "RMS residual" "$LOG_FILE" | awk '{print $3}')

      # 反復回数を抽出（全タイムステップから平均と最大を計算）
      ITERATIONS=$(grep "Iteration=" "$LOG_FILE" | awk '{print $4}' | tr -d ',')
      if [ -n "$ITERATIONS" ]; then
        AVG_ITER=$(echo "$ITERATIONS" | awk '{sum+=$1; n++} END {if(n>0) printf "%.1f", sum/n; else print "-"}')
        MAX_ITER=$(echo "$ITERATIONS" | sort -n | tail -1)
      else
        AVG_ITER="-"
        MAX_ITER="-"
      fi

      # 最終ステップ（t=10/10）の初期残差を抽出
      LAST_RES=$(grep "\[t=10/10\]" "$LOG_FILE" | grep "Res_0=" | awk '{print $7}' || echo "-")

      echo "| $solver | $precond | ${DHCP_TIME}s | $AVG_ITER | $MAX_ITER | $LAST_RES | $RMS_ERROR K |" >> "$SUMMARY_FILE"

      echo "✓ 完了: $solver + $precond (平均反復: $AVG_ITER, 最大反復: $MAX_ITER)"
    else
      echo "| $solver | $precond | エラー | - | - | - | - |" >> "$SUMMARY_FILE"
      echo "✗ エラー: $solver + $precond"
    fi
  done
done

echo ""
echo "========================================================================"
echo "全テスト完了"
echo "========================================================================"
echo "結果サマリー: $SUMMARY_FILE"
echo ""

# サマリー表示
cat "$SUMMARY_FILE"
