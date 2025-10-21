#!/bin/bash
# DHCPソルバー・前処理組み合わせベンチマーク実行スクリプト
#
# 実行方法:
#   bash julia/scripts/run_solver_comparison.sh
#
# 出力先: shared/results/solver_comparison/

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUTPUT_DIR="$PROJECT_ROOT/shared/results/solver_comparison"

# 出力ディレクトリ作成
mkdir -p "$OUTPUT_DIR"

# タイムスタンプ生成
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "========================================================================"
echo "DHCPソルバー・前処理組み合わせベンチマーク"
echo "========================================================================"
echo "実行日時: $(date)"
echo "出力先: $OUTPUT_DIR"
echo ""

# ベンチマーク結果を保存する配列
declare -a results

# テストパターン定義
# ソルバー: pbicgstab, pcg
# 前処理: none, diagonal, gs
solvers=("pbicgstab" "pcg")
preconditioners=("none" "diagonal" "gs")

# 全組み合わせを実行
for solver in "${solvers[@]}"; do
  for precond in "${preconditioners[@]}"; do
    echo "========================================================================"
    echo "テスト: $solver + $precond"
    echo "========================================================================"

    log_file="$OUTPUT_DIR/${solver}_${precond}_${TIMESTAMP}.log"

    # test_dhcp_solver.jlを実行（10ステップ）
    if julia --project=julia julia/scripts/test_dhcp_solver.jl \
        --solver "$solver" \
        --precond "$precond" \
        --nt 10 > "$log_file" 2>&1; then

      echo "✓ 完了: $solver + $precond"

      # 結果を抽出
      elapsed=$(grep "DHCP elapsed time:" "$log_file" | awk '{print $4}')
      iterations=$(grep "converged:" "$log_file" | tail -1 | awk -F'Iteration=' '{print $2}' | awk '{print $1}')
      final_res0=$(grep "converged:" "$log_file" | tail -1 | awk -F'Res_0=' '{print $2}' | awk '{print $1}')

      results+=("| $solver | $precond | $elapsed s | $iterations | $final_res0 |")
    else
      echo "❌ エラー: $solver + $precond"
      results+=("| $solver | $precond | ERROR | - |")
    fi

    echo ""
  done
done

# サマリー生成
summary_file="$OUTPUT_DIR/summary.md"

cat > "$summary_file" <<EOF
# DHCPソルバー・前処理組み合わせベンチマーク結果

**実行日時**: $(date)
**テスト条件**:
- 格子: 80×100×20 (N=160,000)
- 時間ステップ: 10
- 収束条件: rtol=1.0e-6, maxiter=20000
- 熱流束: ゼロ（簡易テスト）

## 結果一覧

| ソルバー | 前処理 | 実行時間 | 最終反復数 | 最終ステップ初期残差 |
|----------|--------|----------|------------|---------------------|
EOF

# 結果を追加
for result in "${results[@]}"; do
  echo "$result" >> "$summary_file"
done

cat >> "$summary_file" <<EOF

## 詳細ログ

各テストの詳細ログは以下のファイルに保存されています:
EOF

for solver in "${solvers[@]}"; do
  for precond in "${preconditioners[@]}"; do
    echo "- \`${solver}_${precond}_${TIMESTAMP}.log\`" >> "$summary_file"
  done
done

cat >> "$summary_file" <<EOF

## 分析

### 実行時間比較

(ベンチマーク完了後、手動で追加)

### 反復回数比較

(ベンチマーク完了後、手動で追加)

### 推奨設定

(ベンチマーク完了後、手動で追加)
EOF

echo "========================================================================"
echo "ベンチマーク完了"
echo "========================================================================"
echo "サマリー: $summary_file"
echo ""
cat "$summary_file"
