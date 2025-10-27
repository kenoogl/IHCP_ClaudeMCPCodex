#!/bin/bash
# nt値とbasesize依存性の体系的測定スクリプト
#
# 目的: nt=(10,30,50,100) × basesize=(400,500,600,700)の組み合わせを測定
#
# 使い方:
#   ./run_nt_basesize_measurements.sh [--resume]
#
# オプション:
#   --resume : 既存のログファイルをスキップして中断から再開
#
# 出力:
#   - 各測定のログ: shared/results/nt_basesize/nt{nt}_bs{basesize}.log
#   - 進捗記録: shared/results/nt_basesize/progress.txt
#   - 測定完了リスト: shared/results/nt_basesize/completed.txt
#
# セッション継続性:
#   - progress.txtで進捗を追跡
#   - completed.txtで完了済み測定を記録
#   - --resumeオプションで中断から再開可能

set -e

# プロジェクトルート
PROJECT_ROOT="/Users/Daily/Development/IHCP/TrialClaudeMCPCodex"
cd "$PROJECT_ROOT"

# 出力ディレクトリ
OUTPUT_DIR="shared/results/nt_basesize"
mkdir -p "$OUTPUT_DIR"

# 進捗ファイル
PROGRESS_FILE="$OUTPUT_DIR/progress.txt"
COMPLETED_FILE="$OUTPUT_DIR/completed.txt"

# 測定パラメータ
NT_VALUES=(30 50 100)
BASESIZE_VALUES=(400 500 600 700)
THREADS=4
SOLVER="pbicgstab"
PRECOND="gs"

# --resumeオプションの処理
RESUME_MODE=false
if [[ "$1" == "--resume" ]]; then
  RESUME_MODE=true
  echo "=== Resume mode: Skipping completed measurements ==="
fi

# 進捗記録の初期化
if [[ ! -f "$PROGRESS_FILE" ]]; then
  echo "0/12" > "$PROGRESS_FILE"
fi

if [[ ! -f "$COMPLETED_FILE" ]]; then
  touch "$COMPLETED_FILE"
fi

# 測定関数
run_measurement() {
  local nt=$1
  local basesize=$2
  local log_file="$OUTPUT_DIR/nt${nt}_bs${basesize}.log"
  local measurement_id="nt${nt}_bs${basesize}"

  # 既に完了している場合はスキップ
  if $RESUME_MODE && grep -q "^${measurement_id}$" "$COMPLETED_FILE"; then
    echo "  [SKIP] nt=$nt, basesize=$basesize (already completed)"
    return 0
  fi

  echo ""
  echo "=================================="
  echo "Measurement: nt=$nt, basesize=$basesize"
  echo "Output: $log_file"
  echo "=================================="

  # 測定実行
  JULIA_NUM_THREADS=$THREADS julia --project=julia julia/scripts/test_dhcp_solver.jl \
    --solver "$SOLVER" \
    --precond "$PRECOND" \
    --nt "$nt" \
    --basesize "$basesize" \
    2>&1 | tee "$log_file"

  # 完了を記録
  echo "$measurement_id" >> "$COMPLETED_FILE"

  # 測定成功確認
  if grep -q "Total runtime:" "$log_file"; then
    echo "  [OK] Measurement completed successfully"
  else
    echo "  [ERROR] Measurement may have failed (check log)"
    return 1
  fi
}

# メイン測定ループ
echo "========================================"
echo "nt-basesize dependency measurements"
echo "========================================"
echo "Configuration:"
echo "  NT values: ${NT_VALUES[*]}"
echo "  basesize values: ${BASESIZE_VALUES[*]}"
echo "  Threads: $THREADS"
echo "  Solver: $SOLVER"
echo "  Preconditioner: $PRECOND"
echo "  Total measurements: 12 (3 nt × 4 basesize)"
echo "  Resume mode: $RESUME_MODE"
echo "========================================"

START_TIME=$(date +%s)
TOTAL_MEASUREMENTS=12
COMPLETED_COUNT=0

# 既に完了した測定数をカウント
if [[ -f "$COMPLETED_FILE" ]]; then
  COMPLETED_COUNT=$(wc -l < "$COMPLETED_FILE" | tr -d ' ')
fi

echo "Progress: $COMPLETED_COUNT/$TOTAL_MEASUREMENTS completed"
echo ""

# 測定実行
for nt in "${NT_VALUES[@]}"; do
  for basesize in "${BASESIZE_VALUES[@]}"; do
    run_measurement "$nt" "$basesize"

    COMPLETED_COUNT=$((COMPLETED_COUNT + 1))
    echo "$COMPLETED_COUNT/$TOTAL_MEASUREMENTS" > "$PROGRESS_FILE"

    echo ""
    echo "Progress: $COMPLETED_COUNT/$TOTAL_MEASUREMENTS"

    # 推定残り時間
    ELAPSED=$(($(date +%s) - START_TIME))
    if [[ $COMPLETED_COUNT -gt 0 ]]; then
      AVG_TIME=$((ELAPSED / COMPLETED_COUNT))
      REMAINING=$((TOTAL_MEASUREMENTS - COMPLETED_COUNT))
      EST_REMAINING=$((AVG_TIME * REMAINING))
      echo "Estimated remaining time: $((EST_REMAINING / 60)) minutes"
    fi
    echo ""
  done
done

# 完了サマリー
END_TIME=$(date +%s)
TOTAL_ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "========================================"
echo "All measurements completed!"
echo "========================================"
echo "Total time: $((TOTAL_ELAPSED / 60)) minutes $((TOTAL_ELAPSED % 60)) seconds"
echo "Results directory: $OUTPUT_DIR"
echo "Log files: $OUTPUT_DIR/nt*_bs*.log"
echo ""
echo "Next steps:"
echo "  1. Extract data: julia scripts/extract_nt_basesize_data.jl"
echo "  2. Generate report: See docs/reports/nt_basesize_optimization_report.md"
echo "========================================"
