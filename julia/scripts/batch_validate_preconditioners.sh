#!/bin/bash
#
# batch_validate_preconditioners.sh
#
# BiCGSTABソルバーの全前処理器（none, diagonal, gs）に対して
# 並列化実装の数値精度検証をバッチ実行する
#
# 使用方法:
#   cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
#   bash julia/scripts/batch_validate_preconditioners.sh
#

set -e  # エラー時に即座に終了

# カラー出力設定
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}================================================================================${NC}"
echo -e "${BLUE}BiCGSTAB前処理器の包括的な並列化精度検証${NC}"
echo -e "${BLUE}================================================================================${NC}"
echo ""

# 前処理器リスト（jacobi は既に実施済み）
PRECONDITIONERS=("none" "diagonal" "gs")

# 検証条件
NT=10
WINDOW=5
OVERLAP=2
CGM_ITER=2

for precond in "${PRECONDITIONERS[@]}"; do
  echo -e "${YELLOW}================================================================================${NC}"
  echo -e "${YELLOW}前処理器: ${precond}${NC}"
  echo -e "${YELLOW}================================================================================${NC}"
  echo ""

  # 1スレッド実行
  echo -e "${GREEN}[1/4] 1スレッド実行中...${NC}"
  JULIA_NUM_THREADS=1 julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt ${NT} --window ${WINDOW} --overlap ${OVERLAP} --cgm-iter ${CGM_ITER} \
    --precond ${precond} \
    2>&1 | tee "shared/results/log_1thread_${precond}.txt"

  # 結果保存
  mv shared/results/julia_sliding_window_cgm2.npz \
     "shared/results/julia_sliding_window_cgm2_1thread_${precond}.npz"
  echo -e "${GREEN}✅ 1スレッド実行完了${NC}"
  echo ""

  # 8スレッド実行
  echo -e "${GREEN}[2/4] 8スレッド実行中...${NC}"
  JULIA_NUM_THREADS=8 julia --project=julia \
    julia/scripts/run_sliding_window.jl \
    --nt ${NT} --window ${WINDOW} --overlap ${OVERLAP} --cgm-iter ${CGM_ITER} \
    --precond ${precond} \
    2>&1 | tee "shared/results/log_8threads_${precond}.txt"

  # 結果保存
  mv shared/results/julia_sliding_window_cgm2.npz \
     "shared/results/julia_sliding_window_cgm2_8threads_${precond}.npz"
  echo -e "${GREEN}✅ 8スレッド実行完了${NC}"
  echo ""

  # 精度検証
  echo -e "${GREEN}[3/4] 精度検証実行中...${NC}"
  julia --project=julia julia/scripts/compare_parallel_results_with_args.jl \
    "shared/results/julia_sliding_window_cgm2_1thread_${precond}.npz" \
    "shared/results/julia_sliding_window_cgm2_8threads_${precond}.npz" \
    > "shared/results/validation_${precond}.txt"

  # 結果表示
  echo -e "${GREEN}[4/4] 検証結果:${NC}"
  cat "shared/results/validation_${precond}.txt" | tail -20
  echo ""

  echo -e "${GREEN}✅ ${precond} 前処理器の検証完了${NC}"
  echo ""
done

echo -e "${BLUE}================================================================================${NC}"
echo -e "${BLUE}全前処理器の検証完了${NC}"
echo -e "${BLUE}================================================================================${NC}"
echo ""
echo -e "${YELLOW}結果ファイル:${NC}"
ls -lh shared/results/julia_sliding_window_cgm2_*.npz
echo ""
echo -e "${YELLOW}検証レポート:${NC}"
ls -lh shared/results/validation_*.txt
