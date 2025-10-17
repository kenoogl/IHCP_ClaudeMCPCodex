#!/bin/bash
# テスト1の完了監視スクリプト

PID=80022
CHECK_INTERVAL=30  # 30秒ごとにチェック

echo "Monitoring PID $PID..."
echo "Press Ctrl+C to stop monitoring"
echo ""

while true; do
  if ps -p $PID > /dev/null 2>&1; then
    ELAPSED=$(ps -p $PID -o etime= | tr -d ' ')
    CPU=$(ps -p $PID -o %cpu= | tr -d ' ')
    MEM=$(ps -p $PID -o %mem= | tr -d ' ')
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running - Elapsed: $ELAPSED, CPU: ${CPU}%, MEM: ${MEM}%"
    
    # ログファイルサイズチェック
    if [ -f test1_pbicgstab.log ]; then
      LOGSIZE=$(wc -c < test1_pbicgstab.log)
      echo "  Log size: $LOGSIZE bytes"
      if [ $LOGSIZE -gt 0 ]; then
        echo "  Last 3 lines:"
        tail -3 test1_pbicgstab.log | sed 's/^/    /'
      fi
    fi
  else
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Process completed!"
    
    # 結果ファイル確認
    if [ -f shared/results/julia_10steps_fullsize.npz ]; then
      echo "  Result file: $(ls -lh shared/results/julia_10steps_fullsize.npz)"
    fi
    
    # ログ最終部分表示
    if [ -f test1_pbicgstab.log ]; then
      echo ""
      echo "=== Last 30 lines of log ==="
      tail -30 test1_pbicgstab.log
    fi
    
    break
  fi
  
  echo ""
  sleep $CHECK_INTERVAL
done
