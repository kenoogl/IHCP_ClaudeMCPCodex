# 次セッションでの作業TODO

**作成日**: 2025年10月27日  
**ブランチ**: basesize-consistency  
**優先タスク**: test_dhcp_solver.jlの可変格子サイズ対応実装

---

## 今セッションで完了した作業 ✅

### 1. DHCP単体テストの再現とスケーラビリティ分析
- ✅ Phase 5.2 Step 1レポートの再現（basesize=1, 1000, 10000）
- ✅ nt依存性スケーラビリティ測定（nt=10, 50, 100）
- ✅ 反復回数履歴の詳細分析
- ✅ ホットスタート効果の定量化（反復回数28%削減）

### 2. ドキュメント作成
- ✅ `docs/scripts/test_dhcp_solver_guide.md` - スクリプト説明書
- ✅ `docs/plans/test_dhcp_solver_grid_size_enhancement.md` - 修正案

### 3. 重要な発見
- ⚠️ test_dhcp_solver.jlはq=0（ゼロ熱流束）で計算
- ✅ 測定データは初期条件と検証にのみ使用
- ✅ ソルバー性能テストとして有効、物理的妥当性テストではない

---

## 次セッションの作業（推定2時間）

### Phase 1: 実装（30分）
**対象**: `julia/scripts/test_dhcp_solver.jl`

```bash
# バックアップ作成
cp julia/scripts/test_dhcp_solver.jl julia/scripts/test_dhcp_solver.jl.bak
```

**実装内容**:
1. `parse_command_line_args()`拡張（--ni, --nj, --nk, --synthetic）
2. `prepare_synthetic_test()`追加
3. `prepare_measurement_test()`抽出
4. `analyze_residuals()`追加
5. `main()`修正

**参照**: `docs/plans/test_dhcp_solver_grid_size_enhancement.md`

### Phase 2-4: テスト（40分）

```bash
# Phase 2: 互換性テスト
julia --project=julia julia/scripts/test_dhcp_solver.jl

# Phase 3: 合成モードテスト
julia --project=julia julia/scripts/test_dhcp_solver.jl --synthetic

# Phase 4: 大規模テスト
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 160 --nj 200 --nk 40 --nt 10 --synthetic --basesize 600
```

### Phase 5: スケーラビリティ測定（30分）

```bash
# 問題サイズを変えて測定
for size in "40 50 10" "60 75 15" "80 100 20" "100 125 25" "120 150 30"
do
  ni=$(echo $size | awk '{print $1}')
  nj=$(echo $size | awk '{print $2}')
  nk=$(echo $size | awk '{print $3}')
  N=$((ni * nj * nk))
  
  JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
    --ni $ni --nj $nj --nk $nk --nt 10 --synthetic --basesize 600 \
    2>&1 | tee shared/results/scaling_N${N}.log
done
```

### Phase 6: 分析とレポート（20分）

**成果物**: `docs/reports/dhcp_problem_size_scaling_report.md`

---

## 重要な留意事項 ⚠️

### test_dhcp_solver.jlの制限
- q=0（ゼロ熱流束）で計算
- 物理的妥当性テストではない
- ソルバー性能テストとしてのみ有効

### 用途の使い分け
- 性能テスト: `test_dhcp_solver.jl` ✅
- 物理的検証: `run_10steps_fullsize_test.jl`
- 逆解析: `run_sliding_window.jl`

---

## 開始時チェックリスト

1. ディレクトリ確認: `pwd` → `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex`
2. ブランチ確認: `git status`
3. Julia確認: `julia --version` (v1.10+)
4. スレッド数: `echo $JULIA_NUM_THREADS` (推奨: 4)
5. 修正案確認: `cat docs/plans/test_dhcp_solver_grid_size_enhancement.md`
6. バックアップ作成: `cp julia/scripts/test_dhcp_solver.jl julia/scripts/test_dhcp_solver.jl.bak`

---

**作成**: 2025年10月27日
**所要時間**: 約2時間
