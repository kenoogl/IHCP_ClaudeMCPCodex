# Phase 5.2 実環境検証計画書

**作成日**: 2025年10月23日
**ブランチ**: parallelization
**前提**: Phase 5.1完了（basesize最適化実装済み）
**目的**: test_basesize_performance.jlで得られた効果が実際のDHCP/CGMソルバーでも再現されるか段階的に検証

---

## 📋 エグゼクティブサマリー

### 検証の背景

Phase 5.1で`test_basesize_performance.jl`により以下の結果を得た：

- **配列サイズ80×100×20（160,000要素）での測定**
- basesize=1: 64.080 ms
- basesize=10000: 0.037 ms
- **高速化率: 1720.3倍**

### 検証の必要性

単純な配列操作（`a = b + c`）での高速化が、実際の以下の計算でも再現されるか確認が必要：

1. DHCPソルバー（PBICGSTAB! + GS前処理）
2. Adjoint/Sensitivityソルバー
3. CGMアルゴリズム全体
4. スライディングウィンドウ計算

### ソルバー設定（全テスト統一）

**全テストで以下に統一**:
- **Linear solver**: PBICGSTAB!
- **Preconditioner**: GS (Gauss-Seidel)
- **収束判定**: rtol=1e-6, atol=1e-10
- **理由**: 最も安定した収束性、Phase 5で実績あり

---

## 🎯 検証目標

### 主要目標

1. **性能向上の実証**: basesize=10000で実用的な高速化を達成
2. **数値精度の保証**: basesizeの変更が結果に影響しないことを確認
3. **最適値の決定**: 実環境での最適basesize値を特定

### 成功基準

| 項目 | 最小要件 | 理想目標 |
|------|----------|----------|
| Step 1 (DHCP単体) | 10倍高速化 | 100倍高速化 |
| Step 2 (10steps full) | 3倍高速化 | 10倍高速化 |
| Step 3 (Sliding window) | 2倍高速化 | 5倍高速化 |
| 数値精度 | 相対誤差 < 1% | 相対誤差 < 0.1% |

---

## 📝 Step 1: test_dhcp_solver.jl の更新とテスト

### 1.1 タスク概要

- **目的**: DHCP単体でのbasesize効果を測定
- **対象**: `julia/test/test_dhcp_solver.jl`
- **更新内容**: backend設定機能追加、性能測定コード追加

### 1.2 必要な更新

#### 更新1: ソルバー設定の統一

```julia
# 全テストケースで統一
solver = "pbicgstab"
precond = "gs"
rtol = 1e-6
atol = 1e-10
```

#### 更新2: backend設定機能の追加

```julia
using IHCP_CGM.Commons: set_backend_config

# basesize可変テスト
@testset "DHCP basesize効果測定" begin
    for basesize in [1, 1000, 10000]
        @testset "basesize=$basesize" begin
            set_backend_config(basesize=basesize)

            # DHCP計算実行
            result = @timed solve_dhcp(...)

            @info "basesize=$basesize" time=result.time
        end
    end
end
```

### 1.3 実行コマンド

```bash
# デフォルト実行（既存テスト）
JULIA_NUM_THREADS=4 julia --project=julia julia/test/test_dhcp_solver.jl

# 性能測定版
JULIA_NUM_THREADS=4 julia --project=julia julia/test/test_dhcp_solver.jl \
  2>&1 | tee shared/results/step1_dhcp_basesize_test.log
```

### 1.4 期待される結果

✅ **動作確認**:
- 全テストケース合格
- ソルバー収束（PBICGSTAB! + GS）
- 数値精度維持

✅ **性能測定**:
- basesize=1: ベースライン（遅い）
- basesize=1000: 中程度の高速化
- basesize=10000: 最大高速化

### 1.5 測定項目

| 項目 | 説明 |
|------|------|
| DHCP solve time | 1ステップあたりの解法時間 |
| CG iterations | 反復回数（変化なし想定） |
| Memory usage | メモリ使用量 |
| Speedup ratio | basesize=1との比較 |

### 1.6 推定所要時間

- コード更新: 20分
- テスト実行: 10分
- **合計: 30分**

---

## 📝 Step 2: run_10steps_fullsize_test.jl での検証

### 2.1 タスク概要

- **目的**: 10ステップフル計算でのbasesize効果測定
- **対象**: DHCP + Adjoint + Sensitivity + CGM の統合動作
- **計算規模**: 80×100×20格子、10時間ステップ

### 2.2 前提条件確認

```bash
# ファイル存在確認
ls -lh julia/scripts/run_10steps_fullsize_test.jl

# --basesizeオプション実装確認
grep -n "basesize" julia/scripts/run_10steps_fullsize_test.jl
```

### 2.3 実行計画

#### Test 2-1: basesize=1（ベースライン）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_10steps_fullsize_test.jl \
  --solver pbicgstab --precond gs \
  --basesize 1 \
  2>&1 | tee shared/results/step2_fullsize_bs1.log
```

#### Test 2-2: basesize=1000（中間値）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_10steps_fullsize_test.jl \
  --solver pbicgstab --precond gs \
  --basesize 1000 \
  2>&1 | tee shared/results/step2_fullsize_bs1000.log
```

#### Test 2-3: basesize=10000（最適値候補）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_10steps_fullsize_test.jl \
  --solver pbicgstab --precond gs \
  --basesize 10000 \
  2>&1 | tee shared/results/step2_fullsize_bs10000.log
```

### 2.4 期待される結果

✅ **動作確認**:
- 10ステップ完全実行
- DHCP/Adjoint/Sensitivity/CGM 全て正常動作
- 温度場・熱流束結果の保存

✅ **性能測定**:
- Total runtime短縮
- 各ソルバーの時間短縮

### 2.5 測定項目

| 項目 | basesize=1 | basesize=1000 | basesize=10000 |
|------|-----------|---------------|----------------|
| Total runtime | ? s | ? s | ? s |
| DHCP avg time | ? ms | ? ms | ? ms |
| Adjoint avg time | ? ms | ? ms | ? ms |
| Sensitivity avg time | ? ms | ? ms | ? ms |
| CGM time | ? s | ? s | ? s |
| Speedup | 1.0x | ?x | ?x |

### 2.6 推定所要時間

- Test 2-1: 10分
- Test 2-2: 5分
- Test 2-3: 3分
- 結果確認: 5分
- **合計: 23分**

---

## 📝 Step 3: run_sliding_window.jl での検証

### 3.1 タスク概要

- **目的**: スライディングウィンドウでの実用性能測定
- **対象**: 実際の運用シナリオに近い計算
- **パターン**: 小ウィンドウ×2 + 大ウィンドウ×2 = 計4回

### 3.2 実行計画

#### Test 3-1: 小ウィンドウ + basesize=1

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 1 \
  --solver pbicgstab --precond gs \
  --basesize 1 \
  2>&1 | tee shared/results/step3_small_bs1.log
```

**特徴**:
- ウィンドウサイズ: 5ステップ
- オーバーラップ: 2ステップ
- CGM反復: 1回（高速）

#### Test 3-2: 小ウィンドウ + basesize=10000

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 1 \
  --solver pbicgstab --precond gs \
  --basesize 10000 \
  2>&1 | tee shared/results/step3_small_bs10000.log
```

#### Test 3-3: 大ウィンドウ + basesize=1

```bash
JULIA_NUM_THREADS=8 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 10 --window 71 --overlap 17 --cgm-iter 3 \
  --solver pbicgstab --precond gs \
  --basesize 1 \
  2>&1 | tee shared/results/step3_large_bs1.log
```

**特徴**:
- ウィンドウサイズ: 71ステップ（大規模）
- オーバーラップ: 17ステップ
- CGM反復: 3回（実用設定）
- スレッド数: 8（並列性向上）

#### Test 3-4: 大ウィンドウ + basesize=10000

```bash
JULIA_NUM_THREADS=8 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 10 --window 71 --overlap 17 --cgm-iter 3 \
  --solver pbicgstab --precond gs \
  --basesize 10000 \
  2>&1 | tee shared/results/step3_large_bs10000.log
```

### 3.3 期待される結果

✅ **動作確認**:
- 全ウィンドウ処理完了
- CGM収束
- 結果ファイル生成

✅ **性能測定**:
- Window処理時間の短縮
- Total runtimeの短縮

### 3.4 測定項目

#### 小ウィンドウ測定

| 項目 | basesize=1 | basesize=10000 | Speedup |
|------|-----------|----------------|---------|
| Window avg time | ? s | ? s | ?x |
| Total runtime | ? s | ? s | ?x |

#### 大ウィンドウ測定

| 項目 | basesize=1 | basesize=10000 | Speedup |
|------|-----------|----------------|---------|
| Window avg time | ? s | ? s | ?x |
| Total runtime | ? s | ? s | ?x |
| Parallel efficiency | ?% | ?% | +?% |

### 3.5 推定所要時間

- Test 3-1: 5分
- Test 3-2: 3分
- Test 3-3: 30分
- Test 3-4: 15分
- 結果確認: 10分
- **合計: 63分**

---

## 📝 Step 4: 結果分析とレポート作成

### 4.1 データ収集

#### 自動抽出スクリプト

```bash
#!/bin/bash
# extract_results.sh

echo "=== Step 1: DHCP単体テスト ==="
grep -E "(basesize|time|iterations)" shared/results/step1_dhcp_basesize_test.log

echo ""
echo "=== Step 2: 10steps fullsize ==="
for bs in 1 1000 10000; do
    echo "--- basesize=$bs ---"
    grep "Total runtime:" shared/results/step2_fullsize_bs${bs}.log
    grep "DHCP.*avg" shared/results/step2_fullsize_bs${bs}.log
    grep "Adjoint.*avg" shared/results/step2_fullsize_bs${bs}.log
done

echo ""
echo "=== Step 3: Sliding window ==="
for test in small_bs1 small_bs10000 large_bs1 large_bs10000; do
    echo "--- $test ---"
    grep "Total runtime:" shared/results/step3_${test}.log
    grep "Window.*avg" shared/results/step3_${test}.log
done
```

### 4.2 比較表作成

#### 総合比較表

| テスト | 計算規模 | basesize=1 | basesize=10000 | 高速化率 |
|--------|----------|-----------|----------------|----------|
| Step 1: DHCP | 1ステップ | ? ms | ? ms | ?x |
| Step 2: 10steps | 10ステップ | ? s | ? s | ?x |
| Step 3: Small window | 5ステップ×2 | ? s | ? s | ?x |
| Step 3: Large window | 71ステップ×1 | ? s | ? s | ?x |

#### test_basesize_performanceとの比較

| 測定環境 | basesize=1 | basesize=10000 | 高速化率 | 実効率 |
|----------|-----------|----------------|----------|--------|
| 理論値（単純計算） | 64.080 ms | 0.037 ms | 1720.3x | 100% |
| Step 1（DHCP） | ? ms | ? ms | ?x | ?% |
| Step 2（Full） | ? s | ? s | ?x | ?% |
| Step 3（Sliding） | ? s | ? s | ?x | ?% |

**実効率の計算**:
```
実効率 = (実測高速化率 / 理論高速化率) × 100%
```

### 4.3 レポート構成

**ファイル**: `docs/reports/phase5_2_real_world_validation_report.md`

#### 目次

1. **エグゼクティブサマリー**
   - 検証結果の概要
   - 主要な発見事項

2. **検証方法**
   - 3段階の検証アプローチ
   - ソルバー設定（bicgstab+gs統一）

3. **Step 1結果: DHCP単体**
   - 性能測定データ
   - 考察

4. **Step 2結果: 10steps fullsize**
   - 性能測定データ
   - 各ソルバーの詳細分析

5. **Step 3結果: Sliding window**
   - 小ウィンドウ結果
   - 大ウィンドウ結果
   - 実用性評価

6. **総合分析**
   - 理論値との比較
   - 実効率の評価
   - ボトルネック分析

7. **結論と推奨事項**
   - 最適basesize値の決定
   - 実装への推奨事項
   - 今後の改善案

8. **付録**
   - 全測定データ
   - 実行コマンド一覧

### 4.4 推定所要時間

- データ収集・整理: 15分
- 比較表作成: 10分
- レポート執筆: 30分
- レビュー・修正: 15分
- **合計: 70分**

---

## 📅 全体スケジュール

### タイムライン

| Step | タスク | 所要時間 | 累積時間 |
|------|--------|----------|----------|
| 1 | test_dhcp_solver.jl更新・実行 | 30分 | 30分 |
| 2 | run_10steps_fullsize_test.jl実行 | 23分 | 53分 |
| 3 | run_sliding_window.jl実行（4回） | 63分 | 116分 |
| 4 | 結果分析・レポート作成 | 70分 | 186分 |
| **合計** | | **186分 (3.1時間)** | |

### マイルストーン

- [x] Phase 5.1完了（basesize実装）
- [ ] Step 1完了（DHCP検証）
- [ ] Step 2完了（Full計算検証）
- [ ] Step 3完了（Sliding window検証）
- [ ] Phase 5.2完了（レポート提出）

---

## 🚨 リスクと対策

### リスク1: 期待した高速化が得られない

**原因候補**:
- ソルバー内部の非並列部分がボトルネック
- メモリアクセスパターンが非効率
- 前処理器のオーバーヘッド

**対策**:
- プロファイリングツールで詳細分析
- 別のbasesize値を試行
- 前処理器の変更を検討

### リスク2: 数値精度の劣化

**症状**:
- CGM収束性の悪化
- 温度場・熱流束の誤差増大

**対策**:
- basesizeを中間値（1000）に調整
- 収束判定基準の見直し
- Python版との詳細比較

### リスク3: メモリ不足

**症状**:
- OOM (Out of Memory)
- スワップ発生による性能低下

**対策**:
- ウィンドウサイズを縮小
- スレッド数を削減
- メモリプロファイリング実施

---

## 📊 成果物

### 実行完了後の成果物

1. **ログファイル**: 7個
   - `step1_dhcp_basesize_test.log`
   - `step2_fullsize_bs1.log`
   - `step2_fullsize_bs1000.log`
   - `step2_fullsize_bs10000.log`
   - `step3_small_bs1.log`
   - `step3_small_bs10000.log`
   - `step3_large_bs1.log`
   - `step3_large_bs10000.log`

2. **レポート**: 1個
   - `docs/reports/phase5_2_real_world_validation_report.md`

3. **更新コード**: 1個
   - `julia/test/test_dhcp_solver.jl`（basesize測定機能追加）

4. **スクリプト**: 1個
   - `extract_results.sh`（結果抽出スクリプト）

---

## 🎯 次のアクション

### 即座に開始可能

```bash
# Step 1開始
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# test_dhcp_solver.jlの確認
head -50 julia/test/test_dhcp_solver.jl

# 必要に応じて更新
# （backend設定機能、性能測定コード追加）
```

### 承認後に実施

- [ ] Step 1実行
- [ ] Step 2実行
- [ ] Step 3実行
- [ ] レポート作成
- [ ] Phase 5.2完了コミット

---

**計画書作成者**: Claude
**レビュー**: 未実施
**承認**: 未承認

---

**準備完了。Step 1から開始しますか？**
