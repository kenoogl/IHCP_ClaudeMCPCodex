# 中規模ウィンドウPython-Julia性能比較計画書（推奨版）

**作成日**: 2025年10月24日
**バージョン**: 1.0（推奨版）
**ステータス**: 実行準備完了
**目的**: 中規模ウィンドウ（nt=100, window=60, overlap=15）でのPython-Julia実装の性能と数値精度を比較検証する
**関連**: Phase 5.2 Step 3 Phase 2の代替計画（現実的な実行時間）

---

## 🎯 この計画の特徴

### ✅ **現実的な実行時間**
- **予想実行時間**: Phase 1（Julia basesize最適化）: 約5-6時間
- **予想実行時間**: Phase 2（Python-Julia比較）: 約6-8時間
- **合計**: 約12-14時間（1日で完了可能）

### ✅ **適切なウィンドウ構成**
- **ウィンドウ数**: 2個（適度なオーバーラップ）
- **カバー率**: 100%（全時間範囲をカバー）
- **Window 1**: [0, 60]（60ステップ）
- **Window 2**: [45, 100]（55ステップ、15ステップオーバーラップ）

### ✅ **小→中のスケーラビリティ評価**
- 小ウィンドウ: window=5 ✅ 完了
- **中ウィンドウ**: window=60 ← 今回
- 大ウィンドウ: window=71 ← 将来

---

## 1. 背景と目的

### 1.1 背景

**完了済み**:
- ✅ 小ウィンドウ（window=5, nt=10）でのPython-Julia比較完了
  - Python版: 3.3~3.6倍高速、数値不安定（発散傾向）
  - Julia版: 数値安定、実用推奨

**未完了**:
- ⏸️ 大ウィンドウ（window=71）での測定
  - 問題: nt=10では不足、nt=100だと実行時間30-45時間

**この計画の位置づけ**:
- 中規模ウィンドウ（window=60, nt=100）で現実的な測定
- 小→中のスケーラビリティ評価
- basesize最適化の完全検証

### 1.2 目的

1. **basesize最適値の特定**: 中規模ウィンドウでの最適basesize特定
2. **性能比較**: Python-Julia実行速度、ソルバー反復回数の比較
3. **数値安定性比較**: 熱流束範囲、CGM収束挙動の比較
4. **スケーラビリティ評価**: 小→中ウィンドウでの性能変化

---

## 2. 実行パラメータ

### 2.1 共通パラメータ

```
格子サイズ:
  - ni × nj × nk = 80 × 100 × 20 (N=160,000セル)
  - dx = 0.12e-3 m
  - dy = 0.12e-3 * sin(80°) / sin(45°) m
  - Lz = 0.5e-3 m (stretch_factor = 3.0)

時間設定:
  - dt = 1.0e-3 s (1ms)
  - nt = 100 ステップ

スライディングウィンドウ（中規模）:
  - window_size = 60 ステップ
  - overlap = 15 ステップ
  - ウィンドウ数: 2個
  - Window 1: [0, 60] (60ステップ)
  - Window 2: [45, 100] (55ステップ)

並列化:
  - スレッド数: 4
  - CPU: Apple Silicon (4コア使用)
```

### 2.2 CGM設定

**Phase 1 & 2共通**:
```
CGM反復数: 3回（標準設定）
ソルバー設定:
  - DHCP: rtol=1e-6, maxiter=20000
  - Adjoint: rtol=1e-8, maxiter=20000
  - Sensitivity: rtol=1e-8, maxiter=20000
```

---

## 3. Phase 1: Julia版basesize最適化

### 3.1 測定パターン（3パターン）

**目的**: 中規模ウィンドウでの最適basesize特定

#### Test 1-1: basesize=500（小ウィンドウ最適値）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 100 \
  --window 60 \
  --overlap 15 \
  --cgm-iter 3 \
  --solver pbicgstab \
  --precond gs \
  --basesize 500 \
  2>&1 | tee shared/results/medium_window_basesize500_cgm3.log
```

**実行時間見込み**: 約1.5-2時間

#### Test 1-2: basesize=1000（Step 2最適値）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 100 \
  --window 60 \
  --overlap 15 \
  --cgm-iter 3 \
  --solver pbicgstab \
  --precond gs \
  --basesize 1000 \
  2>&1 | tee shared/results/medium_window_basesize1000_cgm3.log
```

**実行時間見込み**: 約1.5-2時間

#### Test 1-3: basesize=2000（中規模用推定値）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 100 \
  --window 60 \
  --overlap 15 \
  --cgm-iter 3 \
  --solver pbicgstab \
  --precond gs \
  --basesize 2000 \
  2>&1 | tee shared/results/medium_window_basesize2000_cgm3.log
```

**実行時間見込み**: 約1.5-2時間

**Phase 1合計**: 約5-6時間

---

## 4. Phase 2: Python-Julia完全比較

### 4.1 測定パターン（2パターン、CGM 3回のみ）

#### Test 2-1: Python版 CGM 3回

```bash
export NUMBA_NUM_THREADS=4
python3 python/scripts/run_sliding_window.py \
  --nt 100 \
  --cgm-iter 3 \
  --window 60 \
  --overlap 15 \
  --output python_medium_window_cgm3 \
  2>&1 | tee shared/results/python_medium_window_cgm3.log
```

**実行時間見込み**: 約2-3時間

#### Test 2-2: Julia版 CGM 3回（最適basesize使用）

```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 100 \
  --window 60 \
  --overlap 15 \
  --cgm-iter 3 \
  --solver pbicgstab \
  --precond gs \
  --basesize <最適値> \
  2>&1 | tee shared/results/julia_medium_window_cgm3.log
```

**実行時間見込み**: 約2時間（既にPhase 1で実行済みの可能性）

**Phase 2合計**: 約4-5時間（または2-3時間、Phase 1結果流用の場合）

---

## 5. 比較項目

### 5.1 性能比較

**実行速度**:
- 総実行時間（Python vs Julia）
- ウィンドウあたりの平均時間
- 各ソルバー実行時間（DHCP/Adjoint/Sensitivity）

**ソルバー効率**:
- 各ソルバー反復回数（平均、最小、最大）
- 反復回数の推移（ウィンドウ間、CGM反復間）

### 5.2 数値精度比較

**熱流束**:
- 熱流束範囲（min/max/mean）の比較
- Python-Julia間の相対誤差
- 熱流束振幅の比較（発散傾向の確認）

**CGM収束性**:
- 目的関数Jの推移
- CGM減少率（各反復での減少率）

### 5.3 スケーラビリティ評価

**小→中ウィンドウでの変化**:

| 項目 | 小ウィンドウ<br>(5 steps) | 中ウィンドウ<br>(60 steps) | スケーリング |
|------|-------------------------|--------------------------|------------|
| 実行時間 | 34.34秒（Julia） | 予想: 約2時間 | 約200倍 |
| 反復回数 | 実測値 | 測定予定 | 増加率？ |
| 熱流束振幅 | 1.81e+05 W/m² | 測定予定 | 安定性維持？ |

---

## 6. 期待される成果

### Phase 1完了時（5-6時間後）

- ✅ 中規模ウィンドウでの最適basesize特定
- ✅ basesize依存性の確認（小ウィンドウとの違い）
- ✅ 中規模ウィンドウでのJulia版性能プロファイル

### Phase 2完了時（合計12-14時間後）

- ✅ **中規模ウィンドウでのPython-Julia完全比較**
- ✅ 小→中スケーラビリティの定量評価
- ✅ 実用シナリオでの推奨実装の確定
- ✅ **Phase 5.2の実質完了**（小+中ウィンドウ測定完了）

---

## 7. 成功基準

### Phase 1

- [ ] 3パターンのbasesize測定完了
- [ ] 最速構成の特定（誤差5%以内）
- [ ] 数値精度の完全一致確認（全パターンで同一結果）

### Phase 2

- [ ] Python-Julia両実装がエラーなく完走
- [ ] 実行時間の定量比較完了
- [ ] 熱流束相対誤差の評価（発散傾向の定量化）
- [ ] スケーラビリティの評価完了

---

## 8. 実行スケジュール

### タイムライン（合計約12-14時間）

**Phase 1: Julia版basesize最適化**（5-6時間）
```
00:00 - 02:00  Test 1-1: basesize=500
02:00 - 04:00  Test 1-2: basesize=1000
04:00 - 06:00  Test 1-3: basesize=2000
```

**中間分析**（15分）
```
06:00 - 06:15  Phase 1結果分析、最適basesize判定
```

**Phase 2: Python-Julia比較**（6-8時間）
```
06:15 - 09:00  Test 2-1: Python CGM 3回（2.5-3時間見込み）
09:00 - 11:00  Test 2-2: Julia CGM 3回（2時間見込み、もしくはPhase 1結果流用でスキップ）
```

**最終分析とレポート作成**（2時間）
```
11:00 - 13:00  結果統合、比較分析、レポート作成
```

**⚠️ 注意**: 1日で完了可能（約13時間）

---

## 9. リスクと対策

### リスク1: 実行時間が予想より長い

**対策**:
- Phase 1で実行時間を確認し、Phase 2の見積もりを修正
- 必要に応じてbasesize測定を2パターンに削減（500, 1000のみ）

### リスク2: Python版が発散する

**対策**:
- 小ウィンドウでの発散傾向を踏まえ、中ウィンドウでも同様と予想
- 発散した場合もそれを「結果」として記録
- Julia版の数値安定性を再確認

---

## 10. 成果物

### Phase 1完了時

1. **実行ログ**（3件）:
   ```
   shared/results/medium_window_basesize500_cgm3.log
   shared/results/medium_window_basesize1000_cgm3.log
   shared/results/medium_window_basesize2000_cgm3.log
   ```

2. **中間レポート**:
   ```
   docs/reports/medium_window_basesize_optimization.md
   ```

### Phase 2完了時

3. **実行ログ**（2件）:
   ```
   shared/results/python_medium_window_cgm3.log
   shared/results/julia_medium_window_cgm3.log
   ```

4. **データファイル**（2件）:
   ```
   shared/results/python_medium_window_cgm3.npz
   shared/results/julia_medium_window_cgm3.npz
   ```

5. **最終比較レポート**:
   ```
   docs/reports/medium_window_python_julia_comparison_final.md
   ```

---

## 11. 次のステップ（Phase 2完了後）

### オプション1: Phase 5.2完了宣言（推奨）

- **Phase 1-6**: 505テスト全合格 ✅
- **Phase 5.2**: basesize最適化、Python-Julia比較完了 ✅
  - 小ウィンドウ測定完了 ✅
  - **中ウィンドウ測定完了 ✅**（今回）
- **最終結論**: Julia版実用推奨、Python版は発散傾向
- **次の行動**: 論文・報告書作成、研究成果の公表

### オプション2: 大ウィンドウ測定継続（将来）

- **nt=200-300**: より大規模な測定
- **window=71**: オリジナルPython設定での測定
- **実行時間**: 数日間の計算時間

---

## 付録: コマンド一覧

### Phase 1: Julia版basesize最適化

```bash
# Test 1-1: basesize=500
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 500 \
  2>&1 | tee shared/results/medium_window_basesize500_cgm3.log

# Test 1-2: basesize=1000
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 1000 \
  2>&1 | tee shared/results/medium_window_basesize1000_cgm3.log

# Test 1-3: basesize=2000
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 2000 \
  2>&1 | tee shared/results/medium_window_basesize2000_cgm3.log
```

### Phase 2: Python-Julia比較

```bash
# Test 2-1: Python CGM 3回
export NUMBA_NUM_THREADS=4
python3 python/scripts/run_sliding_window.py \
  --nt 100 --cgm-iter 3 --window 60 --overlap 15 \
  --output python_medium_window_cgm3 \
  2>&1 | tee shared/results/python_medium_window_cgm3.log

# Test 2-2: Julia CGM 3回（最適basesize使用）
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 100 --window 60 --overlap 15 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize <最適値> \
  2>&1 | tee shared/results/julia_medium_window_cgm3.log
```

---

**計画書バージョン**: 1.0（推奨版）
**最終更新**: 2025年10月24日
**作成者**: Claude Code
**承認**: 未承認
**ステータス**: 実行準備完了
**推奨度**: ★★★★★（最推奨）

---

## この計画が推奨される理由

1. **現実的な実行時間**: 1日で完了可能（12-14時間）
2. **適切なウィンドウ構成**: 2個のウィンドウで100%カバー
3. **小→中のスケーラビリティ評価**: 段階的な性能評価が可能
4. **Phase 5.2の実質完了**: 小+中ウィンドウで十分な知見を取得
5. **次フェーズへの円滑な移行**: プロジェクト完了または性能改善へ

**この計画の採用を強く推奨します。**
