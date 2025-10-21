# 次セッション作業ガイド - Julia版CGM反復数増加問題診断

**日付**: 2025年10月21日
**ブランチ**: sliding-window-validation
**状態**: シナリオ2完了、根本原因特定段階

---

## 現状サマリー

### 発見した問題

**Julia版はPython版の2-3倍の線形ソルバー反復数を使用**

| CGM反復 | ソルバー | Python | Julia(pcg+none) | 比率 |
|---------|---------|--------|-----------------|------|
| 0 | DHCP | 231 | 870 | **3.8倍** |
| 0 | Adjoint | 705 | 851 | 1.2倍 |
| 0 | Sensitivity | 491 | 716 | 1.5倍 |
| 1 | DHCP | 334 | 914 | **2.7倍** |

この反復数の差が、CGM 3回目で問題が悪化する根本原因。

---

## 完了した診断 (codexシナリオ1-2)

### ✅ シナリオ1: 現象の再現ログ取得

**Python側**:
- ログ: `shared/results/python_10steps_cgm3_with_iters.log`
- データ: `shared/results/python_10steps_cgm3.npz`
- 各時刻ステップの反復数カウント追加済み (1200-1214, 1355-1370行)

**Julia側**:
- ログ: `shared/results/julia_10steps_cgm3_detailed.log`
- CGMSolver.jlに詳細ログ出力追加 (491-504行目)

### ✅ シナリオ2: 線形ソルバーの内訳比較

**実施した対策と効果**:

1. **`use_previous_solution=false`**:
   - CGM0 DHCP: 870 → 0 (-100%) ← 劇的改善
   - CGM0 Adjoint: 851 → 318 (-62.6%)
   - しかしPython版との差は残る

2. **diagonal前処理**:
   - CGM0合計: 2437 → 2082 (-15%)
   - CGM1合計: 2468 → 1998 (-19%)
   - 改善するが不十分

3. **rtol設定確認**:
   - Python/Julia完全一致 (DHCP: 1e-6, Adjoint/Sens: 1e-8)
   - これは原因ではない ❌

---

## 根本原因（最有力候補）

### 線形ソルバーの実装差

| 項目 | Python | Julia |
|------|--------|-------|
| ソルバー | SciPy CG | 自前PCG |
| 係数行列 | CSR明示的 | マトリクスフリー |

**仮説**:
- 同じrtol設定でも、残差計算や収束判定の実装差で反復数が変わる
- SciPy CGは高度に最適化されており、より少ない反復で収束

---

## 次セッションでの作業候補

### 推奨: codexシナリオ5 (最小ケース検証)

**実行方法**:
```bash
julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 3 --window 2 --cgm-iter 3 \
  --solver pcg --precond diagonal
```

**確認ポイント**:
- 最小ケースでも同様の反復数比率が出るか
- どのソルバー(DHCP/Adjoint/Sensitivity)で顕著か

---

## 重要なファイル

### コード変更
- `julia/src/solvers/CGMSolver.jl` (491-504行: 詳細ログ追加)
- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py` (1200-1214, 1355-1370行: コールバック追加)

### データ
- `shared/results/python_10steps_cgm3_with_iters.log`
- `shared/results/julia_10steps_cgm3_detailed.log`
- `shared/results/julia_10steps_cgm3_pcg_diagonal.log`

---

## codexへの報告（次セッション開始時）

```
【発見した問題】
Julia版はPython版の2-3倍の線形ソルバー反復数を使用（特にDHCPで3.8倍）。

【実施済み診断】
- シナリオ1-2完了
- use_previous_solution=false: 大幅改善するが不十分
- diagonal前処理: 15-20%改善だが不十分
- rtol設定: 完全一致（原因ではない）

【根本原因（仮説）】
線形ソルバーの実装差（SciPy CG vs 自前PCG）

【次のステップ候補】
1. シナリオ5: 最小ケース検証（推奨）
2. シナリオ3: 中間データ数値比較
3. Krylov.jl標準CG使用

指示をお願いします。
```
