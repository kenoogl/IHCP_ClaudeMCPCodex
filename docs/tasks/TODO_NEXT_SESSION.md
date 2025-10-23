# 次セッション: Python版CGM 4桁差問題の解決

**作成日**: 2025年10月24日
**ブランチ**: parallelization
**状態**: 🔴 未解決（Julia版と4桁差）
**前回セッション成果**: 面積スケール・時間ステップ修正完了、10^8倍改善達成

## 📊 現状

### 修正完了項目
1. ✅ **面積スケール問題**（コミット: f2dcc47）
   - 勾配・感度を`cell_area = dx * dy`で割り戻し
   - 効果: 10^8倍改善（熱流束 1e-6 → 1e0 W/m²レベル）

2. ✅ **時間ステップインデックス修正**
   - `T_cal[t]` → `T_cal[t+1]`、`Y_obs[t]` → `Y_obs[t+1]`
   - 効果: わずかに改善（勾配 6.28 → 8.06）

3. ✅ **随伴初期値への残差注入**（外部修正済み）
   - Julia版`:residual`戦略に対応
   - Line 1400-1403追加

### 残る問題: Julia版との4桁差

**CGM 2回、window 5、overlap 2での比較**:

| 項目 | Julia版 | Python版（最新） | 差（倍率） |
|------|---------|-----------------|----------|
| **熱流束最大値** | 1.10e+05 W/m² | 1.56 W/m² | **70,577倍** ❌ |
| **熱流束最小値** | -3.37e+04 W/m² | -9.05 W/m² | **3,724倍** ❌ |
| **勾配最大値** | ? | 8.06 W/m² | ? |

## 🎯 次セッションの対応策

### 優先度1: 残差符号の完全検証

**仮説**: 残差の符号が逆の可能性

#### 確認箇所

**Julia版**: `julia/src/solvers/AdjointSolver.jl` Line 147-154
```bash
grep -B 2 -A 5 "residual =" julia/src/solvers/AdjointSolver.jl
```

**Python版**: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2_linearized.py` Line 1333
```python
rhs += 2 * residual_scale * (T_cal[i, j] - Y_obs[i, j]) * dx * dy
```

#### 修正手順

1. Julia版の符号を確認
2. もし `(Y_obs - T_cal)` なら Python版も修正
3. テスト実行して効果確認

### 優先度2: Julia版勾配値の直接確認

同条件でJulia版の勾配を出力し比較

### 優先度3: 中間値の完全比較

T_cal、lambda_field、gradをnpz保存して直接比較

## 📁 関連ファイル

- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2_linearized.py`
- `julia/src/solvers/AdjointSolver.jl`
- `shared/results/python_timestep_fix_test.log`
- `shared/results/julia_sliding_window_cgm2.npz`

## 🎯 成功基準

- Python版熱流束: 10^4 ~ 10^5 W/m²（Julia版と同じオーダー）
- 相対誤差: 10% 以内
