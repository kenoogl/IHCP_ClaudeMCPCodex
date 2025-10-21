# 次セッション継続ガイド

**日付**: 2025年10月21日
**ブランチ**: `sliding-window-validation`
**最終更新**: Python版とJulia版の3者比較実行完了

---

## 🚨 重大な問題発見

### Python版の熱流束が異常に小さい

**症状**:
- Python版の熱流束: -5.49e-06 ~ 7.79e-07 W/m²（異常に小さい）
- Julia版の熱流束: -4.28e+04 ~ 1.38e+05 W/m²（正常範囲）
- **差**: 約10^10倍の違い 🔴

**実行結果ファイル**:
```
shared/results/python_numba8_cgm3.npz       (Python版 8並列)
shared/results/python_numba2_cgm3.npz       (Python版 2並列)
shared/results/julia_sliding_window_cgm3.npz (Julia版)
```

---

## 📋 前セッションで完了した作業

### ✅ 完了事項

1. **run_sliding_window.pyの修正**
   - パス問題を修正（`temp_filename`をファイル名のみに変更）
   - `python/scripts/run_sliding_window.py`

2. **Python版実行（2ケース）**
   - Numba 8並列: 11.52秒 ✅
   - Numba 2並列: 11.80秒 ✅
   - **数値的に完全一致（誤差=0）**

3. **Julia版実行**
   - pcg+diagonal: 165.71秒 ✅
   - Python版の14.4倍遅い

4. **並列性能分析**
   - 2並列と8並列でほぼ同じ時間（差0.27秒）
   - **理由**: CG法の逐次性、小規模問題、メモリボトルネック
   - 並列化効率: わずか0.8%

---

## 🔍 次セッションの優先タスク

### Task 1: Python版の熱流束異常値の原因調査 🔴 最優先

**仮説**:
1. **CGM反復数不足**: 3回では収束していない可能性
2. **初期値の問題**: `q_init_value=0.0`が不適切
3. **単位変換の問題**: 結果を保存する際の単位ミス
4. **スライディングウィンドウの実装差異**: ウィンドウ結合時のバグ
5. **境界条件の設定ミス**: 熱流束境界条件が正しく設定されていない

**調査手順**:

#### Step 1: CGM収束状況の確認
```bash
# ログから目的関数（J）の値を確認
grep "J = " shared/results/python_numba8_cgm3.log | head -20
grep "J = " shared/results/julia_sliding_window_cgm3.log | head -20

# Jの値が減少しているか確認
# Python版: J = 2.85404e+03 → 大きい（未収束の可能性）
# Julia版: J = 5.753417e+03 → さらに大きい
```

#### Step 2: 中間データの比較
```python
# Python版の中間データを確認
import numpy as np
py_data = np.load('shared/results/python_numba8_cgm3.npz')
print("Python版の全キー:", list(py_data.keys()))
print("T0 range:", py_data['T0'].min(), "~", py_data['T0'].max())
print("Y_obs range:", py_data['Y_obs'].min(), "~", py_data['Y_obs'].max())

# Julia版と比較
jul_data = np.load('shared/results/julia_sliding_window_cgm3.npz')
print("\nJulia版の全キー:", list(jul_data.keys()))
print("T_final range:", jul_data['T_final'].min(), "~", jul_data['T_final'].max())
```

#### Step 3: オリジナルコードの確認
```bash
# 保存処理を確認
grep -A10 "np.save.*filename" python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py

# 単位変換の有無を確認
grep -i "unit\|scale\|convert" python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py
```

#### Step 4: CGM反復数を増やして再実行
```bash
# CGM反復数を20に増やして再実行
python python/scripts/run_sliding_window.py \
  --nt 10 --cgm-iter 20 --window 5 --overlap 2 \
  --output python_numba8_cgm20 \
  | tee shared/results/python_numba8_cgm20.log

# 結果を確認
python -c "
import numpy as np
data = np.load('shared/results/python_numba8_cgm20.npz')
q = data['q_result']
print(f'Heat-flux range: [{q.min():.6e}, {q.max():.6e}] W/m²')
"
```

#### Step 5: Julia版との詳細比較
```python
# ウィンドウごとの結果を比較
import numpy as np

py_data = np.load('shared/results/python_numba8_cgm3.npz')
jul_data = np.load('shared/results/julia_sliding_window_cgm3.npz')

# Juliaのウィンドウ情報
print("Julia版ウィンドウ情報:")
print("  window_q_min:", jul_data['window_q_min'])
print("  window_q_max:", jul_data['window_q_max'])
print("  window_J_final:", jul_data['window_J_final'])
```

---

## 📁 重要ファイル

### 修正済み
- `python/scripts/run_sliding_window.py` ⭐ パス問題修正済み

### 実行結果
- `shared/results/python_numba8_cgm3.npz` - Python 8並列（11.52秒）
- `shared/results/python_numba2_cgm3.npz` - Python 2並列（11.80秒）
- `shared/results/julia_sliding_window_cgm3.npz` - Julia版（165.71秒）
- `shared/results/python_numba8_cgm3.log` ⭐ デバッグ用ログ
- `shared/results/julia_sliding_window_cgm3.log` ⭐ デバッグ用ログ

### オリジナルコード
- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py` ⭐ バグ調査対象

---

## 🔧 未コミットの変更

```bash
# 修正したファイル
modified:   python/scripts/run_sliding_window.py

# 新規ファイル（.npyは中間ファイル、不要）
python/python_numba2_cgm3.npy
python/python_numba8_cgm3.npy
```

---

## 📊 次セッション開始時の手順

1. **状態確認**
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -3
cat TODO_NEXT_SESSION.md
```

2. **Task 1実行: Python版バグ調査**
   - 上記のStep 1~5を順次実行
   - 各ステップで結果を記録

3. **修正とテスト**
   - バグを特定したら修正
   - 再実行して結果を確認
   - Julia版と数値比較

4. **コミット**
   - バグ修正をコミット
   - 検証レポート作成

---

## 📚 参照ドキュメント

### 既存の重要資料
- `.claude/CLAUDE.md` - プロジェクト全体のガイドライン
- `docs/plans/sliding_window_validation_plan.md` - 検証計画
- `shared/results/window_splitting_validation_report.md` - ウィンドウ分割検証

---

**次セッションでの成功を祈ります！**

**重要**: Python版の熱流束が異常に小さい原因を特定し、修正することが最優先課題です。
