# シナリオ3: 中間データ数値比較 - 実行手順書

**目的**: Python版とJulia版の線形ソルバー反復数の差（2-3倍）の原因を、中間データの数値レベルで特定する。

**日付**: 2025年10月21日
**ブランチ**: sliding-window-validation

---

## 📋 概要

### 実施内容

1. **Python側**: CGM反復0〜3の各タイミングで中間データ（T_cal, λ, q, grad, dT）を保存
2. **Julia側**: 同じタイミングで同じデータを保存（配列を転置してPython形状に合わせる）
3. **比較スクリプト**: 要素ごとの差分を計算し、誤差統計を出力

### 期待される発見

- **数値が完全一致** → 線形ソルバーの実装差が原因確定
- **数値に微小な差** → 累積誤差の影響を定量評価
- **配列の対応ミス** → インデックス処理のバグ発見

---

## 🛠️ 実装済み機能

### Python側の変更

**ファイル**: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`

1. **`global_CGM_time`関数**:
   - パラメータ追加: `save_intermediate=False`, `output_dir="shared/results"`
   - CGM反復0〜3で中間データを`numpy.savez`で保存（1527-1540行）

2. **`sliding_window_CGM_q_saving`関数**:
   - パラメータ追加: `save_intermediate=False`, `output_dir="shared/results"`
   - `global_CGM_time`にパラメータ転送（1630行）

**実行スクリプト**: `python/scripts/run_python_10steps_cgm3.py`
- `save_intermediate=True`を指定して実行

### Julia側の変更

**ファイル**: `julia/src/solvers/CGMSolver.jl`

1. **`solve_cgm!`関数**:
   - パラメータ追加: `save_intermediate::Bool=false`, `output_dir::String="shared/results"`
   - CGM反復0〜3で中間データを`NPZ.npzwrite`で保存（493-517行）
   - **重要**: Python形状に転置して保存
     - `(ni, nj, nk, nt)` → `(nt, ni, nj, nk)`
     - `(ni, nj, nt-1)` → `(nt-1, ni, nj)`

**実行スクリプト**: `julia/scripts/run_10steps_fullsize_test.jl`
- 環境変数で制御:
  - `CGM_ITER=3`: CGM反復数
  - `SAVE_INTERMEDIATE=true`: 中間データ保存を有効化

### 比較スクリプト

**ファイル**: `python/validation/compare_intermediate_data.py`

**機能**:
- Python/Juliaの中間データ読み込み
- 配列形状の一致確認
- 要素ごとの絶対誤差・相対誤差計算
- 統計情報の表示
- Markdownレポート自動生成

---

## 🚀 実行手順

### Step 1: Python側実行（約15分）

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

python python/scripts/run_python_10steps_cgm3.py > shared/results/python_scenario3.log 2>&1
```

**生成されるファイル**:
- `shared/results/python_cgm0_data.npz`
- `shared/results/python_cgm1_data.npz`
- `shared/results/python_cgm2_data.npz`
- `shared/results/python_cgm3_data.npz`

**確認**:
```bash
ls -lh shared/results/python_cgm*.npz
```

---

### Step 2: Julia側実行（約5分）

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

CGM_ITER=3 SAVE_INTERMEDIATE=true julia --project=julia julia/scripts/run_10steps_fullsize_test.jl \
  > shared/results/julia_scenario3.log 2>&1
```

**生成されるファイル**:
- `shared/results/julia_cgm0_data.npz`
- `shared/results/julia_cgm1_data.npz`
- `shared/results/julia_cgm2_data.npz`
- `shared/results/julia_cgm3_data.npz`

**確認**:
```bash
ls -lh shared/results/julia_cgm*.npz
```

---

### Step 3: 中間データ比較

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

python python/validation/compare_intermediate_data.py
```

**出力例**:
```
================================================================================
シナリオ3: Python-Julia中間データ数値比較
================================================================================

読み込み: CGM反復0
  Python: python_cgm0_data.npz
  Julia:  julia_cgm0_data.npz

[CGM反復0] 配列形状確認
--------------------------------------------------------------------------------
  ✅ T_cal          : Python (10, 80, 100, 20) vs Julia (10, 80, 100, 20)
  ✅ lambda_field   : Python (10, 80, 100, 20) vs Julia (10, 80, 100, 20)
  ✅ q              : Python (9, 80, 100) vs Julia (9, 80, 100)
  ✅ grad           : Python (9, 80, 100) vs Julia (9, 80, 100)
  ✅ dT             : Python (10, 80, 100, 20) vs Julia (10, 80, 100, 20)

[CGM反復0] 数値誤差比較
================================================================================

T_cal (CGM反復0)
--------------------------------------------------------------------------------
  Python値範囲:  +5.123456e+02 ~ +6.234567e+02 (平均: +5.678901e+02)
  Julia値範囲:   +5.123456e+02 ~ +6.234567e+02 (平均: +5.678901e+02)

  絶対誤差:
    最大: 1.234567e-08
    平均: 2.345678e-10
    中央: 1.234567e-10

  相対誤差:
    最大: 2.345678e-11 (0.0000%)
    平均: 4.567890e-13 (0.000000%)
    中央: 2.345678e-13 (0.000000%)
...
```

**生成されるレポート**:
- `shared/results/scenario3_comparison_report.md`

---

## 📊 期待される結果

### ケース1: 数値が完全一致（相対誤差 < 1e-10）

**解釈**:
- 中間データは一致しているが、線形ソルバー内部の反復過程で差が生じている
- **原因**: 線形ソルバーの実装差（SciPy CG vs 自前PCG）
- **次のアクション**: Krylov.jl標準CGの使用、または収束判定基準の詳細比較

### ケース2: 微小な差（相対誤差 1e-10 〜 1e-6）

**解釈**:
- 浮動小数点演算の累積誤差
- 配列の並び（row-major vs column-major）による計算順序の違い
- **次のアクション**: 誤差の蓄積パターンを分析、許容範囲かを判定

### ケース3: 大きな差（相対誤差 > 1e-6）

**解釈**:
- 配列のインデックス対応ミス
- 境界条件の扱いの違い
- **次のアクション**: 配列の要素を詳細にダンプして対応確認

---

## 🔍 トラブルシューティング

### ファイルが見つからない

**症状**: `FileNotFoundError: Python data not found`

**対処**:
```bash
# Step 1/2が正常に完了しているか確認
ls -lh shared/results/*_cgm*.npz

# ログファイルでエラー確認
tail -50 shared/results/python_scenario3.log
tail -50 shared/results/julia_scenario3.log
```

### 配列形状が一致しない

**症状**: `❌ T_cal: Python (10, 80, 100, 20) vs Julia (80, 100, 20, 10)`

**原因**: Julia側で転置処理が正しく実行されていない

**対処**:
```julia
# CGMSolver.jl の493-517行を確認
# permutedims が正しく呼ばれているか
```

### メモリ不足

**症状**: `MemoryError` または Julia crash

**対処**:
```bash
# CGM反復数を減らす
python python/validation/compare_intermediate_data.py --cgm-iter 0 1
```

---

## 📁 生成されるファイル一覧

```
shared/results/
├── python_cgm0_data.npz          # Python CGM反復0データ
├── python_cgm1_data.npz          # Python CGM反復1データ
├── python_cgm2_data.npz          # Python CGM反復2データ
├── python_cgm3_data.npz          # Python CGM反復3データ
├── julia_cgm0_data.npz           # Julia CGM反復0データ
├── julia_cgm1_data.npz           # Julia CGM反復1データ
├── julia_cgm2_data.npz           # Julia CGM反復2データ
├── julia_cgm3_data.npz           # Julia CGM反復3データ
├── scenario3_comparison_report.md # 比較レポート
├── python_scenario3.log          # Python実行ログ
└── julia_scenario3.log           # Julia実行ログ
```

---

## 🎯 次のステップ

### シナリオ3完了後

1. **レポート確認**: `scenario3_comparison_report.md`を精査
2. **結果に基づく判断**:
   - 数値一致 → シナリオ5（最小ケース）またはKrylov.jl検証
   - 数値不一致 → 配列対応の詳細確認、境界条件の再検証
3. **TODO_NEXT_SESSION.md更新**: 発見事項を記録

---

## 📝 注意事項

### データ品質保証ルール

**必須**:
- 実行完了確認: ログに"Total runtime:"が記録されているか
- ファイル確認: 全CGM反復（0-3）のnpzファイルが生成されているか
- サイズ確認: npzファイルサイズが妥当か（数MB〜数十MB）

**禁止**:
- 推定値・仮定値の使用
- 不完全なデータでのレポート作成
- 未検証データの公開

---

## 🔗 関連ドキュメント

- **診断計画**: `TODO_NEXT_SESSION.md`
- **過去の診断結果**:
  - `shared/results/python_10steps_cgm3_with_iters.log`
  - `shared/results/julia_10steps_cgm3_detailed.log`
- **CLAUDE.md**: データ品質保証ルール（必読）

---

**作成日**: 2025年10月21日
**最終更新**: 2025年10月21日
