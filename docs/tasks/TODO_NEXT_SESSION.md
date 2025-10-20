# 次セッション作業ガイド - スライディングウィンドウ検証

**作成日**: 2025年10月21日
**ブランチ**: `sliding-window-validation`
**最終コミット**: `6d4e5cb` - Julia診断機能とスライディングウィンドウ検証スクリプト作成

---

## 現在の状況

### 完了した作業 ✅

1. **計画書作成** (b4a8e59)
   - `docs/plans/sliding_window_validation_plan.md` (詳細計画、9章構成)
   - `docs/plans/sliding_window_validation_plan_summary.md` (クイックサマリー)

2. **Python側実装** (27c040d)
   - `python/validation/solvers_with_diagnostics.py` (450行)
     - 診断機能付きソルバー（DHCP/Adjoint）
     - 反復回数・残差・実行時間を記録
   - `python/validation/run_sliding_window_validation.py` (380行)
     - スライディングウィンドウ検証スクリプト
     - コマンドライン引数対応

3. **Julia側実装** (6d4e5cb)
   - `julia/src/solvers/SolverDiagnostics.jl` (110行)
     - 診断情報構造体（Python互換）
   - `julia/scripts/run_sliding_window_validation.jl` (570行)
     - スライディングウィンドウ検証スクリプト
     - Python版と同等の機能
   - `julia/src/IHCP_CGM.jl` (SolverDiagnostics統合)

### 残りの作業 ⏳

#### 1. 比較スクリプト作成 🔥 **次の最優先タスク**
- [ ] `python/validation/compare_sliding_window_results.py` 作成
  - Python/Julia両方の結果を読み込み
  - ウィンドウごとの詳細比較
  - Markdownレポート生成

#### 2. Phase 1実行（CGM 3回）
- [ ] Python実行: `python python/validation/run_sliding_window_validation.py --cgm-iter 3`
- [ ] Julia実行: `julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3`
- [ ] 比較実行: `python python/validation/compare_sliding_window_results.py --phase 1`
- [ ] 結果確認・問題修正

#### 3. Phase 2実行（CGM 20000回）
- [ ] Python実行: `python python/validation/run_sliding_window_validation.py --cgm-iter 20000`
- [ ] Julia実行: `julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 20000`
- [ ] 比較実行: `python python/validation/compare_sliding_window_results.py --phase 2`
- [ ] 詳細分析

---

## 次セッション開始手順

### 1. 環境確認

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git branch
# 確認: sliding-window-validation ブランチにいること
```

### 2. 最新状態の取得

```bash
git pull origin sliding-window-validation
git log --oneline -5
# 確認: 6d4e5cb が最新コミット
```

### 3. 計画書の確認

```bash
# クイックサマリー確認
cat docs/plans/sliding_window_validation_plan_summary.md

# 詳細計画確認（必要に応じて）
cat docs/plans/sliding_window_validation_plan.md
```

### 4. 比較スクリプト作成（最優先）

**作成ファイル**: `python/validation/compare_sliding_window_results.py`

**必要機能**:
```python
# 引数パース
--phase 1 or 2
--python-result shared/results/python_sliding_window_cgm{N}.npz
--julia-result shared/results/julia_sliding_window_cgm{N}.npz
--output-report shared/results/sliding_window_comparison_cgm{N}.md

# 比較項目
1. 性能比較
   - 総実行時間
   - ウィンドウあたりの平均時間
   - 各ソルバーの累積実行時間（将来拡張）

2. 精度比較
   - 熱流束誤差: 絶対誤差、相対誤差、RMS誤差
   - 温度場誤差（最終状態）
   - ウィンドウごとの目的関数値比較

3. 収束性比較
   - 各ウィンドウのCGM反復回数
   - J_history比較（目的関数推移）
```

**出力レポート例**:
```markdown
# Python-Julia スライディングウィンドウ比較レポート

## 1. 実行環境
- CGM反復数: 3
- 時間ステップ数: 300
- ウィンドウサイズ: 71, オーバーラップ: 17

## 2. 性能比較
| 項目 | Python | Julia | Julia/Python比 |
|------|--------|-------|----------------|
| 総実行時間 | 123.45s | 67.89s | 54.9% |
| 平均ウィンドウ時間 | 20.58s | 11.32s | 55.0% |

## 3. 精度比較
### 3.1 熱流束誤差
- 絶対誤差: max=1.23e-3, mean=4.56e-5, RMS=7.89e-5
- 相対誤差: max=2.34%, mean=0.12%, RMS=0.45%

### 3.2 ウィンドウ別比較
...
```

---

## 実行パラメータ（計画書より）

### 共通パラメータ
```
格子: 80 × 100 × 20 (フルサイズ)
時間ステップ: 300
dt: 1.0e-3 s
ウィンドウサイズ: 71
オーバーラップ: 17
予想ウィンドウ数: 約6
```

### Phase 1（短時間確認）
```
CGM反復数: 3回固定
目的: スクリプト動作確認、基本一致性検証
推定実行時間: Python約15分、Julia約5分
```

### Phase 2（本格検証）
```
CGM反復数: 20000回（早期停止あり）
目的: オリジナルと同じ設定での完全比較
推定実行時間: 状況により変動
```

---

## 成功基準

### Phase 1（必須）
- [ ] 両スクリプトがエラーなく完走
- [ ] 熱流束相対誤差 < 5%
- [ ] ウィンドウ数の一致
- [ ] 各ウィンドウの目的関数値が近い（相対誤差<10%）

### Phase 2（目標）
- [ ] 温度場相対誤差 < 2%
- [ ] 熱流束相対誤差 < 1%
- [ ] CGM収束挙動の一致（反復回数の差<20%）
- [ ] Julia実行時間がPythonの50%以下

---

## 既知の課題・注意事項

### 1. データファイルの確認
- `shared/data/T_measure_700um_1ms.npy` (1.1GB) が必要
- 存在しない場合は別途取得

### 2. メモリ使用量
- フルサイズ計算: 8GB以上推奨
- メモリ不足時はウィンドウサイズを縮小（71→50）

### 3. Pythonスクリプトの依存関係
- scipy.sparse.linalg.cg使用
- オリジナルコードの関数をインポート
  - thermal_properties_calculator
  - coeffs_and_rhs_building_DHCP
  - assemble_A_DHCP
  - coeffs_and_rhs_building_Adjoint
  - assemble_A_Adjoint

### 4. Juliaスクリプトの依存関係
- ArgParse.jl（コマンドライン引数解析）
- NPZ.jl（Python互換npz形式保存）
- 既存のIHCP_CGM.jlモジュール

---

## 関連ドキュメント

### 必読
- **計画書**: `docs/plans/sliding_window_validation_plan.md`
- **サマリー**: `docs/plans/sliding_window_validation_plan_summary.md`

### 参考
- オリジナルPythonコード: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`
- シングルウィンドウテスト: `julia/scripts/run_10steps_fullsize_test.jl`
- 過去の検証結果: `shared/results/python_julia_10steps_baseline_comparison.md`

---

## クイックスタートコマンド

### 比較スクリプト作成後

```bash
# Phase 1: 動作確認（CGM 3回）
python python/validation/run_sliding_window_validation.py --cgm-iter 3
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
python python/validation/compare_sliding_window_results.py --phase 1

# Phase 2: 本格検証（CGM 20000回） - Phase 1成功後
python python/validation/run_sliding_window_validation.py --cgm-iter 20000
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 20000
python python/validation/compare_sliding_window_results.py --phase 2
```

---

## Todoリスト（優先順位順）

1. 🔥 **比較スクリプト作成** - 最優先
2. Phase 1実行（CGM 3回）
3. Phase 1結果確認・問題修正
4. Phase 2実行（CGM 20000回）
5. Phase 2詳細分析
6. 最終レポート作成
7. プルリクエスト作成（mainブランチへマージ）

---

## 連絡先・リソース

- **プロジェクトルート**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex`
- **ブランチ**: `sliding-window-validation`
- **リモート**: `origin/sliding-window-validation`

---

**次セッションの目標**: 比較スクリプトを作成し、Phase 1（CGM 3回）を実行して基本動作を確認する。

**推定所要時間**:
- 比較スクリプト作成: 30-60分
- Phase 1実行: 20-30分（Python 15分 + Julia 5分 + 比較）
- 問題修正（必要に応じて）: 30-60分

**合計**: 約2-3時間

---

**作成者**: Claude Code
**最終更新**: 2025年10月21日
