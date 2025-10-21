# 次セッション継続ガイド

**日付**: 2025年10月21日
**ブランチ**: `sliding-window-validation`
**最終更新**: ウィンドウ分割ロジック3者完全一致確認完了

---

## 📋 本セッションで完了した作業

### ✅ 完了事項

1. **ウィンドウ分割パラメータテスト（3ケース）**
   - テストケース1: nt=10, window=5, overlap=2 → 5ウィンドウ一致 ✅
   - テストケース2: nt=300, window=71, overlap=17 → 23ウィンドウ一致 ✅
   - テストケース3: nt=10, window=10, overlap=0 → 1ウィンドウ一致 ✅

2. **Julia版スクリプト修正**
   - `julia/scripts/run_sliding_window.jl:581`
   - ドライランモード時の`return 0`修正完了

3. **オリジナルコードへのドライラン機能追加**
   - `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`
   - `sliding_window_CGM_q_saving()`に`dry_run`パラメータ追加
   - ウィンドウ分割のみ表示して早期終了する機能

4. **検証スクリプト作成**
   - `python/scripts/test_original_dryrun.py` - オリジナルコードのドライランテスト
   - `python/scripts/compare_window_splitting.py` - 3者自動比較スクリプト（4テストケース）
   - `python/scripts/run_sliding_window.py` - Python統一実行スクリプト

5. **3者完全一致確認（追加検証4ケース）**
   - テストケース1: nt=10, window=5, overlap=2 → 5ウィンドウ ✅
   - テストケース2: nt=300, window=71, overlap=17 → 23ウィンドウ ✅
   - テストケース3: nt=10, window=10, overlap=0 → 1ウィンドウ ✅
   - テストケース4: nt=20, window=10, overlap=5 → 8ウィンドウ ✅

6. **検証レポート作成・更新**
   - `shared/results/window_splitting_validation_report.md`
   - オリジナル・Python・Julia版の3者比較を中心に簡潔化
   - 検証スクリプトの使用方法を追加

7. **不要ファイル削除**
   - ❌ `python/scripts/run_python_option2.py` (不正確、ウィンドウ分割がJulia版と異なる)
   - ❌ `python/scripts/compare_results_10steps_cgm3.py` (10ステップ専用、今回不要)
   - ❌ `python/scripts/run_python_10steps_cgm3.py` (10ステップ専用、今回不要)

8. **保持した重要ファイル**
   - ✅ `python/scripts/run_python_10steps_cgm3_with_iters.py` (反復数記録ガイド、将来の差異分析用)

---

## 🎯 重要な発見

### ウィンドウ分割ロジックの3者完全一致
**オリジナルコード、Python版（run_sliding_window.py）、Julia版の3者すべてが同一のウィンドウ分割ロジックを使用していることを確認**

| 要素 | 実装 | 一致 |
|------|------|------|
| ループ条件 | `start_idx < (nt - 1)` | ✅ |
| max_L計算 | `min(window_size, (nt-1) - start_idx)` | ✅ |
| step計算 | `max(1, max_L - overlap)` | ✅ |

### Python版の実装構造
- `run_sliding_window.py`はラッパースクリプト
- 実際の計算はオリジナルコードの`sliding_window_CGM_q_saving()`を呼び出し
- ドライラン表示と実際の計算で同じロジックが使用される

---

## 📁 重要ファイル

### 新規作成
- `python/scripts/run_sliding_window.py` ⭐ Python統一実行スクリプト
- `python/scripts/test_original_dryrun.py` - ドライランテスト
- `python/scripts/compare_window_splitting.py` - 3者自動比較
- `shared/results/window_splitting_validation_report.md` ⭐ 検証レポート

### 修正済み
- `julia/scripts/run_sliding_window.jl` - ドライラン機能のバグ修正
- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py` - ドライラン機能追加

### 既存の有効データ
- `shared/results/julia_sliding_window_cgm3.npz` - Julia版結果（nt=10, window=5, overlap=2, cgm_iter=3）
- `shared/results/julia_sliding_window_cgm3_metadata.txt` - メタデータ

---

## 🚀 次セッションの優先タスク

### Task 2: Python版オプション2実行（正しいスクリプト使用）

**重要**: 前セッションの`run_python_option2.py`は削除済み。正しいスクリプト`run_sliding_window.py`を使用すること。

```bash
# Python版実行
python python/scripts/run_sliding_window.py \
  --nt 10 \
  --cgm-iter 3 \
  --window 5 \
  --overlap 2 \
  --output python_option2_correct \
  | tee shared/results/python_option2_correct.log
```

**パラメータ**:
- nt=10 (時間ステップ)
- cgm_iter=3 (CGM最大反復数)
- window=5 (ウィンドウサイズ)
- overlap=2 (オーバーラップ)

**期待される出力**:
- `shared/results/python_option2_correct.npz`
- `shared/results/python_option2_correct_metadata.txt`
- `shared/results/python_option2_correct.log`

**Julia版は既存データ使用可能**:
- `shared/results/julia_sliding_window_cgm3.npz`
- パラメータ: nt=10, window=5, overlap=2, cgm_iter=3, solver=pcg, precond=diagonal

---

### Task 3: Python-Julia結果比較

**比較スクリプト作成が必要**:
- `python/scripts/compare_sliding_window_results.py`

**比較項目**:
1. 熱流束（q_result）の数値誤差
   - 絶対誤差（最大、平均、RMS）
   - 相対誤差（最大、平均）
2. 実行時間の差
3. 各ウィンドウの収束状況
4. メモリ使用量（ログから抽出）

**期待される分析結果**:
- CGM反復数の差（2-3倍）の原因特定
- 数値精度の比較
- 性能差の分析

---

## 🔍 データ品質保証ルール（必須）

**⚠️ 絶対厳守**: レポート・ドキュメント作成時のデータ品質要件

### 禁止事項
1. **推定値・仮定値の使用禁止**
   - 実測データが不完全な場合はレポート作成を中断
2. **不完全データでのレポート作成禁止**
   - 実行が完了していない（"Total runtime:"未記録）
   - 収束していない（発散、エラー終了）
   - データが欠損している
3. **未検証データの公開禁止**
   - ファイルの存在確認必須
   - ファイルサイズ・内容の妥当性確認必須

### 必須手順
```bash
# 実行完了確認
grep "Total runtime:" <logfile>
ls -lh <result_file>
tail -20 <logfile>
```

---

## 📊 次セッション開始時の手順

1. **状態確認**
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5
cat TODO_NEXT_SESSION.md
```

2. **検証レポート確認**
```bash
cat shared/results/window_splitting_validation_report.md
```

3. **Task 2実行（Python版）**
```bash
python python/scripts/run_sliding_window.py \
  --nt 10 --cgm-iter 3 --window 5 --overlap 2 \
  --output python_option2_correct \
  | tee shared/results/python_option2_correct.log
```

4. **実行完了確認**
```bash
grep "Total runtime:" shared/results/python_option2_correct.log
ls -lh shared/results/python_option2_correct.npz
```

5. **Task 3準備（比較スクリプト作成）**

---

## 📚 参照ドキュメント

### 作成済み
- `shared/results/window_splitting_validation_report.md` ⭐ ウィンドウ分割検証レポート
- `docs/plans/sliding_window_validation_plan.md` - 検証計画

### 既存の重要資料
- `.claude/CLAUDE.md` - プロジェクト全体のガイドライン
- `docs/reports/PROJECT_COMPLETION.md` - Phase 1-6完成報告

---

**次セッションでの成功を祈ります！**

**重要**: ウィンドウ分割ロジックの完全一致が確認済みのため、本実行・結果比較に安心して進めます。
