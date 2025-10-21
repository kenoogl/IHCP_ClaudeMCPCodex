# ウィンドウ分割ロジック検証レポート

**日付**: 2025年10月21日
**ブランチ**: `sliding-window-validation`
**目的**: オリジナルコード・Python版・Julia版の3者間でウィンドウ分割ロジックの完全一致を検証

---

## 📋 検証対象

| 実装 | ファイル/関数 |
|------|-------------|
| オリジナル | `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`<br>`sliding_window_CGM_q_saving()` |
| Python版 | `python/scripts/run_sliding_window.py`<br>（オリジナルを呼び出し） |
| Julia版 | `julia/scripts/run_sliding_window.jl` |

---

## ✅ 検証結果サマリー

### 3者間ウィンドウ分割一致確認

| # | テストケース | nt | window | overlap | ウィンドウ数 | 3者一致 |
|---|-------------|-----|---------|---------|------------|---------|
| 1 | 小規模スライディング | 10 | 5 | 2 | 5 | ✅ |
| 2 | 実用規模 | 300 | 71 | 17 | 23 | ✅ |
| 3 | 単一ウィンドウ | 10 | 10 | 0 | 1 | ✅ |
| 4 | エッジケース | 20 | 10 | 5 | 8 | ✅ |

**結論**: **全4テストケースで3者完全一致** ✅

---

## 📊 詳細検証結果

### テストケース1: 小規模スライディング
**パラメータ**: `nt=10, window=5, overlap=2`

**出力（オリジナル・Python・Julia共通）**:
```
[Window 1] range=[0,5] length=5
[Window 2] range=[3,8] length=5
[Window 3] range=[6,9] length=3
[Window 4] range=[7,9] length=2
[Window 5] range=[8,9] length=1
Total windows: 5
```

**判定**: ✅ 完全一致

---

### テストケース2: 実用規模
**パラメータ**: `nt=300, window=71, overlap=17`

**出力（オリジナル・Python・Julia共通）**:
```
[Window 1] range=[0,71] length=71
[Window 2] range=[54,125] length=71
[Window 3] range=[108,179] length=71
[Window 4] range=[162,233] length=71
[Window 5] range=[216,287] length=71
[Window 6] range=[270,299] length=29
[Window 7] range=[282,299] length=17
[Window 8] range=[283,299] length=16
[Window 9] range=[284,299] length=15
[Window 10] range=[285,299] length=14
[Window 11] range=[286,299] length=13
[Window 12] range=[287,299] length=12
[Window 13] range=[288,299] length=11
[Window 14] range=[289,299] length=10
[Window 15] range=[290,299] length=9
[Window 16] range=[291,299] length=8
[Window 17] range=[292,299] length=7
[Window 18] range=[293,299] length=6
[Window 19] range=[294,299] length=5
[Window 20] range=[295,299] length=4
[Window 21] range=[296,299] length=3
[Window 22] range=[297,299] length=2
[Window 23] range=[298,299] length=1
Total windows: 23
```

**判定**: ✅ 完全一致

---

### テストケース3: 単一ウィンドウ
**パラメータ**: `nt=10, window=10, overlap=0`

**出力（オリジナル・Python・Julia共通）**:
```
[Window 1] range=[0,9] length=9
Total windows: 1
```

**判定**: ✅ 完全一致

---

### テストケース4: エッジケース
**パラメータ**: `nt=20, window=10, overlap=5`

**出力（オリジナル・Python・Julia共通）**:
```
[Window 1] range=[0,10] length=10
[Window 2] range=[5,15] length=10
[Window 3] range=[10,19] length=9
[Window 4] range=[14,19] length=5
[Window 5] range=[15,19] length=4
[Window 6] range=[16,19] length=3
[Window 7] range=[17,19] length=2
[Window 8] range=[18,19] length=1
Total windows: 8
```

**判定**: ✅ 完全一致

---

## 🔍 ウィンドウ分割ロジック比較

### コア計算式（3者共通）

```python
# 共通のウィンドウ分割アルゴリズム
start_idx = 0
while start_idx < (nt - 1):  # total_flux_steps
    max_L = min(window_size, (nt - 1) - start_idx)
    end_idx = start_idx + max_L
    # ウィンドウ処理
    step = max(1, max_L - overlap)
    start_idx += step
```

### ロジック要素の一致確認

| 要素 | オリジナル | Python版 | Julia版 | 一致 |
|------|-----------|---------|---------|------|
| ループ条件 | `start_idx < (nt - 1)` | `start_idx < total_flux_steps` | `start_idx < total_flux_steps` | ✅ |
| max_L計算 | `min(window_size, (nt-1) - start_idx)` | `min(args.window, total_flux_steps - start_idx)` | `min(window_size, total_flux_steps - start_idx)` | ✅ |
| end_idx計算 | `start_idx + max_L` | `start_idx + max_L` | `start_idx + max_L` | ✅ |
| step計算 | `max(1, max_L - overlap)` | `max(1, max_L - args.overlap)` | `max(1, max_L - overlap)` | ✅ |
| 更新処理 | `start_idx += step` | `start_idx += step` | `start_idx += step` | ✅ |

**注**: `total_flux_steps = nt - 1`（熱流束ステップ数）

---

## 🎯 結論

### 検証完了事項
1. ✅ **3者完全一致**: オリジナル・Python・Julia版のウィンドウ分割ロジックが同一
2. ✅ **4テストケース合格**: 小規模（5窓）から実用規模（23窓）まで検証
3. ✅ **エッジケース対応**: 単一ウィンドウ、大きいoverlap等も一致

### 重要な確認事項
- **Python版の実装構造**: `run_sliding_window.py`はラッパースクリプト
  - ウィンドウ分割の可視化・検証を提供
  - 実際の計算はオリジナルコードの`sliding_window_CGM_q_saving()`を呼び出し
  - ドライラン表示と実際の計算で同じロジックが使用される

### 次のステップへの準備
✅ **3者のウィンドウ分割ロジック完全一致が確認されたため、Task 2（本実行・結果比較）に進めます。**

---

## 📁 検証スクリプトと使用方法

### 1. オリジナルコード ドライラン機能テスト
**スクリプト**: `python/scripts/test_original_dryrun.py`

**目的**: オリジナルコードのドライラン機能が正しく動作することを確認

**使用方法**:
```bash
python python/scripts/test_original_dryrun.py
```

**期待される出力**:
```
=== Test Case 1: nt=10, window=5, overlap=2 ===
[Dry-run mode] Calculating window configuration...
  [Window 1] range=[0,5] length=5
  [Window 2] range=[3,8] length=5
  [Window 3] range=[6,9] length=3
  [Window 4] range=[7,9] length=2
  [Window 5] range=[8,9] length=1
  Total windows: 5
```

---

### 2. 3者自動比較スクリプト
**スクリプト**: `python/scripts/compare_window_splitting.py`

**目的**: オリジナル・Python・Juliaの3者間でウィンドウ分割を自動比較（4テストケース）

**使用方法**:
```bash
python python/scripts/compare_window_splitting.py
```

**期待される出力**:
```
================================================================================
Test: nt=10, window=5, overlap=2
================================================================================

オリジナルコード:
  [Window 1] range=[0,5] length=5
  ...
  Total windows: 5

run_sliding_window.py:
  [Window 1] range=[0,5] length=5
  ...

✅ ウィンドウ数一致: 5
✅ 全ウィンドウ詳細一致

================================================================================
✅ 全テストケース合格
================================================================================
```

---

### 3. 個別ドライラン実行

**Python版**:
```bash
python python/scripts/run_sliding_window.py \
  --nt 10 --window 5 --overlap 2 --dry-run
```

**Julia版**:
```bash
julia --project=julia julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --dry-run
```

**用途**: 任意のパラメータでウィンドウ分割を確認したい場合

---

**検証者**: Claude Code
**検証日時**: 2025年10月21日
**ステータス**: ✅ 全テストケース合格（4/4）
