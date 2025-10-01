# 計算結果ディレクトリ

このディレクトリには、プロジェクトで生成された計算結果ファイルが格納されています。

## ディレクトリ構成

### ベンチマーク結果 (`benchmarks/`)

性能測定とベンチマークの結果データ。

- タイムステップ推定の計算結果
- 性能プロファイリングデータ
- メモリ使用量の測定結果

### 検証結果 (`validation/`)

Python実装とJulia実装の検証データ。

- **`python_test_1min.npz`** (22MB)
  - Python実装の1分間テスト結果
  - 温度分布、熱流束などの数値データ
  - NumPy圧縮形式

- **`python_summary.json`**
  - 検証結果の要約情報
  - 実行時間、収束性などの統計データ

- **`python_viz_data.json`**
  - 可視化用の抽出データ
  - グラフ描画用の時系列データ

## データファイルの取り扱い

### ファイルサイズ

- 大容量ファイル（1MB以上）は `.gitignore` に追加することを推奨
- 現在、`python_test_1min.npz` (22MB) はgit管理下にあります

### データ形式

- **NumPy形式 (`.npz`)**: 多次元配列データの圧縮保存
- **JSON形式 (`.json`)**: メタデータと小規模な数値データ

### データの再生成

これらの結果ファイルは、以下のコマンドで再生成できます：

```bash
# Python検証データの生成
cd python/original/
python IHCP_CGM_Sliding_Window_Calculation_ver2.py

# Julia検証データの生成（実装後）
cd julia/
julia --project=. src/main.jl
```

## 関連ドキュメント

- 検証レポート: [`../../docs/reports/validation/`](../../docs/reports/validation/)
- ベンチマークレポート: [`../../docs/reports/benchmarks/`](../../docs/reports/benchmarks/)
- 実装ログ: [`../../docs/logs/`](../../docs/logs/)
