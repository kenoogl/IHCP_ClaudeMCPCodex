# 検証レポート

このディレクトリには、実装の正確性を検証したレポートが格納されています。

## レポート一覧

### Python/Julia比較検証

- **`python_julia_comparison_report.md`**
  - Python実装とJulia実装の比較検証レポート
  - 数値計算結果の一致性確認
  - 性能比較
  - 実装の違いと注意点

## 検証の目的

これらの検証レポートは、以下を確認します：

1. **数値精度** - Python版とJulia版の計算結果が一致するか
2. **アルゴリズム正確性** - 物理モデルが正しく実装されているか
3. **境界条件処理** - 境界条件が適切に扱われているか
4. **収束性** - 反復計算が正しく収束するか

## 検証データ

検証に使用されたデータファイルは、以下に格納されています：

- **結果ファイル**: [`../../../shared/results/validation/`](../../../shared/results/validation/)
  - `python_test_1min.npz` - Python実装の検証結果
  - `python_summary.json` - 検証結果の要約
  - `python_viz_data.json` - 可視化用データ

## 関連ドキュメント

- 実装ログ: [`../../logs/`](../../logs/)
- ベンチマークレポート: [`../benchmarks/`](../benchmarks/)
- レポート索引: [`../README.md`](../README.md)
