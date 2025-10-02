# IHCP-CGM: 逆熱伝導問題 Python→Julia移植プロジェクト

[![Tests](https://img.shields.io/badge/tests-505%20passed-brightgreen)](julia/test/)
[![Phase](https://img.shields.io/badge/phase-6%2F6%20complete-blue)](.claude/CLAUDE.md)
[![Validation](https://img.shields.io/badge/validation-exact%20match-success)](shared/results/validation/FINAL_VERIFICATION_SUMMARY.md)

## プロジェクト概要

逆熱伝導問題（IHCP）を共役勾配法（CGM）で解くPythonコードをJuliaに移植したプロジェクトです。IRカメラからの温度測定データを用いてSUS304ステンレス鋼の表面熱流束を逆解析します。

**プロジェクト完成度: 100%** ✅

## 主要な達成事項

### ✅ 完全移植完了（Phase 1-6）
- **Phase 1**: 熱物性値計算（25テスト）
- **Phase 2**: DHCP直接ソルバー（298テスト）
- **Phase 3**: Adjoint随伴ソルバー（13テスト）
- **Phase 4**: CGM最適化（18テスト）
- **Phase 5**: スライディングウィンドウ計算（29テスト）
- **Phase 6**: 検証・実データ読込（122テスト）

### ✅ 完全一致検証成功
- Python版とJulia版で**数値的完全一致**を達成
- 熱流束相対誤差: **最大0.0041%** （基準0.1%を大幅クリア）
- 温度場相対誤差: **最大0.0047%** （基準0.1%を大幅クリア）
- 総テスト数: **505項目 全合格**

### 📊 検証結果詳細
詳細は以下のドキュメントを参照:
- [完全一致検証 最終結果サマリー](shared/results/validation/FINAL_VERIFICATION_SUMMARY.md)
- [完全一致検証 実行手順書](shared/results/validation/README_EXACT_MATCH.md)

## クイックスタート

### 必要環境
- **Python**: 3.12以上（オリジナルコード実行用）
- **Julia**: 1.9以上（移植版実行用）
- **メモリ**: 8GB以上推奨

### Python版実行（オリジナル）
```bash
cd python/original/
python IHCP_CGM_Sliding_Window_Calculation_ver2.py
```

### Julia版実行（移植版）
```bash
cd julia/

# 全テスト実行（505項目）
julia --project=. -e 'using Pkg; Pkg.test()'

# メイン計算実行
julia --project=. src/main.jl --config config/test.toml --output results.jld2
```

### 完全一致検証の再現
```bash
# 1. Python版実行
python python/validation/run_exact_match.py

# 2. Julia版実行
julia julia/examples/run_exact_match.jl

# 3. 結果比較
python python/validation/compare_exact_match.py
```

## ディレクトリ構造

```
TrialClaudeMCPCodex/
├── python/                # Python関連
│   ├── original/          # オリジナルPythonコード
│   └── validation/        # 検証スクリプト
├── julia/                 # Julia移植版
│   ├── src/               # ソースコード
│   ├── test/              # テストスイート（505項目）
│   ├── config/            # 設定ファイル
│   ├── data/              # 参照データ
│   └── examples/          # 使用例
├── shared/                # 共通リソース
│   ├── data/              # 共通データファイル
│   ├── config/            # 共通設定
│   └── results/           # 計算結果・検証データ
├── docs/                  # ドキュメント
│   ├── INDEX.md           # ドキュメント索引
│   ├── logs/              # 実装ログ（Phase 1-6）
│   ├── reports/           # 分析・検証レポート
│   ├── reviews/           # コードレビュー
│   └── plans/             # 計画書
└── .claude/               # Claude Code設定
    └── CLAUDE.md          # プロジェクト指示書
```

## 主要ドキュメント

### スタート地点
- [ドキュメント索引](docs/INDEX.md) - 全ドキュメントへのゲートウェイ
- [プロジェクト指示書](.claude/CLAUDE.md) - 開発方針・実装状況

### 技術ドキュメント
- [Julia移行計画](docs/plans/julia_migration_plan.md) - 14週間移植計画
- [Pythonコード構造レビュー](docs/reviews/python_code_structure.md)
- [Phase 1-6実装ログ](docs/logs/) - 時系列の実装記録

### 検証レポート
- [完全一致検証 最終サマリー](shared/results/validation/FINAL_VERIFICATION_SUMMARY.md) ⭐
- [Python/Julia比較検証](docs/reports/validation/python_julia_comparison_report.md)
- [CGM結果ゼロ問題修正](docs/logs/bug_fix_cgm_zero_result_20250210.md)

## アルゴリズム概要

本プロジェクトは以下の数値解法を実装しています:

1. **熱物性値計算**: 温度依存熱物性値の多項式フィッティング
2. **DHCP直接ソルバー**: 前進時間差分による熱伝導方程式の順解析
3. **Adjoint随伴ソルバー**: 後退時間差分による随伴方程式の解法
4. **CGM最適化**: 共役勾配法による熱流束の逆解析
5. **スライディングウィンドウ**: 長時間計算の効率化手法

### 数値計算設定
- **空間格子**: x,y方向均等格子（dx=0.12mm）、z方向20層非均等格子
- **時間刻み**: 1ms
- **境界条件**: 表面での熱流束境界条件
- **材料**: SUS304ステンレス鋼

## テスト体系

| Phase | テスト数 | 内容 |
|-------|---------|------|
| Phase 1 | 25 | 熱物性値計算の正確性 |
| Phase 2 | 298 | DHCP直接解法の精度・安定性 |
| Phase 3 | 13 | Adjoint随伴解法の収束性 |
| Phase 4 | 18 | CGM最適化の反復収束 |
| Phase 5 | 29 | スライディングウィンドウの連続性 |
| Phase 6 | 122 | 検証器・データ読込・統合機能 |
| **合計** | **505** | **全合格** ✅ |

## パフォーマンス

### 計算時間（356ステップ、window_size=71）
- **Python版**: 106秒（CGM: 92秒、検証: 14秒）
- **Julia版**: 1516秒（CGM: 1179秒、検証: 337秒）

**備考**: Julia版の性能最適化は今後の課題（並列化、型安定性改善）

### 数値精度
- 熱流束範囲: 3.27e+03 ~ 1.90e+06 W/m²
- 温度範囲: 300.0 ~ 452.9 K
- Python-Julia最大相対誤差: **< 0.005%**

## データファイル

### 必須データ
- `shared/data/metal_thermal_properties.csv` (540B): SUS304熱物性値データ
- `shared/data/T_measure_700um_1ms.npy` (1.1GB): 測定温度データ ※gitignore済み

### 検証結果データ
- `shared/results/validation/python_exact.npz` (402MB)
- `shared/results/validation/julia_exact.npz` (481MB)
- `shared/results/validation/exact_match_comparison_stats.json` (2.3KB)

## 今後の展開

プロジェクト完成後の推奨タスク:
- [ ] Julia版の性能最適化（並列化、メモリ最適化）
- [ ] フルスケールデータでの実運用検証
- [ ] GPU対応の検討
- [ ] 自動パラメータチューニング機能の追加

詳細は[今後のTODO](docs/plans/future_todos.md)を参照

## 技術スタック

### Python版
- NumPy 1.26.4 (数値計算)
- SciPy 1.13.1 (sparse行列、最適化)
- Pandas 2.2.2 (データ処理)
- Numba 0.60.0 (高速化、並列処理)

### Julia版
- Julia 1.9以上
- LinearAlgebra (標準ライブラリ)
- SparseArrays (標準ライブラリ)
- NPZ.jl (NumPyデータ読込)
- TOML.jl (設定ファイル)
- JLD2.jl (Julia形式データ保存)
- Plots.jl (可視化)

## ライセンス・貢献

このプロジェクトは研究目的で開発されました。

### 開発環境
- Claude Code + MCP Codex連携による開発
- Test-Driven Development (TDD) による段階的実装
- 完全な日本語ドキュメント化

### 貢献者
- 開発: Claude Code (Anthropic)
- 監修・検証: プロジェクトオーナー

## 問い合わせ

技術的な質問や問題報告は以下のドキュメントを参照:
- [完全一致検証 実行手順書](shared/results/validation/README_EXACT_MATCH.md)
- [プロジェクト指示書](.claude/CLAUDE.md)
- [ドキュメント索引](docs/INDEX.md)

---

**最終更新**: 2025年10月2日
**プロジェクト完成度**: 100% ✅
**総テスト数**: 505項目 全合格 ✅
