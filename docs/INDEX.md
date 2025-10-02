# ドキュメント索引

このファイルは、プロジェクトの全ドキュメントへの索引です。

## 📁 ディレクトリ構成

```
docs/
├── logs/                    # 実装ログ（時系列）
│   ├── phase1/             # Phase 1: 基礎実装
│   ├── phase2/             # Phase 2: JSON処理とテスト
│   ├── phase3/             # Phase 3: バグ修正とテスト
│   ├── phase4/             # Phase 4: DHCP実装
│   ├── phase5/             # Phase 5: Adjoint実装
│   ├── phase6/             # Phase 6: 検証と最適化
│   └── project/            # プロジェクト全体のログ
│
├── reports/                # 分析・検証レポート
│   ├── benchmarks/         # ベンチマーク結果
│   └── validation/         # 検証レポート
│
├── reviews/                # コードレビュー
│
└── plans/                  # 計画書・TODO
```

---

## 📝 実装ログ

### Phase 1: 基礎実装
- [Phase 1 実装ログ](logs/phase1/implementation.md)

### Phase 2: JSON処理とテスト
- [Phase 2 実装ログ](logs/phase2/implementation.md)
- [Phase 2 JSON修正とテスト結果](logs/phase2/json_fix.md)

### Phase 3: バグ修正とテスト
- [Phase 3 実装ログ](logs/phase3/implementation.md)
- [Phase 3 バグ修正とテスト結果](logs/phase3/bug_fix.md)

### Phase 4: 直接熱伝導ソルバー（DHCP）実装
- [Phase 4 実装ログ](logs/phase4/implementation.md)

### Phase 5: 随伴ソルバー（Adjoint）実装
- [Phase 5 実装ログ](logs/phase5/implementation.md)
- [Phase 5 Codex分析](logs/phase5/codex_analysis.md)

### Phase 6: 検証と最適化（完了✅）
- [Phase 6 A1 Python検証データ分析](logs/phase6/a1_validation_analysis.md)
- [Phase 6 C1 MATファイル読み込み分析](logs/phase6/c1_mat_loading_analysis.md)

### プロジェクト全体
- [プロジェクト初期化ログ（2025年9月30日）](logs/project/initialization.md)
- [ドキュメント整理完了レポート（2025年10月2日）](logs/project/documentation_reorganization_report.md)
- [最終ドキュメント整理レポート（2025年10月2日）](logs/project/final_documentation_organization_report.md) ⭐

### バグ修正
- [CGM結果ゼロ問題修正（2025年10月2日）](logs/bug_fix_cgm_zero_result_20250210.md)

---

## 📊 レポート

### ベンチマーク
- [タイムステップ推定レポート（詳細）](reports/benchmarks/timestep_estimation_report.md)

### 検証
- [Python/Julia比較検証レポート](reports/validation/python_julia_comparison_report.md)
- [完全一致検証 最終結果サマリー](../shared/results/validation/FINAL_VERIFICATION_SUMMARY.md) ⭐
- [完全一致検証 詳細レポート](../shared/results/validation/exact_match_comparison_report.md)
- [完全一致検証 実行手順書](../shared/results/validation/README_EXACT_MATCH.md)

### その他
- [Phase A1 メインエントリーポイント実装レポート](reports/Phase_A1_Main_Entry_Point_Implementation_Report.md)
- [完全一致検証セットアップ完了レポート](reports/exact_match_verification_setup_complete.md)

---

## 🔍 レビュー

- [Pythonコード構造レビュー](reviews/python_code_structure.md)

---

## 📋 計画書・TODO

- [Julia移行計画](plans/julia_migration_plan.md)
- [今後のTODO](plans/future_todos.md)

---

## 🎓 プロジェクト完成ドキュメント

- [プロジェクト完成チェックリスト](FINAL_CHECKLIST.md) ⭐
- [プロジェクト完成宣言](PROJECT_COMPLETION.md) ⭐

---

## 💾 計算結果データ

計算結果ファイルは [`../shared/results/`](../shared/results/) に格納されています。

### 検証データ
- [`python_exact.npz`](../shared/results/validation/python_exact.npz) - Python完全一致検証結果（402MB）
- [`julia_exact.npz`](../shared/results/validation/julia_exact.npz) - Julia完全一致検証結果（481MB）
- [`exact_match_comparison_stats.json`](../shared/results/validation/exact_match_comparison_stats.json) - 比較統計データ
- [`comparison_plots.png`](../shared/results/validation/comparison_plots.png) - 比較可視化プロット
- [`python_test_1min.npz`](../shared/results/validation/python_test_1min.npz) - Python実装の1分間テスト結果（22MB）
- [`python_summary.json`](../shared/results/validation/python_summary.json) - 検証結果要約
- [`python_viz_data.json`](../shared/results/validation/python_viz_data.json) - 可視化用データ

---

## 📚 カテゴリ別索引

### 実装関連
- [Phase 1〜6 全実装ログ](#-実装ログ)
- [プロジェクト初期化](logs/project/initialization.md)
- [ドキュメント整理完了レポート](logs/project/documentation_reorganization_report.md)

### 分析・検証関連
- [ベンチマークレポート](#ベンチマーク)
- [検証レポート](#検証)
- [コードレビュー](#-レビュー)

### 計画関連
- [Julia移行計画](plans/julia_migration_plan.md)
- [今後のTODO](plans/future_todos.md)

---

## 🔗 関連READMEファイル

- [実装ログディレクトリ](logs/README.md)
- [レポートディレクトリ](reports/README.md)
- [ベンチマークレポート](reports/benchmarks/README.md)
- [検証レポート](reports/validation/README.md)
- [計算結果ディレクトリ](../shared/results/README.md)

---

## 📖 ドキュメント作成ルール

### 命名規則
- **実装ログ**: `implementation.md` - 各Phaseの主要な実装内容
- **修正ログ**: `*_fix.md` - 特定の問題の修正記録
- **分析ログ**: `*_analysis.md` - コード分析や検証の記録
- **レポート**: `*_report.md` - 分析結果のレポート
- **計画書**: `*_plan.md` - 実装計画や設計書

### ディレクトリルール
- **`logs/`** - 時系列の実装ログ（Phase別に整理）
- **`reports/`** - 分析結果とレポート（種類別に整理）
- **`reviews/`** - コードレビューとアーキテクチャ分析
- **`plans/`** - 計画書とTODOリスト

---

## 🎯 プロジェクト完成度

**Julia移植プロジェクト: 100%完了** ✅

### Phase 1-6 全完了
- ✅ Phase 1: 熱物性値計算（25テスト）
- ✅ Phase 2: DHCP直接ソルバー（298テスト）
- ✅ Phase 3: Adjoint随伴ソルバー（13テスト）
- ✅ Phase 4: CGM最適化（18テスト）
- ✅ Phase 5: スライディングウィンドウ（29テスト）
- ✅ Phase 6: 検証・実データ読込（122テスト）
  - A-1: 検証器関数群（89テスト）
  - A-2: 結果保存・可視化
  - A-3: メインエントリポイント
  - C-1: MATLABファイル読込（33テスト）

### テスト合格状況
**総テスト数: 505項目 全合格** ✅

### 完全一致検証結果
- ✅ Python-Julia完全一致検証成功
- ✅ 熱流束相対誤差: 最大0.0041% （基準0.1%を大幅クリア）
- ✅ 温度場相対誤差: 最大0.0047% （基準0.1%を大幅クリア）

詳細は[完全一致検証 最終結果サマリー](../shared/results/validation/FINAL_VERIFICATION_SUMMARY.md)を参照

---

最終更新: 2025年10月2日
