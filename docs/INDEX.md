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

### Phase 6: 検証と最適化
- [Phase 6 A1 Python検証データ分析](logs/phase6/a1_validation_analysis.md)
- [Phase 6 C1 MATファイル読み込み分析](logs/phase6/c1_mat_loading_analysis.md)

### プロジェクト全体
- [プロジェクト初期化ログ（2025年9月30日）](logs/project/initialization.md)

---

## 📊 レポート

### ベンチマーク
- [タイムステップ推定レポート（詳細）](reports/benchmarks/timestep_estimation_report.md)
- [タイムステップ推定サマリー](reports/benchmarks/timestep_estimation_summary.txt)

### 検証
- [Python/Julia比較検証レポート](reports/validation/python_julia_comparison_report.md)

### その他
- [Phase A1 メインエントリーポイント実装レポート](reports/Phase_A1_Main_Entry_Point_Implementation_Report.md)

---

## 🔍 レビュー

- [Pythonコード構造レビュー](reviews/python_code_structure.md)

---

## 📋 計画書・TODO

- [Julia移行計画](plans/julia_migration_plan.md)
- [今後のTODO](plans/future_todos.md)

---

## 💾 計算結果データ

計算結果ファイルは [`../shared/results/`](../shared/results/) に格納されています。

### 検証データ
- [`python_test_1min.npz`](../shared/results/validation/python_test_1min.npz) - Python実装の1分間テスト結果（22MB）
- [`python_summary.json`](../shared/results/validation/python_summary.json) - 検証結果要約
- [`python_viz_data.json`](../shared/results/validation/python_viz_data.json) - 可視化用データ

---

## 📚 カテゴリ別索引

### 実装関連
- [Phase 1〜6 全実装ログ](#-実装ログ)
- [プロジェクト初期化](logs/project/initialization.md)

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

最終更新: 2025年10月2日
