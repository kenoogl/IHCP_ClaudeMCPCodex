# ワークフローガイド

このディレクトリには、プロジェクト進行に関するワークフロー文書があります。

---

## 📚 主要ドキュメント

### [codex連携ワークフロー](./codex_collaboration_workflow.md)

**目的**: codexとClaude Codeの効率的な連携方法

**内容**:
- 役割分担の明確化
- 標準フロー（準備→実行→報告→検証）
- 並行作業による効率化
- ブランチ戦略
- タスクブリーフィングの書き方
- コミュニケーションプロトコル
- トラブルシューティング

**いつ見る**:
- 🚀 新しいPhase開始時
- 🔄 作業フロー確認時
- 🐛 問題発生時

---

### [セッション管理ガイド](./session_management.md) ⭐ 重要

**目的**: コンテキスト不足時(/compact・/clear)の対処法

**内容**:
- セッション中断が必要な状況の判断
- セッション終了前の必須作業
- SESSION_STATE.md更新方法
- セッション再開時の手順
- 効率的なセッション分割パターン
- 緊急時の対応

**いつ見る**:
- ⚠️ トークン使用量80%超
- 🔄 /compactや/clear実行前
- 🚀 セッション再開時

---

### [セッション状態](./SESSION_STATE.md) 📍 常に最新

**目的**: 現在の作業状態を記録（セッション再開時に必須）

**内容**:
- 完了した作業
- 進行中の作業
- 次のタスク
- 重要な情報（ベースライン性能、最新結果）
- 環境情報
- 参考ドキュメント

**更新タイミング**:
- Phase完了時（必須）
- 重要な決定・発見時
- セッション終了前

---

## 📋 テンプレート

### [タスクブリーフィングテンプレート](../templates/task_template.md)

codexにタスク依頼する際のテンプレート

**使い方**:
1. テンプレートをコピー
2. `docs/tasks/task_<phase名>.md`として保存
3. 各項目を埋める
4. codexに渡す

### [完了チェックリスト](../templates/completion_checklist.md)

作業完了時の確認項目リスト

**使い方**:
1. codex作業完了後に確認
2. Claude Codeが検証時に使用
3. 完了基準の確認

---

## 🚀 クイックスタート

### 新しいPhaseを開始する場合

1. **タスクブリーフィング作成**:
   ```bash
   # テンプレートをコピー
   cp docs/templates/task_template.md docs/tasks/task_phase1a_type_stability.md
   # 内容を編集
   ```

2. **ブランチ作成**:
   ```bash
   git checkout -b tuning5
   ```

3. **codexに依頼**（別ターミナル）:
   ```bash
   cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
   git checkout tuning5
   codex
   # プロンプト: "docs/tasks/task_phase1a_type_stability.mdのタスクを実行してください"
   ```

4. **完了報告**（ユーザー → Claude Code）:
   ```
   「codex完了。レビュー依頼」
   ```

5. **検証**（Claude Code）:
   - 変更確認
   - テスト実行
   - ドキュメント更新
   - コミット・プッシュ

---

## 📁 ディレクトリ構造

```
docs/
├── workflows/
│   ├── README.md                        # このファイル
│   └── codex_collaboration_workflow.md  # メインワークフロー
├── templates/
│   ├── task_template.md                 # タスクテンプレート
│   └── completion_checklist.md          # チェックリスト
├── tasks/
│   ├── task_phase0_allocation.md        # 完了済み（参考用）
│   ├── task_phase1a_type_stability.md   # 次のタスク例
│   └── ...
└── reports/
    └── benchmarks/
        └── ...
```

---

## 🔗 関連ドキュメント

- [性能改善提案書](../performance_improvement_proposals.md)
- [プロジェクト完成報告](../PROJECT_COMPLETION.md)
- [完成度チェックリスト](../FINAL_CHECKLIST.md)

---

## ❓ FAQ

### Q: codexとClaude Codeの使い分けは？

**A**:
- **codex**: 大規模コード分析、複雑な実装（3ファイル以上）
- **Claude Code**: 進行管理、レビュー、テスト、ドキュメント、git操作

### Q: タスクブリーフィングは毎回必要？

**A**: はい。codexに明確な指示を与えるために毎回作成します。テンプレートを使えば15分程度で作成できます。

### Q: ブランチ名の規則は？

**A**: `tuning<番号>`（性能改善）、`feature/<機能名>`（新機能）、`bugfix/<修正内容>`（バグ修正）

### Q: codex作業中にClaude Codeは何をする？

**A**: 次のPhase準備、ドキュメント整備、ベンチマークスクリプト作成など、並行作業で効率化

### Q: コンテキスト不足(/compact必要)になったら？

**A**:
1. 作業中の変更をコミット・プッシュ
2. SESSION_STATE.md更新
3. /compactまたは/clear実行
4. 再開時: 「SESSION_STATE.md確認してください」

### Q: セッション再開時に何を伝えればいい？

**A**: 「SESSION_STATE.md確認してください」だけでOK。必要な情報は全てSESSION_STATE.mdに記録されています。

### Q: SESSION_STATE.mdの更新を忘れたら？

**A**: git logとgit showで最新コミットを確認すれば状況を把握できます。詳細は[セッション管理ガイド](./session_management.md#緊急時の対応)参照。

---

## 🔄 改善履歴

- **2025-10-12**: 初版作成（Phase 0完了後）

---

**最終更新**: 2025年10月12日
