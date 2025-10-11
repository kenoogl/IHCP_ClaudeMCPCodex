# Codex連携ワークフロー

**作成日**: 2025年10月12日
**目的**: codex MCPサーバーの不安定性に対応し、別ターミナルでのcodex直接操作を前提とした効率的な開発フローを確立

---

## 1. 役割分担

### Claude Code（このセッション）
- ✅ プロジェクト全体の進行管理
- ✅ 要件整理・タスク定義
- ✅ コードレビュー
- ✅ テスト実行・検証
- ✅ ドキュメント作成・更新
- ✅ git操作（コミット、プッシュ、PR作成）
- ✅ 簡単なコード修正（1-2ファイル程度）

### codex（別ターミナル操作）
- 🔧 大規模コード分析
- 🔧 複雑な最適化実装（3ファイル以上）
- 🔧 パターン検出・リファクタリング
- 🔧 アルゴリズム改善
- 🔧 性能改善提案の検証実装

---

## 2. 標準フロー

### Phase 1: 準備（Claude Code）

```
1. タスクブリーフィング作成
   場所: docs/tasks/task_<phase名>.md
   テンプレート: docs/templates/task_template.md

2. ブランチ作成
   命名: tuning<番号> (例: tuning5)

3. ユーザーに依頼
   「<ブランチ名>でcodexにタスク実行を依頼してください」
```

### Phase 2: 実行（ユーザー + codex）

```bash
# 別ターミナル
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git checkout <ブランチ名>
codex

# プロンプト例:
"docs/tasks/task_<phase名>.mdのタスクを実行してください"
```

### Phase 3: 報告（ユーザー → Claude Code）

**簡潔な報告形式**:

```
✅ 推奨: 「codex完了。レビュー依頼」
✅ 推奨: 「codex完了。変更: 4ファイル、レポート追加」
✅ 推奨: 「codex完了。エラーあり: <簡単な説明>」
```

### Phase 4: 検証（Claude Code）

```
1. 変更確認
   → git status, git diff

2. コードレビュー
   → 品質、API互換性、安全性

3. テスト実行
   → julia --project=julia julia/test/runtests.jl

4. ベンチマーク実行（必要に応じて）
   → julia --project=julia scripts/run_10steps_fullsize_test.jl

5. ドキュメント統合
   → performance_improvement_proposals.md更新
   → 実装レポート確認

6. コミット・プッシュ
   → git commit && git push
```

---

## 3. 並行作業による効率化

### codex作業中にClaude Codeが実施できること

- 📝 次のPhaseのタスクブリーフィング作成
- 📊 ベンチマークスクリプト準備
- 📚 関連ドキュメントの調査・整理
- 🔍 既存実装の分析レポート作成
- 📖 performance_improvement_proposals.mdの更新準備

### 時間効率の最大化

```
┌─────────────────────┬──────────────────────────┐
│ codex（別ターミナル）  │ Claude Code（このセッション）│
├─────────────────────┼──────────────────────────┤
│ 型安定化実装         │ ・次Phase準備            │
│ （15-30分）         │ ・ドキュメント整備        │
│                     │ ・ベンチマークスクリプト作成│
│                     │ ・既存コードの調査        │
└─────────────────────┴──────────────────────────┘
```

---

## 4. ブランチ戦略

### 命名規則

```
tuning<番号>: 性能改善関連
  例: tuning4 (Phase 0: アロケーション削減)
      tuning5 (Phase 1-A: 型安定化)
      tuning6 (Phase 1-B: 初期推定値改善)

feature/<機能名>: 新機能追加
  例: feature/adaptive-tolerance

bugfix/<修正内容>: バグ修正
  例: bugfix/memory-leak
```

### ブランチ作成タイミング

- ✅ **各Phase開始時**: codex作業前に必ず新ブランチ作成
- ✅ **Phase完了後**: レビュー・テスト完了後にコミット・プッシュ
- ✅ **マージ**: 必要に応じてmainにマージ（通常はPR経由）

---

## 5. タスクブリーフィングの書き方

### テンプレート

`docs/templates/task_template.md`参照

### 必須項目

1. **目的**: 何を達成するか
2. **対象ファイル**: 変更が必要なファイルリスト
3. **要件**: 満たすべき条件
4. **期待される出力**: codexが生成すべきもの
5. **参考資料**: 関連ドキュメント・コード
6. **制約**: 守るべき制約条件

### 良い例・悪い例

```markdown
✅ 良い例:
## 要件
1. @code_warntypeで警告がないこと
2. 動的ディスパッチゼロ
3. 既存テスト505件全通過

❌ 悪い例:
## 要件
1. 型を安定化する
```

---

## 6. コミュニケーションプロトコル

### ユーザーからの報告

**簡潔さを優先**:

```
✅ OK: 「codex完了。レビュー依頼」
✅ OK: 「codex完了。変更: CGMSolver.jl, 3ファイル」
✅ OK: 「codex完了。テスト失敗あり」

❌ NG: 詳細すぎる報告（Claude Codeが確認すればわかる内容）
```

### Claude Codeからの質問

**明確さを優先**:

```
✅ OK: 「codexに型安定化を依頼しますか？ブランチtuning5を作成します」
✅ OK: 「codexの作業完了を確認しますか？テストを実行します」

❌ NG: 「次は何をしますか？」（あいまい）
```

---

## 7. トラブルシューティング

### codex作業中にエラー発生

1. **ユーザー報告**: 「codex完了。エラーあり: <簡単な説明>」
2. **Claude Code対応**:
   - エラー内容確認（git diff, ログ確認）
   - 修正方法提案
   - 必要に応じてClaude Codeが直接修正

### テスト失敗

1. **Claude Code**: テスト実行結果を確認
2. **原因特定**: 失敗したテストの詳細確認
3. **対応方針**:
   - 簡単な修正: Claude Codeが実施
   - 複雑な修正: codexに再依頼（タスクブリーフィング更新）

### 性能改善効果が不十分

1. **ベンチマーク再実行**: 条件を変えて測定
2. **分析**: プロファイリング実施
3. **対応**:
   - 追加最適化検討
   - 別アプローチへの切り替え

---

## 8. ファイル構造

```
docs/
├── workflows/
│   └── codex_collaboration_workflow.md  # このファイル
├── templates/
│   ├── task_template.md                 # タスクブリーフィングテンプレート
│   └── completion_checklist.md          # 完了チェックリスト
├── tasks/
│   ├── task_phase0_allocation.md        # 完了済み（参考用）
│   ├── task_phase1_type_stability.md    # 次のタスク
│   └── ...
└── reports/
    └── benchmarks/
        └── ...
```

---

## 9. 実践例: Phase 1-A実施

### Step 1: 準備（Claude Code）

```markdown
1. タスク作成
   → docs/tasks/task_phase1a_type_stability.md

2. ブランチ作成
   → git checkout -b tuning5

3. 依頼
   → 「tuning5ブランチでcodexにタスク実行を依頼してください」
```

### Step 2: 実行（ユーザー + codex）

```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git checkout tuning5
codex
# タスク実行...
```

### Step 3: 報告（ユーザー）

```
「codex完了。変更: CommonSolver.jl, DHCPSolver.jl, 2ファイル」
```

### Step 4: 検証（Claude Code）

```
1. git diff で変更確認
2. テスト実行（505件）
3. 型安定性検証（@code_warntype）
4. ドキュメント更新
5. コミット・プッシュ
```

---

## 10. チェックリスト

### codex作業前

- [ ] タスクブリーフィング作成完了
- [ ] ブランチ作成完了
- [ ] 作業ディレクトリ確認（git status）

### codex作業後

- [ ] 変更ファイル確認（git status）
- [ ] コードレビュー実施
- [ ] テスト実行・合格
- [ ] ベンチマーク実施（必要に応じて）
- [ ] ドキュメント更新
- [ ] コミット・プッシュ

---

## 11. 改善ポイント

このワークフローは適時改善します。以下のような改善を検討：

- 📝 タスクブリーフィングのフォーマット改善
- 🔄 並行作業の最適化
- 📊 効率測定（各Phase所要時間記録）
- 🛠️ ツール導入（自動化スクリプトなど）

---

## 参考資料

- `docs/templates/task_template.md`: タスクブリーフィングテンプレート
- `docs/templates/completion_checklist.md`: 完了チェックリスト
- `docs/performance_improvement_proposals.md`: 性能改善提案書
