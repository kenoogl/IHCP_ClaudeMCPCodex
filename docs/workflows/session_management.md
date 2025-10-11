# セッション管理ガイド

**目的**: コンテキスト不足による/compact・/clear時の対処法を確立

---

## 1. セッション中断が必要な状況

### よくあるケース

- ✅ トークン使用量が80%超（160k/200k tokens以上）
- ✅ /compactが必要なメッセージ表示
- ✅ レスポンスが遅くなってきた
- ✅ 複雑なタスクで大量のコード読み込み後
- ✅ 長時間作業（2時間以上）

### 対処タイミング

**推奨**: 各Phase完了時（コミット・プッシュ後）にセッションをクリアまたは/compact

---

## 2. セッション終了前の必須作業

### チェックリスト

```markdown
セッション終了前チェックリスト：

□ 作業中の変更をコミット・プッシュ
□ 作業状態を記録（SESSION_STATE.md更新）
□ 未完了タスクを記録（TODO.md更新）
□ 重要な発見・決定事項を記録
□ 次回セッションで必要な情報を確認
```

---

## 3. 状態保存ファイル

### 3-1. SESSION_STATE.md（必須）

**場所**: `docs/workflows/SESSION_STATE.md`

**目的**: セッション再開時に必要な全情報を集約

**内容**:
```markdown
# セッション状態

**最終更新**: 2025-10-12 14:30
**現在のブランチ**: tuning4
**最新コミット**: 2300119

## 現在の作業状態

### 完了した作業
- Phase 0: 配列アロケーション削減（完了、プッシュ済み）
- ワークフロー文書作成（完了、プッシュ済み）

### 進行中の作業
- なし

### 次のタスク
- Phase 1-A: 型安定化
  - ブランチ: tuning5（未作成）
  - タスクブリーフィング: 未作成

## 重要な情報

### ベースライン性能
- コミット: 22fde2d
- 実行時間: 1072.43秒

### Phase 0結果
- 実行時間: 1068.74秒
- 改善: 0.34%
- 採用判断: SlidingWindow実行での効果期待

### 次Phase準備
- performance_improvement_proposals.md参照
- 提案2（型安定化）: 10-20%改善見込み

## 環境情報

### リポジトリ状態
- ブランチ: tuning4（プッシュ済み）
- 未コミット変更: なし
- テスト状態: 505件全通過

### パッケージ
- NPZ: インストール済み
- MAT: インストール済み

## 参考ドキュメント

必ず確認すべきファイル：
1. docs/workflows/codex_collaboration_workflow.md
2. docs/performance_improvement_proposals.md
3. shared/results/performance_tuning4_allocation_reduction.md
```

### 3-2. TODO.md（オプション）

**場所**: プロジェクトルート `TODO.md`

**目的**: 未完了タスクの追跡

---

## 4. セッション再開時の手順

### Step 1: 状態確認（5分）

```bash
# 1. リポジトリ状態確認
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5

# 2. ブランチ確認
git branch

# 3. 最新変更確認
git show --stat
```

### Step 2: コンテキスト復元（Claude Codeに依頼）

**ユーザーからの依頼文（テンプレート）**:

```
セッション再開します。以下を確認してください：
1. docs/workflows/SESSION_STATE.mdを読み込み
2. 現在の作業状態を把握
3. 次のアクションを提案
```

### Step 3: 情報確認（Claude Code実施）

```markdown
□ SESSION_STATE.mdを読み込み
□ git statusで実態確認
□ 最新コミットメッセージ確認
□ performance_improvement_proposals.md確認
□ 次のタスク特定
```

### Step 4: 作業再開

Claude Codeから提案：
- 「現在tuning4ブランチ、Phase 0完了済み」
- 「次はPhase 1-A（型安定化）です。タスクブリーフィングを作成しますか？」

---

## 5. 効率的なセッション管理

### セッション分割の推奨パターン

#### Pattern 1: Phase単位

```
┌─────────────────────────────────┐
│ Session 1: Phase 0              │
│ - タスクブリーフィング作成      │
│ - codex実行（別ターミナル）     │
│ - レビュー・テスト              │
│ - コミット・プッシュ            │
│ - SESSION_STATE.md更新          │
└─────────────────────────────────┘
        ↓ /compact or /clear
┌─────────────────────────────────┐
│ Session 2: Phase 1-A            │
│ - SESSION_STATE.md読み込み      │
│ - タスクブリーフィング作成      │
│ - ...                          │
└─────────────────────────────────┘
```

#### Pattern 2: 作業単位

```
Session 1: 準備フェーズ
  - タスクブリーフィング作成
  - ブランチ作成
  → SESSION_STATE.md更新

Session 2: 実装・検証フェーズ（codex実行後）
  - SESSION_STATE.md読み込み
  - レビュー・テスト
  - ドキュメント統合
  - コミット・プッシュ
```

### トークン使用量モニタリング

**目安**:
- 50k tokens: 余裕あり
- 100k tokens: 注意
- 150k tokens: /compact検討
- 180k tokens: /clear推奨

**確認方法**: セッション中に「/context」コマンド

---

## 6. SESSION_STATE.md更新タイミング

### 必須更新タイミング

1. **Phase完了時**（最重要）:
   - コミット・プッシュ後
   - 次Phaseの準備状況記録

2. **重要な決定・発見時**:
   - 性能測定結果
   - 採用判断
   - 設計変更

3. **セッション終了前**:
   - 作業途中でも状態記録

### オプション更新

- 長時間作業の中間地点
- 大きなマイルストーン達成時

---

## 7. セッション再開時のクイックリファレンス

### 最小限の情報で再開

```bash
# ターミナルで実行
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
cat docs/workflows/SESSION_STATE.md | head -50
git status
git log --oneline -3
```

### Claude Codeへの依頼（最短版）

```
セッション再開。SESSION_STATE.md確認して次のアクション提案してください。
```

---

## 8. 緊急時の対応

### SESSION_STATE.mdがない場合

```bash
# 1. 最新コミットメッセージから状況把握
git log --oneline -10

# 2. 最新コミットの詳細確認
git show

# 3. performance_improvement_proposals.md確認
cat docs/performance_improvement_proposals.md | grep -A 5 "実装状況サマリー"
```

### ブランチが不明な場合

```bash
# 全ブランチ確認
git branch -a

# 最新の作業ブランチ特定
git for-each-ref --sort=-committerdate refs/heads/ --format='%(refname:short) %(committerdate:relative)'
```

---

## 9. テンプレート

### SESSION_STATE.md初期テンプレート

```markdown
# セッション状態

**最終更新**: [日時]
**現在のブランチ**: [ブランチ名]
**最新コミット**: [ハッシュ]

## 現在の作業状態

### 完了した作業
- [Phase名]: [概要]（完了、プッシュ済み）

### 進行中の作業
- [タスク名]: [状態]

### 次のタスク
- [Phase名]
  - ブランチ: [名前]
  - タスクブリーフィング: [状態]

## 重要な情報

### ベースライン性能
- コミット: [ハッシュ]
- 実行時間: [秒]

### 最新結果
- [Phase名]
- 実行時間: [秒]
- 改善: [%]

## 環境情報

### リポジトリ状態
- ブランチ: [名前]
- 未コミット変更: [あり/なし]
- テスト状態: [結果]

## 参考ドキュメント

1. [ドキュメント名]: [パス]
```

---

## 10. ベストプラクティス

### DO（推奨）

✅ **Phase完了毎にセッション区切り**
  - コミット・プッシュ後に/compactまたは/clear

✅ **SESSION_STATE.md更新を習慣化**
  - 重要な作業後は必ず更新

✅ **簡潔な再開依頼**
  - 「SESSION_STATE.md確認」だけで十分

✅ **git statusで実態確認**
  - ドキュメントより実態を優先

### DON'T（避ける）

❌ **長時間セッション継続**
  - トークン不足でコンテキスト喪失リスク

❌ **状態記録なしでセッション終了**
  - 再開時に混乱

❌ **詳細すぎる状態記録**
  - 本質的情報のみ記録

---

## 11. トラブルシューティング

### Q: SESSION_STATE.mdとgit statusが矛盾

**A**: git statusの実態を優先。SESSION_STATE.md更新忘れの可能性。

### Q: セッション再開時に混乱

**A**: 以下の順序で確認：
1. `git log --oneline -5`（最近の作業）
2. `git show`（最新コミットの詳細）
3. `docs/performance_improvement_proposals.md`（全体状況）

### Q: コミット前にセッション終了してしまった

**A**: 次回セッション開始時：
1. `git status`で未コミット変更確認
2. `git diff`で変更内容確認
3. 変更内容をレビューしてコミット判断

---

## 12. 実践例

### 例1: Phase完了後のセッション終了

```markdown
【Session 1終了前】

1. コミット・プッシュ完了
   → git push origin tuning4

2. SESSION_STATE.md更新
   - 完了した作業: Phase 0記録
   - 次のタスク: Phase 1-A記録

3. セッション終了
   → /compact または /clear

【Session 2開始時】

ユーザー: 「SESSION_STATE.md確認してください」

Claude Code:
  1. SESSION_STATE.md読み込み
  2. git status確認
  3. 「Phase 0完了済み。次はPhase 1-A（型安定化）です。
     タスクブリーフィングを作成しますか？」
```

### 例2: 作業途中でのセッション終了

```markdown
【作業途中】

状況: codex実行中、まだコミットしていない

1. 作業状態をSESSION_STATE.mdに記録
   - 進行中の作業: Phase 1-A実装中（codex実行中）
   - 状態: レビュー待ち

2. /clear

【再開時】

ユーザー: 「codex完了しました。SESSION_STATE.md確認してください」

Claude Code:
  1. SESSION_STATE.md確認
  2. 「Phase 1-A実装完了を確認しました。レビューを開始します」
  3. git diff確認
  4. テスト実行
```

---

## 13. まとめ

### セッション管理の3原則

1. **記録する**: SESSION_STATE.mdに作業状態を記録
2. **区切る**: Phase完了毎にセッションクリア
3. **確認する**: 再開時は実態（git status）を優先

### 最小限の手順

```
【終了前】
1. コミット・プッシュ
2. SESSION_STATE.md更新
3. /compact or /clear

【再開時】
1. "SESSION_STATE.md確認してください"
2. 次のアクション提案を受ける
3. 作業再開
```

---

**最終更新**: 2025年10月12日
