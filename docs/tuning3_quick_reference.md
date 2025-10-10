# tuning3作業クイックリファレンス

## 🚀 各セッション開始時の手順

```bash
# 1. ブランチとコミット位置確認
git branch  # tuning3であることを確認
git log --oneline -5

# 2. 計画書確認
cat docs/tuning3_recovery_plan.md | grep "^### Phase" -A 10

# 3. 進捗確認（チェックリストを見る）
cat docs/tuning3_recovery_plan.md | grep -A 100 "## 📊 進捗チェックリスト"
```

---

## 📦 Phase別コミット適用手順

### Phase A: メモリレイアウト最適化

#### コミット一覧
```bash
git log --oneline 99b321f..tuning2 | grep -E "3368d60|c8c99b1|2c00e6c|4c87fcf"
```

#### 適用方法（cherry-pick推奨）
```bash
# コミット3368d60: 配列軸順序統一
git cherry-pick 3368d60
# テスト実行
julia --project=. -e 'using Pkg; Pkg.test()'
python python/validation/run_10steps_fullsize_test.py
julia --project=. julia/scripts/run_10steps_fullsize_test.jl
python python/validation/compare_python_julia_10steps_fullsize.py

# 問題なければ次へ
git cherry-pick c8c99b1
# 同様にテスト...

git cherry-pick 2c00e6c
# テスト...

git cherry-pick 4c87fcf
# テスト...

# Phase A完了コミット
git add -A
git commit -m "[tuning3-PhaseA] メモリレイアウト最適化完了

- 配列軸順序を(ni,nj,nk,nt)に統一
- permutedims対応完了
- 検証: Python-Julia一致確認済み
"
```

---

### Phase B: ガイドセル統合

#### コミット一覧
```bash
git log --oneline 99b321f..tuning2 | grep -E "92d8a40|b1ec95f|233cfc7|5d993f9"
```

#### 適用方法
```bash
git cherry-pick 92d8a40
# テスト実行...

git cherry-pick b1ec95f
# テスト...

git cherry-pick 233cfc7
# テスト...

git cherry-pick 5d993f9
# テスト...

# Phase B完了コミット
git add -A
git commit -m "[tuning3-PhaseB] ガイドセル統合完了

- GridTransformモジュール導入
- 境界条件ユーティリティ追加
- 検証: Python-Julia一致確認済み
"
```

---

### Phase C: マトリクスフリー化

⚠️ **最も大きな変更 - 慎重に進める**

#### コミット一覧
```bash
git log --oneline 99b321f..tuning2 | grep -E "d9799c0|f72e9be|af36461|15719bb|f2e5a70|6af4798|ddf0c3d|4692ae4"
```

#### 適用方法（1コミットずつ慎重に）
```bash
# Step 1: API更新とThermalProperties移行
git cherry-pick d9799c0
julia --project=. -e 'using Pkg; Pkg.test()'
# 失敗した場合は git cherry-pick --abort でロールバック

# Step 2: CommonSolverにPBICGSTAB!実装
git cherry-pick f72e9be
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 3: DHCPから旧cg!撤去
git cherry-pick af36461
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 4: Adjointマトリクスフリー化
git cherry-pick 15719bb
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 5: CGM新API対応
git cherry-pick f2e5a70
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 6: Adjoint底面温度インデックス修正
git cherry-pick 6af4798
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 7: WorkBuffers初期化関数追加
git cherry-pick ddf0c3d
julia --project=. -e 'using Pkg; Pkg.test()'

# Step 8: SensitivitySolver新設
git cherry-pick 4692ae4
julia --project=. -e 'using Pkg; Pkg.test()'

# 全テスト + 10ステップテスト実行
python python/validation/run_10steps_fullsize_test.py
julia --project=. julia/scripts/run_10steps_fullsize_test.jl
python python/validation/compare_python_julia_10steps_fullsize.py

# Phase C完了コミット
git add -A
git commit -m "[tuning3-PhaseC] マトリクスフリー化完了

- PBICGSTAB!実装とWorkBuffers導入
- 全ソルバーマトリクスフリー化
- 検証: Python-Julia一致確認済み
- 性能: 89倍高速化達成
"
```

---

### Phase D: コード整理と安定化

#### コミット一覧
```bash
git log --oneline 99b321f..tuning2 | grep -E "21cbfff|9c68276|3e27f19|49df97e|108bc95|14145c8"
```

#### 適用方法
```bash
git cherry-pick 21cbfff  # FLoopsのlet束縛
julia --project=. -e 'using Pkg; Pkg.test()'

git cherry-pick 9c68276  # レビュー高優先度バグ修正
julia --project=. -e 'using Pkg; Pkg.test()'

git cherry-pick 3e27f19  # レビュー低優先度項目対応
julia --project=. -e 'using Pkg; Pkg.test()'

git cherry-pick 49df97e  # 共通RHSロジック抽出
julia --project=. -e 'using Pkg; Pkg.test()'

git cherry-pick 108bc95  # 未使用import削除
julia --project=. -e 'using Pkg; Pkg.test()'

git cherry-pick 14145c8  # CGMからSparseArrays依存除去
julia --project=. -e 'using Pkg; Pkg.test()'

# 最終検証
python python/validation/run_10steps_fullsize_test.py
julia --project=. julia/scripts/run_10steps_fullsize_test.jl
python python/validation/compare_python_julia_10steps_fullsize.py

# Phase D完了コミット
git add -A
git commit -m "[tuning3-PhaseD] コード整理と安定化完了

- 共通RHSロジック抽出
- 未使用コード削除
- レビュー指摘事項対応
- 検証: Python-Julia一致確認済み
"
```

---

## 🧪 検証コマンド集

### 全テスト実行（505項目）
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

### 10ステップフルサイズテスト
```bash
# Python版実行
python python/validation/run_10steps_fullsize_test.py

# Julia版実行
julia --project=. julia/scripts/run_10steps_fullsize_test.jl

# 結果比較
python python/validation/compare_python_julia_10steps_fullsize.py
```

### 結果確認
```bash
# 最大相対誤差確認（温度場）
cat <(python python/validation/compare_python_julia_10steps_fullsize.py) | grep "Max relative error (T)"

# プロット確認
open shared/results/comparison_10steps_fullsize_heat_flux.png
open shared/results/comparison_10steps_fullsize_temperature.png
open shared/results/comparison_10steps_fullsize_error_hist.png
```

---

## 🚨 トラブルシューティング

### cherry-pick でコンフリクト発生時
```bash
# 1. コンフリクトファイル確認
git status

# 2. 手動修正または中止
# 手動修正する場合:
# - エディタでコンフリクト解消
git add <修正したファイル>
git cherry-pick --continue

# 中止する場合:
git cherry-pick --abort
```

### テスト失敗時
```bash
# 1. 詳細エラーログ確認
julia --project=. -e 'using Pkg; Pkg.test()' 2>&1 | tee test_error.log

# 2. 特定のテストのみ実行
julia --project=. -e 'using Pkg; Pkg.test("IHCP_CGM", test_args=["phase2"])'

# 3. 前のコミットに戻る
git reset --hard HEAD~1
```

### Python-Julia不一致時
```bash
# 1. 詳細な差分確認
python python/validation/compare_python_julia_10steps_fullsize.py

# 2. 中間結果保存・確認
# - 該当コミットを特定
# - git bisectで二分探索

# 3. 問題コミットをスキップ
git cherry-pick --skip
```

---

## 📝 進捗記録方法

### 計画書更新
```bash
# vim/emacsでチェックボックス更新
vim docs/tuning3_recovery_plan.md

# チェックボックス: [ ] → [x]
# トラブルシューティング記録に追記
```

### 性能測定記録
```bash
# Julia実行時間取得
julia --project=. julia/scripts/run_10steps_fullsize_test.jl | grep "CGM elapsed"

# 計画書の「性能測定記録」セクションに追記
```

---

## 🎯 各Phase完了時のチェックリスト

- [ ] 該当コミット全て適用完了
- [ ] Phase 1-6全テスト合格（505項目）
- [ ] 10ステップフルサイズテスト実行
- [ ] Python-Julia一致確認（相対誤差 < 0.01%）
- [ ] 性能測定（Phase Cのみ）
- [ ] Phase完了コミット作成
- [ ] 計画書の進捗チェックリスト更新
- [ ] トラブルがあればトラブルシューティング記録更新
- [ ] 作業ログ更新

---

## 🔗 よく使うコマンド

```bash
# 現在のブランチとコミット確認
git log --oneline -10

# tuning2との差分確認
git diff tuning2

# 特定コミットの内容確認
git show <コミットハッシュ>

# コミット間のファイル変更確認
git diff 99b321f..tuning2 --name-only

# 計画書をMarkdownビューアで開く（macOS）
open -a "Typora" docs/tuning3_recovery_plan.md  # Typoraがある場合
# または
cat docs/tuning3_recovery_plan.md
```

---

## 📅 次回セッション開始テンプレート

```markdown
### YYYY-MM-DD - Phase X 作業開始

**開始時の状態**:
- ブランチ: tuning3
- 最新コミット: <git log --oneline -1>
- 前回完了Phase: <前回完了Phase>

**今回の作業内容**:
- Phase X のコミット適用
- 検証実行

**作業ログ**:
1. HH:MM - コミットXXXXXXX適用
2. HH:MM - テスト実行
3. HH:MM - ...

**結果**:
- テスト結果: 合格/不合格
- Python-Julia一致: ○/×
- 問題点: あれば記載

**次回作業**:
-
```

---

**📌 このファイルと計画書を常に手元に置いて作業してください！**
