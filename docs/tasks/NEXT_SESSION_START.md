# 次セッション開始ガイド

**作成日**: 2025-10-17
**ブランチ**: `tuning7-volume-integral-complete`
**最終コミット**: 4420e41 "feat: Phase 4完了 - HC配列生成と伝播実装"

---

## 🎯 現在の状況

### ✅ Phase 4完了
体積積分形式の対流境界条件実装が完了しました。

**実装内容**:
1. `set_BC_coef`関数でHC配列生成（BoundaryConditions.jl）
2. PBiCGSTAB!、PCG!、solve_linear_system!にHC引数追加
3. Preconditioner!の全バリエーション（:none, :gs, :jacobi）にHC引数追加
4. DHCPSolver、AdjointSolver、SensitivitySolverでHC配列を使用

**テスト結果**:
- ✅ コンパイルテスト成功
- ✅ 全テスト実行: 179 passed, 2 broken（既知の課題）

**コミット状況**:
- ✅ Phase 4の変更をコミット完了（4420e41）
- ✅ リモートへプッシュ完了

---

## 🚀 次セッション開始手順

### Step 1: ブランチ確認
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5
```

**期待結果**:
```
On branch tuning7-volume-integral-complete
Your branch is up to date with 'origin/tuning7-volume-integral-complete'.

4420e41 feat: Phase 4完了 - HC配列生成と伝播実装
344078e docs: セッション継続情報を記録（Phase 2-3完了、Phase 4準備）
4d65584 feat: CommonSolver.jlに対流項h·Aを実装（Phase 3完了）
...
```

### Step 2: 継続情報の確認
```bash
# 詳細な作業内容を確認
cat TODO_NEXT_SESSION.md

# 実装計画を確認
cat docs/VOLUME_INTEGRAL_IMPLEMENTATION_PLAN.md
```

### Step 3: テスト実行（動作確認）
```bash
# コンパイルテスト
julia --project=julia -e 'using Pkg; Pkg.instantiate(); using IHCP_CGM; println("✓ コンパイル成功")'

# 全テスト実行
julia --project=julia julia/test/runtests.jl
```

**期待結果**: 179 passed, 2 broken

---

## 📝 Phase 5の作業内容

### 目的
対流境界条件のスケーリング調整とテスト

### タスク一覧

#### 5.1 スケーリング戦略の決定
- [ ] 係数の相対的な大きさを確認
- [ ] 条件数への影響を評価
- [ ] 必要に応じてスケーリング係数を導入

**検討ポイント**:
- 現在の実装: `HC[i] * area * (1 - m[隣接])`
- 熱伝達係数h: 10〜1000 W/m²·K
- 面積area: 10^-8〜10^-6 m²
- 熱伝導項: λ/dx（約10^-2〜10^2）

#### 5.2 数値テスト
- [ ] 対流境界なしケース（既存テストで確認済み）
- [ ] 対流境界ありケース（新規）
  - X-minus面に対流境界条件を設定
  - h = 10 W/m²·K, T_∞ = 300 K
- [ ] スケール感度テスト
  - h = 1, 10, 100, 1000 W/m²·Kで実行
  - 収束性と精度を評価

#### 5.3 必要に応じた修正
- [ ] スケーリング係数の導入
- [ ] 前処理器の調整
- [ ] 収束判定基準の見直し

---

## 🔧 開発環境

### 必要なツール
- Julia 1.9以降
- Git
- テキストエディタ（VS Code推奨）

### プロジェクト構造
```
TrialClaudeMCPCodex/
├── julia/
│   ├── src/
│   │   ├── solvers/
│   │   │   ├── CommonSolver.jl       # Phase 4で修正
│   │   │   ├── DHCPSolver.jl         # Phase 4で修正
│   │   │   ├── AdjointSolver.jl      # Phase 4で修正
│   │   │   └── SensitivitySolver.jl  # Phase 4で修正
│   │   └── utils/
│   │       └── boundary_conditions.jl # Phase 4で修正
│   └── test/
│       └── runtests.jl
├── docs/
│   └── VOLUME_INTEGRAL_IMPLEMENTATION_PLAN.md
├── TODO_NEXT_SESSION.md      # 詳細な作業ガイド
└── NEXT_SESSION_START.md     # このファイル
```

---

## 📚 重要な参照資料

### 実装計画書
- `docs/VOLUME_INTEGRAL_IMPLEMENTATION_PLAN.md`: 全体計画と実装詳細

### 継続作業ガイド
- `TODO_NEXT_SESSION.md`: Phase 5の詳細タスクと参照情報

### テスト
- `julia/test/runtests.jl`: 全テスト実行
- 179 passed, 2 broken（既知の課題）

### Git履歴
```bash
# Phase 1-4のコミット履歴を確認
git log --oneline --grep="Phase" -10

# 変更ファイルの確認
git show --stat 4420e41
```

---

## ⚠️ 注意事項

### 既知の課題
1. **スケーリング未検証**: Phase 5で対流境界条件のスケーリングを確認必要
2. **対流境界テスト未作成**: 実際の対流境界条件を使ったテストケースが未実装
3. **条件数への影響**: 熱伝達係数hの値によって条件数が変化する可能性

### コーディング規約
- インデント: スペース2文字
- 関数名: snake_case形式
- 配列順序: `(ni, nj, nk)`形状を維持
- コメント: 日本語で詳細に記述

---

## 💡 Claudeへの指示例

次セッション開始時に以下のように指示してください：

```
TODO_NEXT_SESSION.mdを確認して、Phase 5の作業を開始してください。
まず、対流境界条件のスケーリングテストケースを作成します。
```

または：

```
Phase 5のスケーリング戦略を決定するため、
現在の実装の係数の大きさを確認してください。
```

---

## 🔗 関連リンク

- **GitHub**: `tuning7-volume-integral-complete`ブランチ
- **最新コミット**: 4420e41
- **プロジェクトルート**: `/Users/Daily/Development/IHCP/TrialClaudeMCPCodex`

---

**End of Document**
