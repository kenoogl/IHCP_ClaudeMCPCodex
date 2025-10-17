# 次セッション作業ガイド - 体積積分変更後の性能検証

**作成日時**: 2025年10月17日 21:10
**ブランチ**: `tuning7-volume-integral-complete`
**最終コミット**: cf659ae (fix: Jacobi前処理のバグ修正と完全な前処理比較分析)

---

## 🚨 緊急課題：大規模格子での無限ループ問題

### 問題の概要

**症状**:
- `run_10steps_fullsize_test.jl`（80×100×20格子）が最初の時間ステップで無限ループ
- 1分以上経過しても反復が進まない（Phase 1-Eでは841秒で完了していた）

**動作するケース**: ✅
- `test_dhcp_convection.jl`（10×10×5格子）は正常動作
- 各時間ステップ2反復で約0.5ミリ秒で収束

**動作しないケース**: ❌
- `run_10steps_fullsize_test.jl`（80×100×20格子）
- 最初のステップで停止、反復結果が出力されない

### 体積積分変更の影響範囲

#### 主要な変更コミット

```
287c8c7 (2025-10-16) feat: 行列対称性検証とCommonSolver体積積分形式統一
142b94c (2025-10-16) fix: 体積積分形式の境界条件を修正（287c8c7のバグ修正）
27ea92c (2025-10-16) fix: Z方向格子生成の論理エラーを修正、Python版と完全一致を達成
9815b9d fix: RHS対流項をθ非依存に修正（体積積分形式Phase 2完了）
```

#### 変更されたファイル

1. **`julia/src/solvers/CommonSolver.jl`**: 体積積分形式の実装（主要変更）
2. **`julia/src/utils/RHSCore.jl`**: RHS計算の体積積分形式への変更
3. **`julia/src/solvers/DHCPSolver.jl`**: ✅ `test_dhcp_convection.jl`で動作確認済み
4. **🔴 Z方向格子生成**: `build_z_grid()`の実装変更（27ea92c）
   - ループ範囲: `for k in 2:(nk+1)` → `for k in 3:nk`
   - 底面セル（k=2）と表面セル（k=nk+1）を境界として明示的に設定

### Phase 1-E（体積積分変更前）の実績

**参照ファイル**: `shared/results/phase1e_adaptive.npz`

- **実行時間**: 841.20秒（約14分）
- **格子**: 80×100×20（同じサイズ）
- **ソルバー**: PBICGSTAB + GS前処理 + 適応的収束判定
- **RMS誤差**: 7.744 K
- **最大誤差**: 22.615 K
- **コミット**: d9e1b02（Phase 1-E実装完了、2025年10月14日）

---

## 🎯 次セッションの最優先タスク

### Task 1: 原因特定（Codex分析）

**分析対象ファイル（優先度順）**:

1. **🔴 Z方向格子の整合性**:
   ```
   julia/scripts/test_dhcp_convection.jl      # ✅ 動作する（格子生成方法）
   julia/scripts/run_10steps_fullsize_test.jl # ❌ 動作しない（格子生成方法）
   julia/src/main.jl                          # build_z_grid()実装
   ```

2. **体積積分実装でのZ方向格子情報の使用**:
   ```
   julia/src/solvers/CommonSolver.jl  # z_centers, ΔZ, dzの使用箇所
   julia/src/utils/RHSCore.jl         # z_centers, ΔZ, dzの使用箇所
   ```

**重点調査ポイント**:

1. **最優先**: Z方向格子定義（`z_centers`, `ΔZ`, `dz`）の整合性
   - `test_dhcp_convection.jl`と`run_10steps_fullsize_test.jl`の`build_z_grid()`実装の違い
   - 体積積分実装が想定している格子定義
   - 小規模格子（10×10×5）では動作し、大規模格子（80×100×20）で失敗する理由

2. 体積積分実装での境界条件の扱い:
   - `test_dhcp_convection.jl`: 対流境界（h=10.0）明示
   - `run_10steps_fullsize_test.jl`: 分布境界（Distribution）

3. 数値的安定性:
   - 大規模格子での精度問題
   - 反復行列の条件数

### Task 2: 問題の再現と診断

**実行コマンド**:
```bash
# 動作確認（小規模格子）
julia --project=julia julia/scripts/test_dhcp_convection.jl

# 問題の再現（大規模格子）- タイムアウト付き
timeout 120 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs
```

**診断スクリプト作成**:
- 中規模格子（40×50×10）でのテスト
- 格子サイズを徐々に大きくして閾値を特定
- デバッグ出力付きバージョンの作成

### Task 3: 修正と検証

**修正候補**:
1. Z方向格子定義の統一
2. 体積積分実装での格子情報の使い方の修正
3. 数値的安定性の改善（行列スケーリングなど）

**検証手順**:
1. `test_dhcp_convection.jl`で引き続き動作することを確認
2. 中規模格子でのテスト
3. `run_10steps_fullsize_test.jl`での動作確認
4. Phase 1-Eとの性能・精度比較

---

## 📝 本セッションの作業内容

### 完了した作業

1. **`run_10steps_fullsize_test.jl`の機能拡張**:
   - コマンドライン引数でソルバーと前処理を指定可能に
   - `--solver`: pbicgstab, pcg
   - `--precond`: jacobi, gs, none
   - `--help`オプション追加

2. **体積積分変更前の結果特定**:
   - Phase 1-E（適応的収束判定）の結果を確認
   - ファイル: `shared/results/phase1e_adaptive.npz`
   - 実行時間: 841.20秒（GS前処理使用）

3. **問題の特定と整理**:
   - 小規模格子（10×10×5）: ✅ 動作
   - 大規模格子（80×100×20）: ❌ 無限ループ
   - 体積積分変更とZ方向格子修正の影響範囲を特定

### 未完了の作業

1. 無限ループの根本原因特定
2. 体積積分変更前後の性能・精度比較
3. 修正実装とテスト

---

## 🔧 変更されたファイル（本セッション）

**主な変更**:
- `julia/scripts/run_10steps_fullsize_test.jl`: ソルバー・前処理指定機能追加

---

## 📊 参照データ

### Phase 1-E（体積積分変更前のベストタイム）

| 項目 | 値 |
|------|-----|
| ファイル | `shared/results/phase1e_adaptive.npz` |
| 実行時間 | 841.20秒（14.0分） |
| RMS誤差 | 7.744 K |
| 最大誤差 | 22.615 K |
| 格子 | 80×100×20 |
| ソルバー | PBICGSTAB + GS前処理 |
| 適応的収束判定 | ON |
| 体積積分 | 変更前 |

### Phase 1改善履歴

| Phase | 実行時間 | 改善率 | 累積改善 |
|-------|---------|--------|---------|
| ベースライン | 1072.43秒 | - | - |
| Phase 1-B | 890.47秒 | -16.56% | -16.56% |
| Phase 1-E | 841.20秒 | -6.37% | **-21.56%** |

---

## 🚀 次セッション開始時のチェックリスト

1. ✅ ブランチ確認: `tuning7-volume-integral-complete`
2. ✅ git status確認: 未コミット変更の確認
3. ✅ 最新コミット確認: `git log --oneline -5`
4. ⚠️ 問題の再現:
   ```bash
   # 動作するケース
   julia --project=julia julia/scripts/test_dhcp_convection.jl

   # 動作しないケース（タイムアウト付き）
   timeout 120 julia --project=julia julia/scripts/run_10steps_fullsize_test.jl --solver pbicgstab --precond gs
   ```
5. 🔍 Codex分析の準備:
   - `CommonSolver.jl`のZ方向格子使用箇所
   - `RHSCore.jl`のZ方向格子使用箇所
   - 2つのテストスクリプトの格子生成ロジックの違い

---

## 📞 Codex分析用クエリ例

```
体積積分変更後、大規模格子（80×100×20）でrun_10steps_fullsize_test.jlが
無限ループに陥っていますが、小規模格子（10×10×5）のtest_dhcp_convection.jlは
正常に動作します。

以下のファイルを分析して、Z方向格子情報（z_centers, ΔZ, dz）の使用方法に
問題がないか確認してください：

1. julia/src/solvers/CommonSolver.jl
2. julia/src/utils/RHSCore.jl
3. julia/scripts/test_dhcp_convection.jl (動作する)
4. julia/scripts/run_10steps_fullsize_test.jl (動作しない)

特に、格子サイズ依存の数値的問題や配列インデックスの問題に注目してください。
```

---

## 📚 関連ドキュメント

- `shared/results/performance_phase1e_adaptive_tolerance.md`: Phase 1-E詳細レポート
- `docs/performance_improvement_proposals.md`: 性能改善提案
- `.claude/CLAUDE.md`: プロジェクト設定とコーディング規約

---

**次セッションでの成功基準**:
1. 無限ループの原因を特定
2. 修正を実装
3. `run_10steps_fullsize_test.jl`が正常に完了（目標: 約840秒）
4. Phase 1-Eとの精度比較（RMS誤差が同程度であることを確認）
