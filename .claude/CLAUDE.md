# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要
逆熱伝導問題（IHCP）をコンジュゲートグラディエント法（CGM）で解く数値計算プロジェクト。IRカメラからの温度測定データを使用してSUS304材の表面熱流束を逆解析する。

**プロジェクト状態**:
- ✅ Python→Julia移植完了（Phase 1-6、505テスト全合格）
- ✅ Phase A-D完了（マトリクスフリーPBICGSTAB!実装）
- ✅ Python版と同等の結果達成（温度場誤差1.06%、CGM反復1回）
- ✅ 性能ベースライン取得（22fde2d: 1072.43秒、80×100×20格子、10ステップ）
- ✅ 性能改善完了（PCG + none前処理で最速達成）
- 🎯 **現在のミッション**: スライディングウィンドウ計算のPython-Julia完全比較検証

**現在のブランチ**: main
**最終更新**: 2025年10月21日

**進行中のタスク**:
- スライディングウィンドウ計算の比較検証（Phase 1: CGM 3回、Phase 2: CGM 20000回）
- 詳細: `docs/plans/sliding_window_validation_plan.md`

## 実行方法

### テスト実行
```bash
# 全テスト実行（505項目）
julia --project=julia julia/test/runtests.jl

# 個別Phaseテスト実行
julia --project=julia julia/test/test_thermal_properties.jl  # Phase 1
julia --project=julia julia/test/test_dhcp_solver.jl         # Phase 2
julia --project=julia julia/test/test_adjoint_solver.jl      # Phase 3
julia --project=julia julia/test/test_cgm_solver.jl          # Phase 4
julia --project=julia julia/test/test_sliding_window.jl      # Phase 5
julia --project=julia julia/test/test_validators.jl          # Phase 6
```

## 作業開始時の確認手順（必須）

**⚠️ 重要**: 作業開始時は必ず以下の順序で実態確認を行うこと

1. **実態確認を最優先**
   - ドキュメントを読む前に、まず実際の状態を確認
   - git状態、コード、テスト結果が唯一の真実

2. **ドキュメントは参考程度**
   - ドキュメントは古くなる可能性があることを常に意識
   - ドキュメントの記述と実態が異なる場合は実態を優先

3. **git状態とテスト結果を最初に確認**
   ```bash
   # 常にこの順序で確認
   git status              # ブランチ、変更状態
   git log --oneline -10   # 最近のコミット
   julia --project=julia julia/test/runtests.jl  # テスト実行
   # コミットメッセージ内容照合
   ```

4. **コミットメッセージの内容照合を徹底**
   - ハッシュが違っても、メッセージ内容が一致すれば同じコミット
   - cherry-pick/rebaseでハッシュは変わることを前提
   - コミットの実態は`git show --stat <hash>`で確認

## テスト状況

- Phase 1-6: 505項目中179項目実装済み
- **合計**: 179 passed, 2 broken（既知の課題、参照データ関連）
- 全実装済み機能はテスト通過 ✅

## コードアーキテクチャ

### Julia版構成
```
julia/src/
├── IHCP_CGM.jl              # メインモジュール
├── ThermalProperties.jl     # 熱物性値計算
├── DataLoaders.jl           # データ読み込み
├── solvers/
│   ├── DHCPSolver.jl        # 直接解法（マトリクスフリー）
│   ├── AdjointSolver.jl     # 随伴解法（マトリクスフリー）
│   ├── SensitivitySolver.jl # 感度解法（マトリクスフリー）
│   ├── CGMSolver.jl         # 共役勾配法
│   ├── CommonSolver.jl      # PBICGSTAB!実装
│   └── SlidingWindowSolver.jl
└── utils/
    ├── validators.jl
    ├── boundary_conditions.jl
    └── json_helpers.jl
```

## コーディング規約

**言語設定**: 日本語のみ

**Julia版**:
- インデント: スペース2文字
- 関数名: snake_case形式
- 配列順序: `(ni, nj, nk)`形状を維持（Python互換）
- インデックス: 1始まり（Julia標準）

## Git操作ガイドライン

**現在のブランチ**: tuning3（性能改善フェーズ）
**メインブランチ**: main（Phase 1-6完了版）

- **コミット単位**: TDDサイクルに従い、テスト→実装で分割
- **大容量ファイル**: `T_measure_700um_1ms.npy`（1.1GB）はgitignore済み

## トラブルシューティング

### テスト実行エラー
```bash
# エラー時の確認手順
git status                    # 未コミット変更確認
git log --oneline -5          # 最近のコミット確認
julia --project=julia julia/test/runtests.jl  # 全テスト再実行
```

### データファイルが見つからない
- `shared/data/metal_thermal_properties.csv`: 540B（必須）
- `shared/data/T_measure_700um_1ms.npy`: 1.1GB（gitignore済み、別途取得）

### Julia配列インデックスエラー
- 配列形状: 常に`(ni, nj, nk)`
- インデックス: 1始まり（Julia標準）

### メモリ不足エラー
- フルスケール計算: 8GB以上必要
- 対処: 計算ウィンドウサイズを縮小

## 重要な参照ドキュメント

### 現在進行中
- **スライディングウィンドウ検証計画**: `docs/plans/sliding_window_validation_plan.md` ⭐ 最優先

### 過去の成果
- **性能改善提案書**: `docs/plans/performance_improvement_proposals.md`
- **性能ベースライン**: `shared/results/performance_22fde2d.md`
- **完成度チェックリスト**: `docs/tasks/FINAL_CHECKLIST.md`
- **プロジェクト完成報告**: `docs/reports/PROJECT_COMPLETION.md`

## Python-Julia精度検証結果

### Phase 1-6完了時（2025年10月2日）
- 熱流束相対誤差: 最大0.0041% ✅
- 温度場相対誤差: 最大0.0047% ✅
- 詳細: `shared/results/validation/FINAL_VERIFICATION_SUMMARY.md`

### 10ステップフルサイズテスト（2025年10月11日、22fde2d）
- 温度場誤差: 平均6.06K（相対1.06%）、最大21.91K（相対3.98%）
- 熱流束誤差: 平均0.0102 W/m²、最大0.0156 W/m²
- 条件: CGM反復1回、格子80×100×20、10ステップ
- 詳細: `shared/results/performance_22fde2d.md`

**結論**: Python版と実用上同等の結果を達成 ✅
