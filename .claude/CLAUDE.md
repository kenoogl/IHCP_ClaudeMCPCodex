# CLAUDE.md

逆熱伝導問題（IHCP）をCGMで解くJulia実装プロジェクト。

## 状態
✅ Python→Julia移植完了（Phase 1-6、505テスト全合格）
✅ マトリクスフリーPBICGSTAB!実装（89倍高速化）
🔧 **現在作業中**: tuning3ブランチでtuning2修正の段階的適用とデバッグ

**ブランチ**: tuning3 | **ベースコミット**: 99b321f

## 実行

### Julia版
```bash
# 全テスト（505項目）
julia --project=. -e 'using Pkg; Pkg.test()'

# メイン実行
julia --project=. src/main.jl --config config/test.toml
```

### tuning3検証スクリプト（10ステップフルサイズテスト）
```bash
# Python版実行
python python/validation/run_10steps_fullsize_test.py

# Julia版実行
julia --project=. julia/scripts/run_10steps_fullsize_test.jl

# 結果比較
python python/validation/compare_python_julia_10steps_fullsize.py
```

### Phase構成（全6完了）
| Phase | 内容 | テスト |
|-------|------|-------|
| 1 | 熱物性値 | 25 |
| 2 | DHCP（マトリクスフリー版） | 298 |
| 3 | Adjoint | 13 |
| 4 | CGM | 18 |
| 5 | スライディングウィンドウ | 29 |
| 6 | 検証・実データ読込 | 122 |

### 依存
- Julia 1.10+
- コア: LinearAlgebra, FLoops
- I/O: NPZ, JLD2, MAT, CSV
- 詳細: `julia/Project.toml`

## アーキテクチャ

### 主要モジュール
```
julia/src/
├── solvers/
│   ├── CommonSolver.jl      # マトリクスフリーPBICGSTAB!
│   ├── DHCPSolver.jl        # 直接解法（311行、89倍高速）
│   ├── SensitivitySolver.jl # 感度問題
│   ├── AdjointSolver.jl     # 随伴解法
│   ├── CGMSolver.jl         # 共役勾配法
│   └── SlidingWindowSolver.jl
└── utils/
    ├── boundary_conditions.jl
    └── commons.jl           # WorkBuffers
```

### コアアルゴリズム
1. **DHCP**: `solve_dhcp!()` - 前進時間差分、マトリクスフリーPBICGSTAB!
2. **Sensitivity**: `solve_sensitivity!()` - 熱流束微小変化に対する温度応答
3. **Adjoint**: `solve_adjoint!()` - 後退時間差分、勾配計算
4. **CGM**: `solve_cgm!()` - 最適化メインループ

## 規約

- **言語**: 日本語のみ、コメント詳細記述
- **インデント**: 2スペース
- **配列**: `(ni, nj, nk)`順序、1-indexed
- **TDD**: テスト先行開発

## 検証結果

### 完全一致達成 ✅
- 熱流束相対誤差: **最大0.0041%**
- 温度場相対誤差: **最大0.0047%**
- Python完全一致（356ステップ検証）

### 性能
| 版 | 計算時間 | 高速化 |
|----|---------|--------|
| cg!版（旧） | 53.2秒 | 1x |
| **マトリクスフリー版** | **0.6秒** | **89x** |

## tuning3作業（現在進行中）

### 背景
- tuning2ブランチで計算結果不一致問題が発生
- コミット99b321f（安定版）から段階的に修正を適用して原因特定

### ベースライン検証結果（コミット99b321f）
- 温度場T最大相対誤差: **1.36e-05 (0.00136%)**
- 熱流束q最大絶対誤差: **1.98e-07 W/m²**
- ✅ Python-Julia完全一致確認済み

### 段階的適用計画（Phase A-D）
- **Phase A**: メモリレイアウト最適化
- **Phase B**: ガイドセル統合
- **Phase C**: マトリクスフリー化（最重要）
- **Phase D**: コード整理と安定化

### 作業手順
詳細は以下のドキュメントを参照：
- `docs/tuning3_recovery_plan.md`: 詳細計画・進捗チェックリスト
- `docs/tuning3_quick_reference.md`: 実務的な作業手順

## MCP連携とワークフロー

### Codex MCPとの連携方針
- **連携確認**: セッション開始時にCodex MCP接続状態を確認
- **役割分担**:
  - Claude Code: タスク計画・進捗管理・ユーザー対話・レビュー
  - Codex MCP: 詳細調査・コード実装・ファイル分析
- **作業ログ**: Codexへの指示と結果レポートを`docs/logs/`に記録

### コンテキスト管理
- セッション開始時に本CLAUDE.mdを参照して開発方針を確認
- コンテキスト圧縮後も本ファイルから方針を再確認
- **tuning3作業時**: 以下を参照
  - `docs/tuning3_recovery_plan.md`: 全体計画・進捗・トラブルシューティング
  - `docs/tuning3_quick_reference.md`: コマンド・手順クイックリファレンス

## 重要事項
- **データ**: `T_measure_700um_1ms.npy`（1.1GB、gitignore済み）必須
- **メモリ**: 8GB以上推奨
- **数値精度**: 性能改善時は必ず検証

## ドキュメント
- `docs/INDEX.md`: 全ドキュメント索引
- `docs/logs/`: Phase 1-6実装ログ
- `shared/results/validation/`: 検証結果
- **tuning3関連**:
  - `docs/tuning3_recovery_plan.md`: リカバリー計画
  - `docs/tuning3_quick_reference.md`: 作業手順

## 開発方針
- Task toolを用いCodex MCPと連携する
- Claudeの役割：タスク全体の計画と進捗管理、ユーザーとのコミュニケーション、Codexへの指示出しと結果のレビュー
- Codexの役割：詳細な調査作業、コードの実装、ファイル検索や分析
- **tuning3作業**: 各処理前に概要を説明し、ユーザー承認を得る
