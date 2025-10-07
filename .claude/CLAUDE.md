# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要
逆熱伝導問題（IHCP）をコンジュゲートグラディエント法（CGM）で解く数値計算プロジェクト。IRカメラからの温度測定データを使用してSUS304材の表面熱流束を逆解析する。

**プロジェクト状態**:
- ✅ Python→Julia移植完了（Phase 1-6、505テスト全合格）
- 🔄 **現在作業中**: マトリクスフリーPBICGSTAB!ソルバーの実装

**現在のブランチ**: tuning2（マトリクスフリーPBICGSTAB!実装中）
**最終更新**: 2025年10月8日

## 実行方法

### Python版（オリジナル）
```bash
cd python/original/
python IHCP_CGM_Sliding_Window_Calculation_ver2.py
```

### Julia版（移植版、Phase 1-6完了✅）
```bash
cd julia/

# 全テスト実行（505項目）
julia --project=. -e 'using Pkg; Pkg.test()'

# 個別Phaseテスト実行
julia --project=. test/test_thermal_properties.jl  # Phase 1
julia --project=. test/test_dhcp_solver.jl         # Phase 2
julia --project=. test/test_adjoint_solver.jl      # Phase 3
julia --project=. test/test_cgm_solver.jl          # Phase 4
julia --project=. test/test_sliding_window.jl      # Phase 5
julia --project=. test/test_validators.jl          # Phase 6

# メイン実行（設定ファイル指定）
julia --project=. src/main.jl --config config/test.toml --output results.jld2
julia --project=. src/main.jl --config config/example.toml --output results.jld2
```

### 完全一致検証の実行
```bash
# 1. Python版実行
python python/validation/run_exact_match.py

# 2. Julia版実行
julia julia/examples/run_exact_match.jl

# 3. 結果比較
python python/validation/compare_exact_match.py
```

### Julia実装状況

**完了Phase**: Phase 1-6（全6フェーズ完了）✅

| Phase | 内容 | テスト数 | 状態 |
|-------|------|---------|------|
| Phase 1 | 熱物性値計算 | 25 | ✅ |
| Phase 2 | DHCP直接ソルバー | 298 | ✅ |
| Phase 3 | Adjoint随伴ソルバー | 13 | ✅ |
| Phase 4 | CGM最適化 | 18 | ✅ |
| Phase 5 | スライディングウィンドウ | 29 | ✅ |
| Phase 6 | 検証・実データ読込・統合機能 | 122 | ✅ |
| **合計** | | **505** | **全合格** |

**Phase 6詳細**:
- A-1: 検証器関数群（NaN/Inf検出、範囲チェック、勾配検出）- 89テスト
- A-2: 結果保存・可視化（NPZ/JLD2形式）
- A-3: メインエントリポイント（main.jl、設定ファイル、CLI）
- C-1: MATLABファイル読込（IRカメラデータ）- 33テスト

### 必須データファイルの配置
- `shared/data/metal_thermal_properties.csv`: SUS304熱物性値データ（540B）
- `shared/data/T_measure_700um_1ms.npy`: 測定温度データ（1.1GB、gitignore済み）

### 必要な依存関係

**Python版**:
- Python 3.12以上
- NumPy 1.26.4 (数値計算)
- SciPy 1.13.1 (sparse行列、最適化)
- Pandas 2.2.2 (データ処理)
- Numba 0.60.0 (高速化、並列処理)

**Julia版**:
- Julia 1.10以上
- **コアパッケージ**: LinearAlgebra, SparseArrays
- **データI/O**: NPZ, JLD2, JSON, MAT, CSV, DataFrames
- **設定・可視化**: TOML, Plots, ArgParse
- **テスト**: Test, ProgressMeter
- 詳細は`julia/Project.toml`を参照

## コードアーキテクチャ

### Python版構成
- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`: メインプログラム（25KB）
- モノリシック構造、全機能を1ファイルに実装

### Julia版構成（モジュール化）
```
julia/src/
├── IHCP_CGM.jl              # メインモジュール（v0.1.0）
├── main.jl                  # エントリポイント
├── ThermalProperties.jl     # Phase 1: 熱物性値計算
├── DataLoaders.jl           # データ読み込み（CSV, NPY, MAT）
├── solvers/
│   ├── DHCPSolver.jl        # Phase 2: 直接解法
│   ├── AdjointSolver.jl     # Phase 3: 随伴解法
│   ├── StoppingCriteria.jl  # 停止判定ロジック
│   ├── CGMSolver.jl         # Phase 4: 共役勾配法
│   ├── CommonSolver.jl      # 🔄 共通ソルバー（PBICGSTAB!実装）
│   └── SlidingWindowSolver.jl # Phase 5: スライディングウィンドウ
└── utils/
    ├── validators.jl        # 数値検証器（Phase 6）
    ├── io.jl                # 入出力（NPZ, JLD2）
    ├── visualization.jl     # 可視化機能
    ├── json_helpers.jl      # JSON型変換
    └── boundary_conditions.jl # 境界条件処理
```

**重要なモジュール**:
- `CommonSolver.jl`: マトリクスフリー版PBICGSTAB!関数を実装
  - 疎行列の明示的構築を回避
  - 行列-ベクトル積を関数として実装
  - メモリ効率と計算速度の大幅改善

### 核心アルゴリズム
1. **熱物性値計算**: `thermal_properties_calculator()` - 温度依存熱物性値の多項式フィッティング
2. **直接解法（DHCP）**:
   - `coeffs_and_rhs_building_DHCP()` - 係数行列とRHS構築
   - `multiple_time_step_solver_DHCP()` - 前進時間差分解法
3. **随伴解法（Adjoint）**:
   - `coeffs_and_rhs_building_Adjoint()` - 随伴方程式の係数構築
   - `multiple_time_step_solver_Adjoint()` - 後退時間差分解法
4. **最適化**:
   - `global_CGM_time()` - コンジュゲートグラディエント法
   - `sliding_window_CGM_q_saving()` - スライディングウィンドウ計算

### 数値計算設定
- **空間格子**: x,y方向均等格子（dx=0.12mm）、z方向20層非均等格子
- **時間刻み**: 1ms
- **境界条件**: 表面での熱流束境界条件
- **材料**: SUS304ステンレス鋼

## 開発フロー
### Test-Driven Development (TDD)
- 段階的テスト（基本→DHCP→Adjoint→CGM→実データ）で信頼性確保
- テスト実行し、失敗を確認してから実装開始
- すべてのテストが通過するまで繰り返す

## コーディング規約

**言語設定**:
- **回答言語**: 日本語のみ

**共通**:
- **インデント**: スペース2文字
- **変数名**: 英語（物理・数学記法に準拠）
- **コメント**: 日本語で物理的意味を詳細に説明
- **エラー処理**: 数値計算の収束性とメモリ制限を考慮

**Python版**:
- **関数名**: snake_case形式
- **最適化**: Numba `@njit`デコレータで高速化
- **並列処理**: `prange`使用で性能向上
- **テスト**: pytest使用

**Julia版**:
- **関数名**: snake_case形式（Pythonと統一）
- **配列順序**: 常に`(ni, nj, nk)`形状を維持（Python互換）
- **インデックス**: 1始まり（Julia標準）
- **テスト**: `@testset`使用、段階的テスト（Phase 1-6）

## 重要な注意事項
- **大容量データ**: T_measure_700um_1ms.npy（1.1GB）の取り扱いに注意
- **計算時間**: フルスケール計算は数時間を要する場合がある
- **メモリ使用**: 大規模問題では8GB以上のメモリが必要
- **数値精度**: 逆問題の性質上、収束判定とステップサイズ調整が重要

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

## ドキュメント

### 主要ドキュメント索引
- [ドキュメント索引](docs/INDEX.md): 全ドキュメントへのゲートウェイ
- [完成度チェックリスト](docs/FINAL_CHECKLIST.md): プロジェクト完成度の詳細確認
- [プロジェクト完成報告](docs/PROJECT_COMPLETION.md): プロジェクト完了の公式記録

### 技術ドキュメント
- [Julia移行計画](docs/plans/julia_migration_plan.md): 14週間移植計画
- [Pythonコード構造レビュー](docs/reviews/python_code_structure.md)
- [Phase 1-6実装ログ](docs/logs/): 時系列の実装記録

## ディレクトリ構成

```
TrialClaudeMCPCodex/
├── python/                  # Python関連
│   ├── original/            # オリジナルPythonコード（25KB）
│   └── validation/          # 検証スクリプト
├── julia/                   # Julia移植版
│   ├── src/                 # ソースコード（詳細は「コードアーキテクチャ」参照）
│   ├── test/                # テストスイート（505項目、全合格）
│   ├── config/              # 設定ファイル（TOML）
│   ├── data/                # 参照データ生成スクリプト
│   ├── examples/            # 使用例（完全一致検証など）
│   └── Project.toml         # 依存管理
├── shared/                  # 共通リソース
│   ├── data/                # 共通データファイル
│   │   ├── metal_thermal_properties.csv (540B)
│   │   └── T_measure_700um_1ms.npy (1.1GB、gitignore済み)
│   ├── config/              # 共通設定（verification_params.toml）
│   └── results/             # 計算結果・検証データ（NPZ, JSON）
├── docs/                    # ドキュメント
│   ├── INDEX.md             # ドキュメント索引
│   ├── FINAL_CHECKLIST.md   # 完成度チェックリスト
│   ├── logs/                # Phase 1-6実装ログ
│   ├── reports/             # 検証レポート
│   ├── reviews/             # コードレビュー
│   └── plans/               # 計画書・TODO
└── .claude/
    └── CLAUDE.md            # プロジェクト指示書（本ファイル）
```

## Git操作ガイドライン

**現在のブランチ**: tuning2（マトリクスフリーPBICGSTAB!実装中）
**メインブランチ**: main（Phase 1-6完了版）

- **大容量ファイル**: 1MB以上のファイルはステージング前に必ず確認
- **除外ファイル**: `T_measure_700um_1ms.npy`（1.1GB）はgitignore済み
- **コミット単位**: TDDサイクルに従い、テスト作成→実装完了で適切に分割
- **ブランチ運用**: tuning2ブランチでマトリクスフリーPBICGSTAB!実装中、完了後mainにマージ

## トラブルシューティング

### データファイルが見つからない
- `shared/data/metal_thermal_properties.csv`と`T_measure_700um_1ms.npy`が必要
- `T_measure_700um_1ms.npy`は1.1GBのため別途取得が必要（gitignore済み）

### メモリ不足エラー
- フルスケール計算には8GB以上のメモリが必要
- 計算ウィンドウサイズを縮小して段階的に実行

### Numba関連エラー（Python）
- `@njit`デコレータ付き関数では型チェックが厳格
- NumPyの配列形状と型を事前に確認

### Julia配列インデックスエラー
- **重要**: 配列形状は常に`(ni, nj, nk)`とする
- Pythonと同じインデックス順序: `Temperature[i, j, k]`
- 0始まり（Python）→ 1始まり（Julia）のインデックス変換に注意

### Juliaテスト実行エラー
- 個別Phaseのテスト実行コマンドは「実行方法」セクションを参照
- 全テスト実行: `julia --project=. -e 'using Pkg; Pkg.test()'`
- テストが失敗した場合、該当Phaseのログ（`docs/logs/phase*/`）を確認

## Julia移植のベストプラクティス

### 配列順序の統一
- **Python**: `(ni, nj, nk)` 形状、アクセス: `arr[i, j, k]`
- **Julia**: `(ni, nj, nk)` 形状、アクセス: `arr[i, j, k]`（同じ）
- **理由**: Phase 1バグ修正で統一済み

### 多項式係数の扱い
- `polyval_numba(coeffs, x)`: `coeffs = [a, b, c, d]` → `a*x^3 + b*x^2 + c*x + d`
- NumPyの`polyfit`と同じ係数順序

### JSON参照データの使用
- `json_helpers.jl`を使用して型変換
- `reshape_json_array(data, dims, Float64)`でPython-Julia配列順序を変換

## 既知の問題と解決済み事項

### ✅ 解決済み
1. **Phase 1配列インデックスバグ**: `(nk, nj, ni)` → `(ni, nj, nk)`に修正
2. **Phase 2 JSON型変換問題**: `json_helpers.jl`で解決
3. **Phase 3随伴場の異常値**: 配列インデックス修正で解決（精度15桁改善）
4. **Phase 5完了**: スライディングウィンドウ計算実装（全29テスト合格）
5. **Phase 6完了**: 検証器、可視化、メインエントリポイント、実データ読込実装（全122テスト合格）
6. **CGM結果ゼロ問題**: 停止判定タイミング修正で解決（2025年10月2日）

## 完全一致検証結果 ✅

**実行日**: 2025年10月2日
**検証ステップ数**: 356ステップ（window_size=71）

**数値精度達成**:
- 熱流束相対誤差: **最大0.0041%** （基準0.1%を大幅クリア）✅
- 温度場相対誤差: **最大0.0047%** （基準0.1%を大幅クリア）✅
- 検証誤差（DHCP）: Python版と完全一致

**判定**: Python-Julia実用的完全一致を達成 ✅

検証手順と詳細結果は以下を参照:
- `shared/results/validation/FINAL_VERIFICATION_SUMMARY.md`
- `shared/results/validation/README_EXACT_MATCH.md`

## 性能改善タスク

### 現在のベンチマーク（356ステップ、window_size=71）
| 版 | 総時間 | CGM時間 | 検証時間 | 備考 |
|----|--------|---------|----------|------|
| Python | 106秒 | 92秒 | 14秒 | Numba最適化済み |
| Julia（旧版） | 1516秒 | 1179秒 | 337秒 | 行列版CG（14.3倍遅い） |

**改善目標**: Julia版をPython版と同等以下の実行時間に

### 🔄 現在進行中: マトリクスフリーソルバー実装

**目的**:
- 大規模疎行列の明示的構築を回避し、メモリ使用量を削減
- 行列-ベクトル積を関数として実装し、計算効率を向上

**実装方針**:
1. **CommonSolver.jl**に`PBICGSTAB!`関数を自前実装
   - Preconditioned BiConjugate Gradient Stabilized法
   - マトリクスフリー版（疎行列を構築せず、演算のみ実装）
2. DHCP/Adjointソルバーで疎行列構築を削減
   - 行列-ベクトル積関数を直接呼び出し
3. 既存テストとの整合性維持

**期待効果**:
- メモリ使用量: 50%以上削減
- 計算時間: 3-5倍高速化（目標: 300-500秒）
- 外部ライブラリ依存の削減

### 今後の改善候補
1. ✅ **マトリクスフリーPBICGSTAB!**: 実装中（tuning2ブランチ）
2. **並列化**: `@threads`マクロ、マルチスレッド処理
3. **型安定性**: `@code_warntype`による診断と修正
4. **GPU対応**: CUDA.jl/AMDGPU.jlの検討

### プロファイリングコマンド
```bash
# ベンチマーク実行
julia --project=. julia/examples/run_exact_match.jl

# メモリアロケーション追跡
julia --project=. --track-allocation=user src/main.jl --config config/test.toml
```

**重要**: 性能改善時は必ず数値精度を再検証し、完全一致を確認すること

詳細は`docs/plans/future_todos.md`を参照