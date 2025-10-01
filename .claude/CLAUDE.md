# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要
逆熱伝導問題（IHCP）をコンジュゲートグラディエント法（CGM）で解く数値計算プロジェクト。IRカメラからの温度測定データを使用してSUS304材の表面熱流束を逆解析する。

## 実行方法

### Python版（オリジナル）
```bash
cd python/original/
python IHCP_CGM_Sliding_Window_Calculation_ver2.py
```

### Julia版（移植版、Phase 1-5実装完了）
```bash
cd julia/
julia --project=. -e 'using Pkg; Pkg.test()'  # 全テスト実行
```

**Julia実装状況**:
- ✅ Phase 1: 熱物性値計算（全25テストパス）
- ✅ Phase 2: DHCP直接ソルバー（全298テストパス）
- ✅ Phase 3: Adjoint随伴ソルバー（全13テストパス）
- ✅ Phase 4: CGM最適化（18テストパス）
- 🔄 Phase 5: スライディングウィンドウ（実装完了、テスト実行待ち）

### 必須データファイルの配置
- `shared/data/metal_thermal_properties.csv`: SUS304熱物性値データ（540B）
- `shared/data/T_measure_700um_1ms.npy`: 測定温度データ（1.1GB、gitignore済み）

### 必要な依存関係
- Python 3.x
- NumPy (数値計算)
- SciPy (sparse行列、最適化)
- Pandas (データ処理)
- Numba (高速化、並列処理)

## コードアーキテクチャ
### 主要ファイル構成
- `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`: メインプログラム（25KB）
- `shared/data/metal_thermal_properties.csv`: SUS304熱物性値データ
- `shared/data/T_measure_700um_1ms.npy`: 測定温度データ（1.1GB）

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

## 言語設定
- **回答言語**: 日本語のみ
- **コメント**: 日本語で物理的意味を詳細に説明
- **変数名**: 英語（物理・数学記法に準拠）
- **関数名**: 英語（snake_case形式）

## コーディング規約
- **インデント**: スペース2文字（TDD方針に準拠）
- **最適化**: Numba @njitデコレータで高速化
- **並列処理**: prange使用で性能向上
- **エラー処理**: 数値計算の収束性とメモリ制限を考慮
- **テスト**: pytest使用、段階的テスト（基本→DHCP→Adjoint→CGM→実データ）

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

## ディレクトリ構成
```
TrialClaudeMCPCodex/
├── python/           # Python関連
│   └── original/     # オリジナルPythonコード
│       └── IHCP_CGM_Sliding_Window_Calculation_ver2.py (メインプログラム)
├── shared/           # 共通リソース
│   └── data/         # 共通データファイル
│       ├── metal_thermal_properties.csv (熱物性値データ)
│       └── T_measure_700um_1ms.npy (測定温度データ、1.1GB、gitignore済み)
└── .claude/          # Claude Code設定
    └── CLAUDE.md     # プロジェクト指示書
```

### Julia移植実装済みディレクトリ
```
julia/
├── src/
│   ├── IHCP_CGM.jl                    # メインモジュール（v0.5.0）
│   ├── ThermalProperties.jl           # Phase 1: 熱物性値計算 ✅
│   ├── DataLoaders.jl                 # Phase 1: データ読み込み ✅
│   ├── solvers/
│   │   ├── DHCPSolver.jl             # Phase 2: 直接解法 ✅
│   │   ├── AdjointSolver.jl          # Phase 3: 随伴解法 ✅
│   │   ├── StoppingCriteria.jl       # Phase 4: 停止判定 ✅
│   │   ├── CGMSolver.jl              # Phase 4: 共役勾配法 ✅
│   │   └── SlidingWindowSolver.jl    # Phase 5: スライディングウィンドウ 🔄
│   └── utils/
│       └── json_helpers.jl           # JSON型変換ヘルパー ✅
├── test/
│   ├── runtests.jl                   # テストランナー
│   ├── test_thermal_properties.jl    # Phase 1テスト（25項目）✅
│   ├── test_dhcp_solver.jl           # Phase 2テスト（298項目）✅
│   ├── test_adjoint_solver.jl        # Phase 3テスト（13項目）✅
│   ├── test_cgm_solver.jl            # Phase 4テスト（18項目）✅
│   └── test_sliding_window.jl        # Phase 5テスト（実装完了、実行待ち）🔄
├── data/
│   ├── generate_reference_data.py    # Phase 1参照データ
│   ├── generate_phase2_reference.py  # Phase 2参照データ
│   ├── generate_phase3_reference.py  # Phase 3参照データ
│   ├── generate_phase4_reference.py  # Phase 4参照データ
│   ├── generate_phase5_reference.py  # Phase 5参照データ 🔄
│   └── phase*_reference_*.json       # 各Phaseのゴールデンデータ
└── Project.toml                       # 依存管理

docs/
├── plans/
│   └── julia_migration_plan.md       # 14週間移植計画書
├── logs/
│   ├── phase1_implementation.md      # Phase 1完了ログ
│   ├── phase2_implementation.md      # Phase 2完了ログ
│   ├── phase2_json_fix_and_test_results.md
│   ├── phase3_implementation.md      # Phase 3完了ログ
│   ├── phase3_bug_fix_and_test_results.md
│   └── phase4_implementation.md      # Phase 4完了ログ
└── reviews/
    └── python_code_structure.md      # Pythonコード分析
```

## Git操作ガイドライン
- **大容量ファイル**: 1MB以上のファイルはステージング前に必ず確認
- **除外ファイル**: `T_measure_700um_1ms.npy`（1.1GB）はgitignore済み
- **コミット単位**: TDDサイクルに従い、テスト作成→実装完了で適切に分割

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
```bash
# 個別テスト実行
julia --project=. test/test_thermal_properties.jl
julia --project=. test/test_dhcp_solver.jl
julia --project=. test/test_adjoint_solver.jl

# 全テスト実行
julia --project=. -e 'using Pkg; Pkg.test()'
```

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

### 🔄 実装中
1. **Phase 5（スライディングウィンドウ）**: 実装完了、テスト実行待ち（242行）