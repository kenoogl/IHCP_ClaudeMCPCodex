# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要
逆熱伝導問題（IHCP）をコンジュゲートグラディエント法（CGM）で解く数値計算プロジェクト。IRカメラからの温度測定データを使用してSUS304材の表面熱流束を逆解析する。

## 実行方法
### メインプログラム実行
```bash
cd org/
python IHCP_CGM_Sliding_Window_Calculation_ver2.py
```

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
- **インデント**: スペース4文字
- **最適化**: Numba @njitデコレータで高速化
- **並列処理**: prange使用で性能向上
- **エラー処理**: 数値計算の収束性とメモリ制限を考慮

## 重要な注意事項
- **大容量データ**: T_measure_700um_1ms.npy（1.1GB）の取り扱いに注意
- **計算時間**: フルスケール計算は数時間を要する場合がある
- **メモリ使用**: 大規模問題では8GB以上のメモリが必要
- **数値精度**: 逆問題の性質上、収束判定とステップサイズ調整が重要

## MCPによるcodexとの連携
- Task toolを使って、Codexと連携する
- あなたの役割は、タスク全体の計画と進捗管理、ユーザーとのコミュニケーション、Codexへの指示出しと結果のレビューです
- Codexには、詳細な調査作業、コードの実装、ファイル検索や分析を担当してもらってください
- Juliaの実行環境やテスト実行についてcodexでうまくいかない場合には、あなたがサポートしてください
- 私があなたに出した依頼事項を元に、あなたはcodexにプロンプトを渡します。その内容とcodexからのレポートを一つのタスクの記録として、適切な名称をつけて、ログに保存してください。

## Claude codeのコンテキスト圧縮後
- /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/.claude/CLAUDE.mdを読み、開発方針の確認をする

## codexセッション再開時
- セッション毎に、Codex MCPとの連携を確認し、未連携の場合には連携処理をする

## ディレクトリ構成
 IHCP/TrialClaude2MCPcodex/
 ├── docs/          # ドキュメント類
 │  ├── reviews/      # レビュー・報告書
 │  ├── plans/       # 計画・方針文書
 │  └── logs/        # 変換ログ・進捗記録
 ├── python/         # Python関連
 │  ├── original/      # オリジナルPythonコード
 │  ├── improved/      # 改善済みPythonコード
 │  ├── tests/       # Pythonテスト
 │  └── data/        # Python用データ
 ├── julia/         # Julia関連
 │  ├── src/        # Juliaソースコード
 │  ├── test/        # Juliaテスト
 │  ├── data/        # Julia用データ
 │  ├── benchmarks/     # 性能測定
 │  └── examples/      # 使用例
 ├── shared/         # 共通リソース
 │  ├── data/        # 共通データファイル
 │  ├── configs/      # 設定ファイル
 │  └── scripts/      # ユーティリティスクリプト
 └── .claude/        # Claude設定（既存）

 ## git操作
 - 1MB以上のファイルはステージングするかどうか、常に確認すること