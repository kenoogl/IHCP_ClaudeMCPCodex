# PythonとJuliaの完全一致検証環境構築完了レポート

**作成日**: 2025-10-02
**プロジェクト**: 逆熱伝導問題（IHCP）共役勾配法（CGM）ソルバー
**目的**: PythonオリジナルコードとJulia移植版の数値計算完全一致検証

---

## エグゼクティブサマリー

PythonとJuliaで完全に同じパラメータとデータを使用して計算結果の一致性を詳細に検証するための環境を構築しました。

### 達成事項

✅ **完了**: 完全一致検証に必要なすべてのスクリプトとツールの作成

- 共通パラメータファイル（TOML形式）
- Python版実行スクリプト
- Julia版実行スクリプト
- 詳細比較スクリプト
- 実行手順書

### 主要成果物

1. **共通パラメータファイル**: `shared/config/verification_params.toml` (3.2KB)
2. **Python実行スクリプト**: `python/validation/run_exact_match.py` (8.6KB)
3. **Julia実行スクリプト**: `julia/examples/run_exact_match.jl` (8.5KB)
4. **比較スクリプト**: `python/validation/compare_exact_match.py` (15KB)
5. **実行手順書**: `shared/results/validation/README_EXACT_MATCH.md`

---

## 1. プロジェクト背景

### 課題

PythonオリジナルコードをJuliaに移植したが、数値計算の完全一致を確認する必要があった。

### 要件

- PythonとJuliaで**完全に同じパラメータ**を使用
- **同じ測定データ**を読み込み
- **同じアルゴリズム**で計算実行
- 結果の**詳細な数値比較**
- **誤差の可視化と統計分析**

---

## 2. 実装内容

### 2.1 共通パラメータファイル（TOML形式）

**ファイル**: `shared/config/verification_params.toml`

Pythonオリジナルコードから抽出したパラメータを一元管理:

```toml
[problem]
ni = 80
nj = 100
nk = 20
nt = 357  # 1分間テストデータ
dt = 0.001  # [s]
dx = 0.00012  # [m]
dy = 0.00016712741767680456  # [m]

[material]
rho_reference_value = 7823.493962874829  # [kg/m³]
cp_coeffs = [2.009e-10, -3.426e-07, 0.134928, 469.853]
k_coeffs = [4.799e-12, -8.183e-09, 0.016177, 8.118]

[cgm]
window_size = 71
overlap = 17
cgm_iteration = 20  # テスト用（実際は20000）
q_init_value = 0.0  # [W/m²]
```

**特徴**:
- 浮動小数点数を高精度で保存（機械精度レベル）
- PythonとJuliaの両方から読み込み可能
- バージョン管理が容易

### 2.2 パラメータ抽出スクリプト

**ファイル**: `python/validation/extract_params.py`

**機能**:
1. Pythonオリジナルコード (`IHCP_CGM_Sliding_Window_Calculation_ver2.py`) から実際のパラメータを抽出
2. SUS304熱物性値の多項式フィッティング係数を計算
3. z方向非均等格子の生成
4. TOML形式で保存

**実行結果**:
```bash
$ python python/validation/extract_params.py
============================================================
パラメータ抽出開始
============================================================

熱物性値多項式係数:
  密度 (rho @ 498.15K): 7823.493963 kg/m³
  rho_coeffs: [2.45e-09, -0.000122, -0.261, 7983.29]
  cp_coeffs:  [2.01e-10, -3.43e-07, 0.135, 469.85]
  k_coeffs:   [4.80e-12, -8.18e-09, 0.0162, 8.12]

格子パラメータ:
  nz: 20
  dx: 1.2000000000e-04 m
  dy: 1.6712741768e-04 m

パラメータを保存: shared/config/verification_params.toml
============================================================
```

### 2.3 Python版実行スクリプト

**ファイル**: `python/validation/run_exact_match.py`

**処理フロー**:
1. TOMLファイルからパラメータ読み込み
2. 測定データ (`T_measure_test_1min.npy`) 読み込み
3. 初期温度場の厳密な生成（全格子点で300K）
4. スライディングウィンドウCGM計算実行
5. 結果検証（DHCP順解析）
6. NPZ形式で結果保存

**出力**:
- `shared/results/validation/python_exact.npz` (約200MB)
  - `q_result`: 逆解析熱流束 (356, 80, 100)
  - `T_verify`: 検証温度場 (357, 80, 100, 20)
  - `elapsed_time_cgm`: CGM計算時間
  - `rms_error`, `max_error`: 温度誤差統計

### 2.4 Julia版実行スクリプト

**ファイル**: `julia/examples/run_exact_match.jl`

**処理フロー**:
Python版と完全に同じ処理を実行:
1. TOMLファイル読み込み（`TOML.parsefile`）
2. 測定データ読み込み（`npzread`）
3. 初期温度場生成（`fill(300.0, ni, nj, nk)`）
4. スライディングウィンドウCGM計算（`solve_sliding_window_cgm`）
5. DHCP検証計算（`solve_dhcp_multiple_timesteps`）
6. NPZ形式で保存（`npzwrite`）

**Julia側の追加実装**:
- `solve_dhcp_multiple_timesteps`: `solve_dhcp!`のエイリアス関数
- キーワード引数対応: `rtol_dhcp`, `maxiter_dhcp`, `rtol_adjoint`, `maxiter_adjoint`
- IHCP_CGM.jlモジュールへのエクスポート追加

**コード変更箇所**:
1. `julia/src/solvers/DHCPSolver.jl`: エイリアス関数追加（460-482行）
2. `julia/src/solvers/SlidingWindowSolver.jl`: キーワード引数追加（99-102行、165-172行）
3. `julia/src/IHCP_CGM.jl`: エクスポート追加（65行）

### 2.5 詳細比較スクリプト

**ファイル**: `python/validation/compare_exact_match.py`

**機能**:

#### (1) 配列詳細比較
- 形状の一致確認
- 絶対誤差統計: 最大、平均、RMS
- 相対誤差統計: 最大、平均
- 値の範囲比較

#### (2) 時間ステップ毎の誤差分析
サンプリング時刻（t=0, 89, 178, 267, 356）での誤差プロファイル

#### (3) 空間分布の誤差分析
最終ステップでのパーセンタイル分析（0, 25, 50, 75, 95, 100%点）

#### (4) 可視化
6枚組プロット:
- 時間変化（中心点）
- 誤差の時間変化（最大・平均）
- 誤差ヒストグラム
- Python空間分布（最終ステップ）
- Julia空間分布（最終ステップ）
- 誤差マップ（最終ステップ）

#### (5) レポート生成
- Markdown形式: 判定基準、詳細統計、総合判定
- JSON形式: 機械可読な統計データ

**判定基準**:
```python
criteria = {
    '初期条件': max_abs_error < 1e-14,
    '熱流束（最大絶対誤差）': max_abs_error < 1e-6,
    '熱流束（最大相対誤差）': max_rel_error < 1e-8,
}
```

### 2.6 実行手順書

**ファイル**: `shared/results/validation/README_EXACT_MATCH.md`

**内容**:
- 準備（環境確認、ファイル確認）
- 実行手順（4ステップ）
- トラブルシューティング
- 期待される結果の例示
- 数値誤差の解釈ガイド

---

## 3. 技術的詳細

### 3.1 パラメータの完全一致保証

#### 浮動小数点数の扱い
- Python: `float64` (IEEE 754倍精度)
- Julia: `Float64` (IEEE 754倍精度)
- TOML保存時: 最大精度で文字列化（17桁有効数字）

#### 配列順序の統一
- Python: row-major (C順序) → Fortran順序で保存
- Julia: column-major (Fortran順序) → 自然に対応
- NPZ読み書きで順序を保証

#### 無限大の扱い
```python
# Python (TOML保存)
dz_t = [x if not np.isinf(x) else 'inf' for x in dz_t]

# Julia (TOML読込)
dz_t = [x == "inf" ? Inf : Float64(x) for x in dz_t_raw]
```

### 3.2 初期条件の厳密な定義

**Python**:
```python
T_init = np.full((ni, nj, nk), T_init_value, dtype=np.float64)
```

**Julia**:
```julia
T_init = fill(T_init_value, ni, nj, nk)  # Float64型
```

両方とも**機械精度レベルで一致**（誤差 < 1e-15）

### 3.3 アルゴリズムの完全一致

#### 使用ソルバー
- **DHCP**: 共役勾配法（CG）+ 対角プリコンディショナ
- **Adjoint**: 共役勾配法（CG）+ 対角プリコンディショナ
- **CGM**: ポラック・リビエール法

#### 収束判定
- Python: `scipy.sparse.linalg.cg` (rtol, maxiter)
- Julia: `IterativeSolvers.cg!` (reltol, maxiter)
- **完全に同じ収束条件**

#### スライディングウィンドウ
- オーバーラップ平均化: `0.5 * old + new`（Python原典1609行）
- ウィンドウ分割: Python原典1577-1619行に完全準拠

---

## 4. ファイル構成

### 入力ファイル
```
shared/
├── config/
│   └── verification_params.toml          # 共通パラメータ（3.2KB）
└── data/
    └── T_measure_test_1min.npy            # 測定データ（357ステップ、22MB）
```

### 実行スクリプト
```
python/validation/
├── extract_params.py                     # パラメータ抽出（7.5KB）
├── run_exact_match.py                    # Python版実行（8.6KB）
└── compare_exact_match.py                # 比較スクリプト（15KB）

julia/examples/
└── run_exact_match.jl                    # Julia版実行（8.5KB）
```

### 出力ファイル（実行後に生成）
```
shared/results/validation/
├── python_exact.npz                      # Python計算結果（約200MB）
├── julia_exact.npz                       # Julia計算結果（約200MB）
├── exact_match_comparison_report.md      # 詳細レポート
├── exact_match_comparison_stats.json     # 統計データ（JSON）
├── comparison_plots.png                  # 可視化プロット
└── README_EXACT_MATCH.md                 # 実行手順書
```

---

## 5. 実行方法（要約）

### 基本フロー

```bash
# プロジェクトルートへ移動
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# ステップ1: パラメータ抽出（既に完了）
python python/validation/extract_params.py

# ステップ2: Python版実行（5-15分）
python python/validation/run_exact_match.py

# ステップ3: Julia版実行（3-10分）
julia julia/examples/run_exact_match.jl

# ステップ4: 結果比較
python python/validation/compare_exact_match.py
```

### 期待される出力

#### 完全一致の場合
```
熱流束（q_result）の一致性:
  最大絶対誤差: 3.456e-08
  RMS誤差:      1.234e-09
  最大相対誤差: 5.678e-10

✅ 判定: 完全一致確認
```

---

## 6. 次のステップ

### 短期（今後1週間）

1. **実際の実行**:
   - Python版実行（約10分）
   - Julia版実行（約5分）
   - 比較実行（約1分）

2. **結果の確認**:
   - 一致性レポートのレビュー
   - 誤差の原因分析（もしあれば）

### 中期（今後1ヶ月）

1. **フルスケールデータでの検証**:
   - 全データ（357 → 数千ステップ）での実行
   - メモリ最適化の検討

2. **性能ベンチマーク**:
   - Julia版の高速化検証
   - 並列化の検討

### 長期（今後3ヶ月）

1. **実運用への移行**:
   - Julia版をメイン計算エンジンとして採用
   - Python版は検証・可視化用として保持

2. **機能拡張**:
   - リアルタイム計算対応
   - GPU対応の検討

---

## 7. 技術的成果

### 達成した技術要素

1. **クロスランゲージパラメータ管理**:
   - TOML形式による一元管理
   - 高精度浮動小数点数の保存

2. **数値計算の再現性**:
   - 完全に同じアルゴリズム実装
   - 決定的な初期値設定

3. **詳細な数値比較フレームワーク**:
   - 多角的な誤差分析
   - 自動化された判定システム

4. **包括的なドキュメント**:
   - 実行手順書
   - トラブルシューティングガイド
   - 期待される結果の明示

### 技術的課題の解決

| 課題 | 解決策 |
|-----|--------|
| パラメータの不一致 | TOML共通ファイル化 |
| 配列順序の違い | NPZ形式で統一 |
| 無限大の扱い | 文字列 "inf" で統一 |
| 関数名の違い | エイリアス関数追加 |
| デバッグの困難性 | 詳細な統計・可視化 |

---

## 8. まとめ

### 完了事項

✅ 共通パラメータファイルの作成（TOML形式）
✅ Python版実行スクリプトの作成
✅ Julia版実行スクリプトの作成
✅ 詳細比較スクリプトの作成
✅ 実行手順書の作成
✅ Julia側の関数追加・修正

### 残タスク

⏳ 実際の計算実行（Python版）
⏳ 実際の計算実行（Julia版）
⏳ 結果比較の実行
⏳ 完全一致検証レポートの作成

### 推定作業時間

- Python版実行: 5-15分
- Julia版実行: 3-10分（初回は+5分プリコンパイル）
- 比較実行: 1分
- 結果レビュー: 10分

**合計**: 約30-40分

---

## 9. 関連ドキュメント

- **実行手順**: `shared/results/validation/README_EXACT_MATCH.md`
- **Pythonオリジナルコード**: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`
- **Juliaメインモジュール**: `julia/src/IHCP_CGM.jl`
- **共通パラメータ**: `shared/config/verification_params.toml`

---

## 10. 連絡先・サポート

問題が発生した場合は、以下の情報を添えて報告:

1. エラーメッセージ全文
2. 実行環境（OS、Pythonバージョン、Juliaバージョン）
3. 該当するログファイル

---

**作成者**: Claude (Anthropic)
**最終更新**: 2025-10-02
**バージョン**: 1.0
**プロジェクトフェーズ**: 完全一致検証環境構築完了
