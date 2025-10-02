# 完全一致検証 実行手順書

## 概要

PythonオリジナルコードとJulia移植版で完全に同じパラメータとデータを使用して計算を実行し、
結果の数値的一致性を検証します。

## 準備

### 1. 必要なファイルの確認

```bash
# プロジェクトルート
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# 共通パラメータファイル
ls -lh shared/config/verification_params.toml

# 測定データ（1分間テストデータ、357ステップ）
ls -lh shared/data/T_measure_test_1min.npy
```

### 2. 環境確認

#### Python環境
```bash
python --version  # Python 3.12以上
pip list | grep -E "numpy|scipy|pandas|numba|toml"
```

必要なパッケージ:
- numpy >= 1.26.4
- scipy >= 1.13.1
- pandas >= 2.2.2
- numba >= 0.60.0
- toml

#### Julia環境
```bash
julia --version  # Julia 1.9以上推奨
```

必要なパッケージ（Julia REPL内）:
```julia
using Pkg
Pkg.status(["NPZ", "TOML", "LinearAlgebra"])
```

## 実行手順

### ステップ1: パラメータ抽出（既に完了）

共通パラメータファイルが既に生成されています。
再生成する場合:

```bash
python python/validation/extract_params.py
```

生成ファイル:
- `shared/config/verification_params.toml` （約4KB）

### ステップ2: Python版実行

```bash
# プロジェクトルートから実行
python python/validation/run_exact_match.py
```

**予想実行時間**: 約5-15分（マシン性能による）

**出力**:
- `shared/results/validation/python_exact.npz` （約200MB）
- 計算時間、誤差統計が画面に表示

**注意事項**:
- メモリ使用量: 約1-2GB
- CGM計算は20反復（テスト用）で実行

### ステップ3: Julia版実行

```bash
# プロジェクトルートから実行
julia julia/examples/run_exact_match.jl
```

**予想実行時間**: 約3-10分（マシン性能による、Python版より高速）

**出力**:
- `shared/results/validation/julia_exact.npz` （約200MB）
- 計算時間、誤差統計が画面に表示

**初回実行時の注意**:
- Juliaのコンパイル時間（プリコンパイル）がかかるため、初回は数分余分にかかります
- 2回目以降は大幅に高速化されます

### ステップ4: 結果比較

```bash
# プロジェクトルートから実行
python python/validation/compare_exact_match.py
```

**出力ファイル**:
1. `shared/results/validation/exact_match_comparison_report.md` - 詳細レポート
2. `shared/results/validation/exact_match_comparison_stats.json` - 統計データ（JSON形式）
3. `shared/results/validation/comparison_plots.png` - 可視化プロット

**判定基準**:
- ✅ 完全一致: 最大絶対誤差 < 1e-6、最大相対誤差 < 1e-8
- ⚠️  高精度一致: 最大絶対誤差 < 1e-3（実用上問題なし）
- ❌ 要調査: 上記を満たさない場合

## トラブルシューティング

### Python版でメモリエラーが発生する場合

```python
# python/validation/run_exact_match.py の設定を確認
# cgm_iterationを減らす（例: 10反復）
```

### Julia版でエラーが発生する場合

```bash
# パッケージの再インストール
julia -e 'using Pkg; Pkg.update(); Pkg.instantiate()'

# モジュールの動作確認
julia -e 'include("julia/src/IHCP_CGM.jl"); using .IHCP_CGM; version_info()'
```

### 比較で大きな誤差が出る場合

確認項目:
1. 使用データの一致: `shared/data/T_measure_test_1min.npy` が同じか
2. パラメータの一致: `shared/config/verification_params.toml` が同じか
3. 実行順序: Python版 → Julia版 → 比較の順で実行したか

## 期待される結果

### 正常な場合の出力例

#### Python版
```
CGM計算完了！
  計算時間: 342.15秒
  結果形状: (356, 80, 100)
  熱流束範囲: -1.23e+04 ~ 5.67e+05 W/m²
  温度誤差（RMS）: 1.23e-02 K
  温度誤差（最大）: 3.45e-02 K
```

#### Julia版
```
CGM計算完了！
  計算時間: 198.73秒
  結果形状: (356, 80, 100)
  熱流束範囲: -1.23e+04 ~ 5.67e+05 W/m²
  温度誤差（RMS）: 1.23e-02 K
  温度誤差（最大）: 3.45e-02 K
  速度比（Python/Julia）: 1.72x
```

#### 比較結果
```
熱流束（q_result）の一致性:
  最大絶対誤差: 3.456e-08
  RMS誤差:      1.234e-09
  最大相対誤差: 5.678e-10

✅ 判定: 完全一致確認
```

## 結果の解釈

### 数値誤差の原因

完全に同じアルゴリズムでも、以下の理由で微小な誤差が生じる可能性があります:

1. **浮動小数点演算の順序**: 行列演算の順序が異なる場合
2. **線形ソルバーの収束**: CGソルバーの反復回数がわずかに異なる場合
3. **コンパイラ最適化**: Python（Numba）とJuliaで異なる最適化が適用される場合

### 許容範囲

- **機械精度レベル（< 1e-14）**: 完璧な一致
- **高精度（< 1e-8）**: 数値計算上完全一致とみなせる
- **実用精度（< 1e-3）**: 物理的に十分な精度
- **要調査（>= 1e-3）**: アルゴリズムの実装に差異がある可能性

## 次のステップ

完全一致が確認できたら:

1. **フルスケールデータでの検証**:
   ```bash
   # 全データ（357ステップ→数千ステップ）で実行
   # shared/data/T_measure_700um_1ms.npy を使用
   ```

2. **性能ベンチマーク**:
   ```bash
   julia julia/benchmarks/benchmark_full_workflow.jl
   ```

3. **実運用への移行**:
   - Julia版をメイン計算エンジンとして使用
   - Python版は検証・可視化用として保持

## ファイル一覧

### 入力
```
shared/config/verification_params.toml    # 共通パラメータ
shared/data/T_measure_test_1min.npy       # 測定データ（357ステップ）
```

### 実行スクリプト
```
python/validation/extract_params.py       # パラメータ抽出
python/validation/run_exact_match.py      # Python版実行
julia/examples/run_exact_match.jl         # Julia版実行
python/validation/compare_exact_match.py  # 結果比較
```

### 出力
```
shared/results/validation/python_exact.npz                # Python計算結果
shared/results/validation/julia_exact.npz                 # Julia計算結果
shared/results/validation/exact_match_comparison_report.md  # 詳細レポート
shared/results/validation/exact_match_comparison_stats.json # 統計データ
shared/results/validation/comparison_plots.png            # 可視化
```

## 問い合わせ

問題が解決しない場合は、以下の情報を添えて報告してください:

1. エラーメッセージ全文
2. 実行環境（OS、Pythonバージョン、Juliaバージョン）
3. メモリ使用量
4. 該当するログファイル

---

最終更新: 2025-10-02
