# Python-Julia 100ステップ比較テスト パラメータ

## テスト概要
- **目的**: PythonとJuliaの計算結果が一致することを確認
- **データ**: `T_measure_700um_1ms.npy` (1.1GB)
- **計算ステップ数**: 100 (時間方向)
- **Python実行**: Numba並列、8スレッド
- **Julia実行**: 逐次実行

## 計算パラメータ

### グリッド設定
```python
# 空間グリッド
ni = 10          # x方向グリッド数（実データから決定）
nj = 10          # y方向グリッド数（実データから決定）
nk = 20          # z方向グリッド数（固定）

# 時間ステップ
nt = 100         # 時間ステップ数（テスト用）
dt = 1.0e-3      # 時間刻み [s] = 1ms
```

### 空間刻み
```python
# x, y方向（均等格子）
dx = 0.12e-3     # [m] = 0.12mm

# y方向（角度補正）
dy = 0.12e-3 * sin(80°) / sin(45°)  # [m]

# z方向（非均等格子、20層）
# 表面側に密なストレッチングを適用
# dz配列は個別に定義（詳細はコード参照）
```

### 物性値（SUS304）
```python
# 密度 [kg/m³]
rho = 7900.0  # 基準温度での値

# 比熱の多項式係数 [J/(kg·K)]
cp_coeffs = [poly fit from metal_thermal_properties.csv]

# 熱伝導率の多項式係数 [W/(m·K)]
k_coeffs = [poly fit from metal_thermal_properties.csv]

# 温度依存性
cp(T) = cp_coeffs[0] + cp_coeffs[1]*T + cp_coeffs[2]*T² + ...
k(T)  = k_coeffs[0] + k_coeffs[1]*T + k_coeffs[2]*T² + ...
```

### 境界条件
```python
# 上面（z=0, 表面）
# 熱流束 q(x, y, t) を逆解析で推定

# 下面（z=nk-1, 裏面）
# 測定温度 Y_obs(x, y, t) を適用

# 側面（x, y方向）
# 断熱境界条件（熱流束ゼロ）
```

### CGM設定
```python
# 共役勾配法パラメータ
sigma = 1.8           # 測定ノイズ標準偏差 [K]
tolerance = 1.0e-6    # 収束判定閾値
max_iterations = 100  # 最大反復回数

# 初期推定熱流束
q_init = zeros(ni, nj, nt)  # ゼロ初期化
```

### データ読み込み
```python
# Python
T_measure_K = np.load("../shared/data/T_measure_700um_1ms.npy")
Y_obs = T_measure_K[:100, :, :]  # 最初の100ステップを使用
T_init = np.repeat(T_measure_K[0:1, :, :], nk, axis=0)

# Julia
using NPZ
T_measure_K = npzread("../shared/data/T_measure_700um_1ms.npy")
Y_obs = T_measure_K[1:100, :, :]  # 最初の100ステップを使用（1-indexed）
T_init = repeat(T_measure_K[1:1, :, :], nk, 1, 1)
```

## 出力データ

### Python出力
```python
# 保存ファイル
"shared/results/python_100steps.npz"
"shared/results/python_100steps_T.npy"  # 温度場 (ni, nj, nk, nt)
"shared/results/python_100steps_q.npy"  # 熱流束 (ni, nj, nt)
```

### Julia出力
```julia
# 保存ファイル
"shared/results/julia_100steps.npz"
"shared/results/julia_100steps_T.npy"  # 温度場 (ni, nj, nk, nt)
"shared/results/julia_100steps_q.npy"  # 熱流束 (ni, nj, nt)
```

## 比較方法

### 相対誤差計算
```python
# 熱流束の相対誤差
q_diff = abs(q_julia - q_python)
q_rel_err = q_diff / (abs(q_python) + 1e-10)
max_q_rel_err = max(q_rel_err)

# 温度場の相対誤差
T_diff = abs(T_julia - T_python)
T_rel_err = T_diff / (abs(T_python) + 1e-10)
max_T_rel_err = max(T_rel_err)
```

### 合格基準
- **熱流束相対誤差**: < 0.01% (0.0001)
- **温度場相対誤差**: < 0.01% (0.0001)

## 実行コマンド

### Python
```bash
cd python/validation
python run_100steps_test.py
```

### Julia
```bash
cd julia
julia --project=. scripts/run_100steps_test.jl
```

### 比較
```bash
cd python/validation
python compare_python_julia_100steps.py
```
