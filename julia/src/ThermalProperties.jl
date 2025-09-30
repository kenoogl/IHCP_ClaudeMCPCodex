"""
ThermalProperties.jl

Phase 1: 熱物性値計算モジュール
温度依存性を持つSUS304材料の熱物性値（比熱cp、熱伝導率k）を計算する

対応するPython関数:
- polyval_numba() (57-61行)
- thermal_properties_calculator() (63-78行)
"""

module ThermalProperties

export polyval_numba, thermal_properties_calculator

"""
  polyval_numba(coeffs::Vector{Float64}, x::Float64) -> Float64

多項式評価関数（NumPyのpolyvalと同等）

# 引数
- `coeffs`: 多項式係数配列 [a, b, c, d] for y = ax³ + bx² + cx + d
- `x`: 評価点（温度K）

# 戻り値
- 多項式の評価結果

# 実装詳細
Python版のNumba実装と同じアルゴリズム:
result = a*x^3 + b*x^2 + c*x + d

係数の順序はNumPyと同じ（高次から低次）
"""
function polyval_numba(coeffs::Vector{Float64}, x::Float64)::Float64
  result = 0.0
  n = length(coeffs)

  for i in 1:n
    # coeffs[i]はx^(n-i)の係数
    # 例: coeffs = [a, b, c, d] (n=4)
    #   i=1: a * x^3
    #   i=2: b * x^2
    #   i=3: c * x^1
    #   i=4: d * x^0
    result += coeffs[i] * x^(n - i)
  end

  return result
end


"""
  thermal_properties_calculator(
    Temperature::Array{Float64, 3},
    cp_coeffs::Vector{Float64},
    k_coeffs::Vector{Float64}
  ) -> (Array{Float64, 3}, Array{Float64, 3})

3D温度配列から比熱cpと熱伝導率kを計算

# 引数
- `Temperature`: 3D温度配列 (nk, nj, ni) [K]
- `cp_coeffs`: 比熱の多項式係数 [a, b, c, d]
- `k_coeffs`: 熱伝導率の多項式係数 [a, b, c, d]

# 戻り値
- `cp`: 比熱配列 (nk, nj, ni) [J/(kg·K)]
- `k`: 熱伝導率配列 (nk, nj, ni) [W/(m·K)]

# 配列順序の注意
PythonとJuliaで配列の格納順序が異なる:
- Python: 行優先 (i, j, k) → メモリ上でk方向が最も密
- Julia: 列優先 (k, j, i) → メモリ上でk方向が最も密

この実装では、Julia標準の列優先で処理するため、
配列形状を(nk, nj, ni)として扱う。

# 並列処理
Julia標準の型推論による最適化を利用
必要に応じて@threadsマクロで並列化可能
"""
function thermal_properties_calculator(
  Temperature::Array{Float64, 3},
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64}
)
  nk, nj, ni = size(Temperature)

  # 出力配列の事前確保
  cp = Array{Float64, 3}(undef, nk, nj, ni)
  k = Array{Float64, 3}(undef, nk, nj, ni)

  # 3重ループで全格子点の熱物性値を計算
  # Juliaの列優先に合わせてループ順序を最適化
  for i in 1:ni
    for j in 1:nj
      for k_idx in 1:nk
        # 現在の格子点温度
        T_current = Temperature[k_idx, j, i]

        # 多項式評価で比熱と熱伝導率を計算
        cp[k_idx, j, i] = polyval_numba(cp_coeffs, T_current)
        k[k_idx, j, i] = polyval_numba(k_coeffs, T_current)
      end
    end
  end

  return cp, k
end

end # module ThermalProperties