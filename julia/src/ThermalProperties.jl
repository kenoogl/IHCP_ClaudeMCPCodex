"""
ThermalProperties.jl

熱物性値計算モジュール
温度依存性を持つSUS304材料の熱物性値（比熱cp、熱伝導率k）を計算する

対応するPython関数:
- polyval_numba() (57-61行)
- thermal_properties_calculator() (63-78行)
"""

module ThermalProperties

export polyval, thermal_properties!

"""
多項式評価関数

# 引数
- `coeffs`: 多項式係数配列 [a, b, c, d] for y = ax³ + bx² + cx + d
- `x`: 評価点（温度K）

# 戻り値
- 多項式の評価結果
result = a*x^3 + b*x^2 + c*x + d

係数の順序はNumPyと同じ（高次から低次）
"""
function polyval(coeffs::Vector{Float64}, x::Float64)::Float64
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
3D温度から比熱cpと熱伝導率λを計算

# 引数
- `Temperature`: 3D温度配列 (ni, nj, nk) [K]
- `cp`: 比熱配列 (ni, nj, nk) [J/(kg·K)]
- `λ`: 熱伝導率配列 (ni, nj, nk) [W/(m·K)]
- `cp_coeffs`: 比熱の多項式係数 [a, b, c, d]
- `k_coeffs`: 熱伝導率の多項式係数 [a, b, c, d]
"""
function thermal_properties!(
  Temperature::Array{Float64, 3},
  cp::Array{Float64, 3},
  λ::Array{Float64, 3},
  cp_coeffs::Vector{Float64},
  k_coeffs::Vector{Float64}
)
  ni, nj, nk = size(Temperature)

  for k in 1:nk, j in 1:nj, i in 1:ni
    T_current = Temperature[i, j, k]
    cp[i, j, k] = polyval(cp_coeffs, T_current)
    λ[i, j, k]  = polyval(k_coeffs, T_current)
  end
end

end # module ThermalProperties