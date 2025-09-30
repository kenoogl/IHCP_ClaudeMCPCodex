"""
DataLoaders.jl

Phase 1: データ読み込みモジュール
SUS304熱物性値CSVの読み込みと多項式フィッティング

対応するPython処理:
- metal_thermal_properties.csv読み込み (39-52行)
- np.polyfit() による3次多項式フィッティング
"""

module DataLoaders

using DelimitedFiles
using LinearAlgebra

export load_sus304_thermal_properties, polyfit

"""
  load_sus304_thermal_properties(csv_path::String) -> Dict

SUS304熱物性値CSVを読み込み、辞書形式で返す

# 引数
- `csv_path`: CSVファイルのパス

# 戻り値
辞書型で以下のキーを持つ:
- "temperature_K": 温度[K]の配列
- "density": 密度[kg/m³]の配列
- "specific_heat": 比熱[J/(kg·K)]の配列
- "thermal_conductivity": 熱伝導率[W/(m·K)]の配列

# CSVフォーマット
Temperature/C,Temperature/K,Density,Specific_Heat,Thermal_Conductivity
26.85,300,7894,510.37092,12.97
...
"""
function load_sus304_thermal_properties(csv_path::String)
  # CSVファイル読み込み（ヘッダーをスキップ）
  data, header = readdlm(csv_path, ',', Float64, '\n', header=true)

  # UTF-8 BOMを考慮してヘッダーをクリーン化
  header_str = [strip(string(h)) for h in header]

  # 列インデックスの特定
  # ヘッダーに"Temperature/K"が含まれる列を探す
  temp_idx = findfirst(contains("Temperature/K"), header_str)
  density_idx = findfirst(contains("Density"), header_str)
  cp_idx = findfirst(contains("Specific_Heat"), header_str)
  k_idx = findfirst(contains("Thermal_Conductivity"), header_str)

  # データ抽出
  result = Dict{String, Vector{Float64}}(
    "temperature_K" => data[:, temp_idx],
    "density" => data[:, density_idx],
    "specific_heat" => data[:, cp_idx],
    "thermal_conductivity" => data[:, k_idx]
  )

  return result
end


"""
  polyfit(x::Vector{Float64}, y::Vector{Float64}, degree::Int) -> Vector{Float64}

多項式フィッティング（NumPyのpolyfit相当）

# 引数
- `x`: 独立変数配列
- `y`: 従属変数配列
- `degree`: 多項式の次数

# 戻り値
多項式係数配列 [a, b, c, d] for y = ax³ + bx² + cx + d
（高次から低次の順序、NumPyと同じ）

# 実装詳細
最小二乗法により多項式係数を求める
Vandermonde行列を構築して線形方程式を解く

A * coeffs = y
ここで A[i,j] = x[i]^(degree-j)
"""
function polyfit(x::Vector{Float64}, y::Vector{Float64}, degree::Int)::Vector{Float64}
  n = length(x)
  @assert n == length(y) "x と y の長さが一致しません"
  @assert n >= degree + 1 "データ点数が不足しています"

  # Vandermonde行列の構築
  # A[i, j] = x[i]^(degree - j + 1)
  # 例: degree=3 の場合
  #   A = [x[1]^3  x[1]^2  x[1]^1  x[1]^0]
  #       [x[2]^3  x[2]^2  x[2]^1  x[2]^0]
  #       ...
  A = zeros(n, degree + 1)

  for i in 1:n
    for j in 1:(degree + 1)
      A[i, j] = x[i]^(degree - j + 1)
    end
  end

  # 最小二乗法で係数を求める
  # coeffs = (A^T * A)^(-1) * A^T * y
  coeffs = A \ y

  return coeffs
end


"""
  fit_sus304_coefficients(
    sus304_data::Dict{String, Vector{Float64}},
    degree::Int=3
  ) -> Dict{String, Vector{Float64}}

SUS304データから多項式係数を計算

# 引数
- `sus304_data`: load_sus304_thermal_properties()の出力
- `degree`: 多項式次数（デフォルト: 3）

# 戻り値
辞書型で以下のキーを持つ:
- "rho_coeffs": 密度の多項式係数
- "cp_coeffs": 比熱の多項式係数
- "k_coeffs": 熱伝導率の多項式係数
"""
function fit_sus304_coefficients(
  sus304_data::Dict{String, Vector{Float64}},
  degree::Int=3
)
  temp = sus304_data["temperature_K"]

  coeffs = Dict{String, Vector{Float64}}(
    "rho_coeffs" => polyfit(temp, sus304_data["density"], degree),
    "cp_coeffs" => polyfit(temp, sus304_data["specific_heat"], degree),
    "k_coeffs" => polyfit(temp, sus304_data["thermal_conductivity"], degree)
  )

  return coeffs
end

end # module DataLoaders