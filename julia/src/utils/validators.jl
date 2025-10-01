"""
validators.jl

数値計算の健全性チェック機能

機能:
- NaN/Inf値検出
- 温度・熱流束範囲チェック
- 空間勾配急変検出
- 包括的異常検出
"""

module Validators

export check_field_finite,
       check_temperature_range,
       check_flux_range,
       check_gradient_magnitude,
       detect_numerical_anomalies,
       check_temperature_field,
       check_flux_field,
       check_adjoint_field

"""
    check_field_finite(field::Array{Float64}) -> Bool

NaN/Inf値の検出

# 引数
- `field::Array{Float64}`: 検査対象フィールド

# 戻り値
- `Bool`: すべて有限の場合true
"""
function check_field_finite(field::Array{Float64})::Bool
  flat_field = vec(field)
  @inbounds for i in 1:length(flat_field)
    if !isfinite(flat_field[i])
      return false
    end
  end
  return true
end

"""
    check_temperature_range(T, min_temp, max_temp) -> Tuple{Bool, Float64, Float64}

温度場の物理的範囲チェック

# 引数
- `T::Array{Float64}`: 温度場 [K]
- `min_temp::Float64`: 最小許容温度 [K] (デフォルト: 150.0)
- `max_temp::Float64`: 最大許容温度 [K] (デフォルト: 3000.0)

# 戻り値
- `is_valid::Bool`: 範囲内の場合true
- `min_val::Float64`: 実際の最小値 [K]
- `max_val::Float64`: 実際の最大値 [K]
"""
function check_temperature_range(
  T::Array{Float64},
  min_temp::Float64=150.0,
  max_temp::Float64=3000.0
)::Tuple{Bool, Float64, Float64}
  flat_T = vec(T)
  min_val = flat_T[1]
  max_val = flat_T[1]

  @inbounds @simd for i in 1:length(flat_T)
    val = flat_T[i]
    min_val = min(min_val, val)
    max_val = max(max_val, val)
  end

  is_valid = (min_val >= min_temp) && (max_val <= max_temp)
  return is_valid, min_val, max_val
end

"""
    check_flux_range(q, max_abs_flux) -> Tuple{Bool, Float64, Float64}

熱流束の範囲チェック

# 引数
- `q::Array{Float64}`: 熱流束 [W/m²]
- `max_abs_flux::Float64`: 最大許容絶対値 [W/m²] (デフォルト: 1e7)

# 戻り値
- `is_valid::Bool`: 範囲内の場合true
- `min_val::Float64`: 実際の最小値 [W/m²]
- `max_val::Float64`: 実際の最大値 [W/m²]
"""
function check_flux_range(
  q::Array{Float64},
  max_abs_flux::Float64=1e7
)::Tuple{Bool, Float64, Float64}
  flat_q = vec(q)
  min_val = flat_q[1]
  max_val = flat_q[1]

  @inbounds @simd for i in 1:length(flat_q)
    val = flat_q[i]
    min_val = min(min_val, val)
    max_val = max(max_val, val)
  end

  max_abs = max(abs(min_val), abs(max_val))
  is_valid = max_abs <= max_abs_flux
  return is_valid, min_val, max_val
end

"""
    check_gradient_magnitude(field, dx, dy, dz, max_grad_temp, max_grad_flux) -> Tuple{Bool, Float64}

空間勾配の急変検出（2D/3D対応）

# 引数
- `field::Array{Float64}`: フィールド（温度または熱流束）
- `dx::Float64`: x方向格子間隔 [m]
- `dy::Float64`: y方向格子間隔 [m]
- `dz::Vector{Float64}`: z方向格子間隔配列 [m]
- `max_grad_temp::Float64`: 温度勾配の最大許容値 [K/m] (デフォルト: 5000.0)
- `max_grad_flux::Float64`: 熱流束勾配の最大許容値 [W/m³] (デフォルト: 1e8)

# 戻り値
- `is_valid::Bool`: 許容範囲内の場合true
- `max_grad::Float64`: 最大勾配値
"""
function check_gradient_magnitude(
  field::Array{Float64},
  dx::Float64,
  dy::Float64,
  dz::Vector{Float64},
  max_grad_temp::Float64=5000.0,
  max_grad_flux::Float64=1e8
)::Tuple{Bool, Float64}
  max_grad = 0.0
  ndims_field = ndims(field)

  if ndims_field == 2
    # 2次元フィールド（表面熱流束）
    ni, nj = size(field)

    # i方向勾配（物理的にはy方向、注: Pythonコメント参照）
    @inbounds for i in 1:(ni-1)
      for j in 1:nj
        grad = abs(field[i+1, j] - field[i, j]) / dy
        max_grad = max(max_grad, grad)
      end
    end

    # j方向勾配（物理的にはx方向）
    @inbounds for i in 1:ni
      for j in 1:(nj-1)
        grad = abs(field[i, j+1] - field[i, j]) / dx
        max_grad = max(max_grad, grad)
      end
    end

    # 熱流束勾配の閾値で判定
    is_valid = max_grad <= max_grad_flux

  elseif ndims_field == 3
    # 3次元フィールド（温度場）
    ni, nj, nk = size(field)

    # i方向勾配（物理的にはy方向）
    @inbounds for i in 1:(ni-1)
      for j in 1:nj
        for k in 1:nk
          grad = abs(field[i+1, j, k] - field[i, j, k]) / dy
          max_grad = max(max_grad, grad)
        end
      end
    end

    # j方向勾配（物理的にはx方向）
    @inbounds for i in 1:ni
      for j in 1:(nj-1)
        for k in 1:nk
          grad = abs(field[i, j+1, k] - field[i, j, k]) / dx
          max_grad = max(max_grad, grad)
        end
      end
    end

    # k方向勾配（z方向、非均等格子対応）
    @inbounds for i in 1:ni
      for j in 1:nj
        for k in 1:(nk-1)
          grad = abs(field[i, j, k+1] - field[i, j, k]) / dz[k]
          max_grad = max(max_grad, grad)
        end
      end
    end

    # 温度勾配の閾値で判定
    is_valid = max_grad <= max_grad_temp

  else
    error("Unsupported field dimension: $ndims_field (expected 2 or 3)")
  end

  return is_valid, max_grad
end

"""
    detect_numerical_anomalies(field, field_name; kwargs...) -> Tuple{Bool, Vector{String}}

包括的な異常検出メインルーチン

# 引数
- `field::Array{Float64}`: 検査対象の物理場
- `field_name::String`: フィールド名（温度場、熱流束、勾配等）

# キーワード引数
- `iteration::Union{Int, Nothing}`: CGM反復回数
- `timestep::Union{Int, Nothing}`: 時間ステップ
- `temperature_range::Tuple{Float64, Float64}`: 温度の物理的範囲 [K]
- `flux_range::Tuple{Float64, Float64}`: 熱流束の物理的範囲 [W/m²]
- `dx::Float64`: x方向格子間隔 [m]
- `dy::Float64`: y方向格子間隔 [m]
- `dz::Union{Vector{Float64}, Nothing}`: z方向格子間隔配列 [m]

# 戻り値
- `has_anomaly::Bool`: 異常が検出された場合true
- `anomalies::Vector{String}`: 検出された異常のリスト
"""
function detect_numerical_anomalies(
  field::Array{Float64},
  field_name::String;
  iteration::Union{Int, Nothing}=nothing,
  timestep::Union{Int, Nothing}=nothing,
  temperature_range::Tuple{Float64, Float64}=(150.0, 3000.0),
  flux_range::Tuple{Float64, Float64}=(-1e7, 1e7),
  dx::Float64=0.12e-3,
  dy::Float64=0.12e-3*0.866,
  dz::Union{Vector{Float64}, Nothing}=nothing
)::Tuple{Bool, Vector{String}}
  anomalies = String[]

  # 1. NaN/Inf値チェック
  if !check_field_finite(field)
    push!(anomalies, "NaN/Inf値検出")
  end

  # 2. フィールド種類別の詳細チェック
  # 温度場の判定
  is_temperature = occursin("温度", field_name) ||
                   occursin("Temperature", field_name) ||
                   occursin("T_", field_name)

  # 熱流束場の判定
  is_flux = occursin("熱流束", field_name) ||
            occursin("flux", field_name) ||
            occursin("q_", field_name) ||
            field_name == "q"

  # 随伴場・勾配場の判定
  is_adjoint_or_grad = occursin("lambda", field_name) ||
                       occursin("勾配", field_name) ||
                       occursin("gradient", field_name)

  if is_temperature
    # 温度場の範囲チェック
    min_temp, max_temp = temperature_range
    is_valid, min_val, max_val = check_temperature_range(field, min_temp, max_temp)

    if !is_valid
      if min_val < min_temp
        push!(anomalies, "異常低温検出: min=$(min_val)K < $(min_temp)K")
      end
      if max_val > max_temp
        push!(anomalies, "異常高温検出: max=$(max_val)K > $(max_temp)K")
      end
    end

    # 温度勾配チェック（3次元フィールドかつdz指定あり）
    if ndims(field) == 3 && dz !== nothing
      is_grad_valid, max_grad = check_gradient_magnitude(field, dx, dy, dz, 5000.0, 1e8)
      if !is_grad_valid
        push!(anomalies, "異常な温度勾配: $(max_grad)K/m > 5000K/m")
      end
    end

  elseif is_flux
    # 熱流束範囲チェック
    max_abs_flux = max(abs(flux_range[1]), abs(flux_range[2]))
    is_valid, min_val, max_val = check_flux_range(field, max_abs_flux)

    if !is_valid
      max_abs = max(abs(min_val), abs(max_val))
      push!(anomalies, "異常な熱流束: max_abs=$(max_abs)W/m² > $(max_abs_flux)W/m²")
    end

    # 熱流束勾配チェック（2次元フィールド）
    if ndims(field) == 2 && dz !== nothing
      is_grad_valid, max_grad = check_gradient_magnitude(field, dx, dy, dz, 5000.0, 1e8)
      if !is_grad_valid
        push!(anomalies, "異常な熱流束勾配: $(max_grad)W/m³ > 1e8W/m³")
      end
    end

  elseif is_adjoint_or_grad
    # 随伴場・勾配場の範囲チェック（絶対値 ≤ 1e10）
    flat_field = vec(field)
    max_abs = maximum(abs.(flat_field))
    if max_abs > 1e10
      push!(anomalies, "異常な随伴場/勾配場: max_abs=$(max_abs) > 1e10")
    end
  end

  # 3. オーバーフロー検査
  float64_max = floatmax(Float64)
  flat_field = vec(field)
  max_abs_val = maximum(abs.(flat_field))
  if max_abs_val > float64_max * 0.1
    push!(anomalies, "数値オーバーフローの危険: max_abs=$(max_abs_val) > $(float64_max * 0.1)")
  end

  has_anomaly = length(anomalies) > 0
  return has_anomaly, anomalies
end

"""
    check_temperature_field(T; kwargs...) -> Tuple{Bool, String}

温度場の包括的チェック（ラッパー関数）

# 引数
- `T::Array{Float64}`: 温度場 [K]

# キーワード引数
- `timestep::Union{Int, Nothing}`: 時間ステップ
- `dx::Float64`: x方向格子間隔 [m]
- `dy::Float64`: y方向格子間隔 [m]
- `dz::Union{Vector{Float64}, Nothing}`: z方向格子間隔配列 [m]

# 戻り値
- `is_valid::Bool`: 異常検出時false
- `status_msg::String`: 状態メッセージ（"正常" or "異常検出: N件"）
"""
function check_temperature_field(
  T::Array{Float64};
  timestep::Union{Int, Nothing}=nothing,
  dx::Float64=0.12e-3,
  dy::Float64=0.12e-3*0.866,
  dz::Union{Vector{Float64}, Nothing}=nothing
)::Tuple{Bool, String}
  has_anomaly, anomalies = detect_numerical_anomalies(
    T, "温度場",
    timestep=timestep,
    temperature_range=(150.0, 3000.0),
    dx=dx, dy=dy, dz=dz
  )

  if has_anomaly
    return false, "異常検出: $(length(anomalies))件"
  else
    return true, "正常"
  end
end

"""
    check_flux_field(q; kwargs...) -> Tuple{Bool, String}

熱流束場の包括的チェック（ラッパー関数）

# 引数
- `q::Array{Float64}`: 熱流束場 [W/m²]

# キーワード引数
- `iteration::Union{Int, Nothing}`: CGM反復回数
- `timestep::Union{Int, Nothing}`: 時間ステップ
- `dx::Float64`: x方向格子間隔 [m]
- `dy::Float64`: y方向格子間隔 [m]

# 戻り値
- `is_valid::Bool`: 異常検出時false
- `status_msg::String`: 状態メッセージ
"""
function check_flux_field(
  q::Array{Float64};
  iteration::Union{Int, Nothing}=nothing,
  timestep::Union{Int, Nothing}=nothing,
  dx::Float64=0.12e-3,
  dy::Float64=0.12e-3*0.866
)::Tuple{Bool, String}
  # dzはダミー（熱流束勾配チェック用）
  dz_dummy = fill(1.0e-4, 5)

  has_anomaly, anomalies = detect_numerical_anomalies(
    q, "熱流束",
    iteration=iteration,
    timestep=timestep,
    flux_range=(-1e7, 1e7),
    dx=dx, dy=dy, dz=dz_dummy
  )

  if has_anomaly
    return false, "異常検出: $(length(anomalies))件"
  else
    return true, "正常"
  end
end

"""
    check_adjoint_field(lambda_field; kwargs...) -> Tuple{Bool, String}

随伴場の包括的チェック（ラッパー関数）

# 引数
- `lambda_field::Array{Float64}`: 随伴場

# キーワード引数
- `iteration::Union{Int, Nothing}`: CGM反復回数
- `timestep::Union{Int, Nothing}`: 時間ステップ
- `dx::Float64`: x方向格子間隔 [m]
- `dy::Float64`: y方向格子間隔 [m]
- `dz::Union{Vector{Float64}, Nothing}`: z方向格子間隔配列 [m]

# 戻り値
- `is_valid::Bool`: 異常検出時false
- `status_msg::String`: 状態メッセージ
"""
function check_adjoint_field(
  lambda_field::Array{Float64};
  iteration::Union{Int, Nothing}=nothing,
  timestep::Union{Int, Nothing}=nothing,
  dx::Float64=0.12e-3,
  dy::Float64=0.12e-3*0.866,
  dz::Union{Vector{Float64}, Nothing}=nothing
)::Tuple{Bool, String}
  has_anomaly, anomalies = detect_numerical_anomalies(
    lambda_field, "lambda",
    iteration=iteration,
    timestep=timestep,
    dx=dx, dy=dy, dz=dz
  )

  if has_anomaly
    return false, "異常検出: $(length(anomalies))件"
  else
    return true, "正常"
  end
end

end  # module Validators
