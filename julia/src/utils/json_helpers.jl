"""
json_helpers.jl

JSON参照データの型変換ヘルパー関数

JSON.parsefile()が返すネストされたVector{Any}を
Julia型配列（Float64, Int等）に変換する汎用関数群

Phase 2以降のテストで使用
"""

module JSONHelpers

export json_to_array, flatten_nested_array, reshape_json_array


"""
    json_to_array(json_data, target_type::Type{T}) where T -> Array{T}

JSONネスト配列を指定型の多次元配列に変換

Args:
  json_data: JSON.parsefile()が返すVector{Any}
  target_type: 変換先の型（Float64, Int等）

Returns:
  Array{T}: 変換された配列（次元数はデータから自動推定）

Example:
  ```julia
  data = JSON.parsefile("data.json")
  T_initial = json_to_array(data["T_initial"], Float64)
  ```
"""
function json_to_array(json_data, target_type::Type{T}) where T
  # 平坦化
  flat = flatten_nested_array(json_data)

  # 型変換
  return convert(Vector{T}, flat)
end


"""
    flatten_nested_array(arr) -> Vector

再帰的にネスト配列を平坦化

JSONからパースされたネストされたVector{Any}を1次元Vectorに変換
後でreshape()を使用して多次元配列に復元する

Args:
  arr: ネストされた配列（Vector{Any}, Number等）

Returns:
  Vector: 平坦化された1次元配列

Example:
  ```julia
  nested = [[[1, 2], [3, 4]], [[5, 6], [7, 8]]]
  flat = flatten_nested_array(nested)  # [1, 2, 3, 4, 5, 6, 7, 8]
  ```
"""
function flatten_nested_array(arr)
  result = []
  _flatten_recursive!(result, arr)
  return result
end


"""
    _flatten_recursive!(result, arr)

内部再帰関数: 配列を再帰的に平坦化

Args:
  result: 結果を格納するVector（可変）
  arr: 処理対象の配列またはスカラー
"""
function _flatten_recursive!(result, arr)
  if arr isa Number
    # スカラー: そのまま追加
    push!(result, arr)
  elseif arr isa AbstractVector
    # ベクトル: 再帰的に処理
    for elem in arr
      _flatten_recursive!(result, elem)
    end
  else
    error("Unsupported array type: $(typeof(arr))")
  end
end


"""
    reshape_json_array(json_data, dims::Tuple, target_type::Type{T}) where T -> Array{T}

JSON配列を指定次元・型の配列に変換

Pythonとの互換性のため、C順序（行優先）からFortran順序（列優先）に変換
1. JSONデータを平坦化
2. 型変換
3. 次元を逆順でreshape（Pythonの(nt, ni, nj, nk) → Juliaの(nk, nj, ni, nt)）
4. 次元を転置して元の順序に戻す（permutedims）

Args:
  json_data: JSON.parsefile()が返すデータ
  dims: 目標次元 (ni, nj, nk) または (nt, ni, nj, nk)
  target_type: 変換先の型

Returns:
  Array{T}: 指定次元の配列（Julia順序）

Example:
  ```julia
  data = JSON.parsefile("data.json")
  # Python: (5, 5, 10) → Julia: (5, 5, 10)
  T_initial = reshape_json_array(data["T_initial"], (5, 5, 10), Float64)
  ```
"""
function reshape_json_array(json_data, dims::Tuple, target_type::Type{T}) where T
  flat = json_to_array(json_data, target_type)

  # Pythonは行優先（C順序）、Juliaは列優先（Fortran順序）
  # 次元を逆順でreshapeし、元の順序に転置
  ndims_total = length(dims)

  # 次元を逆順にしてreshape（Pythonと同じメモリレイアウトを再現）
  dims_reversed = reverse(dims)
  arr_reversed = reshape(flat, dims_reversed)

  # 次元を転置して元の順序に戻す
  perm = collect(ndims_total:-1:1)  # [ndims, ndims-1, ..., 2, 1]
  arr_transposed = permutedims(arr_reversed, perm)

  return arr_transposed
end


"""
    safe_convert(value, target_type::Type{T}) where T -> T

安全な型変換（エラーハンドリング付き）

Args:
  value: 変換元の値
  target_type: 変換先の型

Returns:
  T: 変換された値
"""
function safe_convert(value, target_type::Type{T}) where T
  try
    if value isa Number
      return T(value)
    elseif value isa AbstractVector
      return T[safe_convert(x, T) for x in value]
    else
      error("Cannot convert $(typeof(value)) to $(target_type)")
    end
  catch e
    error("Type conversion failed: $(e)")
  end
end


end  # module JSONHelpers