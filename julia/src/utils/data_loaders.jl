"""
data_loaders.jl

Phase 6 C-1: 実データ読込機能
IRカメラからのMATLABファイル読み込みとデータ処理

対応するPython処理:
- extract_sorted_mat_files() (58-87行)
- load_region_temperature() (90-135行)
"""

module MATDataLoaders

using MAT
using Logging

export extract_sorted_mat_files, load_region_temperature


"""
    extract_sorted_mat_files(
      folder_path::String;
      prefix::String="SUS",
      extension::String=".MAT"
    ) -> Vector{String}

指定ディレクトリからMATファイルリストを自然順ソートで取得

# 引数
- `folder_path`: MATファイルが格納されているディレクトリパス
- `prefix`: ファイル名の接頭辞（デフォルト: "SUS"）
- `extension`: ファイル拡張子（デフォルト: ".MAT"）

# 戻り値
自然順（数値順）でソートされたファイル名の配列

# エラー処理
- ディレクトリが存在しない → ArgumentError
- 条件に一致するMATファイルが0個 → ArgumentError

# 実装詳細
正規表現パターン: "prefix(数値)extension"
例: "SUS1.MAT", "SUS2.MAT", "SUS10.MAT", "SUS20.MAT"
辞書順ソートでは "SUS10.MAT" < "SUS2.MAT" となるが、
自然順ソートでは数値部分を整数として比較する

# Python版との対応
Python版の extract_sorted_mat_files() と同等
"""
function extract_sorted_mat_files(
  folder_path::String;
  prefix::String="SUS",
  extension::String=".MAT"
)::Vector{String}
  # ディレクトリ存在確認
  if !isdir(folder_path)
    throw(ArgumentError("Directory not found: $folder_path"))
  end

  # ファイル一覧取得
  all_files = readdir(folder_path)

  # 正規表現パターン構築
  # raw文字列を使ってバックスラッシュのエスケープを回避
  pattern = Regex(prefix * raw"(\d+)" * extension)

  # マッチするファイルと数値部分を抽出
  matched_files = Tuple{String, Int}[]
  for file in all_files
    m = match(pattern, file)
    if m !== nothing
      # マッチした場合、ファイル名と数値部分を保存
      file_number = parse(Int, m.captures[1])
      push!(matched_files, (file, file_number))
    end
  end

  # マッチするファイルが存在しない場合はエラー
  if isempty(matched_files)
    throw(ArgumentError("No MAT files found with pattern $(prefix)\\d+$(extension) in: $folder_path"))
  end

  # 数値部分でソート（自然順ソート）
  sort!(matched_files, by = x -> x[2])

  # ファイル名のみを抽出して返す
  return [file for (file, _) in matched_files]
end


"""
    load_region_temperature(
      folder_path::String,
      i_range::Tuple{Int,Int},
      j_range::Tuple{Int,Int}
    ) -> Tuple{Array{Float64,3}, Int}

IRカメラMATLABファイルから指定領域の温度データを読み込む

# 引数
- `folder_path`: MATファイルが格納されているディレクトリパス
- `i_range`: y方向（行方向）の範囲 (start, end) - 1-based index
- `j_range`: x方向（列方向）の範囲 (start, end) - 1-based index

# 戻り値
- `region_data`: 3次元温度配列 [i方向, j方向, 時間フレーム] (Float64)
- `num_frames`: 読み込んだフレーム数 (Int)

# 座標系
- i方向: y軸（行方向、上下）
- j方向: x軸（列方向、左右）
- Juliaは1-based indexing: i_range=(3,7) → arr[3:7]（5要素）

# データ構造
各MATファイルには変数名=ファイル名（拡張子なし）でデータが格納されている
例: SUS1.MAT → mat_data["SUS1"] → 2次元温度配列 [i, j]

# エラー処理
- ディレクトリが存在しない → ArgumentError（extract_sorted_mat_filesで検出）
- 範囲外のインデックス → BoundsError

# データ検証
- NaN/Infが検出された場合は警告を出力（@warn）

# Python版との対応
Python版の load_region_temperature() と同等
- Python: 0-based indexing, Julia: 1-based indexing
- Python: npy配列保存, Julia: 3次元配列返却
"""
function load_region_temperature(
  folder_path::String,
  i_range::Tuple{Int,Int},
  j_range::Tuple{Int,Int}
)::Tuple{Array{Float64,3}, Int}
  # MATファイルリストを自然順ソートで取得
  mat_files = extract_sorted_mat_files(folder_path)
  num_frames = length(mat_files)

  # 領域サイズ計算
  i_size = i_range[2] - i_range[1] + 1
  j_size = j_range[2] - j_range[1] + 1

  # 3次元配列を事前アロケーション
  region_data = zeros(Float64, i_size, j_size, num_frames)

  # 各MATファイルを順次読み込み
  for (frame_idx, mat_file) in enumerate(mat_files)
    # フルパス構築
    mat_path = joinpath(folder_path, mat_file)

    # MATファイル読み込み
    mat_data = matread(mat_path)

    # 変数名 = ファイル名（拡張子なし）
    # 例: "SUS1.MAT" → "SUS1"
    variable_name = splitext(mat_file)[1]

    # 2次元温度配列を取得
    if !haskey(mat_data, variable_name)
      throw(ArgumentError("Variable '$variable_name' not found in $mat_file"))
    end

    frame_data = mat_data[variable_name]

    # 型チェックと次元チェック
    if !(frame_data isa AbstractMatrix)
      throw(ArgumentError("Variable '$variable_name' in $mat_file is not a 2D array"))
    end

    # 領域切り出し
    # Juliaは1-based indexing なので範囲そのまま使用
    try
      region_slice = frame_data[i_range[1]:i_range[2], j_range[1]:j_range[2]]
      region_data[:, :, frame_idx] = region_slice
    catch e
      if e isa BoundsError
        throw(ArgumentError("Index out of bounds: i_range=$i_range, j_range=$j_range for $mat_file with size $(size(frame_data))"))
      else
        rethrow(e)
      end
    end

    # NaN/Inf検出（警告のみ）
    if any(isnan, region_data[:, :, frame_idx])
      @warn "NaN detected in frame $frame_idx ($mat_file)"
    end
    if any(isinf, region_data[:, :, frame_idx])
      @warn "Inf detected in frame $frame_idx ($mat_file)"
    end
  end

  return (region_data, num_frames)
end

end # module MATDataLoaders
