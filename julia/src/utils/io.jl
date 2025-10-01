"""
io.jl

入出力ユーティリティ

機能:
- NumPy配列（.npy）の読み込み
- JLD2形式での結果保存
- MATLAB MAT形式の読み込み
- テキストファイル読み込み
"""

using NPZ
using JLD2
using MAT
using DelimitedFiles
using Printf
using Logging
using Dates
using JSON

"""
    load_array(file_path::String) -> Array

配列データをファイルから読み込む

サポート形式:
- .npy: NumPy配列
- .txt: テキスト形式
- .mat: MATLAB形式

# 引数
- `file_path::String`: ファイルパス

# 戻り値
- `Array`: 読み込んだ配列

# 例外
- ファイルが存在しない場合エラー
- サポートされていない形式の場合エラー
"""
function load_array(file_path::String)
    if !isfile(file_path)
        error("ファイルが見つかりません: $file_path")
    end

    ext = splitext(file_path)[2]

    if ext == ".npy"
        @info "NumPy配列を読み込み中: $file_path"
        data = npzread(file_path)
        @info @sprintf("  形状: %s, 型: %s", size(data), eltype(data))
        return data

    elseif ext == ".txt"
        @info "テキストファイルを読み込み中: $file_path"
        data = readdlm(file_path)
        @info @sprintf("  形状: %s, 型: %s", size(data), eltype(data))
        return data

    elseif ext == ".mat"
        @info "MATLABファイルを読み込み中: $file_path"
        mat_data = matread(file_path)
        # 最初の変数を返す（通常1変数のみ）
        key = first(keys(mat_data))
        data = mat_data[key]
        @info @sprintf("  変数名: %s, 形状: %s, 型: %s", key, size(data), eltype(data))
        return data

    else
        error("サポートされていないファイル形式: $ext")
    end
end

"""
    load_dz_grid(file_path::String) -> Vector{Float64}

z方向格子幅をファイルから読み込む

# 引数
- `file_path::String`: ファイルパス（.txtまたは.npy）

# 戻り値
- `Vector{Float64}`: 格子幅配列
"""
function load_dz_grid(file_path::String)
    data = load_array(file_path)

    # 1次元ベクトルに変換
    if ndims(data) == 2 && size(data, 2) == 1
        data = vec(data)
    elseif ndims(data) == 2 && size(data, 1) == 1
        data = vec(data)
    end

    if ndims(data) != 1
        error("dz格子は1次元配列である必要があります: $(size(data))")
    end

    return Float64.(data)
end

"""
    save_results(output_file::String, q_global, windows_info, config::Dict)

計算結果をJLD2形式で保存

# 引数
- `output_file::String`: 出力ファイルパス
- `q_global`: 全時間領域の熱流束配列
- `windows_info`: ウィンドウ情報配列
- `config::Dict`: 設定情報

# 保存内容
- q_global: 推定された表面熱流束
- windows_info: スライディングウィンドウの詳細情報
- config: 計算に使用した設定
- metadata: メタデータ（バージョン、日時等）
"""
function save_results(output_file::String, q_global, windows_info, config::Dict)
    @info "結果を保存中: $output_file"

    # メタデータ作成
    metadata = Dict(
        "version" => "IHCP_CGM.jl v0.5.0",
        "timestamp" => string(now()),
        "julia_version" => string(VERSION),
        "description" => "逆熱伝導問題 共役勾配法 計算結果"
    )

    try
        jldopen(output_file, "w") do file
            # 主要結果
            file["q_global"] = q_global
            file["windows_info"] = windows_info

            # 設定情報
            file["config"] = config

            # メタデータ
            file["metadata"] = metadata
        end

        @info "✓ 保存完了"
        @info @sprintf("  ファイルサイズ: %.2f MB", filesize(output_file) / 1e6)

    catch e
        error("結果の保存に失敗しました: $e")
    end
end

"""
    save_results_npz(filename::String, q_global, windows_info, config::Dict)

NumPy NPZ形式で結果を保存（Python互換）

Pythonでの読込例:
```python
import numpy as np
data = np.load("results.npz", allow_pickle=True)
q = data["q_global"]
metadata = data["metadata"].item()
```

# 引数
- `filename::String`: 出力ファイルパス
- `q_global`: 全時間領域の熱流束配列
- `windows_info`: ウィンドウ情報配列
- `config::Dict`: 設定情報
"""
function save_results_npz(filename::String, q_global, windows_info, config::Dict)
    @info "NPZ形式で結果を保存中: $filename"

    # メタデータをJSON文字列化（バイト配列として保存）
    metadata = Dict(
        "version" => "IHCP_CGM.jl v0.5.0",
        "timestamp" => string(now()),
        "julia_version" => string(VERSION),
        "shape" => collect(size(q_global)),
        "description" => "逆熱伝導問題 共役勾配法 計算結果"
    )
    metadata_str = JSON.json(metadata)

    # ウィンドウ情報をJSON文字列化
    windows_json = [
        Dict(
            "window_id" => w.window_id,
            "start_idx" => w.start_idx,
            "end_idx" => w.end_idx,
            "n_iterations" => w.n_iterations,
            "J_final" => w.J_final,
            "converged" => w.converged
        ) for w in windows_info
    ]
    windows_str = JSON.json(windows_json)
    config_str = JSON.json(config)

    try
        # NPZ保存（文字列はバイト配列に変換）
        npzwrite(filename, Dict(
            "q_global" => q_global,
            "metadata" => collect(codeunits(metadata_str)),
            "windows_info" => collect(codeunits(windows_str)),
            "config" => collect(codeunits(config_str))
        ))

        @info "✓ NPZ形式で保存完了"
        @info @sprintf("  ファイルサイズ: %.2f MB", filesize(filename) / 1e6)

    catch e
        error("NPZ形式での保存に失敗しました: $e")
    end
end

"""
    load_results_npz(filename::String) -> Tuple

NPZ形式の結果ファイルを読み込む

# 引数
- `filename::String`: 入力ファイルパス

# 戻り値
- `Tuple`: (q_global, metadata, windows_info, config)
"""
function load_results_npz(filename::String)
    if !isfile(filename)
        error("ファイルが見つかりません: $filename")
    end

    @info "NPZ形式の結果ファイルを読み込み中: $filename"

    try
        data = npzread(filename)

        q_global = data["q_global"]

        # バイト配列から文字列に復元
        metadata = JSON.parse(String(UInt8.(data["metadata"])))
        windows_info_json = JSON.parse(String(UInt8.(data["windows_info"])))
        config = JSON.parse(String(UInt8.(data["config"])))

        @info "✓ NPZ形式で読み込み完了"
        @info "  バージョン: $(metadata["version"])"
        @info "  作成日時: $(metadata["timestamp"])"

        return q_global, metadata, windows_info_json, config

    catch e
        error("NPZ形式ファイルの読み込みに失敗しました: $e")
    end
end

"""
    load_results_jld2(filename::String) -> Tuple

JLD2形式の結果ファイルを読み込む

# 引数
- `filename::String`: 入力ファイルパス

# 戻り値
- `Tuple`: (q_global, metadata, windows_info, config)
"""
function load_results_jld2(filename::String)
    if !isfile(filename)
        error("ファイルが見つかりません: $filename")
    end

    @info "JLD2形式の結果ファイルを読み込み中: $filename"

    try
        results = load(filename)

        q_global = results["q_global"]
        metadata = results["metadata"]
        windows_info = results["windows_info"]
        config = results["config"]

        @info "✓ JLD2形式で読み込み完了"
        @info "  バージョン: $(metadata["version"])"
        @info "  作成日時: $(metadata["timestamp"])"

        return q_global, metadata, windows_info, config

    catch e
        error("JLD2形式ファイルの読み込みに失敗しました: $e")
    end
end

"""
    load_results(filename::String) -> Tuple

結果ファイルを読み込む（フォーマット自動判定）

サポート形式:
- .jld2: Julia JLD2形式
- .npz: NumPy NPZ形式

# 引数
- `filename::String`: 入力ファイルパス

# 戻り値
- `Tuple`: (q_global, metadata, windows_info, config)
"""
function load_results(filename::String)
    ext = lowercase(splitext(filename)[2])

    if ext == ".jld2"
        return load_results_jld2(filename)
    elseif ext == ".npz"
        return load_results_npz(filename)
    else
        error("Unsupported format: $ext (対応形式: .jld2, .npz)")
    end
end

"""
    save_array(file_path::String, data::Array)

配列をファイルに保存

サポート形式:
- .npy: NumPy配列
- .jld2: Julia JLD2形式

# 引数
- `file_path::String`: ファイルパス
- `data::Array`: 保存する配列
"""
function save_array(file_path::String, data::Array)
    ext = splitext(file_path)[2]

    @info "配列を保存中: $file_path"

    if ext == ".npy"
        npzwrite(file_path, data)
        @info @sprintf("  形状: %s, サイズ: %.2f MB", size(data), sizeof(data) / 1e6)

    elseif ext == ".jld2"
        jldopen(file_path, "w") do file
            file["data"] = data
        end
        @info @sprintf("  形状: %s, サイズ: %.2f MB", size(data), filesize(file_path) / 1e6)

    else
        error("サポートされていないファイル形式: $ext")
    end

    @info "✓ 保存完了"
end

"""
    print_file_info(file_path::String)

ファイル情報を表示
"""
function print_file_info(file_path::String)
    if !isfile(file_path)
        @warn "ファイルが見つかりません: $file_path"
        return
    end

    file_size = filesize(file_path)
    file_ext = splitext(file_path)[2]

    println(@sprintf("  ファイル: %s", basename(file_path)))
    println(@sprintf("  パス: %s", dirname(file_path)))
    println(@sprintf("  サイズ: %.2f MB", file_size / 1e6))
    println(@sprintf("  形式: %s", file_ext))
end

"""
    check_disk_space(required_mb::Float64) -> Bool

ディスク空き容量をチェック

# 引数
- `required_mb::Float64`: 必要な容量（MB）

# 戻り値
- `Bool`: 十分な空き容量がある場合true
"""
function check_disk_space(required_mb::Float64)
    # Unixシステムでdfコマンドを使用
    try
        df_output = read(`df -m .`, String)
        lines = split(df_output, '\n')
        if length(lines) >= 2
            fields = split(lines[2])
            available_mb = parse(Float64, fields[4])

            if available_mb < required_mb
                @warn @sprintf("ディスク空き容量不足: 必要 %.0f MB, 利用可能 %.0f MB", required_mb, available_mb)
                return false
            else
                @info @sprintf("ディスク空き容量: %.0f MB", available_mb)
                return true
            end
        end
    catch e
        @warn "ディスク空き容量の確認に失敗しました: $e"
    end

    return true  # チェックできない場合は続行
end
