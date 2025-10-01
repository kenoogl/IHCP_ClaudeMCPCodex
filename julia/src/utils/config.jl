"""
config.jl

TOML形式の設定ファイル読み込みと検証

機能:
- TOML設定ファイルのパース
- 必須項目の検証
- 数値範囲チェック
- ファイルパス存在確認
"""

using TOML
using Printf

"""
    load_config(config_file::String) -> Dict

TOML形式の設定ファイルを読み込む

# 引数
- `config_file::String`: 設定ファイルパス

# 戻り値
- `Dict`: 設定情報の辞書

# 例外
- ファイルが存在しない場合エラー
- TOMLパースエラー
"""
function load_config(config_file::String)
    if !isfile(config_file)
        error("設定ファイルが見つかりません: $config_file")
    end

    try
        config = TOML.parsefile(config_file)
        return config
    catch e
        error("設定ファイルの読み込みに失敗しました: $e")
    end
end

"""
    validate_config(config::Dict) -> Bool

設定ファイルの内容を検証

# 引数
- `config::Dict`: 設定情報の辞書

# 戻り値
- `Bool`: 検証成功時true

# 検証項目
- 必須セクションの存在確認
- 数値範囲チェック
- ファイルパス存在確認
"""
function validate_config(config::Dict)
    # 必須セクションの確認
    required_sections = ["problem", "material", "cgm", "input", "output"]
    for section in required_sections
        if !haskey(config, section)
            error("必須セクション '$section' が設定ファイルにありません")
        end
    end

    # problem セクションの検証
    prob = config["problem"]
    validate_problem_config(prob)

    # material セクションの検証
    mat = config["material"]
    validate_material_config(mat)

    # cgm セクションの検証
    cgm = config["cgm"]
    validate_cgm_config(cgm)

    # input セクションの検証
    input = config["input"]
    validate_input_config(input)

    # sliding_window セクションの検証（オプション）
    if haskey(config, "sliding_window")
        sw = config["sliding_window"]
        validate_sliding_window_config(sw)
    end

    return true
end

"""
    validate_problem_config(prob::Dict)

problem セクションの検証
"""
function validate_problem_config(prob::Dict)
    # 必須パラメータ
    required_params = ["ni", "nj", "nk", "dx", "dy", "dt", "nt"]
    for param in required_params
        if !haskey(prob, param)
            error("problem.$param が設定されていません")
        end
    end

    # 格子点数の検証
    for key in ["ni", "nj", "nk", "nt"]
        val = prob[key]
        if val <= 0
            error("problem.$key は正の整数である必要があります: $val")
        end
    end

    # 格子幅の検証
    for key in ["dx", "dy", "dt"]
        val = prob[key]
        if val <= 0.0
            error("problem.$key は正の実数である必要があります: $val")
        end
    end

    # z方向格子の検証
    if haskey(prob, "dz_file")
        if !isfile(prob["dz_file"])
            error("dz_fileが見つかりません: $(prob["dz_file"])")
        end
    elseif haskey(prob, "dz")
        # dz配列の検証
        dz = prob["dz"]
        if length(dz) != prob["nk"]
            error("dz配列の長さがnkと一致しません: $(length(dz)) != $(prob["nk"])")
        end
        if any(dz .<= 0.0)
            error("dz配列の要素は正の実数である必要があります")
        end
    else
        error("problem.dz_file または problem.dz のいずれかが必要です")
    end
end

"""
    validate_material_config(mat::Dict)

material セクションの検証
"""
function validate_material_config(mat::Dict)
    # 必須パラメータ
    required_params = ["rho", "cp_coeffs", "k_coeffs"]
    for param in required_params
        if !haskey(mat, param)
            error("material.$param が設定されていません")
        end
    end

    # 密度の検証
    if mat["rho"] <= 0.0
        error("material.rho は正の実数である必要があります: $(mat["rho"])")
    end

    # 係数配列の検証
    for key in ["cp_coeffs", "k_coeffs"]
        coeffs = mat[key]
        if !isa(coeffs, Vector)
            error("material.$key は配列である必要があります")
        end
        if length(coeffs) != 4
            error("material.$key は4要素の配列である必要があります: $(length(coeffs))")
        end
    end
end

"""
    validate_cgm_config(cgm::Dict)

cgm セクションの検証
"""
function validate_cgm_config(cgm::Dict)
    # 必須パラメータ
    required_params = ["max_iter", "sigma"]
    for param in required_params
        if !haskey(cgm, param)
            error("cgm.$param が設定されていません")
        end
    end

    # 反復回数の検証
    if cgm["max_iter"] <= 0
        error("cgm.max_iter は正の整数である必要があります: $(cgm["max_iter"])")
    end

    # 測定誤差の検証
    if cgm["sigma"] <= 0.0
        error("cgm.sigma は正の実数である必要があります: $(cgm["sigma"])")
    end

    # オプショナルパラメータのデフォルト値設定
    if !haskey(cgm, "rtol_dhcp")
        cgm["rtol_dhcp"] = 1.0e-6
    end
    if !haskey(cgm, "rtol_adjoint")
        cgm["rtol_adjoint"] = 1.0e-8
    end
    if !haskey(cgm, "maxiter_cg")
        cgm["maxiter_cg"] = 20000
    end
end

"""
    validate_input_config(input::Dict)

input セクションの検証
"""
function validate_input_config(input::Dict)
    # T_init_file または T_init 配列が必要
    has_init_file = haskey(input, "T_init_file")
    has_init_array = haskey(input, "T_init")

    if !has_init_file && !has_init_array
        error("input.T_init_file または input.T_init のいずれかが必要です")
    end

    # ファイルパスの検証
    if has_init_file
        if !isfile(input["T_init_file"])
            @warn "初期温度ファイルが見つかりません: $(input["T_init_file"])"
        end
    end

    # Y_obs_file または Y_obs 配列が必要
    has_obs_file = haskey(input, "Y_obs_file")
    has_obs_array = haskey(input, "Y_obs")

    if !has_obs_file && !has_obs_array
        error("input.Y_obs_file または input.Y_obs のいずれかが必要です")
    end

    # ファイルパスの検証
    if has_obs_file
        if !isfile(input["Y_obs_file"])
            @warn "観測温度ファイルが見つかりません: $(input["Y_obs_file"])"
        end
    end

    # MATファイルディレクトリの検証（オプション）
    if haskey(input, "mat_files_dir")
        if !isdir(input["mat_files_dir"])
            @warn "MATファイルディレクトリが見つかりません: $(input["mat_files_dir"])"
        end
    end
end

"""
    validate_sliding_window_config(sw::Dict)

sliding_window セクションの検証
"""
function validate_sliding_window_config(sw::Dict)
    # 必須パラメータ
    required_params = ["window_size", "overlap"]
    for param in required_params
        if !haskey(sw, param)
            error("sliding_window.$param が設定されていません")
        end
    end

    # ウィンドウサイズの検証
    if sw["window_size"] <= 0
        error("sliding_window.window_size は正の整数である必要があります: $(sw["window_size"])")
    end

    # オーバーラップの検証
    if sw["overlap"] < 0
        error("sliding_window.overlap は非負の整数である必要があります: $(sw["overlap"])")
    end

    if sw["overlap"] >= sw["window_size"]
        error("sliding_window.overlap はwindow_sizeより小さい必要があります: $(sw["overlap"]) >= $(sw["window_size"])")
    end

    # 初期熱流束のデフォルト値
    if !haskey(sw, "q_init_value")
        sw["q_init_value"] = 0.0
    end
end

"""
    print_config_summary(config::Dict)

設定情報のサマリーを表示
"""
function print_config_summary(config::Dict)
    println("\n" * "="^60)
    println("設定情報サマリー")
    println("="^60)

    # Problem設定
    prob = config["problem"]
    println("\n[Problem設定]")
    println(@sprintf("  格子点数: ni=%d, nj=%d, nk=%d", prob["ni"], prob["nj"], prob["nk"]))
    println(@sprintf("  格子幅: dx=%.2e m, dy=%.2e m", prob["dx"], prob["dy"]))
    println(@sprintf("  時間設定: nt=%d, dt=%.2e s", prob["nt"], prob["dt"]))

    # Material設定
    mat = config["material"]
    println("\n[Material設定]")
    println(@sprintf("  密度: %.1f kg/m³", mat["rho"]))
    println("  比熱係数: ", mat["cp_coeffs"])
    println("  熱伝導率係数: ", mat["k_coeffs"])

    # CGM設定
    cgm = config["cgm"]
    println("\n[CGM設定]")
    println(@sprintf("  最大反復回数: %d", cgm["max_iter"]))
    println(@sprintf("  測定誤差σ: %.2e K", cgm["sigma"]))

    # Sliding Window設定
    if haskey(config, "sliding_window")
        sw = config["sliding_window"]
        println("\n[Sliding Window設定]")
        println(@sprintf("  ウィンドウサイズ: %d", sw["window_size"]))
        println(@sprintf("  オーバーラップ: %d", sw["overlap"]))
        println(@sprintf("  初期熱流束: %.2e W/m²", sw["q_init_value"]))
    end

    println("\n" * "="^60)
end
