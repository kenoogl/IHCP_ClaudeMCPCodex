"""
main.jl

IHCP CGMソルバー メインエントリポイント

コマンドライン引数処理と実行モード制御
- full: 全時間領域スライディングウィンドウCGM計算
- test: Phase 1-5統合テスト実行
- benchmark: 性能測定

使用方法:
  julia --project=. src/main.jl --config config.toml --output results.jld2 [OPTIONS]
"""

using ArgParse
using ProgressMeter
using Logging
using Printf
using LinearAlgebra
using Dates
using Pkg
using Statistics

# プロジェクトモジュールのインクルード
include("IHCP_CGM.jl")
using .IHCP_CGM

include("utils/config.jl")
include("utils/io.jl")

"""
    parse_commandline() -> Dict

コマンドライン引数をパース
"""
function parse_commandline()
    s = ArgParseSettings(
        description = "IHCP CGM Solver - 逆熱伝導問題 共役勾配法ソルバー",
        version = "0.5.0",
        add_version = true
    )

    @add_arg_table! s begin
        "--config", "-c"
            help = "設定ファイル（TOML形式）[必須]"
            required = true
            arg_type = String

        "--output", "-o"
            help = "出力ファイル（JLD2形式）[必須]"
            required = true
            arg_type = String

        "--mode", "-m"
            help = "実行モード: full/test/benchmark [デフォルト: full]"
            default = "full"
            arg_type = String

        "--threads", "-t"
            help = "スレッド数 [デフォルト: 自動]"
            arg_type = Int
            default = nothing

        "--verbose", "-v"
            help = "詳細ログ出力"
            action = :store_true

        "--quiet", "-q"
            help = "最小限のログのみ"
            action = :store_true

        "--dry-run"
            help = "設定検証のみ（計算は実行しない）"
            action = :store_true
    end

    return parse_args(s)
end

"""
    setup_logging(args::Dict)

ログレベルを設定
"""
function setup_logging(args::Dict)
    if args["verbose"]
        global_logger(ConsoleLogger(stdout, Logging.Debug))
    elseif args["quiet"]
        global_logger(ConsoleLogger(stdout, Logging.Warn))
    else
        global_logger(ConsoleLogger(stdout, Logging.Info))
    end
end

"""
    setup_threads(args::Dict)

スレッド数を設定
"""
function setup_threads(args::Dict)
    if args["threads"] !== nothing
        BLAS.set_num_threads(args["threads"])
        @info "BLASスレッド数: $(args["threads"])"
    else
        @info "BLASスレッド数: $(BLAS.get_num_threads()) (自動)"
    end
end

"""
    run_full_calculation(config::Dict, output_file::String) -> Int

全時間領域スライディングウィンドウCGM計算を実行

# 引数
- `config::Dict`: 設定情報
- `output_file::String`: 出力ファイルパス

# 戻り値
- `Int`: 終了コード（0=成功、1=エラー）
"""
function run_full_calculation(config::Dict, output_file::String)
    println("\n" * "="^60)
    println("全時間領域計算モード")
    println("="^60)

    try
        # パラメータ抽出
        prob = config["problem"]
        mat = config["material"]
        cgm = config["cgm"]
        sw = get(config, "sliding_window", Dict(
            "window_size" => prob["nt"],
            "overlap" => 0,
            "q_init_value" => 0.0
        ))
        input = config["input"]

        # 格子幅配列の準備
        @info "格子設定を読み込み中..."

        # z方向格子幅
        if haskey(prob, "dz_file")
            dz = load_dz_grid(prob["dz_file"])
        else
            dz = Float64.(prob["dz"])
        end

        if length(dz) != prob["nk"]
            error("dz配列の長さがnkと一致しません: $(length(dz)) != $(prob["nk"])")
        end

        # 表面・底面の格子幅（nk要素のベクトル）
        # dz配列（nk要素、各層の厚さ）からセル中心座標を計算し、dz_b、dz_tを構築
        nk = prob["nk"]

        # セル境界座標（z_faces）を計算（底面から表面へ）
        # nk個のセルがあるので、nk+1個の境界
        z_faces = zeros(nk + 1)
        z_faces[1] = 0.0  # 底面
        for k in 1:nk
            z_faces[k+1] = z_faces[k] + dz[k]
        end

        # セル中心座標（z_centers）を計算
        # 底面と表面のセルは境界に重なる（半分のCV）
        # 中間のセルは通常の中心
        z_centers = zeros(nk)
        z_centers[1] = z_faces[1]      # 底面セル中心=底面境界
        z_centers[end] = z_faces[end]  # 表面セル中心=表面境界
        for k in 2:(nk-1)
            z_centers[k] = (z_faces[k] + z_faces[k+1]) / 2.0
        end

        # dz_t: 上側境界からの距離
        dz_t = zeros(nk)
        dz_t[end] = Inf  # 表面境界条件
        for k in 1:(nk-1)
            dz_t[k] = z_centers[k+1] - z_centers[k]
        end

        # dz_b: 下側境界からの距離
        dz_b = zeros(nk)
        dz_b[1] = Inf  # 底面境界条件
        for k in 2:nk
            dz_b[k] = z_centers[k] - z_centers[k-1]
        end

        @info @sprintf("  格子点数: (%d, %d, %d)", prob["ni"], prob["nj"], prob["nk"])
        @info @sprintf("  格子幅: dx=%.2e, dy=%.2e", prob["dx"], prob["dy"])
        @info @sprintf("  z方向格子幅: %d層, 表面=%.2e, 底面=%.2e", length(dz), dz[1], dz[end])

        # 入力データ読み込み
        @info "\n入力データ読み込み中..."

        # 初期温度場
        if haskey(input, "T_init_file")
            T_init = load_array(input["T_init_file"])
        else
            T_init = input["T_init"]
        end

        # 観測温度
        if haskey(input, "Y_obs_file")
            Y_obs = load_array(input["Y_obs_file"])
        else
            Y_obs = input["Y_obs"]
        end

        @info @sprintf("  初期温度場形状（読込時）: %s", size(T_init))
        @info @sprintf("  観測温度形状（読込時）: %s", size(Y_obs))

        # データ形状の調整
        # ソルバーの期待形状: T_init=(ni, nj, nk), Y_obs=(nt, ni, nj)
        # NumPyファイルからの読み込み形状を確認し、必要に応じて転置

        # T_initの検証とreshape
        if size(T_init) != (prob["ni"], prob["nj"], prob["nk"])
            @warn @sprintf("初期温度場の形状が期待と異なります: %s != (%d, %d, %d)",
                          size(T_init), prob["ni"], prob["nj"], prob["nk"])
        end

        # Y_obsの検証とreshape
        # 読み込まれた形状が (ni, nj, nt) の場合、(nt, ni, nj) に変換
        if ndims(Y_obs) == 3
            if size(Y_obs, 1) == prob["ni"] && size(Y_obs, 2) == prob["nj"]
                # (ni, nj, nt) → (nt, ni, nj)
                Y_obs = permutedims(Y_obs, (3, 1, 2))
                @info "  観測温度を (nt, ni, nj) 形式に変換しました"
            end
        end

        @info @sprintf("  初期温度場形状（最終）: %s", size(T_init))
        @info @sprintf("  観測温度形状（最終）: %s", size(Y_obs))

        # スライディングウィンドウCGM実行
        println("\n" * "-"^60)
        println("スライディングウィンドウCGM実行")
        println("-"^60)

        window_size = sw["window_size"]
        overlap = sw["overlap"]
        q_init_value = sw["q_init_value"]
        max_iter = cgm["max_iter"]

        @info @sprintf("  ウィンドウサイズ: %d", window_size)
        @info @sprintf("  オーバーラップ: %d", overlap)
        @info @sprintf("  初期熱流束: %.2e W/m²", q_init_value)
        @info @sprintf("  CGM最大反復回数: %d", max_iter)

        # ウィンドウ数の計算
        nt_total = size(Y_obs, 3)
        n_windows = ceil(Int, (nt_total - overlap) / (window_size - overlap))
        @info @sprintf("  総時間ステップ数: %d", nt_total)
        @info @sprintf("  ウィンドウ数: %d", n_windows)

        # メモリ使用量の推定
        element_size = sizeof(Float64)
        spatial_size = prob["ni"] * prob["nj"]
        q_memory_mb = spatial_size * nt_total * element_size / 1e6
        @info @sprintf("  推定メモリ使用量（q配列）: %.2f MB", q_memory_mb)

        # ディスク空き容量チェック
        required_disk_mb = q_memory_mb * 2  # 余裕を持って2倍
        check_disk_space(required_disk_mb)

        # 計算開始
        println("\n計算開始: $(now())")
        start_time = time()

        q_global, windows_info = solve_sliding_window_cgm(
            Y_obs, T_init,
            prob["dx"], prob["dy"], dz, dz_b, dz_t, prob["dt"],
            mat["rho"], mat["cp_coeffs"], mat["k_coeffs"],
            window_size, overlap, q_init_value, max_iter
        )

        elapsed_time = time() - start_time
        println("計算完了: $(now())")
        @info @sprintf("  計算時間: %.2f 秒 (%.2f 分)", elapsed_time, elapsed_time / 60)

        # 結果の統計情報
        println("\n" * "-"^60)
        println("結果統計")
        println("-"^60)
        @info @sprintf("  熱流束配列形状: %s", size(q_global))
        @info @sprintf("  熱流束範囲: [%.2e, %.2e] W/m²", minimum(q_global), maximum(q_global))
        @info @sprintf("  熱流束平均: %.2e W/m²", mean(q_global))
        @info @sprintf("  ウィンドウ情報数: %d", length(windows_info))

        # 結果保存
        println("\n" * "-"^60)
        save_results(output_file, q_global, windows_info, config)
        println("-"^60)

        println("\n✓ 全時間領域計算が正常に完了しました")
        println("="^60)

        return 0

    catch e
        @error "計算中にエラーが発生しました" exception=(e, catch_backtrace())
        return 1
    end
end

"""
    run_test_mode(config::Dict) -> Int

テストモード: Phase 1-5統合テスト実行

# 引数
- `config::Dict`: 設定情報

# 戻り値
- `Int`: 終了コード（0=成功、1=エラー）
"""
function run_test_mode(config::Dict)
    println("\n" * "="^60)
    println("テストモード: Phase 1-5統合テスト実行")
    println("="^60)

    try
        @info "テストスイートを実行中..."
        Pkg.test("IHCP_CGM")

        println("\n✓ 全テスト合格")
        println("="^60)

        return 0

    catch e
        @error "テスト実行中にエラーが発生しました" exception=(e, catch_backtrace())
        return 1
    end
end

"""
    run_benchmark_mode(config::Dict) -> Int

ベンチマークモード: 性能測定

# 引数
- `config::Dict`: 設定情報

# 戻り値
- `Int`: 終了コード（0=成功、1=エラー）
"""
function run_benchmark_mode(config::Dict)
    println("\n" * "="^60)
    println("ベンチマークモード")
    println("="^60)

    try
        prob = config["problem"]

        # 小規模問題設定
        ni, nj, nk = 5, 5, 10
        nt = 100

        @info "小規模問題でベンチマーク実行中..."
        @info @sprintf("  格子: (%d, %d, %d), 時間ステップ: %d", ni, nj, nk, nt)

        # Phase 1: 熱物性値計算
        println("\nPhase 1: 熱物性値計算")
        T_test = 300.0 .+ 100.0 .* rand(ni, nj, nk)
        cp_coeffs = Float64.(config["material"]["cp_coeffs"])
        k_coeffs = Float64.(config["material"]["k_coeffs"])

        # 簡易ベンチマーク（実行時間計測）
        start_time = time()
        for i in 1:100
            thermal_properties_calculator(T_test, cp_coeffs, k_coeffs)
        end
        elapsed = (time() - start_time) / 100
        println(@sprintf("  平均実行時間: %.6f 秒", elapsed))

        # Phase 2-4: 今後追加
        println("\nPhase 2-4: 統合ベンチマークは今後実装予定")

        println("\n✓ ベンチマーク完了")
        println("="^60)

        return 0

    catch e
        @error "ベンチマーク実行中にエラーが発生しました" exception=(e, catch_backtrace())
        return 1
    end
end

"""
    main() -> Int

メインエントリポイント

# 戻り値
- `Int`: 終了コード（0=成功、1=エラー）
"""
function main()
    # コマンドライン引数パース
    args = parse_commandline()

    # ログ設定
    setup_logging(args)

    # バージョン情報表示
    println("\n" * "="^60)
    println("IHCP CGM Solver - 逆熱伝導問題 共役勾配法ソルバー")
    println("="^60)
    println("Version: $(IHCP_CGM.VERSION)")
    println("Phase: $(IHCP_CGM.PHASE)")
    println("Julia Version: $(VERSION)")
    println("実行日時: $(now())")
    println("="^60)

    # スレッド設定
    setup_threads(args)

    # 設定ファイル読み込み
    @info "\n設定ファイル読み込み中: $(args["config"])"

    local config
    try
        config = load_config(args["config"])
        @info "✓ 設定ファイル読み込み完了"
    catch e
        @error "設定ファイルの読み込みに失敗しました" exception=(e, catch_backtrace())
        return 1
    end

    # 設定検証
    @info "\n設定を検証中..."

    try
        validate_config(config)
        @info "✓ 設定検証完了"
    catch e
        @error "設定検証に失敗しました" exception=(e, catch_backtrace())
        return 1
    end

    # 設定サマリー表示
    print_config_summary(config)

    # Dry-runモード
    if args["dry-run"]
        println("\n" * "="^60)
        println("Dry-runモード: 設定検証のみ完了")
        println("="^60)
        @info "実際の計算は実行されませんでした"
        return 0
    end

    # 実行モード分岐
    mode = args["mode"]
    @info "\n実行モード: $mode"

    if mode == "full"
        return run_full_calculation(config, args["output"])

    elseif mode == "test"
        return run_test_mode(config)

    elseif mode == "benchmark"
        return run_benchmark_mode(config)

    else
        @error "不明な実行モード: $mode"
        @error "有効なモード: full, test, benchmark"
        return 1
    end
end

# エントリポイント
if abspath(PROGRAM_FILE) == @__FILE__
    exit_code = main()
    exit(exit_code)
end
