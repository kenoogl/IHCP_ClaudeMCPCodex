"""
visualization.jl

可視化ユーティリティ

機能:
- 熱流束の時間変化プロット
- 空間分布ヒートマップ
- CGM収束履歴の可視化
"""

using Plots
using Printf
using Logging

"""
    plot_heat_flux(q_global, config; kwargs...)

熱流束の時間変化をプロット

# 引数
- `q_global`: 熱流束配列 (nt-1, ni, nj)
- `config::Dict`: 設定情報

# キーワード引数
- `save_path::String`: 保存先パス（指定時のみ保存）
- `format::Symbol`: 画像形式 (:png, :pdf, :svg)
- `point::Tuple{Int,Int}`: プロット対象の格子点座標 (i, j)（デフォルト: 中心点）

# 戻り値
- プロットオブジェクト
"""
function plot_heat_flux(q_global, config;
                        save_path=nothing,
                        format=:png,
                        point=nothing)
    # データサイズ取得
    nt_minus_1, ni, nj = size(q_global)
    nt = nt_minus_1 + 1

    # デフォルト: 中心点
    if point === nothing
        i, j = div(ni, 2) + 1, div(nj, 2) + 1
    else
        i, j = point
    end

    # 境界チェック
    if i < 1 || i > ni || j < 1 || j > nj
        error("無効な格子点座標: (i=$i, j=$j), 有効範囲: (1:$ni, 1:$nj)")
    end

    # 時間軸（ミリ秒）
    dt = config["problem"]["dt"]
    t_ms = (0:nt-2) .* dt .* 1000  # 0から始まり、nt-1ステップ分

    # 熱流束の時間変化
    q_series = q_global[:, i, j]

    @info @sprintf("熱流束の時間変化プロット: 格子点 (i=%d, j=%d)", i, j)
    @info @sprintf("  時間範囲: %.2f - %.2f ms", t_ms[1], t_ms[end])
    @info @sprintf("  熱流束範囲: %.2e - %.2e W/m²", minimum(q_series), maximum(q_series))

    # プロット作成
    p = plot(t_ms, q_series,
             xlabel="Time [ms]",
             ylabel="Heat Flux [W/m²]",
             title="Heat Flux at (i=$i, j=$j)",
             linewidth=2,
             legend=false,
             grid=true,
             size=(800, 500),
             dpi=150)

    # 保存
    if save_path !== nothing
        savefig(p, save_path)
        @info "✓ プロットを保存: $save_path"
    end

    return p
end

"""
    plot_flux_heatmap(q_global, config, time_index; kwargs...)

指定時刻での熱流束の空間分布をヒートマップ表示

# 引数
- `q_global`: 熱流束配列 (nt-1, ni, nj)
- `config::Dict`: 設定情報
- `time_index::Int`: 時間インデックス（1からnt-1）

# キーワード引数
- `save_path::String`: 保存先パス（指定時のみ保存）
- `format::Symbol`: 画像形式 (:png, :pdf, :svg)

# 戻り値
- プロットオブジェクト
"""
function plot_flux_heatmap(q_global, config, time_index;
                           save_path=nothing,
                           format=:png)
    # データサイズ取得
    nt_minus_1, ni, nj = size(q_global)

    # 境界チェック
    if time_index < 1 || time_index > nt_minus_1
        error("無効な時間インデックス: $time_index, 有効範囲: 1:$nt_minus_1")
    end

    # 指定時刻の熱流束
    q_snapshot = q_global[time_index, :, :]'  # 転置（表示用: y軸が縦、x軸が横）

    # 時刻計算
    dt = config["problem"]["dt"]
    t_ms = (time_index - 1) * dt * 1000

    @info @sprintf("熱流束ヒートマップ: t = %.2f ms (index=%d)", t_ms, time_index)
    @info @sprintf("  空間範囲: (%d, %d)", ni, nj)
    @info @sprintf("  熱流束範囲: %.2e - %.2e W/m²", minimum(q_snapshot), maximum(q_snapshot))

    # ヒートマップ作成
    p = heatmap(q_snapshot,
                xlabel="x grid index",
                ylabel="y grid index",
                title="Heat Flux at t=$(round(t_ms, digits=2)) ms",
                color=:viridis,
                aspect_ratio=:equal,
                colorbar_title="W/m²",
                size=(700, 600),
                dpi=150)

    # 保存
    if save_path !== nothing
        savefig(p, save_path)
        @info "✓ ヒートマップを保存: $save_path"
    end

    return p
end

"""
    plot_cgm_convergence(windows_info; kwargs...)

CGM収束履歴のプロット

# 引数
- `windows_info`: ウィンドウ情報配列

# キーワード引数
- `save_path::String`: 保存先パス（指定時のみ保存）
- `format::Symbol`: 画像形式 (:png, :pdf, :svg)

# 戻り値
- プロットオブジェクト
"""
function plot_cgm_convergence(windows_info;
                               save_path=nothing,
                               format=:png)
    # 各ウィンドウのCGM反復数と目的関数
    n_windows = length(windows_info)
    iterations = [w.n_iterations for w in windows_info]
    J_final = [w.J_final for w in windows_info]

    @info @sprintf("CGM収束履歴プロット: %d ウィンドウ", n_windows)
    @info @sprintf("  反復回数範囲: %d - %d", minimum(iterations), maximum(iterations))
    @info @sprintf("  目的関数範囲: %.2e - %.2e", minimum(J_final), maximum(J_final))

    # プロット1: CGM反復回数
    p1 = bar(1:n_windows, iterations,
             xlabel="Window Index",
             ylabel="CGM Iterations",
             title="CGM Convergence",
             legend=false,
             color=:steelblue,
             grid=true)

    # プロット2: 最終目的関数値
    p2 = plot(1:n_windows, J_final,
              xlabel="Window Index",
              ylabel="Objective Function J",
              title="Final Objective Value",
              marker=:circle,
              legend=false,
              yscale=:log10,
              color=:coral,
              linewidth=2,
              markersize=5,
              grid=true)

    # サブプロットを横並び配置
    p = plot(p1, p2, layout=(1,2), size=(1000, 400), dpi=150)

    # 保存
    if save_path !== nothing
        savefig(p, save_path)
        @info "✓ 収束履歴プロットを保存: $save_path"
    end

    return p
end

"""
    plot_flux_animation(q_global, config; kwargs...)

熱流束の時間発展アニメーション（GIF形式）

# 引数
- `q_global`: 熱流束配列 (nt-1, ni, nj)
- `config::Dict`: 設定情報

# キーワード引数
- `save_path::String`: 保存先パス（必須）
- `fps::Int`: フレームレート（デフォルト: 10）
- `skip::Int`: フレームスキップ数（デフォルト: 1、1フレームごと）

# 戻り値
- アニメーションオブジェクト
"""
function plot_flux_animation(q_global, config;
                             save_path="flux_animation.gif",
                             fps=10,
                             skip=1)
    nt_minus_1, ni, nj = size(q_global)
    dt = config["problem"]["dt"]

    # カラースケールの範囲を固定（全時刻の最小・最大）
    q_min = minimum(q_global)
    q_max = maximum(q_global)

    @info "熱流束アニメーション作成中..."
    @info @sprintf("  フレーム数: %d, FPS: %d, スキップ: %d", nt_minus_1, fps, skip)
    @info @sprintf("  熱流束範囲: %.2e - %.2e W/m²", q_min, q_max)

    # アニメーション作成
    anim = @animate for t_idx in 1:skip:nt_minus_1
        q_snapshot = q_global[t_idx, :, :]'
        t_ms = (t_idx - 1) * dt * 1000

        heatmap(q_snapshot,
                xlabel="x grid index",
                ylabel="y grid index",
                title=@sprintf("Heat Flux at t=%.2f ms", t_ms),
                color=:viridis,
                aspect_ratio=:equal,
                colorbar_title="W/m²",
                clims=(q_min, q_max),  # カラースケール固定
                size=(700, 600))
    end

    # GIF保存
    gif(anim, save_path, fps=fps)
    @info "✓ アニメーションを保存: $save_path"

    return anim
end

"""
    plot_all_diagnostics(q_global, windows_info, config; output_dir="plots/", format=:png)

全ての診断プロットを一括生成

# 引数
- `q_global`: 熱流束配列
- `windows_info`: ウィンドウ情報配列
- `config::Dict`: 設定情報

# キーワード引数
- `output_dir::String`: 出力ディレクトリ
- `format::Symbol`: 画像形式

# 戻り値
- 生成されたファイルパスのベクトル
"""
function plot_all_diagnostics(q_global, windows_info, config;
                               output_dir="plots/",
                               format=:png)
    @info "全診断プロットを生成中..."

    # 出力ディレクトリ作成
    mkpath(output_dir)

    generated_files = String[]

    # 1. 熱流束の時間変化（中心点）
    file_path = joinpath(output_dir, "heat_flux.$format")
    plot_heat_flux(q_global, config; save_path=file_path)
    push!(generated_files, file_path)

    # 2. ヒートマップ（初期、中間、最終）
    nt_minus_1 = size(q_global, 1)
    for (label, t_idx) in [("initial", 1),
                            ("mid", div(nt_minus_1, 2)),
                            ("final", nt_minus_1)]
        file_path = joinpath(output_dir, "heatmap_$label.$format")
        plot_flux_heatmap(q_global, config, t_idx; save_path=file_path)
        push!(generated_files, file_path)
    end

    # 3. CGM収束履歴
    file_path = joinpath(output_dir, "convergence.$format")
    plot_cgm_convergence(windows_info; save_path=file_path)
    push!(generated_files, file_path)

    @info "✓ 全プロット生成完了: $output_dir"
    @info "  生成ファイル数: $(length(generated_files))"

    return generated_files
end
