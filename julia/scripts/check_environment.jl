#!/usr/bin/env julia
"""
環境確認スクリプト

test_dhcp_solver.jlを実行する前に、環境が正しく設定されているか確認します。

Usage:
  julia --project=julia julia/scripts/check_environment.jl
"""

using Printf

const BASE_DIR = @__DIR__
const PROJECT_ROOT = normpath(joinpath(BASE_DIR, "..", ".."))

println("="^80)
println("Julia版DHCP環境確認スクリプト")
println("="^80)
println()

# チェック結果を追跡（グローバル変数として明示的に宣言）
const all_checks_passed = Ref(true)

# ============================================================================
# 1. Juliaバージョン確認
# ============================================================================
println("✓ チェック 1/6: Juliaバージョン")
println("  現在のバージョン: $(VERSION)")
if VERSION >= v"1.10.0"
    println("  ✅ OK (1.10以上)")
else
    println("  ❌ NG: Julia 1.10以上が必要です")
    all_checks_passed[] = false
end
println()

# ============================================================================
# 2. プロジェクトファイル確認
# ============================================================================
println("✓ チェック 2/6: プロジェクトファイル")
project_file = joinpath(PROJECT_ROOT, "julia", "Project.toml")
if isfile(project_file)
    println("  ✅ Project.toml が見つかりました")
    println("     $(project_file)")
else
    println("  ❌ Project.toml が見つかりません")
    println("     期待パス: $(project_file)")
    all_checks_passed[] = false
end
println()

# ============================================================================
# 3. 依存パッケージ確認
# ============================================================================
println("✓ チェック 3/6: 依存パッケージ")
required_packages = [
    "NPZ",
    "CSV",
    "FLoops",
    "IterativeSolvers",
    "LinearAlgebra",
    "Statistics"
]

pkg_errors = []
for pkg_name in required_packages
    try
        @eval using $(Symbol(pkg_name))
        print("  ✅ $(pkg_name)")

        # パッケージバージョン表示
        try
            pkg = Base.PkgId(Base.UUID(Base.loaded_modules_array()[findfirst(m -> m[1].name == pkg_name, Base.loaded_modules_array())][1].uuid), pkg_name)
            ver = string(Pkg.dependencies()[pkg.uuid].version)
            print(" (v$(ver))")
        catch
        end
        println()
    catch e
        println("  ❌ $(pkg_name) が読み込めません")
        push!(pkg_errors, pkg_name)
        all_checks_passed[] = false
    end
end

if !isempty(pkg_errors)
    println()
    println("  以下のコマンドで依存パッケージをインストールしてください:")
    println("  julia --project=julia -e 'using Pkg; Pkg.instantiate()'")
end
println()

# ============================================================================
# 4. ソースコード確認
# ============================================================================
println("✓ チェック 4/6: ソースコード")
required_files = [
    "julia/src/IHCP_CGM.jl",
    "julia/src/ThermalProperties.jl",
    "julia/src/DataLoaders.jl",
    "julia/src/solvers/DHCPSolver.jl",
    "julia/src/solvers/CommonSolver.jl",
    "julia/scripts/test_dhcp_solver.jl"
]

src_errors = Ref(0)
for rel_path in required_files
    full_path = joinpath(PROJECT_ROOT, rel_path)
    if isfile(full_path)
        println("  ✅ $(rel_path)")
    else
        println("  ❌ $(rel_path) が見つかりません")
        src_errors[] += 1
        all_checks_passed[] = false
    end
end

if src_errors[] == 0
    println("  すべてのソースファイルが揃っています")
end
println()

# ============================================================================
# 5. データファイル確認
# ============================================================================
println("✓ チェック 5/6: データファイル")

# 熱物性値ファイル（必須）
thermal_props = joinpath(PROJECT_ROOT, "shared", "data", "metal_thermal_properties.csv")
if isfile(thermal_props)
    file_size = filesize(thermal_props)
    println("  ✅ 熱物性値ファイル: $(thermal_props)")
    println("     サイズ: $(file_size) bytes")
else
    println("  ❌ 熱物性値ファイルが見つかりません（必須）")
    println("     期待パス: $(thermal_props)")
    all_checks_passed[] = false
end

# 測定データファイル（オプション）
measurement_data = joinpath(PROJECT_ROOT, "shared", "data", "T_measure_700um_1ms.npy")
if isfile(measurement_data)
    file_size = filesize(measurement_data)
    file_size_mb = round(file_size / 1024^2, digits=1)
    println("  ✅ 測定データファイル: $(measurement_data)")
    println("     サイズ: $(file_size_mb) MB")
    println("     → 測定データモードで実行可能")
else
    println("  ⚠️  測定データファイルが見つかりません（オプション）")
    println("     期待パス: $(measurement_data)")
    println("     → --synthetic フラグで実行してください")
end
println()

# ============================================================================
# 6. スレッド設定確認
# ============================================================================
println("✓ チェック 6/6: スレッド設定")
num_threads = Threads.nthreads()
println("  利用可能スレッド数: $(num_threads)")
if num_threads == 1
    println("  ⚠️  1スレッドで実行されます")
    println("     並列実行する場合は以下のように実行してください:")
    println("     JULIA_NUM_THREADS=4 julia --project=julia ...")
else
    println("  ✅ マルチスレッド実行が有効です")
end
println()

# ============================================================================
# 総合結果
# ============================================================================
println("="^80)
if all_checks_passed[]
    println("🎉 すべてのチェックが完了しました！")
    println()
    println("実行例:")
    println()
    println("【測定データモード】")
    if isfile(measurement_data)
        println("  # 推奨設定（GS前処理、PCGソルバー）")
        if num_threads > 1
            println("  julia --project=julia julia/scripts/test_dhcp_solver.jl \\")
        else
            println("  JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \\")
        end
        println("    --solver pcg --precond gs --nt 10")
    else
        println("  測定データファイルがないため、syntheticモードで実行してください")
    end
    println()
    println("【Syntheticモード（測定データ不要）】")
    println("  # 小規模テスト（数秒）")
    println("  julia --project=julia julia/scripts/test_dhcp_solver.jl \\")
    println("    --synthetic --ni 40 --nj 50 --nk 10 --nt 5")
    println()
    println("  # フルサイズテスト")
    if num_threads > 1
        println("  julia --project=julia julia/scripts/test_dhcp_solver.jl \\")
    else
        println("  JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \\")
    end
    println("    --synthetic --solver pcg --precond gs --nt 10")
else
    println("❌ いくつかのチェックが失敗しました")
    println()
    println("上記のエラーを修正してから、再度このスクリプトを実行してください。")
    println()
    println("詳細な設定手順は以下を参照してください:")
    println("  julia/SETUP.md")
    println("  julia/QUICKSTART.md")
end
println("="^80)
