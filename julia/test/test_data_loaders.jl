"""
Phase 6 C-1: データ読込機能のテスト

実データ読込機能（extract_sorted_mat_files, load_region_temperature）のテスト
"""

using Test
using MAT
using Logging

# テスト用の一時ディレクトリとMATファイルを作成する関数
function setup_test_mat_files()
  # 一時ディレクトリを作成
  test_dir = mktempdir()

  # テストデータ: 10x8の温度配列（ケルビン単位）
  # 各ファイルで異なる温度パターンを設定

  # SUS1.MAT: 基準温度300K
  mat1_data = fill(300.0, 10, 8)
  matwrite(joinpath(test_dir, "SUS1.MAT"), Dict("SUS1" => mat1_data))

  # SUS2.MAT: 温度勾配あり
  mat2_data = [300.0 + i + j for i in 1:10, j in 1:8]
  matwrite(joinpath(test_dir, "SUS2.MAT"), Dict("SUS2" => mat2_data))

  # SUS10.MAT: 自然順ソートテスト用（10 > 2）
  mat10_data = fill(350.0, 10, 8)
  matwrite(joinpath(test_dir, "SUS10.MAT"), Dict("SUS10" => mat10_data))

  # SUS20.MAT: さらに大きい番号
  mat20_data = fill(400.0, 10, 8)
  matwrite(joinpath(test_dir, "SUS20.MAT"), Dict("SUS20" => mat20_data))

  # ノイズファイル（マッチしないファイル）
  open(joinpath(test_dir, "README.txt"), "w") do f
    write(f, "Test files")
  end

  # 小文字拡張子（マッチしない）
  mat_lower_data = fill(999.0, 10, 8)
  matwrite(joinpath(test_dir, "SUS99.mat"), Dict("SUS99" => mat_lower_data))

  return test_dir
end

# 実装済みデータローダーモジュールを読み込む
# src/utils/data_loaders.jl から MATDataLoaders をインポート
push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
include("../src/utils/data_loaders.jl")
using .MATDataLoaders

@testset "Phase 6 C-1: データ読込機能テスト" begin

  # テスト用MATファイルのセットアップ
  test_dir = setup_test_mat_files()

  @testset "extract_sorted_mat_files" begin

    @testset "正常系: 自然順ソート" begin
      files = MATDataLoaders.extract_sorted_mat_files(test_dir)

      # 期待されるファイル順: SUS1.MAT, SUS2.MAT, SUS10.MAT, SUS20.MAT
      @test length(files) == 4
      @test files[1] == "SUS1.MAT"
      @test files[2] == "SUS2.MAT"
      @test files[3] == "SUS10.MAT"  # 辞書順だとSUS2より前になるが、自然順では後
      @test files[4] == "SUS20.MAT"
    end

    @testset "正常系: カスタムprefix/extension" begin
      # TEMP接頭辞のファイルを作成
      temp_dir = mktempdir()
      matwrite(joinpath(temp_dir, "TEMP1.MAT"), Dict("TEMP1" => fill(300.0, 5, 5)))
      matwrite(joinpath(temp_dir, "TEMP5.MAT"), Dict("TEMP5" => fill(350.0, 5, 5)))

      files = MATDataLoaders.extract_sorted_mat_files(temp_dir; prefix="TEMP")
      @test length(files) == 2
      @test files[1] == "TEMP1.MAT"
      @test files[2] == "TEMP5.MAT"

      rm(temp_dir, recursive=true)
    end

    @testset "異常系: ディレクトリが存在しない" begin
      @test_throws Exception MATDataLoaders.extract_sorted_mat_files("/nonexistent/path")
    end

    @testset "異常系: MATファイルが1つもない" begin
      empty_dir = mktempdir()
      @test_throws Exception MATDataLoaders.extract_sorted_mat_files(empty_dir)
      rm(empty_dir, recursive=true)
    end

    @testset "境界値: 数値部分のゼロパディング" begin
      # ゼロパディング版のテスト
      pad_dir = mktempdir()
      matwrite(joinpath(pad_dir, "SUS001.MAT"), Dict("SUS001" => fill(300.0, 5, 5)))
      matwrite(joinpath(pad_dir, "SUS010.MAT"), Dict("SUS010" => fill(350.0, 5, 5)))
      matwrite(joinpath(pad_dir, "SUS100.MAT"), Dict("SUS100" => fill(400.0, 5, 5)))

      files = MATDataLoaders.extract_sorted_mat_files(pad_dir)
      @test length(files) == 3
      @test files[1] == "SUS001.MAT"
      @test files[2] == "SUS010.MAT"
      @test files[3] == "SUS100.MAT"

      rm(pad_dir, recursive=true)
    end
  end

  @testset "load_region_temperature" begin

    @testset "正常系: 単一フレーム読込" begin
      # SUS1.MATのみのディレクトリ
      single_dir = mktempdir()
      mat_data = fill(300.0, 10, 8)
      matwrite(joinpath(single_dir, "SUS1.MAT"), Dict("SUS1" => mat_data))

      # 全領域読込
      region_data, num_frames = MATDataLoaders.load_region_temperature(
        single_dir,
        (1, 10),
        (1, 8)
      )

      @test num_frames == 1
      @test size(region_data) == (10, 8, 1)
      @test all(region_data[:, :, 1] .≈ 300.0)

      rm(single_dir, recursive=true)
    end

    @testset "正常系: 複数フレーム読込と時系列統合" begin
      # test_dirには4つのMATファイル（SUS1, SUS2, SUS10, SUS20）がある
      region_data, num_frames = MATDataLoaders.load_region_temperature(
        test_dir,
        (1, 10),
        (1, 8)
      )

      @test num_frames == 4
      @test size(region_data) == (10, 8, 4)

      # 各フレームの温度パターンを検証
      @test all(region_data[:, :, 1] .≈ 300.0)  # SUS1: 全て300K
      @test region_data[1, 1, 2] ≈ 302.0       # SUS2: 勾配パターン (i=1, j=1 → 300+1+1)
      @test all(region_data[:, :, 3] .≈ 350.0)  # SUS10: 全て350K
      @test all(region_data[:, :, 4] .≈ 400.0)  # SUS20: 全て400K
    end

    @testset "正常系: 部分領域の切り出し" begin
      # 中央部分を切り出し: i=[3,7], j=[2,6]
      region_data, num_frames = MATDataLoaders.load_region_temperature(
        test_dir,
        (3, 7),
        (2, 6)
      )

      @test num_frames == 4
      @test size(region_data) == (5, 5, 4)  # height=7-3+1=5, width=6-2+1=5

      # SUS2の勾配パターンを部分領域で検証
      # 元データ: mat2_data[i,j] = 300 + i + j
      # 切り出し領域: i=[3,7], j=[2,6]
      # region_data[1,1,2] = mat2_data[3,2] = 300 + 3 + 2 = 305
      @test region_data[1, 1, 2] ≈ 305.0
      @test region_data[5, 5, 2] ≈ 313.0  # mat2_data[7,6] = 300 + 7 + 6
    end

    @testset "異常系: ファイルが存在しない" begin
      @test_throws Exception MATDataLoaders.load_region_temperature(
        "/nonexistent/path",
        (1, 10),
        (1, 8)
      )
    end

    @testset "異常系: 範囲外のインデックス" begin
      # 元データは10x8なので、範囲外を指定
      @test_throws Exception MATDataLoaders.load_region_temperature(
        test_dir,
        (1, 20),  # y方向が範囲外
        (1, 8)
      )
    end

    @testset "境界値: 1x1の最小領域" begin
      region_data, num_frames = MATDataLoaders.load_region_temperature(
        test_dir,
        (5, 5),  # 単一ピクセル
        (4, 4)
      )

      @test size(region_data) == (1, 1, 4)
      @test region_data[1, 1, 1] ≈ 300.0  # SUS1
      @test region_data[1, 1, 2] ≈ 309.0  # SUS2: 300 + 5 + 4
    end

    @testset "データ検証: NaN/Inf検出" begin
      # NaNを含むテストデータ
      nan_dir = mktempdir()
      nan_data = fill(300.0, 10, 8)
      nan_data[5, 5] = NaN
      matwrite(joinpath(nan_dir, "SUS1.MAT"), Dict("SUS1" => nan_data))

      # NaNを含む場合は警告を出すが、エラーは投げない（仕様確認が必要）
      @test_logs (:warn, r"NaN") match_mode=:any begin
        region_data, _ = MATDataLoaders.load_region_temperature(
          nan_dir,
          (1, 10),
          (1, 8)
        )
      end

      rm(nan_dir, recursive=true)
    end
  end

  # テスト後のクリーンアップ
  rm(test_dir, recursive=true)
end
