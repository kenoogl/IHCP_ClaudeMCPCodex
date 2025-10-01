# Phase 6 C-1: IRカメラMATLABファイル読込機能 詳細調査レポート

**作成日**: 2025-10-01
**作成者**: Codex MCP Agent
**目的**: Julia実装のための実データ読込機能の仕様調査

---

## 1. MATLABファイル読込関数の特定

### 1.1 関数一覧

#### 関数1: `extract_sorted_mat_files(folder_path)`
- **行番号**: 1000-1009
- **引数**:
  - `folder_path` (str): MATLABファイルが格納されているディレクトリパス
- **戻り値**:
  - `list[str]`: ソートされたMATファイル名のリスト
- **処理概要**:
  1. 指定フォルダ内のファイルをスキャン
  2. `SUS`で始まり`.MAT`で終わるファイルを抽出
  3. ファイル名から数値を抽出（正規表現: `r'SUS(\d+)\.MAT'`）
  4. 数値順にソート（自然順ソート）

**重要ポイント**:
- **自然順ソート**: 数値部分を整数として比較（例: SUS1.MAT < SUS2.MAT < SUS10.MAT）
- **大文字拡張子**: `.MAT`（大文字）のみをマッチング（`.mat`は非対応）
- **エラー処理なし**: マッチしないファイルは-1を返すが、警告なし

---

#### 関数2: `load_region_temperature(folder_path, i_range, j_range)`
- **行番号**: 1015-1038
- **引数**:
  - `folder_path` (str): MATLABファイルのディレクトリパス
  - `i_range` (tuple): y方向の範囲 `(i_start, i_end)` （両端を含む）
  - `j_range` (tuple): x方向の範囲 `(j_start, j_end)` （両端を含む）
- **戻り値**:
  - `region_data_K` (np.ndarray): 温度データ配列 [K]、shape: `(height, width, num_frames)`
  - `num_frames` (int): 読み込んだフレーム数
- **処理概要**:
  1. `extract_sorted_mat_files()`でMATファイルリストを取得
  2. 各ファイルを`scipy.io.loadmat()`で読み込み
  3. 変数名は**ファイル名から拡張子を除いたもの**（例: `SUS123.MAT` → `SUS123`）
  4. 指定領域を切り出して3次元配列に格納

**重要ポイント**:
- **変数名規則**: MATLABファイル内の変数名 = ファイル名（拡張子なし）
- **座標系**: `frame_data[i, j]` → i:行方向（y物理方向）、j:列方向（x物理方向）
- **インデックス**: 両端を含むスライス（end+1が必要）
- **データ型**: ケルビン単位で保存（コメントに記載）

---

## 2. 読込対象ファイル形式

### 2.1 ファイル名パターン

**正規表現**: `r'SUS(\d+)\.MAT'`

**具体例**:
- `SUS1.MAT`, `SUS2.MAT`, ..., `SUS500.MAT`
- `SUS0001.MAT`, `SUS0002.MAT` (ゼロパディング版も対応)

**非対応パターン**:
- 小文字拡張子: `SUS1.mat` ❌
- 異なる接頭辞: `TEMP1.MAT`, `IR1.MAT` ❌
- 数値なし: `SUS.MAT` ❌

### 2.2 ディレクトリ構造

現在の実装では**ディレクトリパスを直接指定**:
```python
mat_files = extract_sorted_mat_files(folder_path)
```

**想定構造**:
```
folder_path/
├── SUS1.MAT
├── SUS2.MAT
├── SUS3.MAT
...
└── SUS500.MAT
```

**注意**:
- サブディレクトリの再帰探索は**非対応**
- 単一フォルダに全MATファイルが配置される前提

---

## 3. データ構造

### 3.1 MATLABファイル内部構造

**scipy.io.loadmat()の戻り値**:
```python
mat_data = {
    '__header__': b'MATLAB 5.0 MAT-file...',
    '__version__': '1.0',
    '__globals__': [],
    'SUS123': np.ndarray(shape=(ny, nx), dtype=float64)  # 温度データ本体
}
```

**キーとなる変数名**:
- **ファイル名から自動生成**: `file[:-4]`
- 例: `SUS123.MAT` → `mat_data['SUS123']`

### 3.2 配列形状

**元データ（1フレーム）**:
- `frame_data`: shape `(ny, nx)` （2次元配列）
- データ型: `float64`（推定）
- 単位: ケルビン [K]

**切り出し後の領域データ**:
- `sub_region`: shape `(height, width)`
  - `height = i_end - i_start + 1`（y方向ピクセル数）
  - `width = j_end - j_start + 1`（x方向ピクセル数）

**時系列統合後**:
- `region_data_K`: shape `(height, width, num_frames)`
  - 軸0: y方向（i方向）
  - 軸1: x方向（j方向）
  - 軸2: 時間方向（フレーム番号）

### 3.3 座標系の対応

**物理座標 vs 配列インデックス**:
```
frame_data[i, j] ↔ 物理座標(x, y)
  i: 行方向 = y物理方向（鉛直方向）
  j: 列方向 = x物理方向（水平方向）
```

**注意**: Pythonの配列は(row, column) = (y, x)の順序

---

## 4. ファイル選択・ソート機能

### 4.1 ソート順序

**自然順ソート（Natural Sorting）**:
```python
def extract_number(file_name):
    match = re.match(r'SUS(\d+)\.MAT', file_name)
    return int(match.group(1)) if match else -1

sorted(files, key = extract_number)
```

**動作例**:
- 辞書順: `SUS1.MAT, SUS10.MAT, SUS2.MAT` ❌
- 自然順: `SUS1.MAT, SUS2.MAT, SUS10.MAT` ✅

---

## 5. Julia実装設計

### 5.1 必要なJuliaパッケージ

| 機能 | Pythonパッケージ | Juliaパッケージ |
|------|------------------|-----------------|
| MATファイル読込 | `scipy.io` | `MAT.jl` |
| 正規表現 | `re` | `Base.PCRE` (標準) |
| 配列操作 | `numpy` | `Base.Array` (標準) |

### 5.2 実装関数設計

```julia
# julia/src/utils/data_loaders.jl

"""
    extract_sorted_mat_files(folder_path; prefix="SUS", extension=".MAT")

MATファイルリストを自然順ソートで取得

# 引数
- `folder_path::String`: MATファイルのディレクトリパス
- `prefix::String`: ファイル名接頭辞（デフォルト: "SUS"）
- `extension::String`: ファイル拡張子（デフォルト: ".MAT"）

# 戻り値
- `Vector{String}`: ソートされたファイル名リスト

# 例外
- ディレクトリが存在しない場合エラー
- MATファイルが1つも見つからない場合エラー
"""
function extract_sorted_mat_files(
  folder_path::String;
  prefix::String="SUS",
  extension::String=".MAT"
)::Vector{String}
  # 実装...
end

"""
    load_region_temperature(folder_path, i_range, j_range)

IRカメラMATLABファイルから指定領域の温度データを読み込む

# 引数
- `folder_path::String`: MATファイルのディレクトリパス
- `i_range::Tuple{Int,Int}`: y方向の範囲（両端を含む）
- `j_range::Tuple{Int,Int}`: x方向の範囲（両端を含む）

# 戻り値
- `region_data::Array{Float64,3}`: 温度データ [K]、shape: (height, width, num_frames)
- `num_frames::Int`: フレーム数

# 例外
- ファイルが存在しない場合エラー
- 変数名が不一致の場合エラー
- データにNaN/Infが含まれる場合警告
"""
function load_region_temperature(
  folder_path::String,
  i_range::Tuple{Int,Int},
  j_range::Tuple{Int,Int}
)::Tuple{Array{Float64,3}, Int}
  # 実装...
end
```

### 5.3 実装優先順位

**Phase 1（最優先）**:
1. `extract_sorted_mat_files()` の基本実装
2. `load_region_temperature()` の基本実装
3. 小規模テストケース（2-3ファイル）での動作確認

**Phase 2（高優先）**:
4. エラーハンドリング（ファイル存在チェック、データ検証）
5. データ検証（NaN/Inf、温度範囲）
6. ユニットテスト作成

**Phase 3（中優先）**:
7. メモリ効率化（イテレータパターン）
8. 並列読込の実装
9. ベンチマーク（Pythonとの性能比較）

---

## 6. 実装時の注意点

### 6.1 座標系の違い

**Python (row-major)**:
```python
frame_data[i, j]  # i:行, j:列
```

**Julia (column-major)**:
```julia
frame_data[j, i]  # 逆順！または転置が必要
```

### 6.2 配列インデックスの違い

**Python (0-based)**:
```python
sub_region = frame_data[i_start:i_end+1, j_start:j_end+1]  # end+1が必要
```

**Julia (1-based)**:
```julia
sub_region = frame_data[i_start:i_end, j_start:j_end]  # endを含む
```

### 6.3 正規表現の違い

**Python**:
```python
match = re.match(r'SUS(\d+)\.MAT', file_name)
number = int(match.group(1))
```

**Julia**:
```julia
m = match(r"SUS(\d+)\.MAT", file_name)
number = parse(Int, m.captures[1])
```

---

## 7. テスト戦略

### 7.1 ユニットテスト設計

```julia
@testset "MATファイル読込テスト" begin
  @testset "extract_sorted_mat_files" begin
    # Test 1: 自然順ソート
    # Test 2: ファイルが存在しない場合
    # Test 3: MATファイルが1つもない場合
  end

  @testset "load_region_temperature" begin
    # Test 1: 正常読込
    # Test 2: 領域切り出し
    # Test 3: 複数フレームの統合
    # Test 4: データ検証（NaN/Inf）
  end
end
```

### 7.2 統合テスト

- 実際のIRカメラデータ（SUS*.MAT）での読込テスト
- Python版との結果比較（数値一致確認）
- 大容量データでのメモリ使用量測定

---

## 8. まとめ

### 8.1 重要な発見事項

1. **変数名規則の特殊性**: MATLABファイル内の変数名がファイル名（拡張子なし）と一致
2. **自然順ソート**: 数値部分を整数として比較
3. **座標系の複雑さ**: `i`がy方向、`j`がx方向
4. **エラー処理の欠如**: Python実装はシンプルだが堅牢性に欠ける

### 8.2 Julia実装の方針

- TDD方針: テスト先行→実装→検証
- 段階的実装: 基本機能→エラー処理→性能最適化
- Python互換性: 同じインターフェース、同じ結果
