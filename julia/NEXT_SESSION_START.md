# 次回セッション開始手順書

**最終更新**: 2025年10月16日 14:25
**ブランチ**: tuning7
**最新コミット**: 4ffb82a

---

## 📊 現在の状態

### ✅ 完了した作業

#### GridTransform.jl削除とdocstring整理（4ffb82a）

1. **GridTransform.jl完全削除**
   - `src/IHCP_CGM.jl`: include, using, exportを削除

2. **src/main.jl簡素化**
   - `convert_to_guard_cell_grid`呼び出しを削除
   - `dz_b`, `dz_t`計算ロジックを削除（22行削除）
   - `z_centers`を直接計算し、`dz`をそのまま使用

3. **docstring整理**
   - `AdjointSolver.jl`: `dz_b`, `dz_t`パラメータ記載を削除
   - `CGMSolver.jl`: 関数シグネチャから`dz_b`, `dz_t`を削除（2箇所）
   - `SlidingWindowSolver.jl`: 関数シグネチャから`dz_b`, `dz_t`を削除

**変更統計**: 5ファイル修正（+7, -28行）

---

## 🎯 次回の作業: convert_to_guard_cell_grid削除に伴う修正

### 対象ファイル一覧（15箇所）

#### test/（5箇所）
1. **test/test_sliding_window.jl** - 3箇所
   - 24行: `using IHCP_CGM: ..., convert_to_guard_cell_grid` → 削除
   - 261行: `z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)`
   - 364行: `z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)`

2. **test/test_cgm_solver.jl** - 2箇所
   - 18行: `using IHCP_CGM: ..., convert_to_guard_cell_grid` → 削除
   - 109行: `z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)`

#### benchmarks/（7箇所）
3. **benchmarks/benchmark_phase1b.jl** - 310行
4. **benchmarks/benchmark_combined_improvement.jl** - 195行
5. **benchmarks/benchmark_adjoint_scale_02.jl** - 195行
6. **benchmarks/benchmark_phase1e_adaptive.jl** - 131行
7. **benchmarks/benchmark_adjoint_scale_01.jl** - 195行
8. **benchmarks/benchmark_phase1e_tuned.jl** - 137行
9. **benchmarks/benchmark_phase1e_baseline.jl** - 131行

#### scripts/（3箇所）
10. **scripts/run_10steps_fullsize_test.jl** - 133行
11. **scripts/verification/check_matrix_symmetry.jl** - 50行
12. **scripts/verification/test_cg_vs_bicgstab.jl** - 54行

---

## 🔧 修正パターン

### パターン1: using文から削除

```julia
# 修正前
using IHCP_CGM: solve_cgm!, WorkBuffers, convert_to_guard_cell_grid

# 修正後
using IHCP_CGM: solve_cgm!, WorkBuffers
```

### パターン2: 関数呼び出しの削除

```julia
# 修正前
dz_b = convert(Vector{Float64}, prob["dz_b"])
dz_t = convert(Vector{Float64}, prob["dz_t"])
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)
# ... 以降、z_centersとdz_gridを使用

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# z_centersは計算済み、dzをそのまま使用
# dz_gridの代わりにdzを使用
```

### パターン3: 変数名置換

```julia
# 呼び出し部分
# 修正前: dz_grid
# 修正後: dz
```

---

## 📝 作業手順

### ステップ1: 状態確認
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia
git status
git log --oneline -5
cat NEXT_SESSION_START.md
```

### ステップ2: 修正対象の再確認
```bash
# convert_to_guard_cell_grid使用箇所の確認
grep -rn "convert_to_guard_cell_grid" test/ benchmarks/ scripts/ | grep -v deprecated

# 期待値: 15箇所（test:5, benchmarks:7, scripts:3）
```

### ステップ3: 修正作業

#### 優先順位1: テストファイル（重要度高）
1. `test/test_sliding_window.jl`
2. `test/test_cgm_solver.jl`

#### 優先順位2: スクリプトファイル
3. `scripts/run_10steps_fullsize_test.jl`
4. `scripts/verification/check_matrix_symmetry.jl`
5. `scripts/verification/test_cg_vs_bicgstab.jl`

#### 優先順位3: ベンチマークファイル（7ファイル）
6-12. benchmarks/*.jl

### ステップ4: 修正後のテスト実行
```bash
# 全テスト実行
julia --project=. test/runtests.jl

# 個別テスト
julia --project=. test/test_cgm_solver.jl
julia --project=. test/test_sliding_window.jl
```

---

## 🔍 注意点

### 1. dz_b, dz_tの扱い

これらの変数は**テストデータから読み込まれている**場合がある：
```julia
dz_b = convert(Vector{Float64}, prob["dz_b"])  # JSONから読込
dz_t = convert(Vector{Float64}, prob["dz_t"])  # JSONから読込
```

修正方針：
- これらの読み込みコードは**そのまま残す**（テストデータの整合性のため）
- `convert_to_guard_cell_grid`呼び出しのみ削除
- 使用箇所で`dz_grid`を`dz`に置換

### 2. z_centersの計算

main.jlと同様に、各ファイルで`z_centers`が計算されているか確認：
```julia
# z_centersが未定義の場合は追加する必要あり
z_centers = zeros(nk)
z_centers[1] = z_faces[1]
z_centers[end] = z_faces[end]
for k in 2:(nk-1)
    z_centers[k] = (z_faces[k] + z_faces[k+1]) / 2.0
end
```

### 3. deprecatedディレクトリ

`test/deprecated/`内のファイルは**修正不要**（過去のコードとして保存）

---

## 📂 重要なファイル

### ドキュメント
- `NEXT_SESSION_START.md` - このファイル（次回セッション手順書）
- `VARIABLE_RENAME_PROGRESS.md` - 変数名変更作業の進捗記録

### 修正済みファイル
- `src/IHCP_CGM.jl`
- `src/main.jl`
- `src/solvers/AdjointSolver.jl`
- `src/solvers/CGMSolver.jl`
- `src/solvers/SlidingWindowSolver.jl`

---

## 🚀 次回セッション開始コマンド

```bash
# 1. ディレクトリ移動
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex/julia

# 2. 状態確認
git status
git log --oneline -5

# 3. この手順書を表示
cat NEXT_SESSION_START.md

# 4. 修正対象確認
grep -rn "convert_to_guard_cell_grid" test/ benchmarks/ scripts/ | grep -v deprecated | wc -l
# 期待値: 15

# 5. 作業開始
# → test/test_sliding_window.jl から修正開始
```

---

## 📊 進捗管理

- **完了**: src/ディレクトリ（5ファイル）✅
- **残り**: test/, benchmarks/, scripts/（15箇所）
- **推定作業時間**: 30-45分

---

**作成日時**: 2025-10-16 14:25
**ブランチ**: tuning7
**最新コミット**: 4ffb82a
