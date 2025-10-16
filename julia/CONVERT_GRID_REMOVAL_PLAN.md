# convert_to_guard_cell_grid削除作業計画

**作成日**: 2025年10月16日
**ブランチ**: tuning7
**関連コミット**: 4ffb82a

---

## 📋 作業概要

`convert_to_guard_cell_grid`関数の完全削除に伴い、この関数を使用している全てのテスト、ベンチマーク、スクリプトファイルを修正する。

---

## 🎯 修正対象ファイル詳細

### test/test_sliding_window.jl（3箇所）

#### 箇所1: 24行目（using文）
```julia
# 修正前
using IHCP_CGM: solve_sliding_window_cgm, WindowInfo, WorkBuffers, convert_to_guard_cell_grid

# 修正後
using IHCP_CGM: solve_sliding_window_cgm, WindowInfo, WorkBuffers
```

#### 箇所2: 261行目（1Dテスト）
```julia
# 前後のコンテキスト確認が必要
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

#### 箇所3: 364行目（2Dテスト）
```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### test/test_cgm_solver.jl（2箇所）

#### 箇所1: 18行目（using文）
```julia
# 修正前
using IHCP_CGM: solve_cgm!, WorkBuffers, convert_to_guard_cell_grid

# 修正後
using IHCP_CGM: solve_cgm!, WorkBuffers
```

#### 箇所2: 109行目
```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### benchmarks/benchmark_phase1b.jl（310行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# z_centersの計算が必要か確認
# dz_gridの代わりにdzを使用
```

**注意**: `dz_bottom`, `dz_top`という変数名を使用

---

### benchmarks/benchmark_combined_improvement.jl（195行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### benchmarks/benchmark_adjoint_scale_02.jl（195行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### benchmarks/benchmark_phase1e_adaptive.jl（131行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### benchmarks/benchmark_adjoint_scale_01.jl（195行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### benchmarks/benchmark_phase1e_tuned.jl（137行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### benchmarks/benchmark_phase1e_baseline.jl（131行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### scripts/run_10steps_fullsize_test.jl（133行目）

```julia
# 修正前
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# dz_gridの代わりにdzを使用
```

---

### scripts/verification/check_matrix_symmetry.jl（50行目）

```julia
# 修正前
z_centers_grid, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# z_centers_gridの代わりにz_centersを使用
# dz_gridの代わりにdzを使用
```

**注意**: このファイルは`z_centers_grid`という変数名を使用

---

### scripts/verification/test_cg_vs_bicgstab.jl（54行目）

```julia
# 修正前
z_centers_grid, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)

# 修正後
# convert_to_guard_cell_grid呼び出しを削除
# z_centers_gridの代わりにz_centersを使用
# dz_gridの代わりにdzを使用
```

**注意**: このファイルは`z_centers_grid`という変数名を使用

---

## 🔍 各ファイルで確認すべきポイント

### 1. z_centersの計算確認

各ファイルで`z_centers`が計算されているか確認。未定義の場合は追加：

```julia
# z_faces計算
z_faces = zeros(nk + 1)
z_faces[1] = 0.0
for k in 1:nk
    z_faces[k+1] = z_faces[k] + dz[k]
end

# z_centers計算
z_centers = zeros(nk)
z_centers[1] = z_faces[1]
z_centers[end] = z_faces[end]
for k in 2:(nk-1)
    z_centers[k] = (z_faces[k] + z_faces[k+1]) / 2.0
end
```

### 2. 変数名の置換

- `dz_grid` → `dz`
- `z_centers_grid` → `z_centers`（scriptsの一部）

### 3. dz_b, dz_tの扱い

テストデータから読み込んでいる場合は読み込みコードを維持：
```julia
# これは残す（テストデータとの整合性のため）
dz_b = convert(Vector{Float64}, prob["dz_b"])
dz_t = convert(Vector{Float64}, prob["dz_t"])

# これを削除
z_centers, dz_grid = convert_to_guard_cell_grid(nk, dz, dz_b, dz_t)
```

---

## 📝 作業チェックリスト

### Phase 1: テストファイル修正
- [ ] test/test_sliding_window.jl（using文 + 2箇所）
- [ ] test/test_cgm_solver.jl（using文 + 1箇所）
- [ ] テスト実行確認

### Phase 2: スクリプトファイル修正
- [ ] scripts/run_10steps_fullsize_test.jl
- [ ] scripts/verification/check_matrix_symmetry.jl
- [ ] scripts/verification/test_cg_vs_bicgstab.jl
- [ ] 動作確認（可能なら）

### Phase 3: ベンチマークファイル修正
- [ ] benchmarks/benchmark_phase1b.jl
- [ ] benchmarks/benchmark_combined_improvement.jl
- [ ] benchmarks/benchmark_adjoint_scale_02.jl
- [ ] benchmarks/benchmark_phase1e_adaptive.jl
- [ ] benchmarks/benchmark_adjoint_scale_01.jl
- [ ] benchmarks/benchmark_phase1e_tuned.jl
- [ ] benchmarks/benchmark_phase1e_baseline.jl

### Phase 4: 最終確認
- [ ] 全ファイルでconvert_to_guard_cell_grid検索（0件を確認）
- [ ] テスト全実行
- [ ] コミット

---

## 🚀 推奨作業順序

1. **test/test_sliding_window.jl** - 最も重要なテストファイル
2. **test/test_cgm_solver.jl** - CGMソルバーのテスト
3. **テスト実行** - ここで動作確認
4. **scripts/** - 3ファイルを修正
5. **benchmarks/** - 7ファイルを一括修正
6. **最終確認とコミット**

---

## 🎯 成功基準

1. `grep -rn "convert_to_guard_cell_grid" test/ benchmarks/ scripts/ | grep -v deprecated` → 0件
2. `julia --project=. test/runtests.jl` → 全テスト合格
3. GitTransform.jl関連コードが完全に削除されている

---

**作成日時**: 2025-10-16 14:25
**ブランチ**: tuning7
