# スライディングウィンドウ調整システム開発計画

**作成日**: 2025-10-27
**ステータス**: 承認待ち
**担当**: Claude Code
**参照**: `docs/reports/julia_sliding_window_tuning_plan.md`, `docs/reports/basesize_guidelines.md`

---

## 1. 目的・背景

### 1.1 目的
Julia版スライディングウィンドウCGM計算における最適なウィンドウサイズとオーバーラップ値を体系的に探索・評価するシステムを構築する。

### 1.2 背景
- codexにより`basesize_guidelines.md`と`julia_sliding_window_tuning_plan.md`が作成済み
- basesize最適化（600推奨）は既に実施済み
- ウィンドウサイズとオーバーラップの最適値は未確定
- 複数の組み合わせを手動実行するのは非効率

### 1.3 期待される効果
- 計算時間の短縮
- CGM収束性の向上
- ウィンドウ境界での熱流束の滑らかさ向上

---

## 2. 現状分析

### 2.1 既存ファイルの状態

| ファイル | 状態 | 備考 |
|---------|------|------|
| `basesize_guidelines.md` | ✅ 完成 | basesize=600推奨 |
| `julia_sliding_window_tuning_plan.md` | ✅ 完成 | 実装仕様書 |
| `run_10steps_fullsize_test.jl` | ✅ 修正済み | basesize: 1→600 |
| `run_sliding_window.jl` | ✅ 修正済み | basesize: 1→600, --dry-run実装済み |
| `test_dhcp_solver.jl` | ⚠️ 修正必要 | basesize: 1→600に変更必要 |
| `tune_sliding_window.jl` | ❌ 未作成 | 計画書のみ存在 |

### 2.2 参照する計画書の要点

**`julia_sliding_window_tuning_plan.md`の主要機能**:
1. 複数のwindow/overlap組み合わせを自動実行
2. 各組み合わせのログから統計値を抽出（CGM反復数、計算時間、熱流束範囲）
3. 結果をCSV/Markdownで出力
4. `--dry-run`で実験規模を事前確認
5. `--reuse`で既存結果を再利用

---

## 3. 実施する作業

### 3.1 Phase 1: 既存ファイル修正（優先度：高）

**作業内容**: `test_dhcp_solver.jl`のbasesize最適化

**修正箇所**:
```julia
# julia/scripts/test_dhcp_solver.jl:49
basesize = 1              # デフォルト
↓
basesize = 600            # 最適値（スレッド数×basesize最適化実験の結果）
```

**修正理由**:
- 他のスクリプトと一貫性を保つ
- DHCP予備確認時に最適性能を得る

**所要時間**: 5分

---

### 3.2 Phase 2: tune_sliding_window.jl新規作成（優先度：高）

**作業内容**: 計画書仕様に基づくチューニングスクリプト作成

#### 3.2.1 基本機能

**コマンドライン引数**:
```bash
julia tune_sliding_window.jl \
  --nt 120 \
  --windows 50,71,100 \
  --overlap-ratios 0.2,0.3 \
  --cgm-iter 3 \
  --solver pbicgstab \
  --precond gs \
  --basesize 600 \
  --dry-run  # 実験規模プレビュー
```

**または固定値オーバーラップ**:
```bash
julia tune_sliding_window.jl \
  --nt 120 \
  --windows 50,71,100 \
  --overlaps 10,20,30 \
  --cgm-iter 3
```

#### 3.2.2 実装仕様

**1. 組み合わせ生成**
```julia
function generate_combinations(windows, overlaps_or_ratios, is_ratio)
  combinations = []
  for window in windows
    for value in overlaps_or_ratios
      overlap = is_ratio ? max(1, round(Int, window * value)) : value
      if overlap < window
        push!(combinations, (window=window, overlap=overlap))
      end
    end
  end
  return combinations
end
```

**2. ジョブ実行**
```julia
function run_single_job(nt, window, overlap, cgm_iter, solver, precond, basesize, output_dir)
  tag = "julia_sw_w$(window)_o$(overlap)"
  log_path = joinpath(output_dir, "$(tag).log")

  cmd = `julia --project=julia julia/scripts/run_sliding_window.jl
         --nt $nt --window $window --overlap $overlap
         --cgm-iter $cgm_iter --solver $solver --precond $precond
         --basesize $basesize`

  run(pipeline(cmd, stdout=log_path, stderr=log_path))
  return log_path
end
```

**3. ログ解析**
```julia
function extract_metrics(log_path)
  total_time = NaN
  cgm_iters = Int[]
  q_min = NaN
  q_max = NaN

  for line in eachline(log_path)
    # "Total runtime: 123.45 s"
    m = match(r"Total runtime:\s+([\d.]+)\s+s", line)
    if m !== nothing
      total_time = parse(Float64, m.captures[1])
    end

    # "  CGM iterations: 10"
    m = match(r"CGM iterations:\s+(\d+)", line)
    if m !== nothing
      push!(cgm_iters, parse(Int, m.captures[1]))
    end

    # "  q range: [1.23e+05, 4.56e+05] W/m^2"
    m = match(r"q range:\s+\[([\d.e+-]+),\s*([\d.e+-]+)\]", line)
    if m !== nothing
      q_min = parse(Float64, m.captures[1])
      q_max = parse(Float64, m.captures[2])
    end
  end

  return (
    total_time = total_time,
    avg_cgm_iters = isempty(cgm_iters) ? NaN : mean(cgm_iters),
    max_cgm_iters = isempty(cgm_iters) ? 0 : maximum(cgm_iters),
    q_min = q_min,
    q_max = q_max
  )
end
```

**4. 結果出力**
```julia
# CSV形式
using CSV, DataFrames
CSV.write("shared/results/tuning/julia_sliding_window_scan.csv", results_df)

# Markdown形式
function write_markdown_report(results_df, output_path)
  open(output_path, "w") do io
    println(io, "# Sliding Window Tuning Results")
    println(io, "")
    println(io, "| Window | Overlap | Total Time (s) | Avg CGM Iters | Max CGM Iters | q_min | q_max |")
    println(io, "|--------|---------|----------------|---------------|---------------|-------|-------|")
    for row in eachrow(results_df)
      println(io, @sprintf("| %d | %d | %.2f | %.1f | %d | %.3e | %.3e |",
        row.window, row.overlap, row.total_time,
        row.avg_cgm_iters, row.max_cgm_iters,
        row.q_min, row.q_max))
    end
  end
end
```

#### 3.2.3 ファイル構成
```
julia/scripts/tune_sliding_window.jl
├── parse_command_line_args()
├── generate_combinations()
├── run_single_job()
├── extract_metrics()
├── write_markdown_report()
└── main()
```

#### 3.2.4 必須コマンドライン引数

| 引数 | 型 | デフォルト値 | 説明 |
|------|---|-------------|------|
| `--nt` | Int | 300 | タイムステップ数 |
| `--windows` | String | "71" | ウィンドウサイズ（カンマ区切り） |
| `--overlap-ratios` | String | - | オーバーラップ比率（カンマ区切り、0.0-1.0） |
| `--overlaps` | String | - | オーバーラップ固定値（カンマ区切り） |
| `--cgm-iter` | Int | 3 | CGM最大反復数 |
| `--solver` | String | "pbicgstab" | ソルバータイプ |
| `--precond` | String | "gs" | 前処理タイプ |
| `--basesize` | Int | 600 | FLoops basesize |
| `--dry-run` | Flag | false | 実験規模プレビューのみ |
| `--reuse` | Flag | false | 既存結果を再利用 |

**所要時間**: 2-3時間

---

### 3.3 Phase 3: 動作確認テスト（優先度：中）

#### 3.3.1 test_dhcp_solver.jl確認
```bash
# 小規模テスト（nt=10）
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/test_dhcp_solver.jl --nt 10 --basesize 600

# 中規模テスト（nt=50）
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/test_dhcp_solver.jl --nt 50 --basesize 600
```

**期待結果**:
- basesize=600がデフォルト適用されている
- 正常に実行完了する
- RMS residual、Total runtimeが出力される

#### 3.3.2 tune_sliding_window.jl確認

**Step 1: Dry-runテスト**
```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/tune_sliding_window.jl \
  --nt 30 \
  --windows 10,20 \
  --overlap-ratios 0.2,0.3 \
  --cgm-iter 2 \
  --dry-run
```

**期待結果**:
- 4つの組み合わせ（10×0.2, 10×0.3, 20×0.2, 20×0.3）が表示される
- 実際の実行は行われない

**Step 2: 小規模実行テスト**
```bash
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/tune_sliding_window.jl \
  --nt 30 \
  --windows 10,20 \
  --overlap-ratios 0.2,0.3 \
  --cgm-iter 2
```

**期待結果**:
- 4つのジョブが順次実行される
- ログファイルが`shared/results/tuning/`に生成される
- CSVとMarkdownレポートが生成される

**Step 3: 結果確認**
```bash
cat shared/results/tuning/julia_sliding_window_scan.csv
cat shared/results/tuning/julia_sliding_window_scan.md
```

**所要時間**: 1時間

---

### 3.4 Phase 4: ドキュメント更新とコミット（優先度：中）

#### 3.4.1 更新対象ドキュメント

**1. README.md（必要に応じて）**
- スクリプト一覧に`tune_sliding_window.jl`を追加

**2. 新規ドキュメント作成（オプション）**
- `docs/guides/SLIDING_WINDOW_TUNING_GUIDE.md`
  - tune_sliding_window.jlの使い方
  - 推奨パラメータ範囲
  - 結果の解釈方法

#### 3.4.2 Gitコミット

**変更ファイル**:
```
M julia/scripts/test_dhcp_solver.jl          # basesize: 1→600
A julia/scripts/tune_sliding_window.jl       # 新規作成
A docs/plans/sliding_window_tuning_implementation_plan.md  # 本計画書
```

**コミットメッセージ案**:
```
feat: スライディングウィンドウ調整システム実装

**実装内容**:
- tune_sliding_window.jl新規作成（複数window/overlap組み合わせ自動評価）
- test_dhcp_solver.jlのbasesize最適化（1→600）
- 開発計画書作成（sliding_window_tuning_implementation_plan.md）

**機能**:
- 複数のwindow/overlap組み合わせを自動実行
- ログから統計値を抽出（CGM反復数、計算時間、熱流束範囲）
- 結果をCSV/Markdown形式で出力
- --dry-runで実験規模プレビュー
- --reuseで既存結果の再利用

**参照**:
- docs/reports/basesize_guidelines.md
- docs/reports/julia_sliding_window_tuning_plan.md
- docs/plans/sliding_window_tuning_implementation_plan.md

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>
```

**所要時間**: 30分

---

## 4. 成果物

### 4.1 コード
- ✅ `julia/scripts/test_dhcp_solver.jl`（basesize=600修正版）
- ✅ `julia/scripts/tune_sliding_window.jl`（新規作成）

### 4.2 ドキュメント
- ✅ `docs/reports/basesize_guidelines.md`（既存、codex作成）
- ✅ `docs/reports/julia_sliding_window_tuning_plan.md`（既存、codex作成）
- ✅ `docs/plans/sliding_window_tuning_implementation_plan.md`（本計画書）
- ⭕ `docs/guides/SLIDING_WINDOW_TUNING_GUIDE.md`（オプション）

### 4.3 結果ファイル
- `shared/results/tuning/julia_sliding_window_scan.csv`
- `shared/results/tuning/julia_sliding_window_scan.md`
- `shared/results/tuning/julia_sw_w*_o*.log`（各組み合わせのログ）
- `shared/results/julia_sw_w*_o*.npz`（各組み合わせの結果）

---

## 5. スケジュール

| Phase | 作業内容 | 所要時間 | 累積時間 |
|-------|---------|---------|---------|
| Phase 1 | test_dhcp_solver.jl修正 | 5分 | 5分 |
| Phase 2 | tune_sliding_window.jl作成 | 2-3時間 | 2-3時間 |
| Phase 3 | 動作確認テスト | 1時間 | 3-4時間 |
| Phase 4 | ドキュメント更新とコミット | 30分 | 3.5-4.5時間 |

**合計所要時間**: 3.5〜4.5時間

---

## 6. リスクと対応

| リスク | 影響 | 対応策 |
|-------|------|--------|
| tune_sliding_window.jlの実装が計画書と乖離 | 中 | 計画書を厳密に参照、実装前にレビュー |
| ログ解析の正規表現が動作しない | 中 | 実際のログファイルで事前検証 |
| 大規模実験で時間がかかりすぎる | 低 | --dry-runで事前確認、小規模から開始 |
| NPZ読み込みでエラー | 低 | NPZファイルの存在確認を先に実施 |

---

## 7. 次のステップ（本開発計画完了後）

### 7.1 実験実行
- 推奨パラメータ範囲での探索実行
  - windows: 50, 71, 100
  - overlap-ratios: 0.2, 0.25, 0.3
  - nt: 100〜300
  - cgm-iter: 3

### 7.2 性能評価
- Python版との比較
- 最適設定での本番実行

### 7.3 ドキュメント完成
- 最適値の確定と推奨設定の文書化
- 実験結果のレポート作成

---

## 8. 承認・レビュー

### 8.1 確認事項
- [ ] Phase 2の実装仕様は計画書の要求を満たしているか？
- [ ] Phase 3のテスト手順は適切か？
- [ ] 追加で必要な機能はあるか？（例: --reuse機能、NPZ解析機能など）

### 8.2 承認
- **承認者**: ユーザー
- **承認日**: 未定
- **ステータス**: 承認待ち

---

## 変更履歴

| 日付 | 版 | 変更内容 | 変更者 |
|------|---|---------|--------|
| 2025-10-27 | 1.0 | 初版作成 | Claude Code |
