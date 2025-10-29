# Julia版スライディングウィンドウ調整ガイド

**目的**  
Julia 実装の `run_sliding_window.jl` で、ウィンドウサイズ (`--window`) とオーバーラップ (`--overlap`) を体系的に探索し、反復数・計算時間・熱流束の滑らかさを比較できる仕組みを整える。

---

## 1. ツール概要

| 項目 | 内容 |
|------|------|
| ディレクトリ | `julia/scripts` |
| スクリプト名 | `tune_sliding_window.jl`（新規作成） |
| 役割 | `run_sliding_window.jl` を複数回呼び出し、ウィンドウ設定ごとの統計値を収集 |

### 1.1 スクリプトの主機能
1. **組み合わせ生成**  
   - `--windows 50,71,100` のようにカンマ区切りで候補を指定。  
   - `--overlap-ratios 0.2,0.3` で割合を指定するか、`--overlaps 10,20` のように明示値を使う。割合指定時は `round(Int, window * ratio)`→`max(1, …)` で確定。  
2. **ジョブ実行**  
   - 各組み合わせについて `run_sliding_window.jl` を `Cmd(...)/run(...)` で呼び出し、`--output julia_sw_w{window}_o{overlap}` のように識別名を付ける。  
   - `--dry-run` オプションでコマンドのみ表示可能にしておく。  
3. **結果解析**  
   - ログ（標準出力）から以下を正規表現で抽出:  
     - 総計算時間 (`Summary` セクションの `Total runtime`)  
     - 各ウィンドウの `CGM iterations` と `elapsed`  
     - 熱流束レンジ (`Heat-flux range`)  
   - NPZ も存在すれば `NPZ.npzread` で開き、窓境界の熱流束差 (`window_q_max` / `window_q_min`) を確認。  
4. **サマリ出力**  
   - CSV: `shared/results/tuning/julia_sliding_window_scan.csv`  
   - Markdown: `shared/results/tuning/julia_sliding_window_scan.md`  
   - 各行に `window`, `overlap`, `total_time`, `avg_cgm_iters`, `max_cgm_iters`, `boundary_jump` 等を記録。  
5. **再利用モード**  
   - `--reuse` を指定すると、既存の NPZ / ログがある組み合わせは再実行せず、解析のみ行う。  

---

## 2. 実装メモ

- コマンド例  
  ```julia
  cmd = `$(JOINPATH(@__DIR__, "run_sliding_window.jl"))
          --nt $nt --window $win --overlap $ov
          --cgm-iter $iters --output $tag $extra`
  run(pipeline(cmd, stdout=open(log_path, "w"), stderr=devnull))
  ```
- 出力ログ（`log_path`）から `match(r"CGM iterations: (\d+)"...)` で反復数、`match(r"elapsed: ([\d.]+) s")` で時間を取得。  
- NPZ 読み込み例:  
  ```julia
  result = npzread(joinpath(output_dir, "$tag.npz"))
  window_q_max = result["window_q_max"]
  boundary_jump = maximum(diff(window_q_max))
  ```  
- 実行中断に備え、途中結果を都度 CSV に append するか、メモリ上で `Vector{NamedTuple}` に蓄積してから `CSV.write` する。  
- `--mode python` のような分岐が必要であれば、将来の比較用に引数設計だけ先に入れておくと拡張しやすい。  

---

## 3. 典型的な利用手順

1. **DHCPのみで予備確認**  
   - まず `run_sliding_window.jl --dry-run` で窓の分割を理解し、候補となる `window`/`overlap` を選ぶ。  
   - `run_10steps_fullsize_test.jl` といった DHCP テストスクリプトで、候補ウィンドウ長に対応する `nt` を使い、反復数や計算時間が極端に悪化しない範囲を見つける。これにより CGM まで回して失敗するケースを事前に排除できる。  
2. **候補設定**  
   - 上記で残したウィンドウ長（例: 50, 71, 100）と、20〜30% 程度のオーバーラップ比率（例: 0.2, 0.3）を `tune_sliding_window.jl` の引数に渡す。  
3. **探索実行**  
   ```bash
   JULIA_NUM_THREADS=4 julia --project=julia \
     julia/scripts/tune_sliding_window.jl \
     --nt 120 \
     --windows 50,71,100 \
     --overlap-ratios 0.2,0.3 \
     --cgm-iter 3
   ```
   - 実験規模が大きい場合は `--dry-run` でコマンドだけ先に確認する。  
4. **結果の確認**  
   - `shared/results/tuning/julia_sliding_window_scan.csv` をソートし、`total_time` が短く `avg_cgm_iters` が許容範囲で、かつ `boundary_jump` が小さい組み合わせを選ぶ。  
5. **本番候補の再検証**  
   - 最有力候補で改めて `run_sliding_window.jl` を実行し、`window_q_min/max` や `julia_sliding_window_*_metadata.txt` を確認。  
   - Python 版との比較が必要なら、同じ nt/window/overlap で Python を走らせて結果レンジを照合。  
6. **可視化（任意）**  
   - 別途 `plot_sliding_window_scan.jl` を用意し、CSV を読み込んで散布図（横軸: window, 色: overlap, 縦軸: runtime など）を描けば、調整作業が視覚的に楽になる。  

---

## 4. 今後の拡張案

- **ProcessPool/Distributed** 版のバックエンドと basesize を同時に評価できるよう、`--backend thread|sequential|distributed` を加える。  
- **Python 版連携**: `--mode python` を実装し、同スクリプトから Python 版を呼び出すと比較レポートが自動生成できる。  
- **窓境界評価**: NPZ に保存される `window_histories` を用いて、目的関数の収束形状やヒートフラックス連結の滑らかさを数値化。  
- **UI**: Pluto/VSCode から GUI 的にパラメータ選択→実験→結果確認ができるダッシュボードを作る。

---

**まとめ**
スライディングウィンドウの設定調整を Julia で自動化するには、上記のようなチューニングスクリプトを用意し、ウィンドウ長・オーバーラップの候補を一括実行・比較するのが効率的である。最終的な判断は総実行時間と収束特性、窓境界の滑らかさのバランスで行う。***
