# CGM性能改善レポート（配列アロケーション削減）

## 1. レポート概要

- **目的**: Julia版CGM逆解析ループにおける大量アロケーションを抑制し、反復ごとのメモリ再利用を可能にする。
- **対象範囲**: `julia/src/solvers/CGMSolver.jl`, `julia/src/solvers/DHCPSolver.jl`, `julia/src/solvers/SensitivitySolver.jl`, `julia/src/solvers/AdjointSolver.jl`
- **実施者**: Codex（GPT-5 ベース）
- **作業日**: セッション実行日（Codex CLI上）

## 2. 主要変更点

- **CGMループのインプレース化**  
  - `julia/src/solvers/CGMSolver.jl:272` にて残差・勾配・探索方向用バッファを事前確保。  
  - `julia/src/solvers/CGMSolver.jl:325` および `julia/src/solvers/CGMSolver.jl:352` でビュー演算と代入を用い、ドット積以外の一時配列生成を抑制。  
  - `julia/src/solvers/CGMSolver.jl:365` 以降では共役方向更新と熱流束更新をすべて`@.`マクロによるインプレース演算に変更。

- **ソルバーAPIのバッファ再利用対応**  
  - DHCP/感度ソルバー (`DHCPSolver.jl:149`, `SensitivitySolver.jl:145`) と随伴ソルバー (`AdjointSolver.jl:185`) に再利用バッファを受け取るキーワード引数を追加。  
  - 不一致時は明示的な`ArgumentError`を投げ、誤ったサイズのバッファ利用を防止。

- **勾配計算バッファの再活用**  
  - `compute_gradient!`で随伴場・勾配を外部から供給されたバッファに書き戻すよう変更 (`CGMSolver.jl:120`)。

- **テンソル内積の汎用化・高速化**  
  - `tensor_dot`を抽象配列対応＋`@simd`ループに刷新し、`vec`による一時配列生成を回避 (`CGMSolver.jl:65`)。  
  - `compute_step_size`も抽象配列を受け取る形に更新し、ビュー渡し時の型制約を緩和 (`CGMSolver.jl:170`)。

## 3. 期待効果と考察

- CGM外層ループでの大規模アロケーション（温度場・随伴場・感度場など）が一度だけ行われるため、GC負荷と割り当て時間が大幅に削減される見込み。  
- DHCP/Adjoint/Sensitivityソルバーが外部バッファを再利用可能になったことで、スライディングウィンドウや長時間計算でもアロケーション削減の恩恵を受けられる。  
- `tensor_dot`の再実装により、配列ビューをそのまま渡しても余分なコピーが発生しないため、ステップサイズ計算のパスでも不要なメモリ消費を防止。  
- 変更に伴いAPIが拡張されているため、外部から直接ソルバーを呼び出すユースケースでは新キーワード引数を利用しない限り既存挙動を維持（デフォルトで新規バッファを確保する）ことを確認済み。

## 4. テスト結果

- `julia --project=julia -e 'using Pkg; Pkg.test()'`  
  - 505件の既存テストを全て通過。  
  - テスト実行時に作業環境のJuliaバージョンが 1.12.0 であることから、`julia/Manifest.toml` が同バージョン向けに更新されている点に注意。

## 5. 今後の推奨アクション

1. `BenchmarkTools` を用い、代表的なCGMケースでアロケーション数と実行時間がどの程度削減されたかを定量的に測定する。  
2. バッファサイズが可変となる構成（格子数や時間ステップが変化）での安全性を追加テストし、必要に応じて初期化ヘルパーを導入する。  
3. さらなる最適化候補（例: `WorkBuffers` の内部配列共有化、BLASスレッド設定の自動化）について、`docs/performance_improvement_proposals.md` へ追記する。

