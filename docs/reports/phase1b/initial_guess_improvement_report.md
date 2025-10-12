# Phase 1-B – 初期推定値改善レポート

## 1. 実装概要
- DHCP/Sensitivity ソルバーに時間方向ホットスタートと線形/放物型外挿を導入し、`:none`/`:linear`/`:quadratic` をキーワードで切り替え可能にした。
- Adjoint ソルバーに `:previous` / `:residual` 初期値戦略を実装し、最終ステップで残差ベース初期化を行うロジックを追加した。
- SlidingWindow ソルバーにウィンドウ間継承を制御する `use_window_continuation` を追加し、最終温度場・熱流束を次ウィンドウへ受け渡す実装を整備した。
- CGM 上位ループから各ソルバーの初期化戦略をパラメータ指定できるよう拡張し、既存テストが旧参照データと整合するようにデフォルトを上書き可能にした。

## 2. ソルバー別の変更点
### DHCP / Sensitivity ソルバー
- `apply_temporal_initial_guess!` でガイドセルを含む初期値を一括制御。`use_previous_solution=false` で従来のゼロ初期化にフォールバック。
- `:quadratic` では `T(t) ≈ 3T_{t-1} - 3T_{t-2} + T_{t-3}`、` :linear` では `T(t) ≈ 2T_{t-1} - T_{t-2}`。履歴が不足するステップでは自動的に低次へフォールバックする。
- 不正なシンボルを受け取った場合は `ArgumentError` を送出し、明示的な API エラーとした。

### Adjoint ソルバー
- `apply_adjoint_initial_guess!` により、最終ステップ(`t=nt-1`)では残差 `Y_obs - T_cal` を底面に配置、それ以外のステップでは次時刻の随伴場をコピーする。
- 新キーワード `initial_strategy` (既定 `:residual`) を追加し、`compute_gradient!` と `solve_cgm!` から透過的に設定できるようにした。
- テストベンチでは旧 Python 参照値と比較するため `:previous` を指定し、既存のゴールデンデータとの整合性を維持。

### SlidingWindow ソルバー
- `use_window_continuation` 既定値 `true` を導入。オンの場合は直前ウィンドウの最終温度場と熱流束シーケンスを継承し、ウィンドウ境界での初期残差低減を狙う。
- オフの場合は従来どおり定常初期値にリセットされるため、既存参照データと突き合わせるテストでは `false` を指定。
- 既存テストは数値比較の `@test_broken` を維持し、Phase 5 の参照データ更新 TODO を継続記録。

## 3. 外挿法の整理
| モード | 数式 | 有効開始ステップ | 近似次数 | 安定性メモ |
| --- | --- | --- | --- | --- |
| `:none` | \(T^{(t)} = T^{(t-1)}\) | すべて | 0 次 | 最も安全。推定が悪化する場合のフォールバック先。 |
| `:linear` | \(T^{(t)} \approx 2T^{(t-1)} - T^{(t-2)}\) | \(t \ge 3\) | 1 次 | 連続時間解が滑らかな場合に有効。急峻な変化では振動に注意。 |
| `:quadratic` | \(T^{(t)} \approx 3T^{(t-1)} - 3T^{(t-2)} + T^{(t-3)}\) | \(t \ge 4\) | 2 次 | 最も高速だが、ステップ幅が大きい場合は過大予測に注意。コード上は履歴不足時に自動で低次へフォールバック。 |

## 4. 随伴初期値戦略の比較
| 戦略 | 初期化 | 利点 | 留意点 |
| --- | --- | --- | --- |
| `:previous` | 次ステップの随伴場をコピー。最終ステップのみゼロ。 | 既存実装と互換。安定性が高い。 | 大きな残差があるステップでは収束までの反復が増えやすい。 |
| `:residual` (既定) | 最終ステップで \(Y_{obs} - T_{cal}\) を底面に配置。それ以外は `:previous` と同じ。 | 初回反復の初期残差を小さくし、ボトルネック時刻の反復削減が期待できる。 | 小規模問題では `:previous` の方が安定な場合があるため、キーワードで切り替え。 |

## 5. コード変更箇所
- `julia/src/solvers/DHCPSolver.jl:135` – `apply_temporal_initial_guess!` 追加と外挿ハンドリング。
- `julia/src/solvers/SensitivitySolver.jl:131` – DHCP と同等の初期値制御ロジックを導入し、API を拡張。
- `julia/src/solvers/AdjointSolver.jl:44` – `apply_adjoint_initial_guess!` と `initial_strategy` を実装。
- `julia/src/solvers/CGMSolver.jl:258` – CGM パラメータに初期値戦略を受け渡すフックを追加。
- `julia/src/solvers/SlidingWindowSolver.jl:137` – `use_window_continuation` による継承制御とコピーロジックを追加。
- `julia/test/test_cgm_solver.jl:121`, `julia/test/test_sliding_window.jl:269` – 旧参照データとの整合性確保のためのパラメータ指定と `@test_broken` の調整。

## 6. テスト結果
- 実行コマンド: `julia --project=julia -e 'using Pkg; Pkg.test()'`
- 結果: 全 505 テストが成功（Phase 5 数値比較は仕様どおり `@test_broken` を継続）。
- テスト出力に追加されたログ: CGM/SlidingWindow の初期化戦略切り替えに伴う情報出力が増加。

## 7. 初期残差の観測
| ステップ | DHCP res₀ | Adjoint res₀ | Sensitivity res₀ | 備考 |
| --- | --- | --- | --- | --- |
| 1 | **TODO** | **TODO** | **TODO** | 現状 `PBiCGSTAB!` の戻り値を呼び出し側で記録していないため、将来の計測ポイントとして残す。 |

> メモ: 初期残差ログの恒常化は Phase 1-C での観測基盤整備の TODO として連携。

## 8. 既知の課題 / TODO
- 高次外挿は急峻な加熱条件で発散する可能性があるため、今後の検証で安全装置（クリッピング等）を検討。
- SlidingWindow の数値比較テストは参照データが旧アルゴリズム由来のため `@test_broken` を維持。`use_window_continuation` を活かした新ベンチマーク生成が必要。
- 初期残差の自動収集・レポート化は未着手。WorkBuffers へスタブ領域を追加しロギングする案が残課題。
- `:residual` 戦略は小規模問題での挙動が未検証なため、将来的に自動切替ロジック（閾値判定）を検討。

