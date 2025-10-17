適応的収束判定（PBICGSTAB向け）
===============================

背景
----
- 対象: 逆熱伝導問題 (IHCP) の時間発展ソルバー（三種類: DHCP, Adjoint, Sensitivity）。
- 現状: 固定値の相対残差許容値 `tol = 1e-6 (DHCP/Sensitivity)`、`1e-8 (Adjoint)`。
- 課題: 時間ステップ後半で温度場が安定し、固定 `tol` では反復回数が過剰。

理論的根拠
----------
- **時間離散化誤差との整合 (Knoll & Keyes 2004; Saad 2003)**  
  BDF/Crank–Nicolson などで線形ソルバ残差を局所離散化誤差と同程度に抑えるのが理想。p 次精度なら許容誤差を `O(Δt^{p+1})` に設定。
- **Eisenstat–Walker inexact Newton (Eisenstat & Walker 1996)**  
  前ステップの残差情報を用いて許容値を調整する枠組み。時間発展問題でも「求める解に依存して tol を緩めすぎない」理論的裏付けを提供。
- **温度場変化量と収束性 (Elman, Silvester, Wathen 2014)**  
  定常化領域では解の変化量が小さく、変化量に比例した収束基準が合理的。

基準タイプ別の適応式
--------------------

1. **解の相対変化ベース**
   ```
   δₙ = ||θⁿ - θⁿ⁻¹||₂ / (||θⁿ||₂ + ε),  ε = 1e-12
   tolₙ = clamp( κ_θ · δ̄ₙ, tol_min, tol_max )
   ```
   - `κ_θ ∈ [0.1, 0.5]`、`δ̄ₙ` は移動平均 (2–3 ステップ)。

2. **残差履歴ベース**
   ```
   tolₙ = clamp( κ_r · ρ̄_{n-1}^α, tol_min, tol_max )
   ```
   - `ρ̄_{n-1}`: 前ステップ最終残差の移動平均、`κ_r ≈ 0.5`、`α ∈ [0.5,1.0]`。

3. **時間刻み依存**
   ```
   tolₙ = clamp( κ_dt · (Δtₙ / Δt_ref)^{p+1}, tol_min, tol_max )
   ```
   - p: 時間積分の次数（例: Crank–Nicolson で p = 2）。`κ_dt` は離散化誤差目標。`Δt_ref` は代表刻み幅。

4. **複合基準（推奨）**
   ```
   tolₙ_raw = max( c₁ κ_θ δ̄ₙ, c₂ κ_r ρ̄_{n-1}^α, c₃ κ_dt (Δtₙ/Δt_ref)^{p+1} )
   tolₙ = clamp( η tol_{n-1} + (1-η) tolₙ_raw, tol_min, tol_max )
   ```
   - `η ≈ 0.7`（平滑化）、`c_i` は 1.0 で開始し実測に合わせて調整。


既存ソルバー実装例
------------------
- **PETSc**: `KSPSetConvergenceTest` や `SNES` の inexact Newton で Eisenstat–Walker 型の tol 更新が可能。ユーザ定義で残差比に応じて `tol` を動的変更する例が多数。
- **Trilinos (NOX/Belos)**: NOX でインエグザクト戦略 (Eisenstat–Walker II) が利用可能。Belos はステップごとに「目標残差ノルム」を再設定するインタフェースを提供。
- **deal.II**: `SolverControl` + `TimeStepper` を組み合わせ、時間刻みと過去の反復数から tol を変更するチュートリアルが公開されている。

推奨パラメータ範囲
------------------
- `tol_min`: DHCP/Sensitivity で `1e-7`、Adjoint で `1e-9` を下限。
- `tol_max`: 温度誤差許容（例: 0.01 K）から逆算して `1e-5` 前後。
- `κ_θ`: 0.2（初期値）。  
  `κ_r`: 0.5、`α`: 0.7（Suresh & Tamma 1998 の報告値）。  
  `κ_dt`: `0.5 * tol_ref`（既存固定 tol を基準）。  
- 平滑化係数 `η`: 0.6–0.8。

実装時の注意
------------
- **モニタ更新順序**: `θⁿ` または `θⁿ-θⁿ⁻¹` を評価 → `tolₙ` 算出 → 次ステップの `PBiCGSTAB!` 呼び出し前に渡す。
- **ゼロ除算防止**: `||θⁿ||₂ < θ_ref` (例: `1e-8`) の場合は絶対ノルムを用いる。
- **フェイルセーフ**: 反復上限に達したら `tolₙ ← max(tolₙ / 10, tol_min)` で再試行。
- **スレッド安全性**: norm/残差のリダクション後に単一スレッドで `tolₙ` を更新。
- **精度保証**: 最終時刻や観測時刻では強制的に `tol_min` を使用する、または `tolₙ ≤ tol_target` を維持。

参考文献
--------
- Saad, Y. *Iterative Methods for Sparse Linear Systems*, 2nd ed., 2003.
-.Eisenstat, S. C., and Walker, H. F. “Choosing the Forcing Terms in an Inexact Newton Method”, 1996.
- Knoll, D. A., and Keyes, D. E., “Jacobian-Free Newton–Krylov Methods: A Survey”, *J. Comp. Phys.*, 2004.
- Elman, H. C., Silvester, D. J., and Wathen, A. J., *Finite Elements and Fast Iterative Solvers*, 2nd ed., 2014.
- Suresh, S., and Tamma, K. K., “Inexact Krylov methods for transient heat conduction”, *Int. J. Numer. Methods Eng.*, 1998.
