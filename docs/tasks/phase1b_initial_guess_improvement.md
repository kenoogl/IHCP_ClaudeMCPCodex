# Task: Phase 1-B - 初期推定値改善

**作成日**: 2025年10月12日
**担当**: codex
**ブランチ**: tuning6
**優先度**: ★★★★☆

---

## 1. 目的

CG反復における初期推定値を改善することで、収束速度を向上させる。現在は全てのタイムステップで初期値=0から始めているが、前ステップの解や時系列パターンを活用することで、より良い初期推定値を提供し、CG反復回数を削減する。

**期待される効果**:
- CG反復回数15-25%削減
- 総実行時間20-30秒短縮（1,067秒 → 1,040秒程度）
- 特にAdjoint問題での大幅な改善（ステップ9の1348回 → 800回以下）

---

## 2. 対象ファイル

変更が必要なファイルのリスト:

- `julia/src/solvers/DHCPSolver.jl`: 前ステップ解の初期値利用、線形・放物型外挿ロジック実装
- `julia/src/solvers/AdjointSolver.jl`: 後退時間積分における初期値特別対応（従来法/残差ベースの選択可能化）
- `julia/src/solvers/SensitivitySolver.jl`: 前ステップ解の初期値利用（DHCPと同様）
- `julia/src/solvers/SlidingWindowSolver.jl`: スライディングウィンドウ境界での初期値継承ロジック（実装のみ、テストは今後のTODO）

---

## 3. 要件

### 必須要件

1. **機能要件 - DHCP/Sensitivityソルバー**:
   - [ ] 前ステップ解を次ステップの初期推定値として使用
   - [ ] 線形外挿による予測初期値（t>=3で有効化）: `T(t) ≈ 2*T(t-1) - T(t-2)`
   - [ ] 放物型外挿による予測初期値（t>=4で有効化）: `T(t) ≈ 3*T(t-1) - 3*T(t-2) + T(t-3)`
   - [ ] 外挿方法を選択するキーワード引数実装（`:none`, `:linear`, `:quadratic`）

2. **機能要件 - Adjointソルバー**:
   - [ ] 後退時間積分における次ステップ（時間的に後）の解を初期値として使用（従来法）
   - [ ] 最終ステップ（t=nt）では残差ベースの初期値を使用（新規）
   - [ ] 初期値戦略を選択するキーワード引数実装（`:previous`, `:residual`）
   - [ ] `:previous`: 従来の方法（次ステップの解、t=ntのみゼロ初期値）
   - [ ] `:residual`: 残差ベース（t=ntで残差の負、他のステップは次ステップの解）

3. **機能要件 - スライディングウィンドウ**:
   - [ ] ウィンドウ境界での初期値継承ロジック実装
   - [ ] 前ウィンドウ最終ステップの解を次ウィンドウ初期ステップの初期値として使用
   - [ ] **注**: テストは今後のTODO、実装のみ完了させる

4. **品質要件**:
   - [ ] 既存テスト505件全通過
   - [ ] 後方互換性維持（デフォルト: use_previous_solution=true、extrapolation=:none）
   - [ ] API変更は最小限（キーワード引数追加のみ）

5. **性能要件**:
   - [ ] DHCP反復回数: 560 → 450回/ステップ以下（20%削減）
   - [ ] Adjoint反復回数: 715 → 550回/ステップ以下（23%削減）
   - [ ] Adjoint ステップ9: 1348 → 800回以下（40%削減）
   - [ ] 総実行時間: 1067秒 → 1040秒以下（2.5%改善）
   - [ ] 数値計算の正確性維持（収束解の相対誤差1%以内）

### オプション要件

- [ ] 初期推定値戦略の効果を可視化するデバッグモード
- [ ] 各ステップでの初期残差をロギング

---

## 4. 期待される出力

codexが生成すべきもの:

1. **修正済みコード**:
   - `julia/src/solvers/DHCPSolver.jl`: 前ステップ解利用、線形・放物型外挿ロジック
   - `julia/src/solvers/AdjointSolver.jl`: 初期値戦略選択機能（`:previous`, `:residual`）
   - `julia/src/solvers/SensitivitySolver.jl`: 前ステップ解利用、外挿ロジック
   - `julia/src/solvers/SlidingWindowSolver.jl`: ウィンドウ境界初期値継承（実装のみ）

2. **実装レポート**:
   - 場所: `docs/reports/phase1b/initial_guess_improvement_report.md`
   - 内容:
     - 実装の詳細（各ソルバーの変更内容）
     - 外挿方法の説明（線形・放物型の数式、精度、安定性）
     - Adjoint初期値戦略の比較（`:previous` vs `:residual`）
     - 変更箇所の説明（コード差分、ロジック図）
     - 既知の問題・制約（外挿使用時の安定性条件など）
     - テスト結果サマリー（505件全通過確認）
     - 初期残差の変化（ステップ毎の初期残差比較表）
     - スライディングウィンドウのTODO事項（テスト未実装の旨を明記）

3. **その他**（必要に応じて）:
   - ベンチマークスクリプト（初期推定値戦略の効果測定用）

---

## 5. 技術的詳細

### 現状の問題点

1. **全ステップで初期値=0から開始**:
   - 時系列的に連続している温度場を解くにも関わらず、前ステップの情報を活用していない
   - CG法は初期値に依存する反復解法なので、良い初期推定値があれば収束が早い

2. **Adjoint問題での高反復回数**:
   - ステップ9（後退時間積分の最初、t=nt→t=nt-1）: 1348回
   - 初期残差: 6.02e-05（比較的大きい）
   - 残差情報を初期値に反映していない

3. **スライディングウィンドウ境界での非効率**:
   - ウィンドウ間での解の連続性を活用していない
   - 各ウィンドウで初期値=0からリスタート

### 改善アプローチ

#### 5-A: DHCPソルバーの改善

**前ステップ解の利用 + 外挿法（線形・放物型）**:

```julia
# julia/src/solvers/DHCPSolver.jl
function solve_dhcp!(
    T_init::AbstractArray{T,3},
    q::AbstractArray{T,2},
    work::WorkBuffers{T},
    nt::Int,
    # ... 他のパラメータ ...
    ;
    rtol::T = T(1e-6),
    maxiter::Int = 20000,
    smoother::Symbol = :none,
    use_previous_solution::Bool = true,      # 新規: デフォルトON
    extrapolation::Symbol = :none            # 新規: :none, :linear, :quadratic
) where T <: AbstractFloat

    T_all = zeros(T, ni, nj, nk, nt)
    T_all[:, :, :, 1] .= T_init

    iter_counts = zeros(Int, nt)

    for t in 2:nt
        # 初期推定値の設定
        if extrapolation == :quadratic && t >= 4
            # 放物型外挿（2次精度）: T(t) ≈ 3*T(t-1) - 3*T(t-2) + T(t-3)
            @views @. work.Δh = 3 * T_all[:, :, :, t-1] - 3 * T_all[:, :, :, t-2] + T_all[:, :, :, t-3]
        elseif extrapolation == :linear && t >= 3
            # 線形外挿（1次精度）: T(t) ≈ 2*T(t-1) - T(t-2)
            @views @. work.Δh = 2 * T_all[:, :, :, t-1] - T_all[:, :, :, t-2]
        elseif use_previous_solution && t >= 2
            # 前ステップ解: T(t) ≈ T(t-1)
            @views @. work.Δh = T_all[:, :, :, t-1]
        else
            # デフォルト: ゼロ初期値
            fill!(work.Δh, zero(T))
        end

        # PBiCGSTAB!呼び出し（work.Δhが初期値兼解ベクトル）
        isconverged, itr, res0 = PBiCGSTAB!(
            work, work.Δh, dt, Z, ΔZ, z_range, HT, ρ;
            tol=rtol, maxItr=maxiter, smoother=smoother, par=1
        )

        if !isconverged
            @warn "DHCP solver did not converge" t itr res0
        end

        # 解を保存
        T_all[:, :, :, t] .= work.Δh
        iter_counts[t] = itr
    end

    return T_all, iter_counts
end
```

#### 5-B: Adjointソルバーの改善

**初期値戦略の選択機能（`:previous` vs `:residual`）**:

```julia
# julia/src/solvers/AdjointSolver.jl
function compute_gradient!(
    T_cal::AbstractArray{T,4},
    Y_obs::AbstractArray{T,3},
    work::WorkBuffers{T},
    # ... 他のパラメータ ...
    ;
    rtol::T = T(1e-6),
    maxiter::Int = 20000,
    smoother::Symbol = :none,
    initial_guess_strategy::Symbol = :residual   # 新規: :previous または :residual
) where T <: AbstractFloat

    λ = zeros(T, ni, nj, nk, nt)
    grad = zeros(T, ni, nj, nt)
    grad_iters = zeros(Int, nt)

    # 後退時間積分: t=nt → t=2
    for t in nt:-1:2
        # 初期推定値の設定
        if initial_guess_strategy == :residual
            # 残差ベース戦略
            if t == nt
                # 最終ステップ（後退の最初）: 残差の負を初期値
                # ∂J/∂T|_{t=nt} = -(Y_obs - T_cal)
                residual = @view T_cal[:, :, nk, t]  # 表面温度
                fill!(work.Δh, zero(T))
                @views @. work.Δh[:, :, nk] = Y_obs[:, :, t] - residual
            else
                # 通常ステップ: 次ステップ（時間的に後）の解を初期値
                @views @. work.Δh = λ[:, :, :, t+1]
            end
        elseif initial_guess_strategy == :previous
            # 従来の戦略
            if t == nt
                # 最終ステップ: ゼロ初期値
                fill!(work.Δh, zero(T))
            else
                # 通常ステップ: 次ステップの解を初期値
                @views @. work.Δh = λ[:, :, :, t+1]
            end
        else
            error("Unknown initial_guess_strategy: $initial_guess_strategy. Use :previous or :residual")
        end

        # Adjoint方程式を解く
        isconverged, itr, res0 = PBiCGSTAB!(
            work, work.Δh, dt, Z, ΔZ, z_range, HT, ρ;
            tol=rtol, maxItr=maxiter, smoother=smoother, par=1
        )

        if !isconverged
            @warn "Adjoint solver did not converge" t itr res0
        end

        # 解を保存
        λ[:, :, :, t] .= work.Δh
        grad_iters[t] = itr

        # 勾配を計算（表面での熱流束感度）
        @views @. grad[:, :, t] = -k_coeffs[:, :, nk, t] * λ[:, :, nk] / ΔZ[nk]
    end

    return grad, grad_iters
end
```

#### 5-C: Sensitivityソルバーの改善

DHCPソルバーと同様のロジックを適用（線形・放物型外挿サポート）:

```julia
# julia/src/solvers/SensitivitySolver.jl
function compute_sensitivity!(
    # ... パラメータ ...
    ;
    use_previous_solution::Bool = true,      # 新規
    extrapolation::Symbol = :none            # 新規: :none, :linear, :quadratic
) where T <: AbstractFloat

    # DHCPと同様のロジック
    for t in 2:nt
        # 初期推定値の設定（DHCPと同じif-elseif-else構造）
        if extrapolation == :quadratic && t >= 4
            @views @. work.Δh = 3 * S_all[:, :, :, t-1] - 3 * S_all[:, :, :, t-2] + S_all[:, :, :, t-3]
        elseif extrapolation == :linear && t >= 3
            @views @. work.Δh = 2 * S_all[:, :, :, t-1] - S_all[:, :, :, t-2]
        elseif use_previous_solution && t >= 2
            @views @. work.Δh = S_all[:, :, :, t-1]
        else
            fill!(work.Δh, zero(T))
        end

        # Sensitivity方程式を解く
        # ...
    end

    return S_all, sens_iters
end
```

#### 5-D: スライディングウィンドウの改善（実装のみ、テストは今後のTODO）

```julia
# julia/src/solvers/SlidingWindowSolver.jl
function sliding_window_cgm!(
    # ... パラメータ ...
    ;
    window_size::Int = 10,
    use_window_continuation::Bool = true     # 新規: ウィンドウ間初期値継承
) where T <: AbstractFloat

    # ... ウィンドウループ ...

    # 前ウィンドウの最終解を保存
    previous_window_final_T = nothing
    previous_window_final_q = nothing

    for window_idx in 1:num_windows
        # ウィンドウ範囲の計算
        t_start = (window_idx - 1) * window_size + 1
        t_end = min(t_start + window_size - 1, total_steps)

        # 初期値の設定
        if use_window_continuation && window_idx > 1 && previous_window_final_T !== nothing
            # 前ウィンドウの最終解を初期値として使用
            T_init_window = previous_window_final_T
            q_init_guess = previous_window_final_q
        else
            # デフォルト初期値
            T_init_window = T_init
            q_init_guess = zeros(T, ni, nj)
        end

        # ウィンドウ内のCGM計算
        # ...

        # 最終解を保存（次のウィンドウ用）
        previous_window_final_T = T_result[:, :, :, end]
        previous_window_final_q = q_result[:, :, end]
    end

    # TODO: スライディングウィンドウの初期値継承機能のテスト未実装
    # - ウィンドウ境界での初期残差検証
    # - 複数ウィンドウでの収束性確認
    # - 継承ON/OFFでの性能比較

    return results
end
```

### 外挿法の数式と精度

#### 線形外挿（1次精度）
```
T(t) ≈ 2*T(t-1) - T(t-2)

テイラー展開による誤差: O(Δt²)
安定性: 比較的安定、発散リスク低
適用条件: t >= 3
```

#### 放物型外挿（2次精度）
```
T(t) ≈ 3*T(t-1) - 3*T(t-2) + T(t-3)

テイラー展開による誤差: O(Δt³)
安定性: 線形より不安定、滑らかな時間発展で有効
適用条件: t >= 4
```

### 検証方法

1. **初期残差のログ取得**:
   ```julia
   # PBiCGSTAB!内で初期残差を記録
   res0 = sqrt(Fdot2(wk.pcg_r, wk.pcg_r, z_range, par))
   println("Step $t: Initial residual = $res0")
   ```

2. **反復回数の比較**:
   - ベースライン（extrapolation=:none）vs 線形外挿 vs 放物型外挿
   - Adjoint: `:previous` vs `:residual`
   - 各ステップでの反復回数を記録
   - 特にAdjoint問題のステップ9に注目

3. **収束解の精度検証**:
   - 改善版とベースラインの解の相対誤差を計算
   - 数値精度が維持されていることを確認（相対誤差1%以内）

4. **ベンチマーク実行**:
   ```bash
   julia --project=julia julia/scripts/run_10steps_fullsize_test.jl
   ```

---

## 6. 参考資料

### ドキュメント

- `docs/performance_improvement_proposals.md`: 提案4（行207-295）参照
- `shared/results/performance_22fde2d.md`: ベースライン性能
- `shared/results/performance_tuning5_type_stability.md`: Phase 1-A結果

### コード

- `julia/src/solvers/DHCPSolver.jl`: 直接解法実装（~280行）
- `julia/src/solvers/AdjointSolver.jl`: 随伴解法実装（~290行）
- `julia/src/solvers/SensitivitySolver.jl`: 感度解法実装（~270行）
- `julia/src/solvers/CommonSolver.jl`: PBiCGSTAB!実装（~300行）
- `julia/src/solvers/SlidingWindowSolver.jl`: スライディングウィンドウ実装
- `julia/test/runtests.jl`: テストスイート（505項目）

### 外部資料

- Iterative Methods for Linear Systems: 初期値の重要性
- Conjugate Gradient Method: 初期推定値と収束速度の関係
- Time-Stepping Methods: 時系列解における外挿法（線形・放物型）
- Numerical Recipes: Polynomial Extrapolation

---

## 7. 制約条件

### 守るべき制約

1. **後方互換性**:
   - デフォルト動作は改善版（`use_previous_solution=true`, `extrapolation=:none`）
   - 既存テストが変更なしで通過すること
   - キーワード引数でOFF可能

2. **数値精度**:
   - 収束解の相対誤差1%以内
   - 収束判定基準は変更しない（rtol=1e-6維持）
   - CG反復回数が削減されても、最終解の精度は維持

3. **テスト**:
   - 505件のテスト全通過必須
   - スライディングウィンドウの初期値継承機能はテスト未実装でOK（TODOとして明記）

### 避けるべきこと

- [ ] 3次以上の高次外挿（数値的に不安定）
- [ ] 不安定な外挿（発散の可能性がある場合はフォールバック機構検討）
- [ ] 初期値計算での新たなアロケーション（インプレース演算必須）
- [ ] スライディングウィンドウのテスト無しでのマージ（実装のみでOK、TODOを明記）

---

## 8. 完了基準

### Definition of Done

- [ ] 全ての必須要件を満たす
- [ ] テスト505件全通過
- [ ] 実装レポート作成完了（初期残差比較表、反復回数グラフ、外挿法の説明含む）
- [ ] コードレビュー可能な状態
- [ ] ベンチマーク結果記録（CG反復回数削減、総実行時間改善）
- [ ] スライディングウィンドウのTODO事項をレポートに明記

### 確認事項

- [ ] `git status`でuntracked filesなし（意図的なもの除く）
- [ ] コンパイル警告なし
- [ ] 各ソルバーの初期推定値ロジックが正しく動作（ログ確認）
- [ ] Adjoint問題のステップ9で反復回数が大幅に削減されている（1348 → 800回以下）
- [ ] 外挿方法（`:none`, `:linear`, `:quadratic`）が正しく動作
- [ ] Adjoint初期値戦略（`:previous`, `:residual`）が正しく動作

---

## 9. 備考

### 実装上の注意点

1. **外挿法の選択**:
   - デフォルトは`:none`（安全性重視）
   - `:linear`: 1次精度、比較的安定
   - `:quadratic`: 2次精度、高速収束が期待できるが不安定になる可能性
   - ユーザーが明示的に選択

2. **Adjoint問題の初期値戦略**:
   - `:previous`: 従来の方法（安全、実績あり）
   - `:residual`: 残差ベース（ステップ9での高速収束が期待）
   - デフォルトは`:residual`（性能改善が目的）

3. **スライディングウィンドウ**:
   - 実装のみ完了、テストは今後のTODO
   - レポートにTODO事項を明記
   - 実装が正しく動作するかの最小限の確認は実施

4. **メモリ効率**:
   - 初期値設定でのアロケーション発生を避ける
   - `@views`と`@.`マクロでインプレース演算

### デバッグのヒント

- 各ステップの初期残差をロギングして、初期推定値の効果を確認
- `extrapolation=:none`で従来動作と比較
- Adjoint問題のステップ9は必ずチェック（最大のボトルネック）
- 外挿法での発散検出: 反復回数が前ステップより大幅に増加した場合は警告

### スライディングウィンドウのTODO（レポートに記載）

今後実装すべきテスト項目:
- [ ] ウィンドウ境界での初期残差検証テスト
- [ ] 複数ウィンドウでの収束性確認テスト
- [ ] 継承ON/OFFでの性能比較ベンチマーク
- [ ] エッジケース: 最初/最後のウィンドウの動作確認

---

## 10. タイムライン（目安）

- **準備**: 15分（タスクブリーフィング作成、ブランチ作成） ← 完了
- **実装**: 60-75分（codex）
  - DHCPソルバー改善（線形・放物型外挿）: 20分
  - Adjointソルバー改善（初期値戦略選択）: 25分
  - Sensitivityソルバー改善: 15分
  - スライディングウィンドウ（実装のみ）: 10分
  - テスト実行: 15分
- **検証**: 45分（レビュー、ベンチマーク、レポート作成）
- **合計**: 約2-2.5時間

---

**更新履歴**:
- 2025-10-12 21:15: 初版作成（Claude Code）
- 2025-10-12 21:20: ユーザー要件反映（外挿2種類、Adjoint初期値選択、スライディングウィンドウTODO化）
