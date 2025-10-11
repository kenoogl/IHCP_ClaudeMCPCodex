# Task: Phase 1-A - 型安定化

**作成日**: 2025年10月12日
**担当**: codex
**ブランチ**: tuning5
**優先度**: ★★★★★

---

## 1. 目的

Juliaの型推論を最適化し、動的ディスパッチを排除することで、コンパイラ最適化を促進し実行時間を削減する。

**期待される効果**:
- 実行時間10-20%削減
- 動的ディスパッチの排除
- コンパイラによるSIMD最適化の促進
- メモリアクセスパターンの最適化

---

## 2. 対象ファイル

変更が必要なファイルのリスト:

- `julia/src/solvers/CommonSolver.jl`: 型アノテーション追加、Val型による関数特殊化
- `julia/src/solvers/DHCPSolver.jl`: 型アノテーション追加
- `julia/src/solvers/AdjointSolver.jl`: 型アノテーション追加
- `julia/src/solvers/SensitivitySolver.jl`: 型アノテーション追加
- `julia/src/solvers/CGMSolver.jl`: 型アノテーション追加
- `julia/scripts/check_type_stability.jl`: 型安定性チェックスクリプト（新規作成）

---

## 3. 要件

### 必須要件

1. **機能要件**:
   - [ ] 全ての主要関数に適切な型アノテーションを追加
   - [ ] `smoother`引数をVal型による関数特殊化に変更
   - [ ] 型パラメータ`T <: AbstractFloat`を全ソルバーに導入
   - [ ] 内部変数の型推論失敗を排除

2. **品質要件**:
   - [ ] 既存テスト505件全通過
   - [ ] 後方互換性維持（API変更なし）
   - [ ] `@code_warntype`で警告なし

3. **性能要件**:
   - [ ] ベースライン性能を維持または改善（10-20%改善目標）
   - [ ] CG反復回数は変わらない（±1%以内）
   - [ ] 動的ディスパッチゼロ（`@report_opt`で確認）

### オプション要件

- [ ] 型安定性チェックスクリプト作成
- [ ] 実装レポート作成
- [ ] ベンチマーク結果記録

---

## 4. 期待される出力

codexが生成すべきもの:

1. **修正済みコード**:
   - 全ソルバーファイルに型アノテーション追加
   - Val型による関数特殊化実装
   - 型安定性チェックスクリプト

2. **実装レポート**:
   - 場所: `docs/reports/phase1/type_stability_report.md`
   - 内容:
     - 実装の詳細
     - 変更箇所の説明
     - 型不安定箇所の修正内容
     - テスト結果サマリー
     - ベンチマーク結果（実行時間、動的ディスパッチ確認）

3. **その他**:
   - 型安定性チェックスクリプト（`julia/scripts/check_type_stability.jl`）

---

## 5. 技術的詳細

### 現状の問題点

- 型推論失敗による動的ディスパッチの可能性
- 不要な型変換とボクシング
- コンパイラ最適化の阻害
- `smoother`引数が文字列型 → ランタイムディスパッチ

### 改善アプローチ

#### 5-1: 型パラメータの導入

```julia
# 変更前
function PBiCGSTAB!(
    wk::WorkBuffers,
    Δh::AbstractArray,
    dt::Float64,
    Z::AbstractVector,
    ...
)

# 変更後
function PBiCGSTAB!(
    wk::WorkBuffers,
    Δh::AbstractArray{T,3},
    dt::T,
    Z::AbstractVector{T},
    ΔZ::AbstractVector{T},
    z_range::UnitRange{Int},
    HT::AbstractArray{T,3},
    ρ::AbstractArray{T,3};
    tol::T = T(1e-6),
    maxItr::Int = 20000,
    smoother::Symbol = :none,
    par::Int = 1
) where T <: AbstractFloat
    # 内部変数も型注釈
    rho::T = zero(T)
    alpha::T = zero(T)
    omega::T = zero(T)

    # ...
end
```

#### 5-2: Val型による関数特殊化

```julia
# smoother引数のコンパイル時特殊化
@inline function apply_smoother!(
    z::AbstractArray{T,3},
    r::AbstractArray{T,3},
    ::Val{:none},
    args...
) where T
    mycopy!(z, r)
end

@inline function apply_smoother!(
    z::AbstractArray{T,3},
    r::AbstractArray{T,3},
    ::Val{:gs},
    args...
) where T
    GS_smoother_in_place!(z, r, args...)
end

# 呼び出し側
function PBiCGSTAB!(...; smoother::Symbol = :none, ...)
    smoother_val = Val(smoother)
    # ...
    apply_smoother!(wk.pcg_z, wk.pcg_r, smoother_val, ...)
end
```

#### 5-3: 型安定性チェックスクリプト

```julia
# julia/scripts/check_type_stability.jl
using JET

function check_solver_type_stability()
    # テストデータ作成
    ni, nj, nk = 10, 10, 10
    T_init = zeros(Float64, ni, nj, nk)
    q = zeros(Float64, ni, nj, 1)
    work = WorkBuffers(ni, nj, nk)

    # DHCPソルバーの型推論チェック
    @report_opt solve_dhcp!(T_init, q, work, 2, 7900.0, [1.0, 0.0], [1.0, 0.0],
                            0.12e-3, 0.12e-3, [0.0, 0.7e-3], [0.7e-3], 0.001;
                            rtol=1e-6, maxiter=20000)

    # PBiCGSTAB!の型推論チェック
    @report_opt PBiCGSTAB!(work, zeros(ni,nj,nk), 0.001, [0.0, 0.7e-3], [0.7e-3],
                           1:nk, zeros(ni,nj,nk), zeros(ni,nj,nk);
                           tol=1e-6, maxItr=20000, smoother=:gs)
end

check_solver_type_stability()
```

### 検証方法

1. **コンパイル時検証**:
   ```julia
   # 型推論チェック
   julia --project=julia julia/scripts/check_type_stability.jl
   ```

2. **テスト実行**:
   ```bash
   julia --project=julia julia/test/runtests.jl
   ```

3. **ベンチマーク**:
   ```bash
   julia --project=julia julia/scripts/run_10steps_fullsize_test.jl
   ```

---

## 6. 参考資料

### ドキュメント

- `docs/performance_improvement_proposals.md`: 提案2参照（164-262行）
- `shared/results/performance_22fde2d.md`: ベースライン性能（1072.43秒）
- `shared/results/performance_tuning4_allocation_reduction.md`: Phase 0結果（1068.74秒）

### コード

- `julia/test/runtests.jl`: テストスイート（505件）
- `julia/src/solvers/`: 既存ソルバー実装
- `julia/scripts/run_10steps_fullsize_test.jl`: ベンチマークスクリプト

### 外部資料

- Julia Performance Tips: https://docs.julialang.org/en/v1/manual/performance-tips/
- JET.jl Documentation: https://github.com/aviatesk/JET.jl

---

## 7. 制約条件

### 守るべき制約

1. **後方互換性**:
   - API変更禁止（関数シグネチャ維持）
   - デフォルト引数の挙動維持
   - 既存の呼び出しコードが動作すること

2. **数値精度**:
   - CG反復回数が変わらない（±1%以内）
   - 収束挙動の維持
   - Python版との誤差2%以内維持

3. **テスト**:
   - 505件のテスト全通過必須
   - 新規テスト不要（型安定性は実装の内部改善）

### 避けるべきこと

- [ ] 過度な型制約（柔軟性の喪失）
- [ ] 複雑な型階層の導入
- [ ] 実行時の型チェック追加（オーバーヘッド）
- [ ] API変更

---

## 8. 完了基準

### Definition of Done

- [ ] 全ての必須要件を満たす
- [ ] テスト505件全通過
- [ ] 実装レポート作成完了
- [ ] 型安定性チェックスクリプト動作確認
- [ ] ベンチマーク実行（10ステップフルサイズ）
- [ ] 動的ディスパッチゼロ確認（`@report_opt`）

### 確認事項

- [ ] `git status`でuntracked filesなし（意図的なもの除く）
- [ ] コンパイル警告なし
- [ ] `@code_warntype`で警告なし（主要関数）
- [ ] 実行時間がベースライン以下（目標: 10-20%改善）

---

## 9. 備考

### 重要な注意点

1. **smoother引数の変更**:
   - 文字列 `"none"`, `"gs"` → シンボル `:none`, `:gs`
   - テストコードも同様に変更必要
   - 後方互換性のため、文字列も受け付けるラッパー検討

2. **型パラメータT**:
   - `Float64`に固定せず、`T <: AbstractFloat`で抽象化
   - `Float32`への将来的な対応も視野

3. **JET.jlのインストール**:
   ```julia
   using Pkg
   Pkg.add("JET")
   ```

### Phase 0との関連

Phase 0（配列アロケーション削減）の成果を保持しつつ、型安定化を追加実装する。両者は独立しており、競合しない。

---

## 10. タイムライン（目安）

- **準備**: 15分（タスクブリーフィング作成、ブランチ作成） ← 現在ここ
- **実装**: 45-60分（codex）
  - 型アノテーション追加: 30分
  - Val型特殊化: 15分
  - チェックスクリプト: 15分
- **検証**: 30-45分（テスト、型チェック、ベンチマーク、レビュー）
- **合計**: 約1.5-2時間

---

**更新履歴**:
- 2025-10-12: 初版作成
