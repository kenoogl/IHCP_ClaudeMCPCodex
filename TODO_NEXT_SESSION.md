# 次セッション作業ガイド - 体積積分完全実装プロジェクト

**作成日**: 2025-10-17
**ブランチ**: `tuning7-volume-integral-complete`
**最終コミット**: 4d65584 "feat: CommonSolver.jlに対流項h·Aを実装（Phase 3完了）"

---

## 📊 現在の進捗状況

### ✅ 完了したPhase

#### Phase 1: 準備作業 ✅
- ブランチ作成: `tuning7-volume-integral-complete`
- 参照ドキュメント確認完了

#### Phase 2: RHSCore.jl修正 ✅
- **コミット**: 9815b9d
- `RHS_convection!`関数からθ引数を除去
- RHS項を`h·A·T∞`のみに変更（θ依存項除去）
- 全6面（X±, Y±, Z±）で実施
- 動作確認: コンパイル成功

#### Phase 3: CommonSolver.jl修正 ✅
- **コミット**: 4d65584
- 6つの主要関数にHC引数（熱伝達係数配列）を追加:
  1. `CalcAX!` (CommonSolver.jl:373-417)
  2. `CalcRK!` (CommonSolver.jl:310-357)
  3. `jacobi_preconditioner!` (CommonSolver.jl:607-671)
  4. `resSOR` (CommonSolver.jl:794-844)
  5. `rbsor_core!` (CommonSolver.jl:862-915)
  6. `rbsor!` (CommonSolver.jl:933-955)
- **実装方法**: マスク直接利用方式
  ```julia
  axm = λf(...) * area / dx + HC[1] * area * (1 - m[i-1,j,k])
  ```
- 動作確認: コンパイル成功

---

## 🎯 次セッションの作業内容

### Phase 4: DHCPSolver.jl等でのHC配列生成と伝播

**目的**: 各ソルバーでHC配列を生成し、CommonSolver関数に渡す

#### 4.1 DHCPSolver.jl修正

**ファイル**: `julia/src/solvers/DHCPSolver.jl`

**タスク**:
1. `set_BC_coef`関数を拡張してHC配列を生成
   ```julia
   function set_BC_coef(bc_set::BoundaryConditionSet)
       HF = zeros(6)
       HT = zeros(6)
       HC = zeros(6)  # ← 追加

       # X-minus
       if bc_set.x_minus.type == CONVECTION
           HT[1] = bc_set.x_minus.ambient_temperature
           HC[1] = bc_set.x_minus.heat_transfer_coefficient  # ← 追加
       end
       # ... 他の5面も同様

       return HF, HT, HC  # ← HCを返す
   end
   ```

2. `solve_dhcp!`内でHCをPBiCGSTAB!に渡す
   - 現在の呼び出し: `PBiCGSTAB!(wk, Δh, dt, ZC, dz, ρ, ...)`
   - 修正後: `PBiCGSTAB!(wk, Δh, dt, ZC, dz, ρ, HC, ...)`

3. `calRHS!`内でCalcRK!にHCを渡す

**期待結果**:
- HC配列が正しく生成される
- DHCP解法で対流境界条件が機能する

#### 4.2 PBiCGSTAB!とPCG!のシグネチャ修正

**ファイル**: `julia/src/solvers/CommonSolver.jl`

**タスク**:
1. `PBiCGSTAB!`にHC引数を追加 (CommonSolver.jl:194)
   ```julia
   function PBiCGSTAB!(wk::WorkBuffers,
                       Δh::NTuple{3,T},
                       Δt::T,
                       ZC::AbstractVector{T},
                       ΔZ::AbstractVector{T},
                       ρ::T,
                       HC::AbstractVector{T};  # ← 追加（キーワード引数の前）
                       tol::T = T(1e-6),
                       maxItr::Int = 20_000,
                       smoother::Symbol = :none,
                       par::String = "sequential",
                       verbose::Bool=false) where {T <: AbstractFloat}
   ```

2. `PCG!`にも同様にHC引数を追加 (CommonSolver.jl:43)

3. PBiCGSTAB!内でHCを伝播:
   - `CalcRK!`呼び出し (207行目): HC追加
   - `CalcAX!`呼び出し (254, 263行目): HC追加
   - `Preconditioner!`呼び出し (252, 261行目): HC追加

4. PCG!内でも同様にHCを伝播

#### 4.3 Preconditioner!のシグネチャ修正

**ファイル**: `julia/src/solvers/CommonSolver.jl`

**タスク**:
1. `Preconditioner!`にHC引数を追加 (456行目)
2. `_Preconditioner!`の全バリエーションにHC引数を追加:
   - `:none` (486行目)
   - `:gs` (508行目)
   - `:jacobi` (538行目)
   - エラーハンドラ (563行目)

3. `_Preconditioner!(:gs)`内でrbsor!にHCを渡す:
   ```julia
   function _Preconditioner!(xx, bb, λ, cp, mask, ρ, Δh, Δt, ::Val{:gs}, ZC, ΔZ, HC, par, _)
       for _ in 1:PRECONDITIONER_SWEEPS
           rbsor!(xx, λ, cp, bb, mask, ρ, Δh, Δt, one(ρ), ZC, ΔZ, HC, par)  # ← HC追加
       end
       return nothing
   end
   ```

4. `_Preconditioner!(:jacobi)`内でjacobi_preconditioner!にHCを渡す

#### 4.4 AdjointSolver.jl, SensitivitySolver.jl修正

**ファイル**:
- `julia/src/solvers/AdjointSolver.jl`
- `julia/src/solvers/SensitivitySolver.jl`

**タスク**:
- 随伴問題・感度問題では対流境界を通常使用しない
- HC配列をゼロ配列として渡す:
  ```julia
  HC = zeros(6)
  PBiCGSTAB!(wk, Δh, dt, ZC, dz, ρ, HC, ...)
  ```

---

## 🔍 デバッグ手順

### ステップ1: コンパイルテスト
```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate(); using IHCP_CGM; println("✓ コンパイル成功")'
```

**期待結果**: エラーが出た場合、そのエラーメッセージから修正箇所を特定

### ステップ2: 小規模テスト
```bash
# 5×5×5格子での動作確認
julia --project=julia julia/test/test_dhcp_solver.jl
```

**期待結果**: DHCP解法が正常に動作

### ステップ3: 10ステップテスト
```bash
julia --project=julia julia/scripts/run_10steps_fullsize_test.jl
```

**期待結果**:
- DHCP反復回数: < 500回
- DHCP時間（1ステップ）: < 30秒
- 総実行時間: < 900秒

---

## 📚 重要な参照情報

### HC配列の構成
```julia
HC = [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]
```
- 各要素: 熱伝達係数 [W/m²·K]
- 対流境界条件がない面はゼロ

### 境界条件タイプ
```julia
if bc_set.x_minus.type == CONVECTION
    HC[1] = bc_set.x_minus.heat_transfer_coefficient
    HT[1] = bc_set.x_minus.ambient_temperature
end
```

### マスク直接利用方式の仕組み
```julia
# 内部セル: m[隣接]=1 → (1-m)=0 → 対流項=0
# 境界セル: m[隣接]=0 → (1-m)=1 → 対流項=h·A
axm = λf(...) + HC[1] * area * (1 - m[i-1,j,k])
```

---

## ⚠️ 既知の課題

1. **スケーリング未実装**: Phase 4完了後、Phase 5でスケーリング調整が必要
2. **テスト未実行**: Phase 2,3の変更はコンパイルのみ確認、実行テストは未実施
3. **Preconditioner!の大規模修正**: HC引数の伝播が広範囲に及ぶ

---

## 📝 コミット履歴

```
4d65584 feat: CommonSolver.jlに対流項h·Aを実装（Phase 3完了）
9815b9d fix: RHS対流項をθ非依存に修正（体積積分形式Phase 2完了）
05ebe81 docs: セッション継続情報を更新（体積積分完全実装Phase 2開始前）
f024ee2 docs: 体積積分形式完全実装計画書を作成
```

---

## 🚀 次セッション開始時のコマンド

```bash
# 1. ブランチ確認
git status
git log --oneline -5

# 2. このファイルを確認
cat TODO_NEXT_SESSION.md

# 3. 計画書を確認
cat docs/VOLUME_INTEGRAL_IMPLEMENTATION_PLAN.md

# 4. Phase 4開始
# DHCPSolver.jlのset_BC_coef関数から着手
```

---

**End of Document**
