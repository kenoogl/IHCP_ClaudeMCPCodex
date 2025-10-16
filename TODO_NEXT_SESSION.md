# 次セッション作業ガイド - 体積積分完全実装プロジェクト

**作成日**: 2025-10-17
**ブランチ**: `tuning7-volume-integral-complete`
**最終更新**: 2025-10-17（Phase 4完了）

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
  1. `CalcAX!` (CommonSolver.jl:375-417)
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

#### Phase 4: HC配列生成と伝播 ✅
**実装完了**: 2025-10-17

##### 4.1 BoundaryConditions.jlに`set_BC_coef`関数を追加 ✅
- **ファイル**: `julia/src/utils/boundary_conditions.jl`
- 境界条件セットから熱伝達係数配列HC [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]を生成
- 対流境界条件の場合: heat_transfer_coefficient
- それ以外の境界条件: 0.0
- **設計判断**: DHCPSolver専用ではなくBoundaryConditions.jlに配置（他ソルバーでも再利用可能）

##### 4.2 PBiCGSTAB!とPCG!のシグネチャにHC引数を追加 ✅
- **ファイル**: `julia/src/solvers/CommonSolver.jl`
- `PBiCGSTAB!`にHC引数を追加（行195-206）
- `PCG!`にHC引数を追加（行43-54）
- CalcRK!、CalcAX!、Preconditioner!へのHC伝播を実装

##### 4.3 Preconditioner!の全バリエーションにHC引数を追加 ✅
- **ファイル**: `julia/src/solvers/CommonSolver.jl`
- `Preconditioner!`公開API（行458-473）
- `_Preconditioner!`の全実装:
  - `:none` 恒等変換（行489）
  - `:gs` Gauss-Seidel法（行511-515）
  - `:jacobi` Jacobi法（行541-556）
  - エラーハンドラ（行567）
- `rbsor!`と`jacobi_preconditioner!`へのHC伝播

##### 4.4 DHCPSolver.jlでHC配列を生成しsolve_linear_system!に渡す ✅
- **ファイル**: `julia/src/solvers/DHCPSolver.jl`
- `set_BC_coef`をimport（行31）
- 境界条件セットからHC配列を生成（行245）
- `solve_linear_system!`にHCを渡す（行279）

##### 4.5 AdjointSolver.jlとSensitivitySolver.jlでゼロHC配列を渡す ✅
- **ファイル**:
  - `julia/src/solvers/AdjointSolver.jl`（行43, 283, 320）
  - `julia/src/solvers/SensitivitySolver.jl`（行27, 241, 275）
- 両ソルバーで`set_BC_coef`をimport
- HC配列を生成（対流境界なし→ゼロ配列）
- `solve_linear_system!`にHCを渡す

##### 4.6 solve_linear_system!のシグネチャ修正 ✅
- **ファイル**: `julia/src/solvers/CommonSolver.jl`（行1051-1075）
- HC引数を追加（キーワード引数の前）
- PBiCGSTAB!とPCG!にHCを伝播

**テスト結果**:
- ✅ コンパイルテスト成功
- ✅ 全テスト実行: 179 passed, 2 broken（既知の課題）

---

## 🎯 次セッションの作業内容

### Phase 5: スケーリング調整とテスト

**目的**: 対流境界条件の係数スケーリングを調整し、数値的安定性を確保

#### 5.1 スケーリング戦略の決定

**課題**:
- 現在の実装: `HC[i] * area * (1 - m[隣接])`
- 熱伝達係数hのオーダー: 10〜1000 W/m²·K
- 面積areaのオーダー: 10^-8〜10^-6 m²
- 熱伝導項のオーダー: λ/dx（約10^-2〜10^2）

**検討事項**:
1. 係数の相対的な大きさを確認
2. 条件数への影響を評価
3. 必要に応じてスケーリング係数を導入

#### 5.2 数値テスト

**テストケース**:
1. **対流境界なしケース**（既存テストで確認済み）
   - HC = [0, 0, 0, 0, 0, 0]
   - 期待: Phase 1-6と同等の結果

2. **対流境界ありケース**（新規）
   - X-minus面に対流境界条件を設定
   - h = 10 W/m²·K, T_∞ = 300 K
   - 期待: 温度場が物理的に妥当

3. **スケール感度テスト**
   - h = 1, 10, 100, 1000 W/m²·Kで実行
   - 収束性と精度を評価

#### 5.3 必要に応じた修正

**修正候補**:
- スケーリング係数の導入
- 前処理器の調整
- 収束判定基準の見直し

---

## 🔍 デバッグ手順

### ステップ1: 現在の動作確認
```bash
# 既存テストが全て通ることを確認（Phase 4で実施済み）
julia --project=julia julia/test/runtests.jl
```

### ステップ2: 対流境界テストケース作成
```bash
# 新しいテストファイルを作成
julia --project=julia julia/test/test_convection_bc.jl
```

### ステップ3: スケール感度テスト
```bash
# 異なる熱伝達係数でテスト
julia --project=julia julia/scripts/test_convection_scaling.jl
```

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
axm = λf(...) * area / dx + HC[1] * area * (1 - m[i-1,j,k])
```

### 対流項の物理的意味
- 体積積分形式: ∫_Ω ∇·(h·A·T) dV = ∫_∂Ω h·(T - T_∞) dS
- 離散化: h·A·(θ[境界] - T_∞)
- RHS項: -h·A·T_∞
- LHS項（対角）: +h·A（マスク利用で自動的に境界でのみ適用）

---

## ⚠️ 既知の課題

1. **スケーリング未検証**: Phase 5で対流境界条件のスケーリングを確認する必要あり
2. **対流境界テスト未作成**: 実際の対流境界条件を使ったテストケースが未実装
3. **条件数への影響**: 熱伝達係数hの値によって条件数が変化する可能性

---

## 📝 コミット履歴

```
[未コミット] Phase 4完了: HC配列生成と伝播実装
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

# 4. Phase 4の変更をコミット
git add .
git commit -m "feat: Phase 4完了 - HC配列生成と伝播実装

体積積分形式の対流境界条件実装を完了。

**Phase 4.1**: BoundaryConditions.jlにset_BC_coef関数追加
- 境界条件セットからHC配列 [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]を生成
- 対流境界以外はゼロ配列

**Phase 4.2**: PBiCGSTAB!とPCG!のシグネチャにHC引数追加
- CalcRK!, CalcAX!, Preconditioner!へのHC伝播

**Phase 4.3**: Preconditioner!の全バリエーションにHC引数追加
- :none, :gs, :jacobi全てに対応
- rbsor!, jacobi_preconditioner!へのHC伝播

**Phase 4.4**: DHCPSolver.jlでHC配列生成とsolve_linear_system!に渡す
- set_BC_coefで境界条件からHC配列を生成

**Phase 4.5**: AdjointSolver.jl, SensitivitySolver.jlでゼロHC配列を渡す
- 随伴・感度問題では対流境界なし（HC=0配列）

**テスト結果**:
- コンパイルテスト: 成功 ✓
- 全テスト実行: 179 passed, 2 broken（既知の課題）✓

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>"

# 5. Phase 5開始
# 対流境界条件のスケーリングテストケースを作成
```

---

**End of Document**
