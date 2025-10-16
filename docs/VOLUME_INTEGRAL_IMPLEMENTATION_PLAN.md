# 体積積分形式完全実装計画書

**作成日**: 2025年10月16日
**目的**: 287c8c7で導入された体積積分形式を完全に実装し、性能劣化を解消する
**対象ブランチ**: tuning7 → tuning7-volume-integral-complete

---

## 📋 エグゼクティブサマリー

### 問題の本質
コミット287c8c7で左辺（係数行列）を体積積分形式に変更したが、右辺（RHS）と前処理が旧形式のまま残り、以下の問題が発生：

1. **単位の不整合**: 左辺 [W/K] vs 右辺 [W/m³] → **12桁のスケール不整合**
2. **係数行列の条件数悪化**: BiCGSTABの収束が極端に遅化
3. **前処理の機能不全**: Jacobi/GSが新しいスケールに対応していない

### 解決方針
体積積分形式への完全移行 + スケーリング調整により、以下を達成：

- ✅ 左辺と右辺の単位整合（全て [W/K]）
- ✅ 対流境界条件の陰的実装（h·A·θ を左辺に）
- ✅ 係数行列の適切なスケーリング（条件数改善）
- ✅ 前処理の最適化（新しいスケールに対応）

### 開発方針
- codexと連携

#### 役割分担
- claude : タスク全体の計画と進捗管理、ユーザーとのコミュニケーション、Codexへの指示出しと結果のレビュー, Juliaの実行環境やテスト実行のサポート
- codex : 詳細な調査作業、コードの実装、ファイル検索や分析
---

## 🔍 Codex調査報告のサマリー

### 係数スケールの比較（現在の格子）
| 項目 | d9e1b02（旧形式）| 287c8c7（新形式）| スケール変化 |
|------|-----------------|-----------------|-------------|
| 拡散項 λ/dx² | 1.04×10⁹ | 5.22×10⁻⁴ | **10⁻¹³倍** |
| 時間項 1/Δt | 1.0×10¹ | 2.00×10⁻⁵ | **10⁻⁶倍** |
| 熱流束項 q | 4.0×10⁴ | 2.00×10⁻⁸ | **10⁻¹³倍** |

格子パラメータ: dx=1.2e-4 m, dy=1.67e-4 m, dz=2.5e-5 m

### 離散化方程式の形式

**d9e1b02（単位体積あたり形式）**:
```
[(1/Δt) + Σ(λf/(ρcp)Δn²)] θᵢ − Σ(λf/(ρcp)Δn²) θ_nb = −θ_old/Δt − (q + s)/(ρcp)
```
- 単位: 両辺とも [K·s⁻¹]
- 係数: O(10⁹) ~ O(10¹)

**287c8c7（体積積分形式、不完全）**:
```
[ρcp·V/Δt + Σ λf·(A/Δn)] θᵢ − Σ λf·(A/Δn) θ_nb = −ρcp·V·θ_old/Δt − (q/Δn + s)V
```
- 左辺: [W/K], 右辺: [W·m⁻³] ← **不整合！**
- 係数: O(10⁻⁴) ~ O(10⁻⁵)

**目標（体積積分形式、完全）**:
```
[ρcp·V/Δt + Σ λf·(A/Δn) + Σ h·A] θᵢ − Σ λf·(A/Δn) θ_nb
  = −ρcp·V·θ_old/Δt − q·A − h·A·T∞ − s·V
```
- 単位: 両辺とも [W/K] ✓
- 対流項が陰的に実装 ✓
- スケーリング調整が必要 → [K·s⁻¹] 相当に正規化

---

## 🎯 実装計画

### Phase 1: 準備作業

#### 1.1 作業環境のセットアップ
```bash
# 現在: d9e1b02 (detached HEAD)
git checkout tuning7
git stash pop  # 前回の変更を適用（必要に応じて）
git checkout -b tuning7-volume-integral-complete
```

#### 1.2 参照ファイルの確認
- `julia/src/solvers/CommonSolver.jl` (287c8c7時点の状態)
- `julia/src/solvers/RHSCore.jl` (142b94c時点の状態)
- `julia/src/solvers/DHCPSolver.jl`
- `julia/src/solvers/AdjointSolver.jl`
- `julia/src/solvers/SensitivitySolver.jl`

---

### Phase 2: RHSCore.jl の修正

**目的**: 対流境界条件から θ を除去（h·A·T∞ のみをRHSに）

#### 2.1 RHS_convection! 関数の修正

**ファイル**: `julia/src/solvers/RHSCore.jl`
**行番号**: 181-233

**修正内容**:

```julia
# 修正前
function RHS_convection!(b::AbstractArray{T,3},
                        θ::AbstractArray{T,3},  # ← θを使用
                        dx::T,
                        dy::T,
                        dz::AbstractVector{T},
                        bc::BoundaryCondition,
                        face_type::Symbol,
                        par::String)::Nothing where {T <: AbstractFloat}
  # ...
  if face_type == :x_minus
    let i=2, h = bc.heat_transfer_coefficient, ta = bc.ambient_temperature
      @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        area = dy * dz[k]
        b[i,j,k] -= h * area * (ta - θ[i,j,k])  # ← θを含む
      end
    end
  end
  # ... 他の5面も同様
end

# 修正後
function RHS_convection!(b::AbstractArray{T,3},
                        # θ引数を削除
                        dx::T,
                        dy::T,
                        dz::AbstractVector{T},
                        bc::BoundaryCondition,
                        face_type::Symbol,
                        par::String)::Nothing where {T <: AbstractFloat}
  # ...
  if face_type == :x_minus
    let i=2, h = bc.heat_transfer_coefficient, ta = bc.ambient_temperature
      @floop backend for k in 2:SZ[3]-1, j in 2:SZ[2]-1
        area = dy * dz[k]
        b[i,j,k] -= h * area * ta  # ← θを除去、T∞のみ
      end
    end
  end
  # ... 他の5面も同様（全てθを除去）
end
```

**影響範囲**: 6面全ての対流境界条件（X±, Y±, Z±）

#### 2.2 呼び出し側の修正

`RHS_convection!` を呼び出している箇所で θ 引数を削除：

**BoundaryConditions.jl** (apply_boundary_conditions! 内):
```julia
# 修正前
RHS_convection!(b, θ, dx, dy, dz, bc, face_type, par)

# 修正後
RHS_convection!(b, dx, dy, dz, bc, face_type, par)
```

---

### Phase 3: CommonSolver.jl の修正

**目的**:
1. 対流項 h·A を係数行列の対角成分に追加
2. スケーリング調整を実装

#### 3.1 境界条件情報の伝達機構

**問題**: 現在、CommonSolver は境界条件の情報を持っていない

**解決策**: HC配列（Heat transfer Coefficient）を追加

```julia
# 現在の関数シグネチャ
function CalcAX!(ax, θ, Δh, Δt, λ, cp, m, ρ, Z, ΔZ, z_range, HT, par)
  # HT = [HT_xm, HT_xp, HT_ym, HT_yp, HT_zm, HT_zp]  # 熱流束密度

# 修正後
function CalcAX!(ax, θ, Δh, Δt, λ, cp, m, ρ, Z, ΔZ, z_range, HT, HC, par)
  # HT = [HT_xm, HT_xp, HT_ym, HT_yp, HT_zm, HT_zp]  # 熱流束密度
  # HC = [h_xm, h_xp, h_ym, h_yp, h_zm, h_zp]        # 熱伝達係数（新規）
```

#### 3.2 CalcAX! の修正

**ファイル**: `julia/src/solvers/CommonSolver.jl`
**行番号**: 385-434

**修正内容**:

```julia
function CalcAX!(ax::AbstractArray{T,3},
                  θ::AbstractArray{T,3},
                  Δh::NTuple{3,T},
                  Δt::T,
                  λ::AbstractArray{T,3},
                  cp::AbstractArray{T,3},
                  m::AbstractArray{T,3},
                  ρ::T,
                  Z::AbstractVector{T},
                  ΔZ::AbstractVector{T},
                  z_range::AbstractVector{<:Integer},
                  HT::AbstractVector{T},
                  HC::AbstractVector{T},  # ← 追加
                  par::String) where {T <: AbstractFloat}
    backend = get_backend(par)
    SZ = size(θ)
    dx0 = Δh[1]
    dy0 = Δh[2]
    # ... (既存のコード)

    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        dz_k = ΔZ[k]
        λ0 = λ[i,j,k]
        m0 = m[i,j,k]

        # 拡散項（既存）
        axm = λf(λ[i-1,j,k], λ0, m[i-1,j,k], m0) * dy0 * dz_k / dx0
        axp = λf(λ[i+1,j,k], λ0, m[i+1,j,k], m0) * dy0 * dz_k / dx0
        aym = λf(λ[i,j-1,k], λ0, m[i,j-1,k], m0) * dx0 * dz_k / dy0
        ayp = λf(λ[i,j+1,k], λ0, m[i,j+1,k], m0) * dx0 * dz_k / dy0
        azm = λf(λ[i,j,k-1], λ0, m[i,j,k-1], m0) * dx0 * dy0 / (Z[k]-Z[k-1])
        azp = λf(λ[i,j,k+1], λ0, m[i,j,k+1], m0) * dx0 * dy0 / (Z[k+1]-Z[k])

        # 時間項（既存）
        a_p_0 = ρ * cp[i,j,k] * dx0 * dy0 * dz_k / Δt

        # 対流項（新規追加）★
        h_conv = zero(T)

        # X-minus面（i=2）
        if i == 2 && HC[1] != zero(T)
            mw = one(T) - m[i-1,j,k]
            h_conv += HC[1] * dy0 * dz_k * mw
        end

        # X-plus面（i=SZ[1]-1）
        if i == SZ[1]-1 && HC[2] != zero(T)
            me = one(T) - m[i+1,j,k]
            h_conv += HC[2] * dy0 * dz_k * me
        end

        # Y-minus面（j=2）
        if j == 2 && HC[3] != zero(T)
            ms = one(T) - m[i,j-1,k]
            h_conv += HC[3] * dx0 * dz_k * ms
        end

        # Y-plus面（j=SZ[2]-1）
        if j == SZ[2]-1 && HC[4] != zero(T)
            mn = one(T) - m[i,j+1,k]
            h_conv += HC[4] * dx0 * dz_k * mn
        end

        # Z-minus面（k=z_st）
        if k == z_st && HC[5] != zero(T)
            mb = m[i,j,k-1]
            h_conv += HC[5] * dx0 * dy0 * (one(T) - mb)
        end

        # Z-plus面（k=z_ed）
        if k == z_ed && HC[6] != zero(T)
            mt = m[i,j,k+1]
            h_conv += HC[6] * dx0 * dy0 * (one(T) - mt)
        end

        # 対角成分（対流項を追加）★
        dd = (one(T)-m0) + (axp + axm + ayp + aym + azp + azm + a_p_0 + h_conv)*m0

        ss = ( axp * θ[i+1,j  ,k  ] + axm * θ[i-1,j  ,k  ]
             + ayp * θ[i  ,j+1,k  ] + aym * θ[i  ,j-1,k  ]
             + azp * θ[i  ,j  ,k+1] + azm * θ[i  ,j  ,k-1] )
        ax[i,j,k] = (ss - dd*θ[i,j,k]) * m0
    end
end
```

**重要**: CalcRK!, jacobi_preconditioner!, rbsor_core! にも同様の修正が必要

#### 3.3 CalcRK! の修正

**ファイル**: `julia/src/solvers/CommonSolver.jl`
**行番号**: 311-366

CalcAX! と同じロジックで対流項を追加

#### 3.4 jacobi_preconditioner! の修正

**ファイル**: `julia/src/solvers/CommonSolver.jl`
**行番号**: 605-674

Jacobi前処理の対角項に対流項を追加

#### 3.5 rbsor_core! の修正

**ファイル**: `julia/src/solvers/CommonSolver.jl`
**行番号**: 856-914

SOR前処理の対角項に対流項を追加

---

### Phase 4: スケーリング調整の実装

**目的**: 係数行列の条件数を改善（10⁻¹³ → 1 程度に正規化）

#### 4.1 スケーリング戦略

**方法**: 体積積分形式の係数を単位体積あたり形式に正規化

```julia
scale = 1.0 / (ρ * cp_avg * volume_ref)
```

これにより:
- 係数: O(10⁻⁴) → O(10⁹) に戻る
- 前処理の効果が復活

#### 4.2 PBiCGSTAB! 内でのスケーリング

**ファイル**: `julia/src/solvers/CommonSolver.jl`
**関数**: `PBiCGSTAB!`

**実装方針**:
1. RHSを事前にスケーリング
2. CalcAX!, CalcRK! 内で同じスケールを適用
3. 解は元のスケールに戻す必要なし（温度は無次元化されていない）

**実装例**:

```julia
function PBiCGSTAB!(wk::WorkBuffers,
                    Δh::NTuple{3,T},
                    Δt::T,
                    Z::AbstractVector{T},
                    ΔZ::AbstractVector{T},
                    z_range::AbstractVector{<:Integer},
                    HT::AbstractVector{T},
                    HC::AbstractVector{T},  # ← 追加
                    ρ::T;
                    tol::T = T(1e-6),
                    maxItr::Int = 20_000,
                    smoother::Symbol = :none,
                    par::String = "sequential",
                    verbose::Bool=false) where {T <: AbstractFloat}

    # スケーリング係数を計算
    dx0, dy0 = Δh[1], Δh[2]
    z_st, z_ed = Int(first(z_range)), Int(last(z_range))
    dz_avg = mean(ΔZ[z_st:z_ed])
    volume_ref = dx0 * dy0 * dz_avg

    # 平均比熱を計算（内点のみ）
    SZ = size(wk.cp)
    cp_sum = zero(T)
    cp_count = 0
    for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        cp_sum += wk.cp[i,j,k]
        cp_count += 1
    end
    cp_avg = cp_sum / T(cp_count)

    # スケーリング係数（単位体積あたり形式に戻す）
    scale = inv(ρ * cp_avg * volume_ref)

    if verbose
        println("Scaling factor: $scale (cp_avg=$cp_avg, vol=$volume_ref)")
    end

    # RHSをスケーリング
    myfill!(wk.pcg_q, zero(T), par)
    SZ = size(wk.θ)
    backend = get_backend(par)
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1, i in 2:SZ[1]-1
        wk.b[i,j,k] *= scale
    end

    # 残差計算（スケーリング適用）
    res0 = CalcRK!(wk.pcg_r, wk.θ, wk.b, wk.λ, wk.cp, wk.mask, ρ, Δh, Δt, Z, ΔZ, z_range, HT, HC, scale, par)

    # ... (BiCGSTABループは既存のまま)

    # CalcAX!, CalcRK! にscale引数を追加
end
```

#### 4.3 CalcAX!, CalcRK! でのスケーリング適用

各関数に `scale` 引数を追加し、最終結果に乗算：

```julia
function CalcAX!(ax, θ, Δh, Δt, λ, cp, m, ρ, Z, ΔZ, z_range, HT, HC, scale, par)
    # ... (既存の計算)
    ax[i,j,k] = (ss - dd*θ[i,j,k]) * m0 * scale  # ← scale適用
end
```

---

### Phase 5: 各ソルバーでの統合

#### 5.1 DHCPSolver.jl

**修正箇所**:
1. `set_BC_coef` の拡張 → HC配列の生成
2. `solve_dhcp!` 内で HC を PBiCGSTAB! に渡す

```julia
function set_BC_coef(bc_set::BoundaryConditionSet)
    HF = zeros(6)
    HT = zeros(6)
    HC = zeros(6)  # ← 追加

    # X-minus
    if bc_set.x_minus.type == HEAT_FLUX
        HF[1] = bc_set.x_minus.value
    elseif bc_set.x_minus.type == CONVECTION
        HT[1] = bc_set.x_minus.ambient_temperature
        HC[1] = bc_set.x_minus.heat_transfer_coefficient  # ← 追加
    end

    # ... 他の5面も同様

    return HF, HT, HC  # ← HC を返す
end
```

#### 5.2 AdjointSolver.jl

随伴問題では対流境界は通常使用しない → HC = zeros(6) で対応

#### 5.3 SensitivitySolver.jl

感度問題も対流境界は通常使用しない → HC = zeros(6) で対応

---

### Phase 6: テストと検証

#### 6.1 単体テスト（新規作成）

**ファイル**: `julia/test/test_volume_integral_scaling.jl`

**テスト内容**:
1. 1セルでの対流境界条件テスト
2. スケーリング前後での解の一致確認
3. 係数行列の条件数測定

#### 6.2 小規模テスト

```bash
# 5×5×5格子での動作確認
julia --project=julia -e 'include("julia/test/test_small_grid.jl")'
```

期待結果:
- BiCGSTAB反復回数が大幅減少
- 収束時間が短縮

#### 6.3 10ステップテスト

```bash
# 80×100×20格子での性能確認
julia --project=julia julia/scripts/run_10steps_fullsize_test.jl
```

期待結果:
- DHCP: 各ステップ20-30秒で収束（d9e1b02と同等）
- 総実行時間: 900秒以内

---

## 📊 成功基準

### 定量的指標

| 指標 | 現状（142b94c）| 目標 | 測定方法 |
|------|----------------|------|---------|
| DHCP反復回数（平均） | > 10000 | < 500 | 10ステップテスト |
| DHCP時間（1ステップ） | > 90秒 | < 30秒 | 10ステップテスト |
| 総実行時間 | タイムアウト | < 900秒 | 10ステップテスト |
| 温度場誤差（RMS） | N/A | < 10 K | Python比較 |

### 定性的指標

- ✅ 左辺と右辺の単位が整合している
- ✅ 対流境界条件が陰的に実装されている
- ✅ 係数行列のスケールが適切である
- ✅ コードの可読性が保たれている

---

## ⚠️ リスクと対策

### リスク1: 実装の複雑さ

**影響**: バグ混入の可能性、デバッグに時間
**対策**:
- 段階的実装（Phase 2 → 3 → 4）
- 各Phaseでテスト実行
- コード自動生成を活用（対流項の6面分）

### リスク2: スケーリングの副作用

**影響**: 別の数値不安定性が発生する可能性
**対策**:
- スケーリング係数をログ出力
- verbose=true でBiCGSTAB反復を監視
- 異常検知時は即座にロールバック

### リスク3: 既存テストの破壊

**影響**: Phase 1-6の505テストが失敗する可能性
**対策**:
- 修正前に全テストを実行してベースライン確認
- 各Phase完了後に回帰テスト実行
- 失敗時は該当Phaseのみロールバック

---

## 📅 実装スケジュール

### Day 1（約4時間） ✅ 完了
- [x] Phase 1: 準備作業（30分）
- [x] Phase 2: RHSCore.jl修正（1時間）
  - RHS_convection! 関数修正 ✅
  - 呼び出し側修正 ✅
  - 動作確認 ✅
  - コミット: 9815b9d
- [x] Phase 3.1-3.2: CalcAX!修正（1.5時間）
  - HC配列の追加 ✅
  - 対流項の実装（マスク直接利用方式） ✅
  - テスト ✅

### Day 2（約4時間） ✅ 完了
- [x] Phase 3.3-3.5: 残りの修正（2時間）
  - CalcRK!修正 ✅
  - jacobi_preconditioner!修正 ✅
  - rbsor_core!修正 ✅
  - resSOR修正 ✅
  - rbsor!修正 ✅
  - コミット: 4d65584
- [ ] Phase 4: DHCPSolver等での統合（次セッション）
  - set_BC_coef拡張（HC配列生成）
  - PBiCGSTAB!シグネチャ修正
  - Preconditioner!へのHC伝播

### Day 3（約3時間）
- [ ] Phase 5: ソルバー統合（1.5時間）
  - DHCPSolver修正
  - AdjointSolver確認
  - SensitivitySolver確認
- [ ] Phase 6: テスト（1.5時間）
  - 小規模テスト
  - 10ステップテスト
  - 結果検証

**合計見積もり**: 11時間（実作業）

---

## 🔄 ロールバック計画

問題発生時の対応：

### レベル1: 当該Phaseのみロールバック
```bash
git diff HEAD julia/src/solvers/CommonSolver.jl > phase3_backup.patch
git checkout HEAD julia/src/solvers/CommonSolver.jl
```

### レベル2: ブランチ全体を破棄
```bash
git checkout tuning7
git branch -D tuning7-volume-integral-complete
```

### レベル3: 287c8c7を撤回（最終手段）
```bash
git checkout tuning7
git revert 287c8c7
git commit -m "Revert volume integral form due to unresolved issues"
```

---

## 📚 参考資料

### 関連コミット
- `d9e1b02`: Phase 1-E（正常動作版）
- `287c8c7`: 体積積分形式統一（問題発生）
- `142b94c`: RHS修正試行（不完全）

### Codex調査レポート
- 詳細分析結果は本プロジェクト直前に取得
- 係数スケール、離散化方程式、推奨修正方法を含む

### 技術文献
- 有限体積法における保存則の離散化
- BiCGSTABソルバーの前処理技術
- マトリクスフリー実装のスケーリング手法

---

## ✅ チェックリスト

実装完了時の確認項目：

### コード品質
- [ ] 全ファイルがコンパイル可能
- [ ] 警告メッセージがない
- [ ] コメントが適切に記述されている
- [ ] コードスタイルが統一されている

### 機能性
- [ ] 対流境界条件が陰的に実装されている
- [ ] スケーリング調整が適用されている
- [ ] 全てのソルバー（DHCP/Adjoint/Sensitivity）で動作する
- [ ] エッジケース（全断熱境界など）が正しく処理される

### 性能
- [ ] 10ステップテストが900秒以内に完了
- [ ] BiCGSTAB反復回数が500回以下
- [ ] Phase 1-6の505テストが全て合格

### ドキュメント
- [ ] CURRENT_SESSION_STATE.mdを更新
- [ ] コミットメッセージが明確
- [ ] 本実装計画書が最新状態

---

## 📝 実装ログ

### 2025-10-16 22:00 - 計画書作成完了
- Codex調査報告に基づき詳細計画を策定
- tuning7-volume-integral-complete ブランチ作成
- コミット: f024ee2, 05ebe81

### 2025-10-17 - Phase 2完了
- RHSCore.jl修正完了
- `RHS_convection!`からθ引数除去、全6面でRHS項を`h·A·T∞`のみに
- コンパイル成功確認
- コミット: 9815b9d

### 2025-10-17 - Phase 3完了
- CommonSolver.jl修正完了
- 6つの主要関数にHC引数追加:
  - CalcAX!, CalcRK!, jacobi_preconditioner!
  - resSOR, rbsor_core!, rbsor!
- マスク直接利用方式で対流項実装
- コンパイル成功確認
- コミット: 4d65584

### 次セッション: Phase 4開始
- DHCPSolver.jl等でのHC配列生成と伝播
- PBiCGSTAB!, PCG!へのHC引数追加
- Preconditioner!システム全体へのHC伝播
- 詳細: `TODO_NEXT_SESSION.md`参照

---

**End of Document**
