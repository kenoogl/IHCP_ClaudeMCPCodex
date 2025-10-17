# 次セッション作業ガイド - Phase 6対流境界条件修正

**作成日**: 2025-10-17
**ブランチ**: `tuning7-volume-integral-complete`
**セッション状況**: Phase 6実装中（NaN/Inf問題の修正途中）

---

## 🎯 現在の状況

### 問題の発見
Phase 6のDHCP対流境界条件テストで**NaN/Inf発生**を検出しました。

**根本原因**（Codex分析済み）:
```
誤った実装: 対流項を隣接セル温度にも掛けていた
結果: 境界セルで h·A·(θ - 2T∞) = 0 という誤った式になり、発散
症状: PBICGSTAB!が20000反復で収束せず、残差がNaN
```

### 修正方針
**対流項は対角項のみに追加、隣接項には熱伝導項のみを使用**

**修正前（誤り）**:
```julia
axm = λf(...) * area / dx + HC[1] * area * (1 - m[i-1,j,k])
ss = axm * θ[i-1,j,k] + ...  # 対流項も掛かる（誤り）
```

**修正後（正しい）**:
```julia
cond_xm = λf(...) * area / dx  # 熱伝導項のみ
conv_xm = HC[1] * area * (1 - m[i-1,j,k])  # 対流項
dd = ... + (cond_xm + conv_xm + ...) * m0  # 対角項に両方
ss = cond_xm * θ[i-1,j,k] + ...  # 隣接項には熱伝導項のみ
```

---

## ✅ 完了した修正

### CommonSolver.jl（julia/src/solvers/）

1. ✅ **CalcAX!** (line 412-446): 修正完了
2. ✅ **CalcRK!** (line 335-372): 修正完了

---

## ⏳ 残りの修正作業

### 修正が必要な関数（CommonSolver.jl）

#### 1. `jacobi_preconditioner!` (line 660-686付近)
**対象ループ**: line 669-693の`@floop`ブロック内

**修正箇所**:
```julia
# 現在（line 661-666付近）
axm = λf(...) * dy0 * dz_k * ddx + HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
axp = λf(...) * dy0 * dz_k * ddx + HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
...

# 修正後
cond_xm = λf(...) * dy0 * dz_k * ddx
cond_xp = λf(...) * dy0 * dz_k * ddx
...
conv_xm = HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
conv_xp = HC[2] * dy0 * dz_k * (oneT - m[i+1,j,k])
...
dd = (oneT-m0) + (cond_xp + cond_xm + ... + conv_xp + conv_xm + ... + a_p_0)*m0
ss = ( cond_xp * scratch[i+1,j,k] + cond_xm * scratch[i-1,j,k] + ... )
```

#### 2. `resSOR` (line 231-260付近)
**対象ループ**: line 231-259の`@floop`ブロック内

**修正箇所**:
```julia
# 現在（line 238-243付近）
axm = λf(...) * dy0 * dz_k * ddx + HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
...

# CalcAX!と同じパターンで修正
```

#### 3. `rbsor_core!` (line 301-329付近)
**対象ループ**: line 301-329の`@floop`と`@simd`ブロック内

**修正箇所**:
```julia
# 現在（line 309-314付近）
axm = λf(...) * dy0 * dz_k * ddx + HC[1] * dy0 * dz_k * (oneT - m[i-1,j,k])
...

# CalcAX!と同じパターンで修正
```

---

## 🔧 次セッションの作業手順

### Step 1: ブランチとコード確認
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -5
```

**期待される状態**:
- ブランチ: `tuning7-volume-integral-complete`
- 未コミットの変更: `julia/src/solvers/CommonSolver.jl` (CalcAX!, CalcRK!修正済み)

### Step 2: 残りの3関数を修正
```bash
# エディタでCommonSolver.jlを開く
# 修正対象: jacobi_preconditioner!, resSOR, rbsor_core!
```

**修正パターン**（CalcAX!と同じ）:
1. 熱伝導項と対流項を分離（`cond_*`と`conv_*`）
2. 対角項`dd`には両方を含める
3. 隣接項`ss`には熱伝導項のみ

### Step 3: テスト実行
```bash
# コンパイルテスト
julia --project=julia -e 'using IHCP_CGM; println("✓ コンパイル成功")'

# 対流境界条件テスト
julia --project=julia julia/scripts/test_dhcp_convection.jl
```

**期待される結果**:
- PBICGSTAB!が収束（反復回数 < 20000）
- NaN/Infが発生しない
- 温度が300K〜550Kの範囲内
- 対流による冷却効果を確認

### Step 4: コミット
```bash
git add julia/src/solvers/CommonSolver.jl
git add julia/src/solvers/DHCPSolver.jl
git add julia/scripts/test_dhcp_convection.jl

git commit -m "$(cat <<'EOF'
fix: Phase 6 - 対流境界条件のステンシル計算を修正

**問題**:
- 対流項を隣接セル温度にも掛けていた（誤り）
- 境界セルで h·A·(θ - 2T∞) = 0 という誤った式になり発散
- PBICGSTAB!が収束せず、残差がNaNになる

**修正内容**（CommonSolver.jl）:
1. CalcAX! - 対流項を隣接項から除外
2. CalcRK! - 同上
3. jacobi_preconditioner! - 同上
4. resSOR - 同上
5. rbsor_core! - 同上

**修正方針**:
- 熱伝導項と対流項を分離（cond_*, conv_*）
- 対角項（dd）: 熱伝導項 + 対流項 + 時間項
- 隣接項（ss）: 熱伝導項のみ（対流項は除外）

**根本原因**（Codex分析）:
対流項がゴーストセル温度（T_∞）に掛かることで、
境界セルの離散式が h·A·(θ - 2T∞) = 0 となり、
θが2T∞に向かって発散していた。

**DHCPSolver.jl**:
- set_dhcp_bc_parameters_convection()関数を追加
- solve_dhcp!にbc_setキーワード引数を追加

**テストスクリプト**:
- julia/scripts/test_dhcp_convection.jl作成
- Z-plus面に対流境界条件（h=1 W/m²·K, T_∞=300K）
- 10×10×5格子、5ステップ

**テスト結果**:
✓ PBICGSTAB!が収束
✓ NaN/Inf検出なし
✓ 温度が物理的に妥当（300K〜550K）
✓ 対流による冷却効果を確認

🤖 Generated with [Claude Code](https://claude.com/claude-code)

Co-Authored-By: Claude <noreply@anthropic.com>
EOF
)"
```

### Step 5: Phase 6完了確認
- TODO_NEXT_SESSION.mdを更新
- Phase 6完了報告を作成

---

## 📝 関連ファイル

### 修正対象
- `julia/src/solvers/CommonSolver.jl`: ステンシル計算修正（main）
- `julia/src/solvers/DHCPSolver.jl`: 境界条件オプション追加

### テストファイル
- `julia/scripts/test_dhcp_convection.jl`: Phase 6テストスクリプト

### 参照ドキュメント
- `docs/VOLUME_INTEGRAL_IMPLEMENTATION_PLAN.md`: 実装計画
- Codex分析結果（セッション内で取得済み）

---

## 🔍 デバッグ情報

### Codex分析サマリー
- **発見**: ゴーストレイヤーのmask=0、λ=0、θ=T∞設定は正しい
- **問題**: CommonSolver内で対流項を隣接温度にも掛けている
- **結果**: `h·A·(θ - 2T∞) = 0`という誤った式
- **修正**: 隣接項から対流項を除外、対角項のみに含める

### テスト条件
- 格子: 10×10×5
- 時間ステップ: 5
- 初期温度: 550K
- Z-plus面: 対流（h=1 W/m²·K, T_∞=300K）
- 他の面: 断熱

---

## ⚠️ 注意事項

1. **CalcAX!とCalcRK!は修正済み** - 再修正不要
2. **残り3関数の修正パターンは同じ** - CalcAX!を参考に
3. **テスト実行前に必ずコンパイル確認**
4. **NaN/Inf検出は数値異常の兆候** - テストで必ず確認

---

## 📊 進捗状況

- [x] Phase 6作業計画立案
- [x] テストケース設計
- [x] テストスクリプト作成
- [x] NaN/Inf問題の根本原因特定（Codex分析）
- [x] CalcAX!修正
- [x] CalcRK!修正
- [ ] jacobi_preconditioner!修正
- [ ] resSOR修正
- [ ] rbsor_core!修正
- [ ] テスト再実行
- [ ] テスト結果ドキュメント化
- [ ] Phase 6完了

---

**End of Document**
