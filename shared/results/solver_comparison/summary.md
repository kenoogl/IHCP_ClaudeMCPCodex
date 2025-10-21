# DHCPソルバー・前処理組み合わせベンチマーク結果

**実行日時**: 2025年 10月21日 火曜日 09時49分28秒 JST

**使用スクリプト**:
- ベンチマーク実行: `julia/scripts/run_solver_comparison.sh`
- 単体テストスクリプト: `julia/scripts/test_dhcp_solver.jl`

**テスト条件**:
- 格子: 80×100×20 (N=160,000)
- 時間ステップ: 10
- 収束条件: rtol=1.0e-6, maxiter=20000
- 熱流束: ゼロ（簡易テスト）

**実行方法**:
```bash
# 全6パターン自動実行
bash julia/scripts/run_solver_comparison.sh

# 個別実行例
julia --project=julia julia/scripts/test_dhcp_solver.jl --solver pcg --precond none --nt 10
```

## 結果一覧

| ソルバー | 前処理 | 実行時間 | 最終反復数 |
|----------|--------|----------|------------|
| pbicgstab | none | 10.43 s | 61 |
| pbicgstab | diagonal | 14.21 s | 47 |
| pbicgstab | gs | 11.08 s | 11 |
| pcg | none | 8.96 s | 95 |
| pcg | diagonal | 12.17 s | 81 |
| pcg | gs | 9.52 s | 20 |

## 詳細ログ

各テストの詳細ログは以下のファイルに保存されています:
- `pbicgstab_none_20251021_094731.log`
- `pbicgstab_diagonal_20251021_094731.log`
- `pbicgstab_gs_20251021_094731.log`
- `pcg_none_20251021_094731.log`
- `pcg_diagonal_20251021_094731.log`
- `pcg_gs_20251021_094731.log`

## 分析

### 実行時間比較

**ソルバー別**:
- PCG: 8.96s (none) < 9.52s (gs) < 12.17s (diagonal)
- PBICGSTAB: 10.43s (none) < 11.08s (gs) < 14.21s (diagonal)

**前処理別**:
- none: PCG 8.96s < PBICGSTAB 10.43s
- gs: PCG 9.52s < PBICGSTAB 11.08s
- diagonal: PCG 12.17s < PBICGSTAB 14.21s

**最速**: PCG + none (8.96s)

### 反復回数比較

**前処理の効果**:
- PBICGSTAB: 61 (none) → 47 (diagonal) → 11 (gs) - GS前処理で大幅削減
- PCG: 95 (none) → 81 (diagonal) → 20 (gs) - GS前処理で約1/5に削減

**ソルバー別**:
- none: PBICGSTAB 61回 < PCG 95回
- diagonal: PBICGSTAB 47回 < PCG 81回
- gs: PBICGSTAB 11回 < PCG 20回

**反復回数削減率** (none基準):
- PBICGSTAB+gs: 61→11 (82%削減)
- PCG+gs: 95→20 (79%削減)

### 推奨設定

**1. 最速を優先する場合**:
```julia
solver = :pcg
smoother = :none
```
- 実行時間: 8.96s
- 反復回数: 95回
- **理由**: 前処理オーバーヘッドが反復削減効果を上回る

**2. 反復回数を抑えたい場合**:
```julia
solver = :pbicgstab
smoother = :gs
```
- 実行時間: 11.08s (最速比+23.7%)
- 反復回数: 11回（全パターン中最少）
- **理由**: GS前処理で反復回数を大幅削減

**3. バランス型**:
```julia
solver = :pcg
smoother = :gs
```
- 実行時間: 9.52s (最速比+6.3%)
- 反復回数: 20回
- **理由**: 実行時間と反復回数のバランスが良い

### 重要な知見

1. **対角前処理の非効率性**: diagonal前処理は反復回数を削減するが、オーバーヘッドが大きく実行時間は最も遅い（PCG 12.17s, PBICGSTAB 14.21s）

2. **前処理なし (none) の優位性**: 小規模問題では前処理オーバーヘッドが支配的で、none が最速

3. **GS前処理の効果**: 反復回数は大幅削減（82%減）するが、実行時間では none に対し+6~23%程度

4. **ソルバー選択**: PCGは全ての前処理パターンでPBICGSTABより高速

### 今後の検討事項

- より大規模な問題（格子数増加）での前処理効果の再評価
- 非ゼロ熱流束条件での性能比較
- より多くのタイムステップでの長時間計算での性能比較
