# ベンチマーク実行状況

**日時**: 2025-10-16 19:35
**ブランチ**: tuning7
**最新コミット**: 0dc54f7 (fix: Jacobi前処理のバグ修正)

## 目的
PBICGSTAB vs PCG の性能比較測定

## 完了した作業

### 1. ソルバー選択機能の実装 ✅
- **コミット**: d1562ce
- CGMSolver.jlに以下のパラメータを追加：
  - `dhcp_solver` (デフォルト: `:pbicgstab`)
  - `adjoint_solver` (デフォルト: `:pbicgstab`)
  - `sensitivity_solver` (デフォルト: `:pbicgstab`)
- AdjointSolver.jl、SensitivitySolver.jlは既に対応済み

### 2. Jacobi前処理のバグ修正 ✅
- **コミット**: 0dc54f7
- 修正内容：
  - 635行目: `zeroT` → `zero(T)` （未定義変数）
  - 652-656行目: `θ` → `scratch` （変数名の誤り）

## 次のステップ

### テスト1: デフォルト設定（PBICGSTAB + GS）
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status  # tuning7ブランチであることを確認
julia --project=julia julia/scripts/run_10steps_fullsize_test.jl
```

**期待される結果**:
- CGM実行時間: 約1000-1200秒（80×100×20格子、10ステップ）
- 数値異常なく完了
- 結果ファイル: `shared/results/julia_10steps_fullsize.npz`

### テスト2: PCG設定（PCG + GS）

`run_10steps_fullsize_test.jl`を編集してPCGを使用：

```julia
# 183-200行目のcgm_paramsを以下に変更：
cgm_params = (
    max_iter = 1,
    sigma = 1.8,
    rtol_dhcp = 1.0e-6,
    rtol_adjoint = 1.0e-8,
    maxiter_cg = 20000,
    dire_reset_every = 5,
    eps = 1.0e-12,
    min_iter = 10,
    P = 10,
    eta = 1.0e-4,
    beta_max = 1.0e8,
    verbose = true,
    dhcp_extrapolation = :quadratic,
    adjoint_residual_scale = 0.5,
    # PCG設定を追加
    dhcp_solver = :pcg,
    adjoint_solver = :pcg,
    sensitivity_solver = :pcg
)
```

その後、再実行：
```bash
julia --project=julia julia/scripts/run_10steps_fullsize_test.jl
```

## 比較基準（22fde2d ベースライン）

**10ステップフルサイズ（80×100×20格子）**:
- 総実行時間: 1072.43秒（約17.9分）
- CGM時間: 1069.53秒（99.7%）
- DHCP時間: 321.406秒（30.0%）
- Adjoint時間: 447.191秒（41.7%）
- Sensitivity時間: 300.695秒（28.0%）

## 注意事項

1. **simple_benchmark.jlは使用しない**
   - ランダムデータで数値的に不安定（t=2でNaN/Inf発生）
   - 過去のコミット（287c8c7、230745b）でも同じ問題

2. **実行環境**
   - Julia 1.12.0
   - スレッド数: 1
   - BLASスレッド数: 1

3. **データファイル**
   - `shared/data/T_measure_700um_1ms.npy` (1.1GB) が必要
   - 存在しない場合はエラーで停止

## 結果の記録方法

各テスト完了後、以下を記録：
1. 総実行時間（Total runtime）
2. CGM実行時間（CGM elapsed time）
3. DHCP/Adjoint/Sensitivity各時間
4. 収束状況（反復回数、残差）
5. npzファイルの保存確認

## トラブルシューティング

**数値異常（NaN/Inf）が発生した場合**:
- データファイルの確認
- 格子生成の確認（z_centers, ΔZの値）
- ソルバー設定の確認（rtol, maxiter）

**実行が遅い場合**:
- 正常（1000秒以上かかる）
- バックグラウンド実行推奨: `nohup julia ... > log.txt 2>&1 &`
