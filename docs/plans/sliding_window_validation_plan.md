# Python-Juliaスライディングウィンドウ計算比較検証計画

**作成日**: 2025年10月21日
**ステータス**: 計画中
**目的**: オリジナルPythonコードとJuliaコードでスライディングウィンドウ計算を実行し、性能と精度を完全に比較検証する

---

## 1. 背景と目的

### 1.1 背景

これまでに以下の検証が完了している:
- ✅ Phase 1-6: 505テスト全合格（Julia移植完了）
- ✅ シングルウィンドウ10ステップテスト: Python-Julia同等性確認済み
  - 温度場誤差: 平均6.06K（相対1.06%）
  - CGM反復1回での検証完了

### 1.2 目的

**スライディングウィンドウ適用時のPython-Julia完全互換性検証**

1. **精度検証**: 各ウィンドウでの計算結果の一致確認
2. **性能検証**: 各ソルバー（DHCP/Adjoint/Sensitivity）の反復回数・実行時間比較
3. **安定性検証**: 長時間計算での収束性・数値安定性確認
4. **詳細プロファイル取得**: 性能改善のボトルネック特定

---

## 2. 検証計画の全体構成

### 2段階アプローチ

#### **Phase 1: オプション2（短時間動作確認）**
- **CGM反復数**: 3回固定
- **目的**: 両コードの基本的な一致確認、スクリプト動作検証
- **所要時間**: Python約15分、Julia約5分（見積もり）

#### **Phase 2: オプション1（本格検証）**
- **CGM反復数**: 20000回（早期停止あり）
- **目的**: オリジナルと同じ設定での完全比較
- **所要時間**: Python約5分×反復回数、Julia約1-2分×反復回数（見積もり）

---

## 3. 実行パラメータ

### 3.1 共通パラメータ

```
格子サイズ:
  - ni × nj × nk = 80 × 100 × 20 (フルサイズ)
  - dx = 0.12e-3 m
  - dy = 0.12e-3 * sin(80°) / sin(45°) m
  - Lz = 0.5e-3 m (stretch_factor = 3.0)

時間設定:
  - dt = 1.0e-3 s (1ms)
  - nt = 300 ステップ (約5分のPython numba処理時間を想定)

スライディングウィンドウ:
  - window_size = 71 ステップ
  - overlap = 17 ステップ
  - 予想ウィンドウ数: 約6ウィンドウ

熱物性値:
  - rho = 7823.493962874829 kg/m³
  - cp_coeffs = [2.00929659e-10, -3.42605571e-07, 1.34927936e-01, 4.69852860e+02]
  - k_coeffs = [4.79912245e-12, -8.18299348e-09, 1.61765445e-02, 8.11751748e+00]
```

### 3.2 CGM設定

#### Phase 1 (オプション2):
```
CGM反復数: 3回固定
早期停止: なし
ソルバー設定:
  - DHCP: rtol=1e-6, maxiter=20000
  - Adjoint: rtol=1e-8, maxiter=20000
  - Sensitivity: rtol=1e-8, maxiter=20000
```

#### Phase 2 (オプション1):
```
CGM反復数: 最大20000回
早期停止: あり
  - min_iter = 10
  - P = 10 (プラトー検出ウィンドウ)
  - eta = 1e-4 (平均相対減少閾値)
ソルバー設定: Phase 1と同じ
```

---

## 4. 実装タスク

### 4.1 出力フォーマット整備

#### 4.1.1 Python側（`global_CGM_time()`関数拡張）

**現在の出力**（1535-1541行）:
```python
print(f"@ ___ Iter {it:3d} ___ @ wall_s = {wall_s:.3f}s")
print(f"J = {J:.5e}, beta = {beta:.4e}, rel_drop = {rel_drop:.3e}")
print(f"|T - Y|: max={delta_T.max():.3e}, min={delta_T.min():.3e}, mean={delta_T.mean():.3e}")
print(f"grad: min={grad.min():.4e}, max={grad.max():.4e}, mean={grad.mean():.4e}")
print(f"dT:   min={dT[1:].min():.4e}, max={dT[1:].max():.4e}, mean={dT[1:].mean():.4e}")
print(f"q:    min={q.min():.4e}, max={q.max():.4e}, mean={q.mean():.4e}")
```

**追加する情報**:
- DHCPソルバー: 反復回数、実行時間、初期残差、最終残差
- Adjointソルバー: 反復回数、実行時間、初期残差、最終残差
- Sensitivityソルバー: 反復回数、実行時間、初期残差、最終残差
- 温度場範囲: T_cal全体のmin/max/mean

#### 4.1.2 Julia側（各ソルバー拡張）

**DHCPSolver.jl**:
- `solve_dhcp_pbicgstab!()`: 反復回数、残差履歴を返す
- 実行時間計測追加

**AdjointSolver.jl**:
- `solve_adjoint_pbicgstab!()`: 反復回数、残差履歴を返す
- 実行時間計測追加

**SensitivitySolver.jl**:
- `solve_sensitivity_pbicgstab!()`: 反復回数、残差履歴を返す
- 実行時間計測追加

**CGMSolver.jl**:
- Python出力フォーマットと同等の詳細ログ追加
- `verbose=true`時に各ソルバー情報を出力

### 4.2 テストスクリプト作成

#### 4.2.1 Python版: `python/validation/run_sliding_window_validation.py`

**機能**:
- フルサイズ格子（80×100×20）、300ステップ
- スライディングウィンドウ適用
- 各ウィンドウごとの詳細ログ出力
- CGM反復数を引数で指定可能（デフォルト3回）
- 結果をnpz形式で保存

**出力ファイル**:
- `shared/results/python_sliding_window_cgm3.npz` (Phase 1)
- `shared/results/python_sliding_window_cgm20000.npz` (Phase 2)

**保存内容**:
```python
{
  'q_global': 熱流束全体 (nt-1, ni, nj),
  'T_final': 最終温度場 (ni, nj, nk),
  'elapsed_total': 総実行時間,
  'windows_info': [
    {
      'window_id': ウィンドウID,
      'start_idx': 開始インデックス,
      'end_idx': 終了インデックス,
      'cgm_iterations': CGM反復数,
      'J_final': 最終目的関数値,
      'dhcp_iters': DHCP反復回数リスト,
      'adjoint_iters': Adjoint反復回数リスト,
      'sensitivity_iters': Sensitivity反復回数リスト,
      'dhcp_time': DHCP実行時間,
      'adjoint_time': Adjoint実行時間,
      'sensitivity_time': Sensitivity実行時間,
      'elapsed_window': ウィンドウ実行時間
    },
    ...
  ],
  # パラメータ
  'nt': 300,
  'window_size': 71,
  'overlap': 17,
  'cgm_max_iter': 3 or 20000,
  ...
}
```

#### 4.2.2 Julia版: `julia/scripts/run_sliding_window_validation.jl`

**機能**:
- Python版と同じパラメータ設定
- 同じ出力フォーマット
- CGM反復数をコマンドライン引数で指定

**使用方法**:
```bash
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 20000
```

**出力ファイル**:
- `shared/results/julia_sliding_window_cgm3.npz` (Phase 1)
- `shared/results/julia_sliding_window_cgm20000.npz` (Phase 2)

#### 4.2.3 比較スクリプト: `python/validation/compare_sliding_window_results.py`

**機能**:
- Python/Juliaの両結果を読み込み
- ウィンドウごとの詳細比較
- レポート生成（Markdown形式）

**比較項目**:

1. **性能比較**:
   - 総実行時間
   - ウィンドウあたりの平均時間
   - 各ソルバーの累積実行時間
   - ソルバー反復回数の統計

2. **精度比較**:
   - 熱流束誤差: 絶対誤差、相対誤差、RMS誤差
   - 温度場誤差（最終状態）
   - ウィンドウごとの目的関数値比較

3. **収束性比較**:
   - 各ウィンドウのCGM反復回数
   - 各ソルバーの平均反復回数
   - 残差推移（Phase 2のみ）

**出力レポート**:
- `shared/results/sliding_window_comparison_cgm3.md` (Phase 1)
- `shared/results/sliding_window_comparison_cgm20000.md` (Phase 2)

---

## 5. 実行手順

### Phase 1: オプション2（CGM 3回固定）

#### Step 1: スクリプト作成
```bash
# 出力フォーマット整備、テストスクリプト作成
# タスク4.1, 4.2を実施
```

#### Step 2: Python実行
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
python python/validation/run_sliding_window_validation.py --cgm-iter 3
# 推定実行時間: 約15分
```

#### Step 3: Julia実行
```bash
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
# 推定実行時間: 約5分
```

#### Step 4: 結果比較
```bash
python python/validation/compare_sliding_window_results.py --phase 1
```

#### Step 5: 結果確認・問題修正
- 誤差が許容範囲内（相対誤差<5%）か確認
- 問題があればデバッグ・修正

### Phase 2: オプション1（CGM 20000回、早期停止あり）

#### Step 1: Phase 1の問題修正完了確認

#### Step 2: Python実行
```bash
python python/validation/run_sliding_window_validation.py --cgm-iter 20000
# 推定実行時間: 状況により変動（数時間の可能性）
```

#### Step 3: Julia実行
```bash
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 20000
# 推定実行時間: Pythonより高速（見積もり30-50%削減）
```

#### Step 4: 結果比較
```bash
python python/validation/compare_sliding_window_results.py --phase 2
```

#### Step 5: 詳細分析
- 各ソルバーの性能プロファイル
- ボトルネック特定
- 性能改善提案書へのフィードバック

---

## 6. 期待される成果

### Phase 1完了時:
- ✅ スライディングウィンドウ機能の基本動作確認
- ✅ Python-Julia間の計算一致性確認（短時間実行）
- ✅ 出力フォーマット整備完了

### Phase 2完了時:
- ✅ オリジナルPythonと同等設定での完全検証
- ✅ 長時間計算での安定性・精度確認
- ✅ 各ソルバーの詳細性能プロファイル取得
- ✅ 性能改善の優先順位特定

---

## 7. 成功基準

### Phase 1（必須）:
- [ ] 両スクリプトがエラーなく完走
- [ ] 熱流束相対誤差 < 5%
- [ ] ウィンドウ数の一致
- [ ] 各ウィンドウの目的関数値が近い（相対誤差<10%）

### Phase 2（目標）:
- [ ] 温度場相対誤差 < 2%
- [ ] 熱流束相対誤差 < 1%
- [ ] CGM収束挙動の一致（反復回数の差<20%）
- [ ] Julia実行時間がPythonの50%以下

---

## 8. リスクと対策

### リスク1: 実行時間が長すぎる（Phase 2）
**対策**:
- まずPhase 1で動作確認
- Phase 2は時間ステップ数を段階的に増やす（100→200→300）

### リスク2: メモリ不足
**対策**:
- ガベージコレクション強化
- ウィンドウサイズを小さくする（71→50）

### リスク3: 精度が一致しない
**対策**:
- シングルウィンドウ単位で比較
- 差分が発生するウィンドウを特定
- 該当ウィンドウのみ詳細デバッグ

---

## 9. 次のステップ（Phase 2完了後）

1. **性能改善実施**: 特定されたボトルネックの最適化
2. **フルスケールテスト**: 全測定データ（数千ステップ）での検証
3. **論文・報告書作成**: 検証結果のまとめ

---

## 付録A: ファイル構成

```
TrialClaudeMCPCodex/
├── python/
│   └── validation/
│       ├── run_sliding_window_validation.py      # [新規作成]
│       └── compare_sliding_window_results.py     # [新規作成]
├── julia/
│   ├── scripts/
│   │   └── run_sliding_window_validation.jl      # [新規作成]
│   └── src/solvers/
│       ├── DHCPSolver.jl                          # [出力追加]
│       ├── AdjointSolver.jl                       # [出力追加]
│       ├── SensitivitySolver.jl                   # [出力追加]
│       └── CGMSolver.jl                           # [出力追加]
├── shared/results/
│   ├── python_sliding_window_cgm3.npz            # Phase 1 Python結果
│   ├── julia_sliding_window_cgm3.npz             # Phase 1 Julia結果
│   ├── sliding_window_comparison_cgm3.md         # Phase 1 比較レポート
│   ├── python_sliding_window_cgm20000.npz        # Phase 2 Python結果
│   ├── julia_sliding_window_cgm20000.npz         # Phase 2 Julia結果
│   └── sliding_window_comparison_cgm20000.md     # Phase 2 比較レポート
└── docs/plans/
    └── sliding_window_validation_plan.md         # 本ドキュメント
```

---

**計画書バージョン**: 1.0
**最終更新**: 2025年10月21日
