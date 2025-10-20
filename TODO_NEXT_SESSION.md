# 次セッションの作業内容

## 現在の状況（2025-10-21）

### ✅ 完了した作業

1. **Python版Adjoint診断ソルバーのバグ修正**
   - `python/validation/solvers_with_diagnostics.py:211`
   - `Y_obs[n + 1]` → `Y_obs[n]` に修正
   - コミット: `70dbf92`

2. **test_single_window.jl修正**
   - 存在しない関数呼び出しを削除
   - 熱物性値をハードコードに変更
   - コミット: `70dbf92`

3. **Python版Phase 1実行成功**
   - 実行時間: 251.38秒（約4.2分）
   - 23ウィンドウ完了、CGM 3回
   - 結果ファイル: `shared/results/python_sliding_window_cgm3.npz` (19MB)
   - ログ: `shared/results/python_phase1_fixed.log` (2801行)

### ⚠️ 重大な問題発見

**Julia版がt=64以降で発散**

#### 残差の比較（ウィンドウ1、CGM反復0）

**Python版（正常）**:
```
t=61: 16回収束, Res_0=5.63e-3 (0.00563)
t=62: 16回収束, Res_0=5.60e-3 (0.00560)
t=63: 16回収束, Res_0=5.57e-3 (0.00557)
t=64: 16回収束, Res_0=5.54e-3 (0.00554)
```

**Julia版（発散）**:
```
t=61: 179回収束, Res_0=51.93
t=62: 1260回収束, Res_0=63.54
t=63: 17256回収束, Res_0=483.41
t=64: 20000回発散, Res_0=659.49
t=65: 20000回発散, Res_0=7998.62
t=66: 20000回発散, Res_0=7.30e6
...
t=71: 20000回発散, Res_0=1.88e21
```

#### 重要な観察

1. **残差が100倍以上違う**: PythonとJuliaで初期残差が全く異なる
2. **t=22から兆候**: Julia版はt=22頃から残差が増加傾向
3. **指数関数的発散**: t=64以降、残差が指数関数的に増加

---

## 🎯 次セッションの最優先タスク

### タスク1: ウィンドウ1の最初のステップを比較

**仮説**: ウィンドウ1では両方とも初期値が同じはずなので、最初のステップ（t=1-10）で既に違いがあるかを確認。

**実行手順**:
```bash
# Python版のウィンドウ1最初の10ステップ
grep "^\[t=" shared/results/python_phase1_fixed.log | head -10

# Julia版のウィンドウ1最初の10ステップ（過去のログから）
# ログはタイムアウトで失われているため、再実行が必要
```

### タスク2: Julia版を小さなウィンドウサイズで再実行

**目的**: どのステップから問題が始まるかを特定

**実行コマンド**:
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# ウィンドウサイズ10、CGM 1回で実行
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --window 10 --overlap 5 --nt 50 \
  > shared/results/julia_small_window_test.log 2>&1
```

### タスク3: 最初のステップの詳細比較

**確認項目**:
1. 初期温度場（T_init）の値
2. 観測データ（Y_obs）の値
3. 初期熱流束（q_init）の値
4. DHCPソルバーの係数行列の値
5. 前処理の有無

### タスク4: インデックスの確認

**疑われる原因**:
- Python: 0-indexed、配列形状 (nt, ni, nj)
- Julia: 1-indexed、配列形状 (ni, nj, nt)

**確認ファイル**:
- `julia/scripts/run_sliding_window_validation.jl:283`
  ```julia
  Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]
  ```
- Python版:
  ```python
  Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]
  ```

---

## 📁 重要なファイル

### 実装ファイル
- **Python診断ソルバー**: `python/validation/solvers_with_diagnostics.py`
- **Python検証スクリプト**: `python/validation/run_sliding_window_validation.py`
- **Julia検証スクリプト**: `julia/scripts/run_sliding_window_validation.jl`

### 結果ファイル
- **Python結果**: `shared/results/python_sliding_window_cgm3.npz` (19MB)
- **Pythonログ**: `shared/results/python_phase1_fixed.log` (2801行)
- **Julia結果**: なし（発散のため未完成）

### ドキュメント
- **検証計画**: `docs/plans/sliding_window_validation_plan.md`
- **前回のサマリー**: `docs/tasks/SESSION_SUMMARY_2025-10-21.md`

---

## 💡 推測される原因

### 仮説1: 配列インデックスのズレ
**可能性**: 高
- PythonとJuliaで配列の次元順序が異なる
- スライディングウィンドウでのインデックス計算にバグ

### 仮説2: ウィンドウ間継承の問題
**可能性**: 中
- 前のウィンドウの最終温度場が正しく継承されていない
- オーバーラップ部分の処理が不正確

### 仮説3: 数値精度の問題
**可能性**: 低
- Julia版のPCGソルバーの精度不足
- しかし、シングルウィンドウテストは成功しているため可能性は低い

---

## 🔧 デバッグ戦略

### ステップ1: 最小限のテストケース
```bash
# ウィンドウ1つ、ステップ10で実行
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --window 10 --nt 10 > test_10steps.log 2>&1
```

### ステップ2: 詳細ログ出力
- Julia版のDHCPソルバーに初期残差のログを追加
- Python版と同じ形式で出力

### ステップ3: 中間データ保存
- 各ウィンドウの終了時に温度場を保存
- 次のウィンドウの開始時と比較

---

## 📊 統計データ

### Python版（成功）
- 総実行時間: 251.38秒
- ウィンドウ数: 23個
- 平均DHCP反復数: 16-17回/ステップ
- 残差範囲: 3.6e-1 → 5.3e-3

### Julia版（失敗）
- 実行時間: タイムアウト
- ウィンドウ数: 1個未完成
- DHCP反復数: t=64で発散（20000回）
- 残差: t=64で659.49 → t=71で1.88e21

---

## 🚀 次セッション開始手順

1. **環境確認**
   ```bash
   cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
   git status
   git log --oneline -5
   ```

2. **ブランチ確認**
   ```bash
   git branch  # → sliding-window-validation
   ```

3. **最新コミット**
   ```bash
   git log -1 --stat
   # → 70dbf92: Python版Adjoint修正、test_single_window.jl追加
   ```

4. **バックグラウンドプロセス確認**
   ```bash
   ps aux | grep "python\|julia" | grep validation
   # 必要に応じてkill
   ```

5. **小規模テスト実行**
   ```bash
   julia --project=julia julia/scripts/run_sliding_window_validation.jl \
     --cgm-iter 1 --window 10 --nt 20 \
     > shared/results/julia_debug_small.log 2>&1 &

   tail -f shared/results/julia_debug_small.log
   ```

---

## 📝 メモ

- Python版は正常動作（251秒で完了）
- Julia版はスライディングウィンドウ特有の問題
- シングルウィンドウテストは両方とも成功
- 問題はウィンドウ継承またはインデックス計算にある可能性が高い

---

**作成日時**: 2025年10月21日
**最終コミット**: 70dbf92
**ブランチ**: sliding-window-validation
