# 次セッション クイックスタート

## 📌 現在地

**ブランチ**: `sliding-window-validation`
**最新コミット**: `3836f0c` - セッション2サマリー追加
**日時**: 2025年10月21日 03:30 JST

---

## ⚡ 即座に実行すべきコマンド

### 1. 環境確認（5秒）
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
git log --oneline -3
```

### 2. 前セッションの成果確認（5秒）
```bash
# Python版Phase 1成功確認
ls -lh shared/results/python_sliding_window_cgm3.npz
# → 19MB、23ウィンドウ完了

# Julia版小規模テスト成功確認
ls -lh shared/results/julia_debug_small.log
# → window=10で成功
```

### 3. 中規模テスト実行（開始）（5分）
```bash
# Julia版window=40でテスト（二分探索の開始）
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --window 40 --nt 70 \
  > shared/results/julia_window40.log 2>&1 &

# ログ監視
tail -f shared/results/julia_window40.log
# Ctrl+Cで停止

# 5分後に結果確認
tail -50 shared/results/julia_window40.log
```

---

## 🎯 最優先タスク（30分以内）

### タスク1: 発散境界の特定（20分）

**目的**: どのウィンドウサイズから発散するかを二分探索で特定

**現状**:
- window=10: ✅ 成功（確認済み）
- window=71: ❌ t=61から発散（確認済み）

**二分探索手順**:
1. window=40を試す（中間地点）
2. 成功なら → window=55を試す
3. 失敗なら → window=25を試す
4. 境界が見つかるまで繰り返す

**期待結果**: 例「window=35まで成功、window=40から発散」

### タスク2: ソルバー設定の調整（10分）

発散境界が見つかったら、以下の設定を試す：

**1. 許容誤差を厳しくする**
```julia
# julia/scripts/run_sliding_window_validation.jl
# 修正候補:
rtol_dhcp: 1.0e-6 → 1.0e-8
rtol_adjoint: 1.0e-8 → 1.0e-10
```

**2. 前処理を変更**
```julia
dhcp_smoother: :jacobi → :ssor
adjoint_smoother: :jacobi → :ssor
```

**3. テスト実行**
```bash
# 修正後、発散したウィンドウサイズで再テスト
julia --project=julia julia/scripts/run_sliding_window_validation.jl \
  --cgm-iter 1 --window XX --nt YY \
  > shared/results/julia_windowXX_fixed.log 2>&1
```

---

## 🔍 問題の詳細

### セッション2で判明した事実

1. **Python版は完全成功** ✅
   - 23ウィンドウ完了、約4分
   - 残差は安定（5e-3レベル）
   - 結果: `python_sliding_window_cgm3.npz` (19MB)

2. **Julia版の問題はウィンドウサイズ依存** ❌
   - window=10: ✅ 成功
   - window=71: ❌ t=61から発散
   - 残差がPython版の**9000倍**（51.93 vs 5.63e-3）

3. **配列インデックス・データ形状は正しい** ✅
   - Python版: `(300, 80, 100)` 正しくスライス
   - Julia版: `(80, 100, 300)` 正しく変換済み
   - ウィンドウスライスも正しい（両方とも72個の時刻）

### 推測される原因（優先度順）

1. **ウィンドウサイズ依存の数値安定性問題（80%）**
   - 大きなウィンドウで数値誤差が蓄積
   - 対策: 許容誤差を厳しく、前処理を改善

2. **ソルバー実装の微妙な違い（15%）**
   - PCG前処理の違い（jacobi vs ssor）
   - 対策: Python版と同じ設定に変更

3. **初期化の問題（5%）**
   - t=61までは正常に動作
   - 対策: 詳細ログで確認

---

## 📂 重要なファイル

### 読むべきファイル
1. `TODO_NEXT_SESSION.md` - 詳細な状況説明
2. `julia/scripts/run_sliding_window_validation.jl:273-320` - スライディングウィンドウループ
3. `python/validation/run_sliding_window_validation.py:102-211` - Python版の同じ部分

### 結果ファイル
- `shared/results/python_sliding_window_cgm3.npz` - Python成功結果（19MB）
- `shared/results/python_phase1_fixed.log` - Pythonログ（2801行）

---

## 🚨 デバッグのヒント

### 最初に確認すべきこと

1. **Y_obsのスライス**
   ```julia
   # Julia版: これが正しいか？
   Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]

   # start_idx=0, end_idx=70の場合
   # → Y_obs[:, :, 1:71] (71ステップ)
   ```

2. **ウィンドウの範囲**
   ```julia
   # ウィンドウ1: start_idx=0, max_L=71, end_idx=71
   # → Y_obsは1:72が必要？（nt+1個）
   ```

3. **温度場の時間次元**
   - T_cal: `(ni, nj, nk, nt)` - 71ステップなら72個の時刻（0を含む）
   - Y_obs: `(ni, nj, nt)` - 71ステップなら72個の時刻

### 簡単なテスト

```julia
# julia REPLで確認
using NPZ
Y_obs = npzread("shared/data/T_measure_700um_1ms.npy")
println(size(Y_obs))  # (18143, 80, 100)

# ウィンドウ1のスライス
start_idx = 0
end_idx = 70
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]
println(size(Y_obs_win))  # (80, 100, 71) のはず
```

---

## 📝 次回への引き継ぎ

**もし今回のセッションで解決できなかった場合**:

1. 小規模テストのログを保存
2. Python版とJulia版の最初の10ステップを比較した結果をメモ
3. 疑わしいコード箇所を特定
4. 次回は特定箇所のデバッグから開始

**もし解決した場合**:

1. 修正内容をコミット
2. Julia版Phase 1を完全実行
3. Python版との結果比較
4. Phase 2（CGM 20000回）へ進む

---

**作成日時**: 2025年10月21日
**所要時間**: このセッションで約2時間作業
**次回予想時間**: 30分-1時間（問題特定と修正）
