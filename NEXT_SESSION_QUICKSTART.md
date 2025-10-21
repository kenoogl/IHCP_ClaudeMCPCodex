# 次セッション クイックスタート

**📌 5分で状況把握 + 調査開始**

---

## 🔥 現在の状況（重要）

### 問題
Julia版スライディングウィンドウ計算で**初期ステップから発散**
→ Python版と比較して反復数が3-5倍、残差が10倍以上

### 既に完了した作業
- ✅ Y_obsスライス修正（71ステップ）
- ✅ Python版完全動作確認（23ウィンドウ完了）
- ✅ 詳細ログ比較完了

### 未解決の問題
❌ **Julia版はDHCPソルバーの最初のステップから収束が悪い**

---

## 📊 決定的な証拠

### Python版（正常）
```
[t=2/71] Iteration= 32  Res_0= 0.169
[t=10/71] Iteration= 17 Res_0= 0.013  ← 順調に減少
```

### Julia版（異常）
```
[t=2/71] Iteration= 96  Res_0= 0.365  ← 反復3倍、残差2倍
[t=10/71] Iteration= 96 Res_0= 0.162  ← 残差13倍！
[t=62/71] Iteration= 1260 Res_0= 63.5 ← 完全発散
```

---

## 🎯 次の調査ステップ（優先順）

### ステップ1: 時間ステップインデックス確認（最優先）

**仮説**: 時間ループのインデックスが1つずれている

**確認コマンド**:
```bash
# Python版
grep -A 20 "def solve_dhcp" python/validation/run_sliding_window_validation.py

# Julia版
grep -A 20 "function solve_dhcp!" julia/src/solvers/DHCPSolver.jl
```

**チェックポイント**:
- [ ] 時間ループ: `for t in range(1, nt)` vs `for t in 2:nt`
- [ ] T_oldの初期化タイミング
- [ ] 境界条件適用タイミング（t=1 or t=2から）

---

### ステップ2: q_init配列形状確認

**仮説**: q_initの次元順序が逆（Python: (nt,ni,nj) vs Julia: (ni,nj,nt)）

**確認コマンド**:
```bash
# Python版
grep -B 3 -A 3 "q_init.*zeros\|q_init.*shape" python/validation/run_sliding_window_validation.py

# Julia版
grep -B 3 -A 3 "q_init.*zeros" julia/scripts/run_sliding_window_validation.jl
```

**チェックポイント**:
- [ ] Python: `q_init = np.zeros((nt-1, ni, nj))` → (70, 80, 100)
- [ ] Julia: `q_init = zeros(ni, nj, nt-1)` → (80, 100, 70)
- [ ] **次元順序の違いが問題か？**

---

### ステップ3: 初期温度場T_initの値確認

**確認コマンド**:
```bash
# Python版とJulia版でT_initの値が一致するか
python -c "import numpy as np; T=np.load('shared/data/T_measure_700um_1ms.npy'); print(T[0,:5,:5,0])"

julia --project=julia -e 'using NPZ; T=npzread("shared/data/T_measure_700um_1ms.npy"); println(T[1:5,1:5,1])'
```

---

## 📂 重要ファイル

### 必読ドキュメント
1. **SESSION_SUMMARY_SLIDING_WINDOW.md** ← 詳細な調査結果
2. このファイル（クイックスタート）

### ログファイル
- `shared/results/python_phase1_fixed.log` - Python版成功ログ
- `shared/results/julia_phase1_FINAL.log` - Julia版失敗ログ

### 修正済みコード
- `julia/scripts/run_sliding_window_validation.jl:283` - Y_obsスライス

---

## 🚀 作業開始手順

### 1. 状況確認（30秒）
```bash
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex
git status
cat SESSION_SUMMARY_SLIDING_WINDOW.md | head -50
```

### 2. ステップ1実行（2分）
時間ステップインデックスの比較確認

### 3. ステップ2実行（2分）
q_init配列形状の確認

### 4. 問題箇所特定後、修正実施

---

## 💡 重要な気づき

**t=3の異常パターン**:
```
Julia版のt=3だけ残差が0.0688（他のステップは0.09-0.19）
→ インデックスエラーまたは境界条件適用ミスの可能性大
```

この異常を手がかりに調査を進めてください！

---

**最終更新**: 2025-10-21 08:00
**推定作業時間**: 30分～1時間
