# 次セッションで実施するタスク

**作成日**: 2025年10月24日
**ブランチ**: parallelization
**現在の状態**: Python版CGM反復回数検証完了、Julia版比較検証準備中

---

## 📊 このセッションで完了したこと

### ✅ Python版CGM反復回数検証（完了）

**実行済みテスト**:
1. **CGM反復10回**: 42.49秒完了
2. **CGM反復20回**: 43.74秒完了

**主要な発見**:
- **Discrepancy条件による早期停止**: CGM反復20回でも実質10回で自動停止
- **計算時間がほぼ同じ理由**: 両方とも実質10回で収束判定により停止
- **最適CGM反復回数**: 2-3回が実用的（10回以上は発散傾向）

**詳細データ**:
| 項目 | CGM 2回 | CGM 10回 | CGM 20回 |
|-----|---------|----------|----------|
| 実行時間 | 8.39秒 | 42.49秒 | 43.74秒 |
| 熱流束範囲 | -9.46e+07 ~ 3.91e+08 | -1.82e+08 ~ 4.22e+08 | -2.47e+08 ~ 4.40e+08 W/m² |
| 実際の反復回数 | 2回 | 10回 | 10回（自動停止） |
| 数値安定性 | 良好 ✅ | 発散傾向 ⚠️ | さらに発散 ⚠️ |

**結果ファイル**:
- `shared/results/python_cgm_iter10.log` - CGM反復10回の詳細ログ
- `shared/results/python_cgm_iter20.log` - CGM反復20回の詳細ログ
- `shared/results/python_cgm_iter10.npz` - CGM反復10回の熱流束データ
- `shared/results/python_cgm_iter20.npz` - CGM反復20回の熱流束データ

---

## 🎯 次セッションで実施すべきタスク

### 優先度A: 重要・緊急

#### 1. Julia版との比較検証（CGM反復2-3回）

**実行コマンド**:
```bash
# Julia版CGM反復2回
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 2 \
  --solver pbicgstab --precond gs --basesize 500 \
  2>&1 | tee shared/results/julia_cgm_iter2_basesize500.log

# Julia版CGM反復3回
JULIA_NUM_THREADS=4 julia --project=julia \
  julia/scripts/run_sliding_window.jl \
  --nt 10 --window 5 --overlap 2 --cgm-iter 3 \
  --solver pbicgstab --precond gs --basesize 500 \
  2>&1 | tee shared/results/julia_cgm_iter3_basesize500.log
```

**比較ポイント**:
- 実行時間の比較
- 熱流束範囲の比較
- CGM収束性の比較
- 数値安定性の比較

**期待される結果**:
- Julia版の方が高速（Phase 5.2 Step 3では34.34秒）
- 熱流束範囲は同オーダー（10^7-10^8 W/m²）
- 両版ともCGM 2-3回で安定した結果

#### 2. Python-Julia比較レポート作成

**レポートファイル**: `docs/reports/python_julia_cgm_iteration_comparison.md`

**含めるべき内容**:
1. CGM反復2回での比較
2. CGM反復3回での比較
3. CGM反復回数依存性の分析
4. 早期停止機能の有無の比較
5. 推奨CGM反復回数の結論

### 優先度B: 重要・非緊急

#### 3. Phase 5.2完了報告書の更新

**ファイル**: `docs/reports/PHASE5_2_COMPLETION_SUMMARY.md`

**追加すべき内容**:
- Python版CGM反復回数検証結果（セクション追加）
- Discrepancy条件による早期停止の解説
- 最適CGM反復回数の推奨（2-3回）
- Python-Julia比較結果の統合

#### 4. Step 3 Phase 2（大ウィンドウ）測定

**実行予定**:
- basesize=500: window=71, overlap=17, CGM反復3回
- basesize=1000: window=71, overlap=17, CGM反復3回
- basesize=2000: window=71, overlap=17, CGM反復3回

**実行時間見込み**: 各30分×3 = 90分

### 優先度C: 改善提案

#### 5. basesize自動調整機能の実装
#### 6. 性能プロファイリング詳細分析

---

## 📁 重要なファイルの場所

### 実行結果ログ
```
shared/results/python_cgm_iter10.log        # Python CGM 10回
shared/results/python_cgm_iter20.log        # Python CGM 20回
shared/results/python_cgm_iter10.npz        # 熱流束データ
shared/results/python_cgm_iter20.npz        # 熱流束データ
```

### レポート
```
docs/reports/PHASE5_2_COMPLETION_SUMMARY.md # Phase 5.2完了報告書
docs/plans/phase5_2_real_world_validation_plan.md # Phase 5.2計画書
```

### スクリプト
```
julia/scripts/run_sliding_window.jl        # Julia版スライディングウィンドウ
python/scripts/run_sliding_window.py       # Python版スライディングウィンドウ
```

---

## 🔧 注意事項

### バックグラウンドプロセス

**以下のプロセスが実行中の可能性があります（確認してkillが必要）**:
- Julia版basesize測定テスト（複数）
- Python版テスト（複数）

**確認コマンド**:
```bash
ps aux | grep julia
ps aux | grep python
```

**killが必要な場合**:
```bash
pkill -f "julia.*run_sliding_window"
pkill -f "python.*run_sliding_window"
```

### Phase 5.2の進捗状況

**達成度**: 約85%完了
- ✅ Step 1: DHCP単体basesize検証完了
- ✅ Step 2: 10ステップCGM basesize検証完了
- ✅ Step 3 Phase 1: 小ウィンドウbasesize検証完了
- ✅ Python版CGM反復回数検証完了（追加作業）
- ⏸️ Step 3 Phase 2: 大ウィンドウbasesize検証未完
- ⏸️ Python-Julia完全比較未完

---

## 🚀 次セッション開始時の手順

### 1. 実態確認（最優先）

```bash
# 現在のディレクトリ確認
pwd

# git状態確認
git status
git log --oneline -5

# ブランチ確認
git branch

# 最新の結果ファイル確認
ls -lht shared/results/*.log | head -10
```

### 2. バックグラウンドプロセス確認

```bash
# 実行中プロセス確認
ps aux | grep -E "(julia|python).*sliding_window"

# 必要に応じてkill
pkill -f "julia.*run_sliding_window"
pkill -f "python.*run_sliding_window"
```

### 3. Julia版CGM反復2-3回の実行

Phase 5.2完了報告書の内容を参照し、basesize=500で実行。

### 4. 比較レポート作成

Python版とJulia版のCGM反復回数依存性を比較分析。

---

## 📝 重要な技術的知見

### Python版CGMの特徴

**Discrepancy条件による早期停止**:
```python
# 停止条件
J < J_threshold かつ max|ΔT| ≤ σ (標準偏差閾値 = 1.8K)
```

**利点**:
- 過度な反復による発散を自動防止
- 無駄な計算を回避
- 数値安定性の向上

**CGM反復20回でも実質10回で停止した例**:
- Window 1: J=2.29e+03 < 1.30e+05 かつ max|ΔT|=1.59K ≤ 1.8K → 停止
- Window 2: J=3.26e+03 < 1.30e+05 かつ max|ΔT|=1.51K ≤ 1.8K → 停止
- Window 3-5: 同様に10回で自動停止

### 最適CGM反復回数

**推奨**: **2-3回**

**理由**:
1. 数値安定性が高い
2. 計算時間が短い（8-15秒程度）
3. 熱流束が物理的に妥当な範囲
4. CGM 10回以上では発散傾向

---

**次セッション開始時は、このファイルを必ず最初に読んでください。**
