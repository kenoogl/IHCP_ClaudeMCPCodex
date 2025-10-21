# スライディングウィンドウ実装差異分析レポート

**日付**: 2025年10月21日
**作成者**: Claude Code
**目的**: Python版とJulia版のスライディングウィンドウ実装差異を特定し、CGM反復数の違いの原因を解明

---

## 📋 エグゼクティブサマリー

### 重大な発見

**Python版とJulia版のスライディングウィンドウアルゴリズムに本質的な差異が存在**

| 項目 | Python版 | Julia版 |
|------|---------|---------|
| **ウィンドウ数** | 9個 | 1個 |
| **ウィンドウサイズ** | 縮小（9→8→7→...→1） | 固定（9） |
| **処理方法** | 1ステップずつスライド | 1回で終了 |

### 影響

この差異により、Python版では**9回のCGM最適化**が実行されるのに対し、Julia版では**1回のみ**となります。

---

## 🔍 詳細分析

### テストケース

- **総時間ステップ**: 10 (nt=10)
- **熱流束ステップ**: 9 (nt-1)
- **ウィンドウサイズ**: 71
- **オーバーラップ**: 17

### Python版の動作

**ファイル**: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py`

#### アルゴリズム (1603-1653行)

```python
start_idx = 0
while start_idx < nt - 1:
    max_L = min(window_size, (nt - 1) - start_idx)  # 1610行
    end_idx = start_idx + max_L
    Y_obs_win = Y_obs[start_idx: end_idx + 1, :, :]

    # CGM計算実行
    q_win, T_win_last, J_hist = global_CGM_time(...)

    step = max(1, max_L - overlap)  # 1652行
    start_idx += step  # 1653行
```

#### 実行結果

| ウィンドウID | 開始 | 終了 | 長さ | ステップ | 時間範囲 (ms) |
|------------|------|------|------|---------|-------------|
| 0 | 0 | 9 | 9 | 1 | 0.0 - 9.0 |
| 1 | 1 | 9 | 8 | 1 | 1.0 - 9.0 |
| 2 | 2 | 9 | 7 | 1 | 2.0 - 9.0 |
| 3 | 3 | 9 | 6 | 1 | 3.0 - 9.0 |
| 4 | 4 | 9 | 5 | 1 | 4.0 - 9.0 |
| 5 | 5 | 9 | 4 | 1 | 5.0 - 9.0 |
| 6 | 6 | 9 | 3 | 1 | 6.0 - 9.0 |
| 7 | 7 | 9 | 2 | 1 | 7.0 - 9.0 |
| 8 | 8 | 9 | 1 | 1 | 8.0 - 9.0 |

**総ウィンドウ数**: 9個
**総実行時間**: 28.60秒

#### キーポイント

```python
step = max(1, max_L - overlap)
# 第1ウィンドウ: max_L=9, overlap=17
# step = max(1, 9 - 17) = max(1, -8) = 1
```

**`max_L < overlap`の場合、強制的に1ステップずつスライド**

---

### Julia版の動作

**ファイル**: `julia/scripts/run_sliding_window.jl`

#### アルゴリズム (286-398行)

```julia
while prev_flux_end < total_flux_steps
    start_idx = if window_id == 1
        0
    else
        last_start = window_starts[end]
        last_length = window_lengths[end]
        last_start + max(1, last_length - overlap)  # 293行
    end

    max_L = min(window_size, total_flux_steps - start_idx)  # 300行
    end_idx = start_idx + max_L

    # CGM計算実行
    q_win, T_win, J_hist = solve_cgm!(...)

    # 次のウィンドウへ
    flux_end = start_idx + max_L
    prev_flux_end = max(prev_flux_end, flux_end)  # 372行
    window_id += 1
end
```

#### 実行結果

| ウィンドウID | 開始 | 終了 | 長さ | カバー範囲 |
|------------|------|------|------|----------|
| 1 | 0 | 9 | 9 | 全範囲 |

**総ウィンドウ数**: 1個
**状態**: 最初のウィンドウで全熱流束ステップ（0-8）をカバーしたため、ループ終了

#### ループ終了条件

```julia
while prev_flux_end < total_flux_steps  # 286行
```

- 第1ウィンドウ後: `prev_flux_end = 9`
- `total_flux_steps = 9`
- 条件: `9 < 9` → **False** → ループ終了

---

## 🔬 根本原因分析

### Python版の終了条件

```python
while start_idx < nt - 1:  # 1603行
```

- **基準**: ウィンドウの開始位置
- **動作**: `start_idx`が熱流束範囲内である限り継続
- **結果**: 縮小するウィンドウでも処理継続

### Julia版の終了条件

```julia
while prev_flux_end < total_flux_steps  # 286行
```

- **基準**: カバーした熱流束の終端位置
- **動作**: 全範囲をカバーしたら終了
- **結果**: 1ウィンドウで全範囲カバー → 即座に終了

---

## ⚠️ 実装差異の本質

### 設計思想の違い

| 観点 | Python版 | Julia版 |
|------|---------|---------|
| **終了判定** | 開始位置ベース | カバー範囲ベース |
| **小ウィンドウ** | 処理する | 処理しない |
| **効率性** | 低（重複計算多） | 高（最小計算） |
| **精度** | 段階的改善 | 一括最適化 |

### Python版の「バグ」または「仕様」？

コードコメント (1601行):
```python
safety_limit = nt * 5  # 经验值
```

中国語コメント「経験値」→ 無限ループ防止のための安全措置

**推測**:
- Python版は意図的に1ステップずつスライドする設計
- 小さいウィンドウで段階的に熱流束を改善
- 計算コストは高いが、各時刻で局所最適化

---

## 📊 CGM反復数への影響

### Python版（10ステップ、CGM最大3回）

| ウィンドウ | サイズ | CGM反復 | 目的関数J |
|----------|--------|---------|----------|
| 0 | 9 | 3 | 6926.093 |
| 1 | 8 | 3 | 9246.580 |
| 2 | 7 | 3 | 9539.347 |
| 3 | 6 | 3 | 9295.892 |
| 4 | 5 | 3 | 8582.703 |
| 5 | 4 | 3 | 7436.788 |
| 6 | 3 | 3 | 5935.158 |
| 7 | 2 | 3 | 4140.566 |
| 8 | 1 | 3 | 2126.959 |

**総CGM実行回数**: 9ウィンドウ × 3反復 = **27回**

### Julia版（期待動作）

| ウィンドウ | サイズ | CGM反復 | 状態 |
|----------|--------|---------|------|
| 1 | 9 | ? | 実行中（停止） |

**総CGM実行回数**: 1ウィンドウ × ? 回

---

## 🎯 診断結果と結論

### 線形ソルバー反復数の差（2-3倍）の原因

**結論**: スライディングウィンドウ実装の差異ではなく、**各ウィンドウ内のCGM計算**における線形ソルバー呼び出し回数の違い

### 証拠

1. **Python版**: 9個の小ウィンドウで各3回CGMを実行
   - 各ウィンドウは異なるサイズ（9→1ステップ）
   - 小ウィンドウほど計算時間短縮（7.21s → 0.61s）

2. **Julia版**: 1個の大ウィンドウで計算中
   - 3分41秒経過でも完了せず
   - おそらくCGM内で大量の線形ソルバー反復

### 次のアクション

**推奨事項**:

1. **Juliaのループ終了条件を修正** → Python版と同じ動作に統一
2. **小ウィンドウでの比較** → 同じサイズのウィンドウで線形ソルバー反復数を比較
3. **詳細ログ追加** → Julia版CGMSolver内の線形ソルバー呼び出し回数を記録

---

## 📝 修正提案

### Julia版 `run_sliding_window.jl` の修正

**変更箇所**: 286行

```julia
# 修正前
while prev_flux_end < total_flux_steps

# 修正後（Python版に合わせる）
start_idx = 0
while start_idx < total_flux_steps
    max_L = min(window_size, total_flux_steps - start_idx)
    # ... (CGM計算) ...
    step = max(1, max_L - overlap)
    start_idx += step
end
```

この修正により、Julia版もPython版と同じ9ウィンドウを処理するようになります。

---

## 🔗 関連ファイル

- **Python実装**: `python/original/IHCP_CGM_Sliding_Window_Calculation_ver2.py` (1587-1659行)
- **Julia実装**: `julia/scripts/run_sliding_window.jl` (286-398行)
- **Python実行ログ**: `shared/results/python_scenario3.log`
- **Julia実行ログ**: `shared/results/julia_sliding_window_scenario3.log`

---

**次セッション優先タスク**:

1. ✅ スライディングウィンドウ実装差異特定 → **完了**
2. ⏭️ Julia版のループロジック修正
3. ⏭️ 同一ウィンドウでの線形ソルバー反復数比較
4. ⏭️ シナリオ3の再実行（修正後）
