# 次セッション クイックスタート

**📌 5分で始められる修正タスク**

## 🔥 最優先修正（5分）

### 問題
Julia版スライディングウィンドウでt=63から発散
→ **原因特定完了**: Y_obsスライスのオフバイワンエラー

### 修正箇所
**ファイル**: `julia/scripts/run_sliding_window_validation.jl:283`

**修正内容**:
```julia
# 修正前（❌ 72ステップ取得）
Y_obs_win = Y_obs[:, :, start_idx+1:end_idx+1]

# 修正後（✅ 71ステップ）
Y_obs_win = Y_obs[:, :, start_idx+1:start_idx+max_L+1]
```

### 実行コマンド
```bash
# 1. 修正実施（上記の1行を変更）

# 2. 再実行
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3

# 3. 成功確認（発散しないことを確認）
```

## 📋 詳細情報
→ `TODO_NEXT_SESSION.md`を参照

---
**最終更新**: 2025-10-21 03:44
