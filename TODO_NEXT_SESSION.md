# 次セッション作業TODO

**作成日時**: 2025年10月21日
**ブランチ**: `sliding-window-validation`
**最新コミット**: `4b6e259` - バグのあるスクリプト削除完了

---

## 🎯 現在の状態

### ✅ 完了した作業

1. **バグのあるスクリプト削除（3ファイル）**
2. **関連ログファイル削除（6ファイル）**
3. **ドキュメント更新** - NEXT_SESSION_QUICKSTART.md全面改訂
4. **Gitコミット & プッシュ完了** - コミット4b6e259
5. **性能比較完了** - Python版45.57秒、Julia版105秒（2.3倍遅い）

### ❌ 未完了: 計算結果の数値的一致性が未検証 ⚠️

**理由**: Julia版スクリプトに結果保存機能がない

---

## 📋 次セッションの最優先タスク

### Phase 1: Julia版に--outputオプション追加（30分）

**実装内容**:
- julia/scripts/run_sliding_window.jlに結果保存機能追加
- 保存データ: q_result, T_result, J_history, 実行時間等
- NPZパッケージ使用: npzwrite(path, Dict(...))

**テスト**:
```bash
julia --project=julia julia/scripts/run_sliding_window.jl \
  --cgm-iter 3 --nt 10 \
  --output shared/results/julia_10steps_cgm3.npz
```

### Phase 2: 計算結果の数値比較（10分）

**実装済みスクリプト**: python/scripts/compare_results_10steps_cgm3.py ✅

```bash
python python/scripts/compare_results_10steps_cgm3.py
```

### Phase 3: 検証結果をレポートに追記（10分）

**更新対象**: PYTHON_JULIA_10STEPS_CGM3_COMPARISON.md

---

## 🗂️ 重要ファイル

- **要修正**: julia/scripts/run_sliding_window.jl
- **完成済み**: python/scripts/compare_results_10steps_cgm3.py
- **データ**: shared/results/python_10steps_cgm3.npz（存在）
- **未作成**: shared/results/julia_10steps_cgm3.npz ⚠️

---

## ⚡ 次セッション開始時の確認

```bash
ps aux | grep -E "(python|julia)" | grep -v grep
git status
git log --oneline -5
```

---

## ✅ チェックリスト

- [ ] Julia版に--outputオプション追加
- [ ] テスト実行して結果保存確認
- [ ] 比較スクリプト実行
- [ ] 数値的一致性確認
- [ ] レポート更新
- [ ] 未追跡ファイルコミット

---

**準備完了！** 次セッションではJulia版改善と計算精度検証を優先実施。
