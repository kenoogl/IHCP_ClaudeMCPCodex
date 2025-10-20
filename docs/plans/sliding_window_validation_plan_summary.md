# スライディングウィンドウ検証計画 - クイックサマリー

**作成日**: 2025年10月21日
**詳細計画**: `sliding_window_validation_plan.md`

---

## 目的

オリジナルPythonコードとJuliaコードでスライディングウィンドウ計算を実行し、性能と精度を完全比較検証する。

---

## 2段階アプローチ

### Phase 1: オプション2（短時間動作確認）
- **CGM反復数**: 3回固定
- **目的**: 基本動作確認、スクリプト検証
- **所要時間**: Python約15分、Julia約5分

### Phase 2: オプション1（本格検証）
- **CGM反復数**: 20000回（早期停止あり）
- **目的**: オリジナルと同じ設定での完全比較
- **所要時間**: 状況により変動

---

## 実行パラメータ

```
格子: 80 × 100 × 20 (フルサイズ)
時間ステップ: 300ステップ
ウィンドウサイズ: 71
オーバーラップ: 17
予想ウィンドウ数: 約6
```

---

## 作業タスク

### 1. 出力フォーマット整備
- [ ] Python側: DHCP/Adjoint/Sensitivityソルバーの詳細ログ追加
- [ ] Julia側: 各ソルバーに反復回数・残差・時間計測追加

### 2. テストスクリプト作成
- [ ] `python/validation/run_sliding_window_validation.py`
- [ ] `julia/scripts/run_sliding_window_validation.jl`
- [ ] `python/validation/compare_sliding_window_results.py`

### 3. 実行と検証
- [ ] Phase 1実行（CGM 3回）
- [ ] Phase 1結果比較・問題修正
- [ ] Phase 2実行（CGM 20000回）
- [ ] Phase 2詳細分析

---

## 成功基準

### Phase 1（必須）
- [ ] 両スクリプトがエラーなく完走
- [ ] 熱流束相対誤差 < 5%
- [ ] ウィンドウ数の一致

### Phase 2（目標）
- [ ] 温度場相対誤差 < 2%
- [ ] 熱流束相対誤差 < 1%
- [ ] Julia実行時間がPythonの50%以下

---

## クイックスタート

```bash
# Phase 1: 動作確認
python python/validation/run_sliding_window_validation.py --cgm-iter 3
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 3
python python/validation/compare_sliding_window_results.py --phase 1

# Phase 2: 本格検証（Phase 1成功後）
python python/validation/run_sliding_window_validation.py --cgm-iter 20000
julia --project=julia julia/scripts/run_sliding_window_validation.jl --cgm-iter 20000
python python/validation/compare_sliding_window_results.py --phase 2
```

---

**次のステップ**: タスク3「Python側出力フォーマット整備」から開始
