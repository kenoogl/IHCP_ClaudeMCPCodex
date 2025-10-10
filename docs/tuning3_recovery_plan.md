# tuning3リカバリー計画書

## 📋 背景と目的

### 問題
- tuning2ブランチでPython-Julia間の計算結果が一致しない
- 複数の最適化が同時に適用されており、原因特定が困難

### 解決方針
1. コミット99b321f（安定版）から新ブランチtuning3を作成
2. tuning2の修正を段階的に適用し、各段階で検証
3. 問題発生箇所を特定し、修正

### 現状確認 ✅
- **日時**: 2025-10-10
- **ブランチ**: tuning3 (コミット99b321f)
- **検証**: フルサイズnt=10ステップ、シングルウインドウ
- **結果**: Python-Julia完全一致
  - 温度場T最大相対誤差: **1.36e-05 (0.00136%)**
  - 熱流束q最大絶対誤差: **1.98e-07 W/m²**

---

## 🎯 段階的適用計画

### Phase A: メモリレイアウト最適化
**コミット範囲**: 3368d60, c8c99b1, 2c00e6c, 4c87fcf

**変更内容**:
- 温度・観測・熱流束配列の軸順序を(ni,nj,nk,nt)系へ統一
- 主要ソルバーとテストのアクセスロジック更新（permutedims対応）
- 数値一致検証スクリプト追加
- ベンチマーク・プロファイリングスクリプト追加

**検証方法**:
1. Phase 1-6全テスト（505項目）実行
2. 10ステップフルサイズテスト実行（Julia版のみ）
3. Python-Julia結果比較（既存Python結果と比較）

**💡 検証効率化のヒント**:
- **Phase AではPython版の再実行は不要**（Julia側のみの変更のため）
- Python結果は準備フェーズで既に保存済み（`shared/results/python_10steps_fullsize.npz`）
- Julia版のみ実行して既存Python結果と比較すればOK

**期待される結果**:
- 全テスト合格
- Python-Julia一致維持（相対誤差 < 0.01%）
- 性能への影響確認

---

### Phase B: ガイドセル統合と計画立案
**コミット範囲**: 92d8a40, b1ec95f, 233cfc7, 5d993f9

**変更内容**:
- GridTransformモジュール導入（Z/ΔZ変換と境界管理）
- IHCP_CGMのinclude/export構造再構築
- Phase2.3b移行計画とドキュメント集約
- 共通境界条件ユーティリティ追加

**検証方法**:
1. Phase 1-6全テスト（505項目）実行
2. 10ステップフルサイズテスト実行
3. Python-Julia結果比較

**期待される結果**:
- 全テスト合格
- Python-Julia一致維持
- ガイドセル変換の正確性確認

---

### Phase C: マトリクスフリー化
**コミット範囲**: d9799c0, f72e9be, af36461, 15719bb, f2e5a70, 6af4798, ddf0c3d, 4692ae4

**変更内容**:
- CommonSolverにマトリクスフリーPBICGSTAB!実装
- DHCPソルバーからcg!実装撤去、WorkBuffers構造軽量化
- Adjointソルバーマトリクスフリー化、境界条件ユーティリティ拡張
- SensitivitySolver新設
- CGMソルバー新API対応、テストWorkBuffers対応
- 作業バッファ初期化関数追加（reset_work_buffers!, myfill!, mycopy!）

**検証方法**:
1. Phase 1-6全テスト（505項目）実行
2. 10ステップフルサイズテスト実行
3. Python-Julia結果比較
4. 性能測定（期待: 89倍高速化）

**期待される結果**:
- 全テスト合格
- Python-Julia一致維持
- 大幅な性能向上確認

**⚠️ 重要**: この段階で最も大きな変更があり、問題が発生する可能性が高い

---

### Phase D: コード整理と安定化
**コミット範囲**: 21cbfff, 9c68276, 3e27f19, 49df97e, 108bc95, 14145c8

**変更内容**:
- FLoopsループでlet束縛導入（BoxedVariable警告解消）
- 共通RHSロジック抽出（RHSCoreモジュール）
- 未使用依存・import削除
- 境界条件やログ制御の統一
- レビュー指摘事項の修正（高優先度バグ、可読性・保守性向上）

**検証方法**:
1. Phase 1-6全テスト（505項目）実行
2. 10ステップフルサイズテスト実行
3. Python-Julia結果比較

**期待される結果**:
- 全テスト合格
- Python-Julia一致維持
- 型安定性向上確認

---

## ⚠️ 作業上の注意点

### 一般的な注意事項
1. **各Phase完了後に必ずコミット**
   - コミットメッセージに `[tuning3-PhaseX]` プレフィックスを付ける
   - 問題発生時のロールバックを容易にする

2. **検証を必ず実施**
   - Phase 1-6全テスト（505項目）
   - 10ステップフルサイズテスト
   - Python-Julia結果比較

3. **問題発生時の対応**
   - 該当Phaseの適用を一時中断
   - 前のPhaseにロールバック
   - 問題箇所を特定し、個別に修正
   - 再度検証

4. **進捗記録**
   - このドキュメントの「進捗チェックリスト」を更新
   - 問題が発生した場合は「トラブルシューティング記録」に記載

5. **ユーザーインタラクション**
   - 処理を行う前に適用する処理の概要を伝え、実行して良いか必ずユーザに確認する。

### コード変更時の注意
1. **配列の軸順序に注意**
   - tuning3（Phase A前）: (nt, ni, nj, nk)
   - tuning3（Phase A後）: (ni, nj, nk, nt)
   - Phase A適用時は特に注意

2. **API変更に注意**
   - Phase C前: `solve_cgm!(T_init, Y_obs, q_init, dx, dy, dz, dz_b, dz_t, ...)`
   - Phase C後: `solve_cgm!(T_init, Y_obs, q_init, work, dx, dy, Z, dZ_guard, ...)`

3. **ガイドセル変換**
   - Phase B以降: `convert_to_guard_cell_grid(nk, dz, dz_bottom, dz_top)`が利用可能
   - WorkBuffersサイズ: `(ni+2, nj+2, nk+2)` に注意

### 検証スクリプト
- **全テスト実行**: `julia --project=. -e 'using Pkg; Pkg.test()'`
- **10ステップテスト**:
  ```bash
  # Python版（Phase A-Dでは通常不要、準備フェーズで実行済み）
  python python/validation/run_10steps_fullsize_test.py

  # Julia版（これだけ実行すればOK）
  julia --project=. julia/scripts/run_10steps_fullsize_test.jl

  # 比較（既存のPython結果と比較）
  python python/validation/compare_python_julia_10steps_fullsize.py
  ```

**💡 ヒント**: Phase A-DではJulia側のみ変更するため、Python版の再実行は不要です。Julia版実行後に既存のPython結果（`shared/results/python_10steps_fullsize.npz`）と比較してください。

---

## 📊 進捗チェックリスト

### 準備フェーズ ✅
- [x] tuning3ブランチ作成（コミット99b321fから）
- [x] 10ステップフルサイズ検証テストコード作成
- [x] Python-Julia一致確認（ベースライン）
- [x] 計画書作成（tuning3_recovery_plan.md, tuning3_quick_reference.md）
- [x] CLAUDE.md更新（tuning3作業対応）
- [x] run_10steps_fullsize_test.jl修正（rho_coeffs削除）
- [x] 最終動作確認（Python-Julia一致維持）

### Phase A: メモリレイアウト最適化 ✅
- [x] コミット3368d60適用
- [x] コミットc8c99b1適用
- [x] コミット2c00e6c適用
- [x] コミット4c87fcf適用
- [x] 全テスト（505項目）実行・合格確認
- [x] 10ステップフルサイズテスト実行
- [x] Python-Julia一致確認
- [x] Phase Aコミット作成（72370b2）

### Phase B: ガイドセル統合
- [ ] コミット92d8a40適用
- [ ] コミットb1ec95f適用
- [ ] コミット233cfc7適用
- [ ] コミット5d993f9適用
- [ ] 全テスト（505項目）実行・合格確認
- [ ] 10ステップフルサイズテスト実行
- [ ] Python-Julia一致確認
- [ ] Phase Bコミット作成

### Phase C: マトリクスフリー化
- [ ] コミットd9799c0適用
- [ ] コミットf72e9be適用
- [ ] コミットaf36461適用
- [ ] コミット15719bb適用
- [ ] コミットf2e5a70適用
- [ ] コミット6af4798適用
- [ ] コミットddf0c3d適用
- [ ] コミット4692ae4適用
- [ ] 全テスト（505項目）実行・合格確認
- [ ] 10ステップフルサイズテスト実行
- [ ] Python-Julia一致確認
- [ ] 性能測定（89倍高速化確認）
- [ ] Phase Cコミット作成

### Phase D: コード整理と安定化
- [ ] コミット21cbfff適用
- [ ] コミット9c68276適用
- [ ] コミット3e27f19適用
- [ ] コミット49df97e適用
- [ ] コミット108bc95適用
- [ ] コミット14145c8適用
- [ ] 全テスト（505項目）実行・合格確認
- [ ] 10ステップフルサイズテスト実行
- [ ] Python-Julia一致確認
- [ ] Phase Dコミット作成

### 最終確認
- [ ] 全Phase完了確認
- [ ] 最終的な全テスト実行
- [ ] 最終的なPython-Julia一致確認
- [ ] tuning3ブランチをmainにマージ（または新しいブランチ名でプッシュ）

---

## 🔧 トラブルシューティング記録

### 問題発生時の記録フォーマット
```markdown
### [YYYY-MM-DD] Phase X: 問題タイトル

**症状**:
-

**原因**:
-

**対応**:
-

**結果**:
-
```

### 既知の問題

#### [2025-10-10] 準備: Codex作成スクリプトのAPI不一致
**症状**:
- Codexが作成したJuliaスクリプトがtuning2以降のAPI（convert_to_guard_cell_grid, WorkBuffers）を使用
- コミット99b321fには存在しない関数のため実行エラー

**原因**:
- Codexが最新のコードベースを参照してスクリプトを作成

**対応**:
- 手動でコミット99b321fのAPI（dz, dz_bottom, dz_top）に修正
- WorkBuffers呼び出しをコメントアウト
- 配列順序を(nt, ni, nj)に修正

**結果**:
- 実行成功、Python-Julia完全一致確認

---

#### [2025-10-10] 準備: run_10steps_fullsize_test.jlの不要なrho_coeffs
**症状**:
- `load_material_properties()`が不要な`rho_coeffs`を戻り値として返していた
- `rho`は定数として扱うため、多項式係数は不要

**原因**:
- Codexが作成時に誤って含めた（cp_coeffsと混同した可能性）

**対応**（コミット21742a1）:
- `load_material_properties()`から`rho_coeffs`定義と戻り値を削除
- 呼び出し側の受け取り変数を修正
- メタデータ保存から`rho_coeffs`を削除

**結果**:
- 正常動作、Python-Julia一致維持（温度T相対誤差1.36e-05、熱流束q絶対誤差1.98e-07）

---

## 📈 性能測定記録

### ベースライン（tuning3初期状態、コミット21742a1）
- **CGM実行時間**: 22.72秒（10ステップ、フルサイズ）
- **Total実行時間**: 25.62秒
- **メモリ使用量**: TBD
- **検証結果**: 温度T最大相対誤差 1.36e-05、熱流束q最大絶対誤差 1.98e-07 W/m²

### Phase A完了時 ✅
- **CGM実行時間**: 22.65秒（10ステップ、フルサイズ）
- **Total実行時間**: 25.49秒
- **メモリ使用量**: TBD
- **性能変化**: ベースラインとほぼ同等（25.62秒 → 25.49秒）
- **検証結果**: 温度T最大相対誤差 1.36e-05、熱流束q最大絶対誤差 1.98e-07 W/m²

### Phase B完了時
- **CGM実行時間**: TBD
- **Total実行時間**: TBD
- **メモリ使用量**: TBD
- **性能変化**: TBD

### Phase C完了時（期待: 89倍高速化）
- **CGM実行時間**: TBD（期待: ~0.26秒）
- **Total実行時間**: TBD
- **メモリ使用量**: TBD
- **性能変化**: TBD

### Phase D完了時
- **CGM実行時間**: TBD
- **Total実行時間**: TBD
- **メモリ使用量**: TBD
- **性能変化**: TBD

---

## 📝 作業ログ

### 2025-10-10 - 準備フェーズ完了 ✅

#### セッション1（午前）
**コミット**: 89c6a2e, 90240d2, 21742a1

1. **tuning3ブランチ作成**（コミット99b321fから）
2. **10ステップフルサイズ検証テストコード作成**
   - `python/validation/run_10steps_fullsize_test.py`
   - `julia/scripts/run_10steps_fullsize_test.jl`
   - `python/validation/compare_python_julia_10steps_fullsize.py`
   - `config/test_10steps_fullsize.toml`
   - Julia版スクリプトをコミット99b321fのAPIに手動修正
3. **ベースライン検証完了**
   - Python版実行: 4.64秒
   - Julia版実行: 25.85秒
   - 温度場T最大相対誤差: 1.36e-05 (0.00136%)
   - 熱流束q最大絶対誤差: 1.98e-07 W/m²
   - ✅ Python-Julia完全一致確認
4. **ドキュメント作成**
   - `docs/tuning3_recovery_plan.md`: 詳細計画・進捗チェックリスト
   - `docs/tuning3_quick_reference.md`: 実務手順・コマンド集
   - コミット89c6a2e
5. **CLAUDE.md更新**（コミット90240d2）
   - tuning3作業セクション追加
   - コンテキスト管理強化（tuning3ドキュメント参照追加）
   - 不要な詳細削除、現在の状態に焦点
6. **run_10steps_fullsize_test.jl修正**（コミット21742a1）
   - 不要なrho_coeffsを削除（rhoは定数として扱う）
   - 動作確認完了（25.62秒、Python-Julia一致維持）

**準備完了**: Phase A開始準備完了、全ドキュメント整備完了

---

## 🔗 参考リンク

- **tuning2コミット履歴**: `git log 99b321f..tuning2 --oneline`
- **Codex調査結果**: 本セッション冒頭のCodex出力参照
- **検証スクリプト**:
  - `python/validation/run_10steps_fullsize_test.py`
  - `julia/scripts/run_10steps_fullsize_test.jl`
  - `python/validation/compare_python_julia_10steps_fullsize.py`
- **ベースライン結果**: `shared/results/python_10steps_fullsize.npz`, `julia_10steps_fullsize.npz`

---

## 次回セッションでの作業開始方法

1. **ブランチ確認**
   ```bash
   git branch  # tuning3であることを確認
   git log --oneline -5  # 現在位置確認（21742a1が最新）
   ```

2. **このドキュメント確認**
   ```bash
   cat docs/tuning3_recovery_plan.md
   ```

3. **進捗チェックリスト確認**
   - **準備フェーズ**: ✅ 完了
   - **次のステップ**: Phase A（メモリレイアウト最適化）開始

4. **作業開始**
   - Phase Aのコミット3368d60から適用開始
   - 各コミット後に検証実行
   - 結果記録

5. **クイックリファレンス参照**
   ```bash
   cat docs/tuning3_quick_reference.md
   ```

**⚠️ 必ず各Phase完了後にこのドキュメントを更新してください！**

---

## 📌 現在の状態（2025-10-10）

- **ブランチ**: tuning3
- **最新コミット**: 72370b2（[tuning3-PhaseA] メモリレイアウト最適化完了）
- **準備フェーズ**: ✅ 完了
- **Phase A**: ✅ 完了
- **次のステップ**: Phase B（ガイドセル統合）開始
- **Phase A性能**: CGM 22.65秒、Total 25.49秒
- **Phase A精度**: 温度T相対誤差1.36e-05、熱流束q絶対誤差1.98e-07 W/m²
