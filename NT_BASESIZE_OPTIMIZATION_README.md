# nt-basesize最適化測定ガイド

このガイドは、nt値（タイムステップ数）とbasesize（並列化粒度）の最適な組み合わせを測定・分析するための手順書です。

---

## 概要

### 測定目的
- nt値(10/30/50/100) × basesize(400/500/600/700)の16パターンを測定
- nt値に応じた最適basesizeを特定
- 計算時間のスケーリング特性を把握

### 測定対象
- **スクリプト**: `julia/scripts/test_dhcp_solver.jl`
- **ソルバー**: DHCP単体（直接熱伝導問題）
- **格子サイズ**: 80×100×20 (N=160,000セル)
- **スレッド数**: 4固定

---

## クイックスタート

### 1. 測定実行（約1時間）

```bash
# プロジェクトルートに移動
cd /Users/Daily/Development/IHCP/TrialClaudeMCPCodex

# 測定スクリプトに実行権限付与
chmod +x run_nt_basesize_measurements.sh

# 測定開始
./run_nt_basesize_measurements.sh
```

**出力先**: `shared/results/nt_basesize/*.log` (12ファイル)

### 2. データ抽出（数秒）

```bash
# CSVファイル生成（既存nt=10データも含む）
julia julia/scripts/extract_nt_basesize_data.jl
```

**出力**: `shared/results/nt_basesize/measurement_results.csv`

### 3. レポート完成

```bash
# CSVデータを確認
cat shared/results/nt_basesize/measurement_results.csv

# レポートを開いて [TBD] 箇所を実データで更新
# ファイル: docs/reports/nt_basesize_optimization_report.md
```

---

## 詳細手順

### 測定の進行状況確認

#### 進捗ファイル
```bash
# 現在の進捗を確認
cat shared/results/nt_basesize/progress.txt
# 出力例: 5/12 → 12測定中5個完了

# 完了済み測定のリスト
cat shared/results/nt_basesize/completed.txt
# 出力例:
# nt30_bs400
# nt30_bs500
# ...
```

#### リアルタイムログ監視
```bash
# 最新のログファイルを監視
tail -f shared/results/nt_basesize/nt30_bs400.log
```

### 中断からの再開

測定が途中で中断した場合、再開可能です：

```bash
# --resumeオプションで再開
./run_nt_basesize_measurements.sh --resume
```

完了済み測定は自動的にスキップされます。

### 個別測定の実行

特定のnt/basesize組み合わせのみ再実行する場合：

```bash
# 例: nt=50, basesize=600のみ測定
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --solver pbicgstab \
  --precond gs \
  --nt 50 \
  --basesize 600 \
  2>&1 | tee shared/results/nt_basesize/nt50_bs600.log
```

---

## 測定時間の見積もり

### 1測定あたりの推定時間

| nt値 | 推定時間 | 根拠 |
|------|---------|------|
| 10 | 約5秒 | 既存データより |
| 30 | 約15秒 | nt値に線形スケール |
| 50 | 約25秒 | nt値に線形スケール |
| 100 | 約50秒 | nt値に線形スケール |

### 全測定の所要時間

```
合計: 12測定 (nt=30,50,100 × basesize=4種)
  nt=30: 4測定 × 15秒 = 1分
  nt=50: 4測定 × 25秒 = 1分40秒
  nt=100: 4測定 × 50秒 = 3分20秒

推定合計: 約6分 + オーバーヘッド = 10-15分
```

---

## データ構造

### ディレクトリ構成

```
shared/results/
├── nt_basesize/              # 新規測定データ
│   ├── nt30_bs400.log
│   ├── nt30_bs500.log
│   ├── nt30_bs600.log
│   ├── nt30_bs700.log
│   ├── nt50_bs400.log
│   ├── nt50_bs500.log
│   ├── nt50_bs600.log
│   ├── nt50_bs700.log
│   ├── nt100_bs400.log
│   ├── nt100_bs500.log
│   ├── nt100_bs600.log
│   ├── nt100_bs700.log
│   ├── progress.txt          # 進捗記録
│   ├── completed.txt         # 完了リスト
│   └── measurement_results.csv  # 集約データ
└── threads_basesize/         # 既存nt=10データ
    ├── step1_t4_bs400.log
    ├── step1_t4_bs500.log
    ├── step1_t4_bs600.log
    └── step1_t4_bs700.log
```

### CSVデータ形式

`measurement_results.csv`:
```csv
nt,basesize,threads,dhcp_time,total_time,solver,precond
10,400,4,5.06,6.85,pbicgstab,gs
10,500,4,4.91,6.71,pbicgstab,gs
...
```

---

## トラブルシューティング

### 測定が失敗する

#### 症状
```bash
[ERROR] Measurement may have failed (check log)
```

#### 対処法
1. ログファイルを確認:
   ```bash
   tail -50 shared/results/nt_basesize/nt30_bs400.log
   ```

2. エラーメッセージを確認して原因特定:
   - メモリ不足 → Juliaプロセスを再起動
   - ファイル未発見 → `T_measure_700um_1ms.npy`の存在確認
   - パッケージエラー → `julia --project=julia`で環境確認

3. 個別に再実行:
   ```bash
   JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
     --nt 30 --basesize 400
   ```

### データ抽出でエラー

#### 症状
```
WARN: DHCP elapsed time not found in: ...
```

#### 対処法
1. ログファイルに"Total runtime:"が含まれているか確認:
   ```bash
   grep "Total runtime:" shared/results/nt_basesize/nt30_bs400.log
   ```

2. ない場合は測定が完了していないため、再実行

### 進捗が更新されない

#### 対処法
```bash
# 進捗ファイルを手動確認
ls shared/results/nt_basesize/*.log | wc -l
# 期待値: 12（全測定完了時）

# 完了リストを手動更新
grep "Total runtime:" shared/results/nt_basesize/*.log | wc -l
```

---

## セッション継続ガイド

測定作業を別セッションで継続する場合の手順:

### 状態確認
```bash
# 1. 現在の進捗を確認
cat shared/results/nt_basesize/progress.txt

# 2. 完了済み測定を確認
cat shared/results/nt_basesize/completed.txt

# 3. ログファイル数を確認
ls shared/results/nt_basesize/nt*.log | wc -l
```

### 再開
```bash
# --resumeで中断箇所から自動再開
./run_nt_basesize_measurements.sh --resume
```

### データ抽出（測定完了後）
```bash
# 既存データ(nt=10)も含めてCSV化
julia julia/scripts/extract_nt_basesize_data.jl
```

### レポート完成
1. CSVデータを確認
2. `docs/reports/nt_basesize_optimization_report.md`を開く
3. [TBD]箇所を実データで更新
4. グラフ・考察を追加

---

## 期待される結果

### 既知の知見（nt=10での測定より）

- **最適basesize**: 600付近（4スレッド）
- **性能低下**: basesize≥1000で顕著
- **最適範囲**: basesize=400-600で安定

### 今回の測定で明らかにすること

1. **nt依存性**: nt値が変わると最適basesizeも変わるか?
2. **スケーリング**: nt値に対して実行時間は線形か?
3. **推奨設定**: 各nt値に対する推奨basesize

---

## 次のステップ

測定完了後の推奨作業:

1. **レポート完成**: 実データでテンプレート更新
2. **グラフ作成**: nt-basesize依存性の可視化
3. **理論分析**: 最適basesizeのnt依存性の考察
4. **コード更新**: 自動basesize選択関数の実装

---

## 参考資料

- **既存レポート**: `docs/reports/threads_basesize_optimization_report.md`
- **測定スクリプト**: `julia/scripts/test_dhcp_solver.jl`
- **抽出スクリプト**: `julia/scripts/extract_nt_basesize_data.jl`

---

## お問い合わせ

質問・不具合報告:
- このREADMEを参照
- ログファイルを確認
- 必要に応じてissue登録

---

**作成日**: 2025年10月25日
**バージョン**: 1.0
**作成者**: Claude Code
