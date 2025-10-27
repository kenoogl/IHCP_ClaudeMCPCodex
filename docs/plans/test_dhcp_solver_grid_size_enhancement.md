# test_dhcp_solver.jl 可変格子サイズ対応 修正案

**作成日**: 2025年10月27日
**目的**: 格子解像度を可変にして問題サイズに対するスケーラビリティ測定を可能にする
**対象ファイル**: `julia/scripts/test_dhcp_solver.jl`

---

## 修正の動機

### 現状の課題

- 格子サイズが80×100×20に固定されている
- 問題サイズに対するスケーラビリティ測定ができない
- 測定データ（1.1GB）に依存するため、大規模問題のテストが困難

### 期待される効果

1. **柔軟な性能測定**
   - 問題サイズ（N）に対する計算時間の依存性評価（O(N), O(N log N)等）
   - メモリ使用量のスケーリング測定
   - キャッシュ効率の分析

2. **測定データ依存の解消**
   - 任意サイズでテスト可能
   - 大規模問題のテスト（測定データサイズを超える）
   - weak scaling / strong scalingテストの実施

3. **互換性の維持**
   - デフォルト動作は変更なし
   - 既存のテストスクリプトは影響なし

---

## 修正内容の概要

### 1. 新規コマンドライン引数

| オプション | 説明 | デフォルト値 | 型 |
|-----------|------|------------|-----|
| `--ni` | X方向格子数 | 80 | Int |
| `--nj` | Y方向格子数 | 100 | Int |
| `--nk` | Z方向格子数 | 20 | Int |
| `--synthetic` | 合成テストモード | false | Bool |

### 2. 動作モードの追加

#### 測定データモード（デフォルト）
- 既存の動作を維持
- `T_measure_700um_1ms.npy`から測定データを読み込み
- 格子サイズは測定データと一致する必要あり（80×100）

#### 合成テストモード（新規）
- `--synthetic` フラグで有効化
- 測定データを使わない
- 任意の格子サイズで実行可能
- 初期温度場は一様温度または勾配温度

---

## 詳細設計

### 1. parse_command_line_args()の拡張

```julia
function parse_command_line_args()
  solver_type = :pbicgstab
  precond_type = :diagonal
  nt = 10
  basesize = 600

  # 新規追加
  ni = 80
  nj = 100
  nk = 20
  synthetic = false

  # ... オプション解析 ...

  return solver_type, precond_type, nt, basesize, ni, nj, nk, synthetic
end
```

**変更点**:
- 戻り値に`ni, nj, nk, synthetic`を追加
- `--ni`, `--nj`, `--nk`, `--synthetic`オプションを解析

### 2. prepare_synthetic_test()関数の追加

```julia
"""
    prepare_synthetic_test(ni, nj, nk, nt) -> (T_init, Y_obs)

合成テストモード用の初期化

# 引数
- ni, nj, nk: 格子数
- nt: タイムステップ数

# 戻り値
- T_init: 初期温度場 (ni, nj, nk) [K]
- Y_obs: 観測データ（合成） (ni, nj, nt) [K]

# 温度設定
オプション1: 一様温度（T_base = 573.15K）
オプション2: 温度勾配（z方向に10K/層ずつ上昇）
"""
function prepare_synthetic_test(ni::Int, nj::Int, nk::Int, nt::Int)
  println("\n[2/5] Preparing synthetic test data")
  flush(stdout)

  # 一様温度
  T_base = 573.15  # 300°C in Kelvin
  T_init = fill(T_base, ni, nj, nk)
  Y_obs = fill(T_base, ni, nj, nt)

  println("  Mode: Synthetic test")
  println(@sprintf("  Initial temperature: %.2f K (uniform)", T_base))
  println(@sprintf("  Grid size: %d × %d × %d (N=%d)", ni, nj, nk, ni*nj*nk))
  flush(stdout)

  return T_init, Y_obs
end
```

**実装オプション**:

#### オプション1: 一様温度（推奨）
```julia
T_base = 573.15  # 300°C
T_init = fill(T_base, ni, nj, nk)
Y_obs = fill(T_base, ni, nj, nt)
```
- 最もシンプル
- 性能測定に集中
- q=0での挙動が明確

#### オプション2: 温度勾配
```julia
T_base = 573.15
T_init = zeros(ni, nj, nk)
for k in 1:nk
  T_init[:, :, k] .= T_base + (k - 1) * 10.0  # 10K/層
end
Y_obs = repeat(T_init[:, :, 1:1], 1, 1, nt)
```
- より現実的
- 熱拡散の影響を確認可能

### 3. prepare_measurement_test()関数の抽出

```julia
"""
    prepare_measurement_test(ni, nj, nk, nt) -> (T_init, Y_obs)

測定データモード用の初期化（既存の処理を関数化）

# サイズ検証
測定データのサイズ（80×100）と一致する必要あり
"""
function prepare_measurement_test(ni::Int, nj::Int, nk::Int, nt::Int)
  println("\n[2/5] Loading measurement data")
  flush(stdout)

  Y_obs_python, ni_file, nj_file = load_measurement_subset(nt)

  # サイズ検証
  if ni_file != ni || nj_file != nj
    error("Measurement grid mismatch: file has $(ni_file)×$(nj_file), expected $(ni)×$(nj)")
  end

  # メモリレイアウト最適化
  Y_obs = permutedims(Y_obs_python, (2, 3, 1))
  Y_obs_python = nothing

  T_init = build_initial_temperature(@view(Y_obs[:, :, 1]), nk)
  println(@sprintf("  initial temperature range: %.2f ~ %.2f K", minimum(T_init), maximum(T_init)))
  flush(stdout)

  return T_init, Y_obs
end
```

### 4. analyze_residuals()関数の追加

```julia
"""
    analyze_residuals(T_result, Y_obs, synthetic)

残差分析（モード依存）

# 合成モード
- 温度統計のみ表示
- 残差計算なし（Y_obsが意味を持たない）

# 測定データモード
- RMS残差、最大残差を計算
- 測定データとの詳細比較
"""
function analyze_residuals(T_result::AbstractArray{<:Real,4},
                          Y_obs::AbstractArray{<:Real,3},
                          synthetic::Bool)
  T_bottom = @view T_result[:, :, 1, :]  # (ni, nj, nt)

  if synthetic
    # 合成モード: 簡易分析
    println("  Mode: Synthetic test")
    println(@sprintf("  Mean temperature: %.2f K", mean(T_result)))
    println(@sprintf("  Temperature range: %.2f ~ %.2f K", minimum(T_result), maximum(T_result)))
    println(@sprintf("  Temperature change: %.4e K (max-min)", maximum(T_result) - minimum(T_result)))
    println("  Note: No residual analysis (synthetic mode)")
  else
    # 測定データモード: 詳細分析
    residual = T_bottom .- Y_obs
    rms_error = sqrt(mean(residual .^ 2))
    max_error = maximum(abs.(residual))
    mean_temp = mean(T_result)

    println("  Mode: Measurement data")
    println(@sprintf("  RMS residual: %.4e K", rms_error))
    println(@sprintf("  Max residual: %.4e K", max_error))
    println(@sprintf("  Mean temperature: %.2f K", mean_temp))
    println(@sprintf("  Temperature range: %.2f ~ %.2f K", minimum(T_result), maximum(T_result)))
  end
  flush(stdout)
end
```

### 5. main()関数の修正

```julia
function main()
  println("="^80)
  println("DHCP（直接熱伝導問題）ソルバー単体テスト")
  println("="^80)

  # コマンドライン引数パース（拡張版）
  solver_type, precond_type, nt, basesize, ni, nj, nk, synthetic = parse_command_line_args()

  # Backend設定
  set_backend_config(basesize=basesize)

  println("\n[Configuration]")
  println("  Solver: $solver_type")
  println("  Preconditioner: $precond_type")
  println("  Time steps: $nt")
  println("  Grid size: $ni × $nj × $nk (N=$(ni*nj*nk))")  # 追加
  println("  FLoops basesize: $basesize")
  println("  Julia threads: $(Threads.nthreads())")
  println("  Test mode: $(synthetic ? "Synthetic" : "Measurement data")")  # 追加

  # ... 問題定義 ...

  # モード分岐（重要な変更点）
  if synthetic
    T_init, Y_obs = prepare_synthetic_test(ni, nj, nk, nt)
  else
    T_init, Y_obs = prepare_measurement_test(ni, nj, nk, nt)
  end

  # ... DHCP solve（共通処理）...

  # 残差分析（モード依存）
  analyze_residuals(T_result, Y_obs, synthetic)

  # ... Summary ...
end
```

---

## 使用例

### 1. デフォルト（互換性テスト）

```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl
```

**期待される動作**:
- 測定データモード
- 格子サイズ: 80×100×20
- 既存の動作と完全に同一

### 2. 合成テストモード（小規模）

```bash
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 40 --nj 50 --nk 10 --nt 10 --synthetic
```

**期待される結果**:
- 格子数: 40×50×10 (N=20,000)
- 実行時間: 約1-2秒
- メモリ使用量: 約100MB

### 3. 合成テストモード（中規模）

```bash
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 80 --nj 100 --nk 20 --nt 10 --synthetic --basesize 600
```

**期待される結果**:
- 格子数: 80×100×20 (N=160,000)
- 実行時間: 約3-5秒
- 測定データモードと同等の性能

### 4. 合成テストモード（大規模）

```bash
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 160 --nj 200 --nk 40 --nt 10 --synthetic --basesize 600
```

**期待される結果**:
- 格子数: 160×200×40 (N=1,280,000、8倍サイズ）
- 実行時間: 約20-40秒（理想的には8倍）
- メモリ使用量: 約8GB

### 5. スケーラビリティテスト

#### 問題サイズを変えて測定

```bash
# N=20,000
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 40 --nj 50 --nk 10 --nt 10 --synthetic --basesize 600 \
  2>&1 | tee results/scaling_N20k.log

# N=80,000
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 60 --nj 80 --nk 15 --nt 10 --synthetic --basesize 600 \
  2>&1 | tee results/scaling_N80k.log

# N=160,000
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 80 --nj 100 --nk 20 --nt 10 --synthetic --basesize 600 \
  2>&1 | tee results/scaling_N160k.log

# N=320,000
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 100 --nj 128 --nk 25 --nt 10 --synthetic --basesize 600 \
  2>&1 | tee results/scaling_N320k.log
```

#### Weak scalingテスト

```bash
# 1スレッド、N=20,000
JULIA_NUM_THREADS=1 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 40 --nj 50 --nk 10 --synthetic --basesize 600

# 2スレッド、N=40,000
JULIA_NUM_THREADS=2 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 50 --nj 64 --nk 12 --synthetic --basesize 600

# 4スレッド、N=80,000
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 60 --nj 80 --nk 15 --synthetic --basesize 600
```

**期待**: スレッド数を増やしても実行時間がほぼ一定

---

## テスト計画

### Phase 1: 互換性テスト

**目的**: 既存の動作が維持されることを確認

```bash
# 既存のデフォルト動作
julia --project=julia julia/scripts/test_dhcp_solver.jl

# 明示的に同じパラメータ指定
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 80 --nj 100 --nk 20 --nt 10
```

**期待される結果**:
- 両方の実行結果が完全に一致
- DHCP時間: 約4.86s
- RMS residual: 2.9516e-01 K

### Phase 2: 合成モード基本テスト

**目的**: 合成テストモードが正常に動作することを確認

```bash
# 同じサイズで合成モード
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 80 --nj 100 --nk 20 --nt 10 --synthetic
```

**期待される結果**:
- DHCP時間: 約4-5s（測定データモードと同等）
- Temperature change: ほぼゼロ（q=0のため）
- No residual analysis メッセージ

### Phase 3: 小規模問題テスト

**目的**: 小さい問題サイズでも正常動作することを確認

```bash
# 最小サイズ
julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 10 --nj 10 --nk 5 --nt 5 --synthetic
```

**期待される結果**:
- 実行時間: 1秒以内
- 正常に収束

### Phase 4: 大規模問題テスト

**目的**: メモリ管理が適切に機能することを確認

```bash
# 8倍サイズ
JULIA_NUM_THREADS=4 julia --project=julia julia/scripts/test_dhcp_solver.jl \
  --ni 160 --nj 200 --nk 40 --nt 10 --synthetic --basesize 600
```

**期待される結果**:
- メモリエラーなし
- 実行時間: 約30-50秒
- 収束性が維持される

### Phase 5: スケーラビリティ測定

**目的**: N vs 計算時間の関係を定量評価

```bash
# 一連のサイズで測定
for size in "40 50 10" "60 75 15" "80 100 20" "100 125 25" "120 150 30"
do
  ni=$(echo $size | awk '{print $1}')
  nj=$(echo $size | awk '{print $2}')
  nk=$(echo $size | awk '{print $3}')
  N=$((ni * nj * nk))

  echo "=== Testing N=$N (${ni}×${nj}×${nk}) ==="
  julia --project=julia julia/scripts/test_dhcp_solver.jl \
    --ni $ni --nj $nj --nk $nk --nt 10 --synthetic --basesize 600 \
    2>&1 | tee results/scaling_N${N}.log
done
```

**分析項目**:
- log(実行時間) vs log(N) のプロット → 傾きからO(N^α)を推定
- 反復回数 vs N の関係
- メモリ使用量 vs N の関係

---

## 期待される成果

### 1. 性能特性の定量評価

- **計算複雑度**: O(N)、O(N log N)等を実測で確認
- **スケーラビリティ限界**: どのサイズから性能劣化が始まるか
- **最適basesize**: 問題サイズに応じた最適値

### 2. メモリスケーリング

- **メモリ使用量**: N vs メモリのプロット
- **キャッシュ効率**: サイズ変化に伴うキャッシュミス率
- **メモリ帯域幅**: 実効帯域幅の測定

### 3. アーキテクチャ依存性

- **Apple Silicon特性**: 統合メモリの効果
- **スレッド数依存性**: 1,2,4スレッドでの比較
- **NUMA効果**: （該当する場合）

---

## リスクと対策

### リスク1: メモリ不足

**問題**: 大規模問題でメモリ不足が発生

**対策**:
- 段階的にサイズを増やす
- メモリ使用量を監視
- スワップを避ける

### リスク2: 収束性の悪化

**問題**: 大規模問題で収束しない

**対策**:
- maxiterを増やす
- rtolを緩める（一時的）
- 前処理器を変更（gs → diagonal）

### リスク3: 計算時間の爆発的増加

**問題**: 予想以上に時間がかかる

**対策**:
- タイムアウトを設定
- 段階的に増やす
- 必要に応じてntを減らす

---

## 実装チェックリスト

- [ ] `parse_command_line_args()`の拡張
  - [ ] `--ni`, `--nj`, `--nk`オプション追加
  - [ ] `--synthetic`オプション追加
  - [ ] ヘルプメッセージ更新
  - [ ] 戻り値に`ni, nj, nk, synthetic`追加

- [ ] `prepare_synthetic_test()`関数の実装
  - [ ] 一様温度初期化
  - [ ] ログ出力
  - [ ] 戻り値`(T_init, Y_obs)`

- [ ] `prepare_measurement_test()`関数の抽出
  - [ ] 既存コードを関数化
  - [ ] サイズ検証
  - [ ] 戻り値`(T_init, Y_obs)`

- [ ] `analyze_residuals()`関数の実装
  - [ ] 合成モード分岐
  - [ ] 測定データモード分岐
  - [ ] 適切なログ出力

- [ ] `main()`関数の修正
  - [ ] 引数パースの更新
  - [ ] モード分岐の追加
  - [ ] ログ出力の更新

- [ ] テストの実施
  - [ ] Phase 1: 互換性テスト
  - [ ] Phase 2: 合成モード基本テスト
  - [ ] Phase 3: 小規模問題テスト
  - [ ] Phase 4: 大規模問題テスト
  - [ ] Phase 5: スケーラビリティ測定

- [ ] ドキュメント更新
  - [ ] `test_dhcp_solver_guide.md`更新
  - [ ] READMEに新オプション追記
  - [ ] 使用例の追加

---

## 次セッションでの作業手順

1. **修正実装**（30分）
   ```bash
   # バックアップ作成
   cp julia/scripts/test_dhcp_solver.jl julia/scripts/test_dhcp_solver.jl.bak

   # 修正実施
   # - parse_command_line_args()
   # - prepare_synthetic_test()
   # - prepare_measurement_test()
   # - analyze_residuals()
   # - main()
   ```

2. **Phase 1テスト**（10分）
   ```bash
   # 互換性確認
   julia --project=julia julia/scripts/test_dhcp_solver.jl
   ```

3. **Phase 2テスト**（10分）
   ```bash
   # 合成モード確認
   julia --project=julia julia/scripts/test_dhcp_solver.jl --synthetic
   ```

4. **Phase 3-4テスト**（20分）
   ```bash
   # 小規模・大規模テスト
   ```

5. **Phase 5スケーラビリティ測定**（30分）
   ```bash
   # 複数サイズで測定
   ```

6. **結果分析とレポート作成**（20分）

**合計所要時間**: 約2時間

---

**作成者**: Claude Code
**最終更新**: 2025年10月27日
