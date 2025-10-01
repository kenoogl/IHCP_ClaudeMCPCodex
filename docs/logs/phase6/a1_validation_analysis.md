# Phase 6 A-1: Python検証機能詳細調査レポート

**作成日**: 2025-10-01
**作成者**: Codex MCP Agent
**目的**: Julia実装のための検証器関数群の仕様調査

---

## 1. 検証機能の概要

Python原典コード（1100-1260行付近）および関連する検証システム全体を調査した結果、**包括的な数値異常検出システム**と**メモリ監視システム**が実装されていることを確認しました。

---

## 2. 発見した検証機能の一覧

### 2.1 数値異常検出機能（Numerical Anomaly Detection）

| 関数名 | 行番号 | 概要 |
|--------|--------|------|
| `check_field_finite()` | 87-104 | NaN/Inf値の検出（高速版） |
| `check_temperature_range()` | 106-138 | 温度場の物理的範囲チェック（Numba高速化） |
| `check_flux_range()` | 140-171 | 熱流束の範囲チェック（Numba高速化） |
| `check_gradient_magnitude()` | 173-251 | 空間勾配の急変検出（2D/3D対応、Numba高速化） |
| `detect_numerical_anomalies()` | 253-332 | 包括的な異常検出メインルーチン |
| `handle_numerical_anomaly()` | 334-373 | 異常検出時の警告表示 |
| `check_temperature_field()` | 375-404 | 温度場の包括的チェック（ラッパー関数） |
| `check_flux_field()` | 406-435 | 熱流束場の包括的チェック（ラッパー関数） |
| `check_adjoint_field()` | 437-468 | 随伴場の包括的チェック（ラッパー関数） |
| `handle_critical_anomaly()` | 470-518 | 重大異常時の緊急対応 |
| `enhanced_anomaly_response()` | 520-602 | 拡張された異常対応システム |

### 2.2 メモリ監視機能（Memory Monitoring）

| 関数名 | 行番号 | 概要 |
|--------|--------|------|
| `get_memory_info()` | 608-638 | システムメモリ情報の取得 |
| `safe_file_size_check()` | 640-662 | ファイルサイズの安全性チェック |
| `monitor_memory_usage()` | 664-707 | メモリ使用量の監視と警告 |
| `safe_load_large_data()` | 709-790 | 大容量データの安全読み込み |
| `optimize_memory_usage()` | 792-822 | メモリ最適化（GC実行） |
| `memory_warning_system()` | 824-851 | メモリ警告システム |
| `check_memory_critical()` | 853-890 | 重大なメモリ不足への対応 |

### 2.3 数値安定性チェック

| 関数名 | 行番号 | 概要 |
|--------|--------|------|
| `check_diffusion_stability()` | 894-963 | 拡散数による数値安定性チェック |

---

## 3. 各検証機能の詳細仕様

### 3.1 基本検証関数

#### 3.1.1 `check_field_finite(field)` [行87-104]

**目的**: NaN/Inf値の高速検出

**引数**:
- `field`: numpy.ndarray - 検査対象フィールド

**戻り値**:
- `bool`: すべて有限の場合True

**ロジック**:
```python
flat_field = field.ravel()
for i in range(flat_field.size):
    if not np.isfinite(flat_field[i]):
        return False
return True
```

**特徴**:
- NumPyのブロードキャスト演算を使わず、明示的ループで高速化
- すべての要素を線形走査

---

#### 3.1.2 `check_temperature_range(T, min_temp=150.0, max_temp=3000.0)` [行106-138]

**目的**: 温度場の物理的範囲チェック（Numba高速化）

**引数**:
- `T`: numpy.ndarray - 温度場 [K]
- `min_temp`: float - 最小許容温度 [K] (デフォルト: 150.0K)
- `max_temp`: float - 最大許容温度 [K] (デフォルト: 3000.0K)

**戻り値**:
- `is_valid`: bool - 範囲内の場合True
- `min_val`: float - 実際の最小値 [K]
- `max_val`: float - 実際の最大値 [K]

**ロジック**:
```python
flat_T = T.ravel()
min_val = flat_T[0]
max_val = flat_T[0]
for i in range(flat_T.size):
    val = flat_T[i]
    if val < min_val: min_val = val
    if val > max_val: max_val = val
is_valid = (min_val >= min_temp) and (max_val <= max_temp)
```

**デコレータ**: `@njit` - Numbaによる高速化

---

#### 3.1.3 `check_flux_range(q, max_abs_flux=1e7)` [行140-171]

**目的**: 熱流束の範囲チェック（Numba高速化）

**引数**:
- `q`: numpy.ndarray - 熱流束 [W/m²]
- `max_abs_flux`: float - 最大許容絶対値 [W/m²] (デフォルト: 1e7)

**戻り値**:
- `is_valid`: bool - 範囲内の場合True
- `min_val`: float - 実際の最小値 [W/m²]
- `max_val`: float - 実際の最大値 [W/m²]

**ロジック**:
```python
flat_q = q.ravel()
min_val = flat_q[0]
max_val = flat_q[0]
for i in range(flat_q.size):
    val = flat_q[i]
    if val < min_val: min_val = val
    if val > max_val: max_val = val
max_abs = max(abs(min_val), abs(max_val))
is_valid = max_abs <= max_abs_flux
```

**デコレータ**: `@njit` - Numbaによる高速化

---

#### 3.1.4 `check_gradient_magnitude(field, dx, dy, dz, max_grad_temp=5000.0, max_grad_flux=1e8)` [行173-251]

**目的**: 空間勾配の急変検出（2D/3D対応、Numba高速化）

**引数**:
- `field`: numpy.ndarray - フィールド（温度または熱流束）
- `dx`, `dy`: float - x, y方向の格子間隔 [m]
- `dz`: numpy.ndarray - z方向の格子間隔配列 [m]
- `max_grad_temp`: float - 温度勾配の最大許容値 [K/m] (デフォルト: 5000.0)
- `max_grad_flux`: float - 熱流束勾配の最大許容値 [W/m³] (デフォルト: 1e8)

**戻り値**:
- `is_valid`: bool - 許容範囲内の場合True
- `max_grad`: float - 最大勾配値

**ロジック**:
- **2次元フィールド（表面熱流束）の場合**:
  ```python
  # x方向勾配
  for i in range(ni-1):
      for j in range(nj):
          grad = abs(field[i+1, j] - field[i, j]) / dy
          if grad > max_grad: max_grad = grad

  # y方向勾配
  for i in range(ni):
      for j in range(nj-1):
          grad = abs(field[i, j+1] - field[i, j]) / dx
          if grad > max_grad: max_grad = grad

  is_valid = max_grad <= max_grad_flux
  ```

- **3次元フィールド（温度場）の場合**:
  ```python
  # x, y, z方向すべての勾配をチェック
  # z方向は非均等格子に対応（dz[k]を使用）
  for i in range(ni):
      for j in range(nj):
          for k in range(nk-1):
              grad = abs(field[i, j, k+1] - field[i, j, k]) / dz[k]
              if grad > max_grad: max_grad = grad

  is_valid = max_grad <= max_grad_temp
  ```

**デコレータ**: `@njit` - Numbaによる高速化

**注意**: コメントに「注：i方向がy物理方向」とあり、インデックスと物理座標の対応に留意

---

### 3.2 包括的異常検出システム

#### 3.2.1 `detect_numerical_anomalies(field, field_name, ...)` [行253-332]

**目的**: 数値計算における異常値を包括的に検出するメインルーチン

**引数**:
- `field`: numpy.ndarray - 検査対象の物理場
- `field_name`: str - フィールド名（温度場、熱流束、勾配等）
- `iteration`: int, optional - CGM反復回数
- `timestep`: int, optional - 時間ステップ
- `temperature_range`: tuple - 温度の物理的範囲 [K] (デフォルト: (150.0, 3000.0))
- `flux_range`: tuple - 熱流束の物理的範囲 [W/m²] (デフォルト: (-1e7, 1e7))
- `dx`, `dy`: float - x, y方向の格子間隔 [m]
- `dz`: numpy.ndarray, optional - z方向の格子間隔配列 [m]

**戻り値**:
- `has_anomaly`: bool - 異常が検出された場合True
- `anomalies`: list - 検出された異常のリスト（文字列）

**検出項目**:

1. **NaN/Inf値チェック**:
   ```python
   if not check_field_finite(field):
       anomalies.append(f"NaN/Inf値検出")
   ```

2. **フィールド種類別の詳細チェック**:
   - **温度場**（"温度" or "Temperature" or "T_" in field_name）:
     - 物理的範囲チェック（150K-3000K）
     - 温度勾配チェック（最大5000K/m）

   - **熱流束場**（"熱流束" or "flux" or "q_" or "q" == field_name）:
     - 熱流束範囲チェック（絶対値≤1e7 W/m²）
     - 熱流束勾配チェック（最大1e8 W/m³）

   - **随伴場・勾配場**（"lambda" or "勾配" or "gradient" in field_name）:
     - 随伴場範囲チェック（絶対値≤1e10）

3. **オーバーフロー検査**:
   ```python
   float64_max = np.finfo(np.float64).max
   max_abs_val = np.max(np.abs(field.ravel()))
   if max_abs_val > float64_max * 0.1:  # float64最大値の10%を超える場合
       anomalies.append(f"数値オーバーフローの危険: max_abs={max_abs_val:.2e}")
   ```

**異常メッセージ例**:
- "異常低温検出: min=120.5K < 150K"
- "異常高温検出: max=3500.2K > 3000K"
- "異常な温度勾配: 6500.1K/m > 5000K/m"
- "異常な熱流束: max_abs=1.50e+08W/m² > 1.00e+07W/m²"
- "数値オーバーフローの危険: max_abs=1.80e+307 > 1.80e+307"

---

## 4. 検証機能の呼び出し箇所

### 4.1 DHCPソルバーでの呼び出し

**場所**: `multiple_time_step_solver_DHCP()` 関数内 [行1156-1219]

**呼び出し内容**:

1. **メモリ監視（開始時）** [行1171-1177]:
   ```python
   if N > 100000:  # 約10万要素以上の場合
       current_memory, is_critical = monitor_memory_usage(
           f"DHCP開始(N={N})",
           warning_threshold_gb=6.0,
           critical_threshold_gb=8.0
       )
       if is_critical:
           raise MemoryError("DHCP: メモリ不足により計算を停止")
   ```

2. **温度場異常検出（各時間ステップ）** [行1211-1216]:
   ```python
   is_valid, status_msg = check_temperature_field(T_all[t], timestep=t, dx=dx, dy=dy, dz=dz)
   if not is_valid:
       print(f"[DHCP警告][t={t}] 温度場異常: {status_msg}")
       # 重大な異常の場合は計算を停止する選択肢もある
       # 現在は警告のみで継続
   ```

3. **Krylov法収束チェック** [行1202-1205]:
   ```python
   x, info = cg(A_csr, b, M=M, rtol=rtol, maxiter=maxiter, x0=x0)
   if info > 0:
       print(f"[警告][t={t}] Krylov法が収束しませんでした (info={info})")
   elif info < 0:
       print(f"[エラー][t={t}] Krylov法の入力が不正です (info={info})")
   ```

---

### 4.2 Adjointソルバーでの呼び出し

**場所**: `multiple_time_step_solver_Adjoint()` 関数内 [行1289-1368]

**呼び出し内容**:

1. **メモリ監視（開始時）** [行1313-1319]:
   ```python
   if N > 100000:
       current_memory, is_critical = monitor_memory_usage(
           f"Adjoint開始(N={N})",
           warning_threshold_gb=6.0,
           critical_threshold_gb=8.0
       )
       if is_critical:
           raise MemoryError("Adjoint: メモリ不足により計算を停止")
   ```

2. **随伴場異常検出（各時間ステップ）** [行1359-1365]:
   ```python
   is_valid, status_msg = check_adjoint_field(lambda_field[t], timestep=t, dx=dx, dy=dy, dz=dz)
   if not is_valid:
       print(f"[Adjoint警告][t={t}] 随伴場異常: {status_msg}")
       # 重大な異常の場合は計算を停止する選択肢もある
       # 現在は警告のみで継続
   ```

3. **Krylov法収束チェック**（DHCPと同様）

---

### 4.3 CGMソルバーでの呼び出し

**場所**: `global_CGM_time()` 関数内 [行1373-1585]

**呼び出し内容**:

1. **メモリ監視（CGMループ開始前）** [行1403-1408]:
   ```python
   current_memory, is_critical = monitor_memory_usage(
       "CGMループ開始",
       warning_threshold_gb=6.0,
       critical_threshold_gb=8.0
   )
   if is_critical:
       raise MemoryError("CGM: メモリ不足により計算を停止")
   ```

2. **定期的なメモリ監視（10反復ごと）** [行1414-1420]:
   ```python
   if it % 10 == 0 or it < 5:  # 最初の5回は毎回、その後は10回ごと
       current_memory, is_critical = monitor_memory_usage(
           f"CGM反復{it}",
           warning_threshold_gb=6.0,
           critical_threshold_gb=8.0
       )
       if is_critical:
           raise MemoryError(f"CGM反復{it}: メモリ不足により計算を停止")
   ```

3. **随伴場の異常検出** [行1463-1468]:
   ```python
   is_lambda_valid, lambda_status = check_adjoint_field(
       lambda_field, iteration=it, dx=dx, dy=dy, dz=dz
   )
   if not is_lambda_valid:
       print(f"[CGM警告][反復{it}] 随伴場全体異常: {lambda_status}")
       # 深刻な場合は計算を停止する選択肢もある
   ```

4. **勾配場の異常検出** [行1475-1480]:
   ```python
   is_grad_valid, grad_status = check_flux_field(grad, iteration=it, dx=dx, dy=dy)
   if not is_grad_valid:
       print(f"[CGM警告][反復{it}] 勾配場異常: {grad_status}")
       # 深刻な場合はステップサイズの調整を検討
   ```

5. **感度場dTの異常検出** [行1506-1511]:
   ```python
   is_dT_valid, dT_status = check_temperature_field(dT, dx=dx, dy=dy, dz=dz)
   if not is_dT_valid:
       print(f"[CGM警告][反復{it}] 感度場dT異常: {dT_status}")
       # 深刻な場合はラインサーチの停止を検討
   ```

6. **熱流束qの異常検出** [行1546-1551]:
   ```python
   is_q_valid, q_status = check_flux_field(q, iteration=it, dx=dx, dy=dy)
   if not is_q_valid:
       print(f"[CGM警告][反復{it}] 熱流束q異常: {q_status}")
       # 深刻な場合は反復の停止または巻き戻しを検討
   ```

---

## 5. Julia実装の設計方針

### 5.1 実装優先順位

**Phase 1（最優先）**:
1. `check_field_finite()` - NaN/Inf検出
2. `check_temperature_range()` - 温度範囲チェック
3. `check_flux_range()` - 熱流束範囲チェック

**Phase 2（高優先）**:
4. `check_gradient_magnitude()` - 勾配急変検出
5. `detect_numerical_anomalies()` - 包括的異常検出

**Phase 3（中優先）**:
6. `check_temperature_field()` - 温度場ラッパー
7. `check_flux_field()` - 熱流束ラッパー
8. `check_adjoint_field()` - 随伴場ラッパー

**Phase 4（低優先、将来拡張）**:
9. メモリ監視機能
10. 数値安定性チェック

### 5.2 モジュール構成

```julia
# julia/src/utils/validators.jl

module Validators

export check_field_finite,
       check_temperature_range,
       check_flux_range,
       check_gradient_magnitude,
       detect_numerical_anomalies,
       check_temperature_field,
       check_flux_field,
       check_adjoint_field

# 実装...

end
```

### 5.3 型システムの活用

```julia
# 検証結果を表す構造体
struct ValidationResult
  is_valid::Bool
  min_val::Float64
  max_val::Float64
  anomalies::Vector{String}
end

# 温度範囲チェック（型安定性確保）
function check_temperature_range(
  T::Array{Float64},
  min_temp::Float64=150.0,
  max_temp::Float64=3000.0
)::Tuple{Bool, Float64, Float64}
  # 実装...
end
```

### 5.4 パフォーマンス最適化

```julia
# @inboundsで境界チェック省略
function check_field_finite(field::Array{Float64})::Bool
  flat_field = vec(field)
  @inbounds for i in 1:length(flat_field)
    if !isfinite(flat_field[i])
      return false
    end
  end
  return true
end

# @simdでSIMD最適化
function check_temperature_range(T::Array{Float64}, min_temp::Float64, max_temp::Float64)
  flat_T = vec(T)
  min_val = flat_T[1]
  max_val = flat_T[1]

  @inbounds @simd for i in 1:length(flat_T)
    val = flat_T[i]
    min_val = min(min_val, val)
    max_val = max(max_val, val)
  end

  is_valid = (min_val >= min_temp) && (max_val <= max_temp)
  return is_valid, min_val, max_val
end
```

### 5.5 テスト方針（TDD準拠）

```julia
# julia/test/test_validators.jl

@testset "Validators Tests" begin
  @testset "check_field_finite" begin
    # 正常値
    @test check_field_finite([1.0, 2.0, 3.0])

    # NaN検出
    @test !check_field_finite([1.0, NaN, 3.0])

    # Inf検出
    @test !check_field_finite([1.0, Inf, 3.0])
  end

  @testset "check_temperature_range" begin
    # 正常範囲
    is_valid, min_val, max_val = check_temperature_range([300.0, 500.0, 800.0], 150.0, 3000.0)
    @test is_valid
    @test min_val ≈ 300.0
    @test max_val ≈ 800.0

    # 異常低温
    is_valid, min_val, max_val = check_temperature_range([100.0, 500.0, 800.0], 150.0, 3000.0)
    @test !is_valid
    @test min_val ≈ 100.0

    # 異常高温
    is_valid, min_val, max_val = check_temperature_range([300.0, 500.0, 3500.0], 150.0, 3000.0)
    @test !is_valid
    @test max_val ≈ 3500.0
  end
end
```

---

## 6. Julia実装時の注意点

### 6.1 インデックスの違い

**Python**: 0-based
```python
for i in range(ni-1):  # 0 ~ ni-2
    grad = abs(field[i+1, j] - field[i, j]) / dy
```

**Julia**: 1-based
```julia
for i in 1:(ni-1)  # 1 ~ ni-1
    grad = abs(field[i+1, j] - field[i, j]) / dy
end
```

### 6.2 配列の平坦化

**Python**: `ravel()`
```python
flat_field = field.ravel()
```

**Julia**: `vec()`
```julia
flat_field = vec(field)
```

### 6.3 文字列フォーマット

**Python**: f-string
```python
print(f"[警告][t={t}] 温度場異常: {status_msg}")
```

**Julia**: 文字列補間
```julia
println("[警告][t=$t] 温度場異常: $status_msg")
```

### 6.4 メモリ監視

**Python**: `psutil`ライブラリ
```python
import psutil
process = psutil.Process()
mem_info = process.memory_info()
```

**Julia**: `Sys`モジュール
```julia
free_gb = Sys.free_memory() / (1024^3)
total_gb = Sys.total_memory() / (1024^3)
```

---

## 7. まとめ

### 7.1 検証システムの全体像

Python原典には以下の包括的な検証システムが実装されている:
- **数値異常検出**: NaN/Inf、範囲チェック、勾配急変（11関数）
- **メモリ監視**: 使用量監視、警告、最適化（7関数）
- **数値安定性**: 拡散数チェック（1関数）

### 7.2 呼び出し箇所

- **DHCPソルバー**: メモリ監視（開始時）、温度場チェック（毎ステップ）
- **Adjointソルバー**: メモリ監視（開始時）、随伴場チェック（毎ステップ）
- **CGMソルバー**: メモリ監視（定期的）、5種類のフィールドチェック（毎反復）

### 7.3 Julia実装方針

1. **Phase 1**: 基本検証関数（NaN/Inf、範囲チェック）
2. **Phase 2**: 高度な検証（勾配検出、包括的異常検出）
3. **Phase 3**: ラッパー関数（ソルバー統合）
4. **Phase 4**: メモリ監視（将来拡張）

TDD方針に従い、**テスト先行→実装→検証**サイクルで進める。
