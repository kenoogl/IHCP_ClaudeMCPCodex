# Backend設定のグローバル化提案

**作成日**: 2025年10月23日
**目的**: basesizeパラメータを全関数に逐一渡す必要をなくす

---

## 現状の問題点

### 1. パラメータ伝播の複雑さ

```julia
# 現在の実装（get_backendが12箇所で呼ばれる）
function CalcRK!(...)
    backend = get_backend(par)  # ← basesize/strideを渡せない
    ...
end

# basesize/strideを渡すには全関数の引数を変更する必要がある
function CalcRK!(..., basesize::Int=1, stride::Int=1)
    backend = get_backend(par, basesize=basesize, stride=stride)
    ...
end
```

**影響範囲**:
- `CommonSolver.jl`: 12関数（PCG!, CalcRK!, CalcAX!, Fdot1, Fdot2, BiCG1!, BICG2!, Triad!, Preconditioner!, 等）
- `DHCPSolver.jl`: solve_linear_system!
- `AdjointSolver.jl`: solve_linear_system!
- `SensitivitySolver.jl`: solve_linear_system!
- `CGMSolver.jl`: solve_dhcp!, compute_gradient!, solve_sensitivity!

**合計**: 約20箇所の関数シグネチャを変更する必要がある

### 2. メンテナンス性の低下

- パラメータ追加のたびに全関数を修正
- 関数呼び出しが冗長になる
- デフォルト値の管理が分散

---

## 解決策: Task-local storage方式（推奨）

### 設計方針

1. **グローバル設定をTask-local storageで管理**
2. **スレッドセーフ**: 各タスク（≒各スレッド）で独立した設定
3. **後方互換性**: 既存のコードは修正不要

### 実装案

#### 1. `commons.jl`の修正

```julia
"""
Backend設定をTask-local storageに保存
"""
function set_backend_config(; basesize::Int=1, stride::Int=1)
    task_local_storage(:backend_basesize, basesize)
    task_local_storage(:backend_stride, stride)
    return nothing
end

"""
Backend設定をTask-local storageから取得
"""
function get_backend_config()
    basesize = get(task_local_storage(), :backend_basesize, 1)
    stride = get(task_local_storage(), :backend_stride, 1)
    return (basesize=basesize, stride=stride)
end

"""
並列動作バックエンドを返す（Task-local storageから自動取得）
"""
function get_backend(par::String)
    if par == "thread"
        config = get_backend_config()
        return ThreadedEx(basesize=config.basesize, stride=config.stride)
    else
        return SequentialEx()
    end
end
```

#### 2. `run_sliding_window.jl`での使用

```julia
using ..Commons: set_backend_config

# メイン処理の最初で一度だけ設定
function main()
    # コマンドライン引数解析
    basesize = parse(Int, get(args, "basesize", "1"))
    stride = parse(Int, get(args, "stride", "1"))

    # グローバル設定（Task-local storageに保存）
    set_backend_config(basesize=basesize, stride=stride)

    # 以降、全ての get_backend(par) 呼び出しで自動的に設定が適用される
    result = sliding_window_cgm_q_saving(...)

    return result
end
```

#### 3. 既存コードは変更不要

```julia
# CommonSolver.jl - 変更不要！
function CalcRK!(...)
    backend = get_backend(par)  # ← Task-local storageから自動取得
    ...
end
```

---

## メリット

### 1. 実装の簡潔性
- ✅ 既存の関数シグネチャを変更不要
- ✅ パラメータ伝播の複雑さを解消
- ✅ `commons.jl`のみの修正で完結（約20行追加）

### 2. スレッドセーフ性
- ✅ 各タスク（スレッド）で独立した設定
- ✅ 並列実行時も安全

### 3. 柔軟性
- ✅ 実行時に動的に設定変更可能
- ✅ テストごとに異なる設定が可能
- ✅ デフォルト値の一元管理

### 4. 後方互換性
- ✅ 既存のコードは一切変更不要
- ✅ set_backend_config()を呼ばなければデフォルト動作（basesize=1）

---

## 代替案との比較

### 案1: Task-local storage（推奨）⭐

**メリット**:
- スレッドセーフ
- 実装が簡潔
- Juliaの標準機能

**デメリット**:
- Task-local storageの概念理解が必要

### 案2: Ref型グローバル変数

```julia
# commons.jl
const BACKEND_BASESIZE = Ref{Int}(1)
const BACKEND_STRIDE = Ref{Int}(1)

function set_backend_config(; basesize::Int=1, stride::Int=1)
    BACKEND_BASESIZE[] = basesize
    BACKEND_STRIDE[] = stride
end

function get_backend(par::String)
    if par == "thread"
        return ThreadedEx(basesize=BACKEND_BASESIZE[], stride=BACKEND_STRIDE[])
    else
        return SequentialEx()
    end
end
```

**メリット**:
- シンプル
- 高速アクセス

**デメリット**:
- ⚠️ マルチスレッド時に競合の可能性
- ⚠️ グローバル変数の変更は推奨されない

### 案3: 環境変数

```julia
function get_backend(par::String)
    if par == "thread"
        basesize = parse(Int, get(ENV, "JULIA_BACKEND_BASESIZE", "1"))
        stride = parse(Int, get(ENV, "JULIA_BACKEND_STRIDE", "1"))
        return ThreadedEx(basesize=basesize, stride=stride)
    else
        return SequentialEx()
    end
end
```

**メリット**:
- 外部から設定可能

**デメリット**:
- ⚠️ 環境変数のパース毎回発生（遅い）
- ⚠️ 実行時に動的変更できない
- ⚠️ テストごとに変更しにくい

---

## 実装計画

### Phase 1: Task-local storage方式の実装（推奨）

#### Step 1: commons.jlの修正（10分）
- `set_backend_config()` 追加
- `get_backend_config()` 追加
- `get_backend()` をTask-local storage対応に修正

#### Step 2: run_sliding_window.jlの修正（5分）
- `set_backend_config()` 呼び出しを追加
- コマンドライン引数からbasesize/stride取得

#### Step 3: テスト実行（5分）
- デフォルト動作確認（basesize=1）
- basesize指定動作確認（basesize=10000）

#### Step 4: 既存の未完成修正の削除（5分）
- `CommonSolver.jl`のPCG!/PBiCGSTAB!のbasesize/stride引数削除
- `CGMSolver.jl`のparams取得削除

**合計**: 約25分

### Phase 2: ドキュメント更新（10分）
- `TODO_NEXT_SESSION.md` 更新
- 実装方針ドキュメント作成

---

## コード例

### commons.jl（追加部分のみ）

```julia
export set_backend_config, get_backend_config

"""
Backend設定をTask-local storageに保存

# 引数
- `basesize::Int`: チャンクサイズ（デフォルト: 1）
- `stride::Int`: メモリアクセスストライド（デフォルト: 1）

# 使用例
```julia
# メイン処理の開始時に一度だけ設定
set_backend_config(basesize=10000, stride=1)

# 以降、全ての get_backend(par) で自動的に適用される
backend = get_backend("thread")  # basesize=10000が適用される
```
"""
function set_backend_config(; basesize::Int=1, stride::Int=1)
    task_local_storage(:backend_basesize, basesize)
    task_local_storage(:backend_stride, stride)
    return nothing
end

"""
Backend設定をTask-local storageから取得

設定が存在しない場合はデフォルト値（basesize=1, stride=1）を返す。
"""
function get_backend_config()
    basesize = get(task_local_storage(), :backend_basesize, 1)
    stride = get(task_local_storage(), :backend_stride, 1)
    return (basesize=basesize, stride=stride)
end

"""
並列動作バックエンドを返す

Task-local storageから自動的に設定を取得する。
set_backend_config()で事前に設定されていない場合はデフォルト値を使用。
"""
function get_backend(par::String)
    if par == "thread"
        config = get_backend_config()
        return ThreadedEx(basesize=config.basesize, stride=config.stride)
    else
        return SequentialEx()
    end
end
```

### run_sliding_window.jl（修正例）

```julia
using ..Commons: set_backend_config

function main()
    # コマンドライン引数解析
    args = parse_commandline()
    basesize = args["basesize"]
    stride = args["stride"]

    # グローバル設定（一度だけ）
    set_backend_config(basesize=basesize, stride=stride)

    println("Backend設定: basesize=$basesize, stride=$stride")

    # 以降、全てのget_backend()で自動適用
    result = sliding_window_cgm_q_saving(...)

    return result
end
```

---

## 期待される効果

### 1. コード量削減
- 既存の20関数の修正が不要
- 新規追加コード: 約30行（commons.jl）のみ

### 2. メンテナンス性向上
- パラメータ設定が一箇所に集約
- 関数シグネチャが簡潔に保たれる

### 3. 柔軟性向上
- 実行時に設定変更可能
- テストごとに異なる設定が容易

### 4. 性能への影響
- Task-local storageのアクセスは非常に高速（ハッシュテーブル1回ルックアップ）
- ホットループ外での呼び出しなので影響皆無

---

## 結論

**Task-local storage方式を推奨**します。

理由:
1. ✅ 既存コード変更不要（後方互換性）
2. ✅ スレッドセーフ
3. ✅ Juliaの標準機能を活用
4. ✅ 実装が簡潔（30行の追加のみ）
5. ✅ メンテナンス性向上

次のセッションで実装を進めます。
