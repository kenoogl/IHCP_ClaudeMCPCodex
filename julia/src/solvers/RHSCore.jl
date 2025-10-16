"""
RHSCore.jl

calRHS!関数のコア処理（初期化 + 6面一様境界条件）を提供

3つのソルバー（DHCP, Sensitivity, Adjoint）で共通する
RHSベクトル計算のコア部分を一元管理する。

共通化範囲:
- wk.bの初期化（ゼロクリア）
- 6面の一様熱流束境界条件（X±, Y±, Z±）

各ソルバー固有処理:
- 分布境界条件（calRHS!内で個別実装）
- 最終RHS計算（内部熱源項の有無）
"""

module RHSCore

using FLoops

import ..Commons
using ..Commons: WorkBuffers, get_backend

export calRHS_core!

"""
    calRHS_core!(wk, HF, dx, dy, Δt, ΔZ, par)
      -> (SZ, dx1, dy1, z_st, z_ed, ddt)

RHSベクトルのコア計算（3ソルバー共通部分）

初期化とX/Y/Z方向の一様熱流束境界条件を処理する。
分布境界条件と最終RHS計算は呼び出し側で実装する。

# 引数
- `wk::WorkBuffers`: ワーク配列群（出力先: wk.b）
- `HF::Vector{Float64}`: 6面の一様熱流束 [X-, X+, Y-, Y+, Z-, Z+] [W/m²]
- `dx::Float64`: x方向格子幅 [m]
- `dy::Float64`: y方向格子幅 [m]
- `Δt::Float64`: 時間刻み [s]
- `ΔZ::Vector{Float64}`: z方向格子幅配列 [m]
- `par::String`: 並列化バックエンド（"sequential" or "thread"）

# 戻り値
- `(SZ, dx1, dy1, z_st, z_ed, ddt)`: 後続処理で使用する計算済み値
  - `SZ`: wk.bのサイズ (ni+2, nj+2, nk+2)
  - `dx1`: 1.0 / dx
  - `dy1`: 1.0 / dy
  - `z_st`: z方向開始インデックス（常に2）
  - `z_ed`: z方向終了インデックス（常にSZ[3]-1）
  - `ddt`: 1.0 / Δt

# 処理内容
1. wk.bをゼロクリア（全要素を0.0に初期化）
2. 6面の一様熱流束境界条件を適用:
   - X_minus (i=2): wk.b += HF[1] / dx
   - X_plus (i=SZ[1]-1): wk.b -= HF[2] / dx
   - Y_minus (j=2): wk.b += HF[3] / dy
   - Y_plus (j=SZ[2]-1): wk.b -= HF[4] / dy
   - Z_minus (k=z_st): wk.b += HF[5] / ΔZ[z_st]
   - Z_plus (k=z_ed): wk.b -= HF[6] / ΔZ[z_ed]

# 注意
- HF[i] == 0.0 の場合、その面の処理はスキップ（性能最適化）
- ガイドセル（境界外）は処理対象外
- 戻り値は呼び出し元で分布境界条件・最終RHS計算に使用
"""
function calRHS_core!(
  wk::WorkBuffers,
  HF::Vector{Float64},
  dx::Float64,
  dy::Float64,
  Δt::Float64,
  ΔZ::Vector{Float64},
  par::String
)
  backend = get_backend(par)
  SZ = size(wk.b)
  dx1 = 1.0 / dx
  dy1 = 1.0 / dy
  z_st = 2
  z_ed = SZ[3] - 1
  ddt = 1.0 / Δt

  # wk.bをゼロクリア
  @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
    wk.b[i,j,k] = 0.0
  end

  # 6面の一様熱流束境界条件

  # X_minus (i=2面)
  if HF[1] != 0.0
    i = 2
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1
      wk.b[i,j,k] += HF[1] * dx1
    end
  end

  # X_plus (i=SZ[1]-1面)
  if HF[2] != 0.0
    i = SZ[1]-1
    @floop backend for k in z_st:z_ed, j in 2:SZ[2]-1
      wk.b[i,j,k] -= HF[2] * dx1
    end
  end

  # Y_minus (j=2面)
  if HF[3] != 0.0
    j = 2
    @floop backend for k in z_st:z_ed, i in 2:SZ[1]-1
      wk.b[i,j,k] += HF[3] * dy1
    end
  end

  # Y_plus (j=SZ[2]-1面)
  if HF[4] != 0.0
    j = SZ[2]-1
    @floop backend for k in z_st:z_ed, i in 2:SZ[1]-1
      wk.b[i,j,k] -= HF[4] * dy1
    end
  end

  # Z_minus (k=z_st面)
  if HF[5] != 0.0
    k = z_st
    a = HF[5] / ΔZ[k]
    @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
      wk.b[i,j,k] += a
    end
  end

  # Z_plus (k=z_ed面)
  if HF[6] != 0.0
    k = z_ed
    a = HF[6] / ΔZ[k]
    @floop backend for j in 2:SZ[2]-1, i in 2:SZ[1]-1
      wk.b[i,j,k] -= a
    end
  end

  return (SZ, dx1, dy1, z_st, z_ed, ddt)
end

end # module RHSCore
