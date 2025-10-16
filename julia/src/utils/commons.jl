"""
commons.jl

共通ユーティリティ

"""
module Commons

using FLoops

export initialize_cells!, λf, WorkBuffers, get_backend,
       BoundaryType, ISOTHERMAL, HEAT_FLUX, CONVECTION, FloatMin, reset_work_buffers!,
       myfill!, mycopy!

# ゼロ除算回避用の小さな数値（BiCGSTAB法などで使用）
const FloatMin = 1.0e-37
      
"""
Harmonic mean with mask correction (Heat3ds互換)
@param a left value
@param b right value
@param ma mask for left (1.0=interior, 0.0=boundary)
@param mb mask for right (1.0=interior, 0.0=boundary)
"""
λf(a, b, ma, mb) = 2.0*a*b / (a+b) * (2.0-div(ma+mb,2))


# 境界条件タイプの列挙型
@enum BoundaryType begin
    ISOTHERMAL   # 等温条件 (Dirichlet)
    HEAT_FLUX    # 熱流束条件 (Neumann)
    CONVECTION   # 熱伝達条件 (Robin)
end


"""
並列動作バックエンドを返す
"""
function get_backend(par::String)
    return (par == "thread") ? ThreadedEx() : SequentialEx()
end


"""
 Heat3D 配列
"""
struct WorkBuffers
    θ      ::Array{Float64,3}
    b      ::Array{Float64,3}
    mask   ::Array{Float64,3}
    cp     ::Array{Float64,3}
    λ      ::Array{Float64,3}
    pcg_p  ::Array{Float64,3}
    pcg_p_ ::Array{Float64,3}
    pcg_r  ::Array{Float64,3}
    pcg_r0 ::Array{Float64,3}
    pcg_q  ::Array{Float64,3}
    pcg_s  ::Array{Float64,3}
    pcg_s_ ::Array{Float64,3}
    pcg_t_ ::Array{Float64,3}
    hsrc   ::Array{Float64,3}
    tmp    ::Array{Float64,3}  # Jacobi前処理用ワーク配列
end

"""
 Heat3D 配列確保
"""
function WorkBuffers(mx::Int64, my::Int64, mz::Int64)
  WorkBuffers(
    zeros(Float64, mx, my, mz), # θ
    zeros(Float64, mx, my, mz), # b
     ones(Float64, mx, my, mz), # mask
    zeros(Float64, mx, my, mz), # cp
    zeros(Float64, mx, my, mz), # λ
    zeros(Float64, mx, my, mz), # pcg_p
    zeros(Float64, mx, my, mz), # pcg_p_
    zeros(Float64, mx, my, mz), # pcg_r
    zeros(Float64, mx, my, mz), # pcg_r0
    zeros(Float64, mx, my, mz), # pcg_q
    zeros(Float64, mx, my, mz), # pcg_s
    zeros(Float64, mx, my, mz), # pcg_s_
    zeros(Float64, mx, my, mz), # pcg_t_
    zeros(Float64, mx, my, mz), # hsrc
    zeros(Float64, mx, my, mz)  # tmp
  )
end

"""
 WorkBuffersのリセット（CGM反復間でクリーンな状態を保証）
"""
function reset_work_buffers!(wk::WorkBuffers)
  # 反復ソルバー用配列をゼロクリア
  wk.pcg_p  .= 0.0
  wk.pcg_p_ .= 0.0
  wk.pcg_r  .= 0.0
  wk.pcg_r0 .= 0.0
  wk.pcg_q  .= 0.0
  wk.pcg_s  .= 0.0
  wk.pcg_s_ .= 0.0
  wk.pcg_t_ .= 0.0

  # θとbもリセット（solve内で上書きされるが、念のため）
  wk.θ    .= 0.0
  wk.b    .= 0.0
  wk.hsrc .= 0.0
  wk.tmp  .= 0.0

  # maskは1.0に保持（境界条件で上書きされる）
  wk.mask .= 1.0

  # cp、λは物性値計算で毎回設定されるのでリセット不要
end


"""
並列対応のcopy!関数

3次元配列のコピーを並列実行可能な形で行う。
シングルスレッド時は組み込みのcopyto!を使用し、
マルチスレッド時はFLoopsで並列化する。

# 引数
- `a::AbstractArray{T,3}`: コピー先配列
- `b::AbstractArray{T,3}`: コピー元配列
- `par::String`: 並列化バックエンド（"sequential" or "thread"）
"""
function mycopy!(a::AbstractArray{T,3}, b::AbstractArray{T,3}, par::String)::Nothing where {T <: AbstractFloat}
  if par == "sequential" || Threads.nthreads() == 1
    # シングルスレッド: 高速な組み込み関数を使用
    copyto!(a, b)
  else
    # マルチスレッド: 並列ループ
    backend = get_backend(par)
    SZ = size(a)
    @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
      a[i,j,k] = b[i,j,k]
    end
  end
  return nothing
end

"""
並列対応のfill!関数

3次元配列を指定値で埋める処理を並列実行可能な形で行う。
シングルスレッド時は組み込みのfill!を使用し、
マルチスレッド時はFLoopsで並列化する。

# 引数
- `a::AbstractArray{T,3}`: 埋める配列
- `val::T`: 埋める値
- `par::String`: 並列化バックエンド（"sequential" or "thread"）
"""
function myfill!(a::AbstractArray{T,3}, val::T, par::String)::Nothing where {T <: AbstractFloat}
  if par == "sequential" || Threads.nthreads() == 1
    # シングルスレッド: 高速な組み込み関数を使用
    fill!(a, val)
  else
    # マルチスレッド: 並列ループ
    backend = get_backend(par)
    SZ = size(a)
    @floop backend for k in 1:SZ[3], j in 1:SZ[2], i in 1:SZ[1]
      a[i,j,k] = val
    end
  end
  return nothing
end


end # module GridTransform
