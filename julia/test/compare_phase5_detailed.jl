# Phase 5詳細比較スクリプト
# Julia実装とPython参照データの完全な数値比較

using JSON

# モジュール読み込み
include("../src/IHCP_CGM.jl")
using .IHCP_CGM

# 参照データ読み込み
function load_reference_data(filename::String)
  filepath = joinpath(@__DIR__, "..", "data", filename)
  return JSON.parsefile(filepath)
end

println("="^70)
println("1D問題の詳細比較")
println("="^70)

# 1D参照データ
ref_data = load_reference_data("phase5_reference_sliding_window_1D.json")

prob = ref_data["problem"]
sw = ref_data["sliding_window"]
input = ref_data["input"]
output = ref_data["output"]

# パラメータ抽出
ni, nj, nk = prob["ni"], prob["nj"], prob["nk"]
nt = prob["nt"]
dx = prob["dx"]
dy = prob["dy"]
dz = Float64.(prob["dz"])
dz_b = Float64.(prob["dz_b"])
dz_t = Float64.(prob["dz_t"])
dt = prob["dt"]
rho = prob["rho"]
cp_coeffs = Float64.(prob["cp_coeffs"])
k_coeffs = Float64.(prob["k_coeffs"])

window_size = sw["window_size"]
overlap = sw["overlap"]
q_init_value = sw["q_init_value"]
cgm_iteration = sw["CGM_iteration"]

# 入力データ
T_init = permutedims(cat([[[input["T_init"][i+1][j+1][k+1]
                              for k in 0:nk-1] for j in 0:nj-1] for i in 0:ni-1]..., dims=3), (3,2,1))
Y_obs = permutedims(cat([[[input["Y_obs"][t+1][i+1][j+1]
                            for j in 0:nj-1] for i in 0:ni-1] for t in 0:nt-1]..., dims=3), (3,2,1))

# Julia実行
q_global_julia, windows_info = IHCP_CGM.solve_sliding_window_cgm(
  Y_obs, T_init, dx, dy, dz, dz_b, dz_t, dt, rho, cp_coeffs, k_coeffs,
  window_size, overlap, q_init_value, cgm_iteration
)

# Python参照データ
q_global_ref = permutedims(cat([[[output["q_global"][t+1][i+1][j+1]
                                   for j in 0:nj-1] for i in 0:ni-1] for t in 0:nt-2]..., dims=3), (3,2,1))

# 詳細比較
println("\n全15ステップの数値比較:")
println("t  |        Julia          |       Python          |      差分       | 相対誤差")
println("-"^70)
for t in 1:15
  j_val = q_global_julia[t, 1, 1]
  p_val = q_global_ref[t, 1, 1]
  diff = j_val - p_val
  rel_err = abs(diff) / (abs(p_val) + 1e-10)
  println("$(t-1)  | $(j_val) | $(p_val) | $(diff) | $(rel_err)")
end

# ウィンドウ情報比較
println("\nウィンドウ情報比較:")
windows_info_ref = output["windows_info"]
for (i, (jw, pw)) in enumerate(zip(windows_info, windows_info_ref))
  println("\nウィンドウ $i:")
  println("  Julia: start=$(jw[:start_idx]), end=$(jw[:end_idx]), J=$(jw[:J_final])")
  println("  Python: start=$(pw["start_idx"]), end=$(pw["end_idx"]), J=$(pw["J_final"])")
  println("  overlap_steps: Julia=$(jw[:overlap_steps]), Python=$(pw["overlap_steps"])")
end

println("\n" * "="^70)
println("2D問題の詳細比較")
println("="^70)

# 2D参照データ
ref_data_2d = load_reference_data("phase5_reference_sliding_window_2D.json")

prob2 = ref_data_2d["problem"]
sw2 = ref_data_2d["sliding_window"]
input2 = ref_data_2d["input"]
output2 = ref_data_2d["output"]

ni2, nj2, nk2 = prob2["ni"], prob2["nj"], prob2["nk"]
nt2 = prob2["nt"]
dx2 = prob2["dx"]
dy2 = prob2["dy"]
dz2 = Float64.(prob2["dz"])
dz_b2 = Float64.(prob2["dz_b"])
dz_t2 = Float64.(prob2["dz_t"])
dt2 = prob2["dt"]
rho2 = prob2["rho"]
cp_coeffs2 = Float64.(prob2["cp_coeffs"])
k_coeffs2 = Float64.(prob2["k_coeffs"])

window_size2 = sw2["window_size"]
overlap2 = sw2["overlap"]
q_init_value2 = sw2["q_init_value"]
cgm_iteration2 = sw2["CGM_iteration"]

# 入力データ
T_init2 = permutedims(cat([[[input2["T_init"][i+1][j+1][k+1]
                               for k in 0:nk2-1] for j in 0:nj2-1] for i in 0:ni2-1]..., dims=3), (3,2,1))
Y_obs2 = permutedims(cat([[[input2["Y_obs"][t+1][i+1][j+1]
                             for j in 0:nj2-1] for i in 0:ni2-1] for t in 0:nt2-1]..., dims=3), (3,2,1))

# Julia実行
q_global_julia2, windows_info2 = IHCP_CGM.solve_sliding_window_cgm(
  Y_obs2, T_init2, dx2, dy2, dz2, dz_b2, dz_t2, dt2, rho2, cp_coeffs2, k_coeffs2,
  window_size2, overlap2, q_init_value2, cgm_iteration2
)

# Python参照データ
q_global_ref2 = permutedims(cat([[[output2["q_global"][t+1][i+1][j+1]
                                    for j in 0:nj2-1] for i in 0:ni2-1] for t in 0:nt2-2]..., dims=3), (3,2,1))

# 詳細比較（中心点）
println("\n全13ステップの数値比較（中心点 [1,1]）:")
println("t  |        Julia          |       Python          |      差分       | 相対誤差")
println("-"^70)
for t in 1:13
  j_val = q_global_julia2[t, 1, 1]
  p_val = q_global_ref2[t, 1, 1]
  diff = j_val - p_val
  rel_err = abs(diff) / (abs(p_val) + 1e-10)
  println("$(t-1)  | $(j_val) | $(p_val) | $(diff) | $(rel_err)")
end

# ウィンドウ情報比較
println("\nウィンドウ情報比較:")
windows_info_ref2 = output2["windows_info"]
for (i, (jw, pw)) in enumerate(zip(windows_info2, windows_info_ref2))
  println("\nウィンドウ $i:")
  println("  Julia: start=$(jw[:start_idx]), end=$(jw[:end_idx]), J=$(jw[:J_final])")
  println("  Python: start=$(pw["start_idx"]), end=$(pw["end_idx"]), J=$(pw["J_final"])")
  println("  overlap_steps: Julia=$(jw[:overlap_steps]), Python=$(pw["overlap_steps"])")
end
