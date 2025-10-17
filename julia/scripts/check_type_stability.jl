#!/usr/bin/env julia

"""
Type stability sanity checks for the IHCP CGM solvers.

Run with:
  julia --project=julia julia/scripts/check_type_stability.jl
"""

using IHCP_CGM
using InteractiveUtils

const Commons = IHCP_CGM.Commons
const BC = IHCP_CGM.BoundaryConditions
const DHCP = IHCP_CGM.DHCPSolver
const Sens = IHCP_CGM.SensitivitySolver
const Adj = IHCP_CGM.AdjointSolver
const CGM = IHCP_CGM.CGMSolver

struct SampleProblem
  T_init::Array{Float64,3}
  q_surface::Array{Float64,3}
  Y_obs::Array{Float64,3}
  T_series::Array{Float64,4}
  q_init::Array{Float64,3}
  dx::Float64
  dy::Float64
  dt::Float64
  z_centers::Vector{Float64}
  dz::Vector{Float64}
  ρ::Float64
  cp_coeffs::Vector{Float64}
  k_coeffs::Vector{Float64}
  ni::Int
  nj::Int
  nk::Int
  nt::Int
end

function build_sample_problem()
  ni, nj, nk, nt = 3, 3, 3, 4
  T_init = fill(300.0, ni, nj, nk)
  q_surface = fill(5.0, ni, nj, nt - 1)
  Y_obs = fill(295.0, ni, nj, nt)
  T_series = fill(300.0, ni, nj, nk, nt)
  q_init = fill(1.0, ni, nj, nt - 1)
  dx = 0.05
  dy = 0.05
  dt = 0.01
  z_centers = collect(range(0.0, step=0.02, length=nk + 2))
  dz = fill(0.02, nk + 2)
  ρ = 1.0
  cp_coeffs = Float64[0.0, 0.0, 0.0, 1.0]
  k_coeffs = Float64[0.0, 0.0, 0.0, 1.0]
  return SampleProblem(T_init, q_surface, Y_obs, T_series, q_init, dx, dy, dt, z_centers, dz, ρ,
                       cp_coeffs, k_coeffs, ni, nj, nk, nt)
end

function prepare_workbuffers(prob::SampleProblem)
  wk = Commons.WorkBuffers(prob.ni + 2, prob.nj + 2, prob.nk + 2)
  Commons.reset_work_buffers!(wk)
  IHCP_CGM.set_properties!(prob.T_init, wk.cp, wk.λ, prob.cp_coeffs, prob.k_coeffs)
  return wk
end

function check_pbicgstab(prob::SampleProblem)
  println("== @code_warntype PBiCGSTAB! ==")
  wk = prepare_workbuffers(prob)
  bc_set = DHCP.set_dhcp_bc_parameters()
  BC.apply_boundary_conditions!(wk.θ, wk.λ, wk.cp, wk.mask, bc_set)
  @code_warntype IHCP_CGM.PBiCGSTAB!(
    wk,
    (prob.dx, prob.dy, 1.0),
    prob.dt,
    prob.z_centers,
    prob.dz,
    prob.ρ;
    tol=1e-6,
    maxItr=5,
    smoother=:none,
    par="sequential"
  )
  println()
end

function check_solve_dhcp(prob::SampleProblem)
  println("== @code_warntype solve_dhcp! ==")
  wk = prepare_workbuffers(prob)
  @code_warntype IHCP_CGM.solve_dhcp!(
    prob.T_init,
    prob.q_surface,
    wk,
    prob.nt,
    prob.ρ,
    prob.cp_coeffs,
    prob.k_coeffs,
    prob.dx,
    prob.dy,
    prob.z_centers,
    prob.dz,
    prob.dt;
    rtol=1e-6,
    maxiter=10,
    verbose=false,
    par="sequential"
  )
  println()
end

function check_solve_sensitivity(prob::SampleProblem)
  println("== @code_warntype solve_sensitivity! ==")
  wk = prepare_workbuffers(prob)
  @code_warntype IHCP_CGM.solve_sensitivity!(
    prob.T_init,
    prob.q_surface,
    wk,
    prob.nt,
    prob.ρ,
    prob.cp_coeffs,
    prob.k_coeffs,
    prob.dx,
    prob.dy,
    prob.z_centers,
    prob.dz,
    prob.dt;
    rtol=1e-6,
    maxiter=10,
    verbose=false,
    par="sequential"
  )
  println()
end

function check_solve_adjoint(prob::SampleProblem)
  println("== @code_warntype solve_adjoint_mf! ==")
  wk = prepare_workbuffers(prob)
  @code_warntype IHCP_CGM.solve_adjoint_mf!(
    prob.T_series,
    prob.Y_obs,
    wk,
    prob.nt,
    prob.ρ,
    prob.cp_coeffs,
    prob.k_coeffs,
    prob.dx,
    prob.dy,
    prob.z_centers,
    prob.dz,
    prob.dt;
    rtol=1e-8,
    maxiter=10,
    verbose=false,
    par="sequential"
  )
  println()
end

function check_compute_gradient(prob::SampleProblem)
  println("== @code_warntype compute_gradient! ==")
  wk = prepare_workbuffers(prob)
  @code_warntype IHCP_CGM.compute_gradient!(
    prob.T_series,
    prob.Y_obs,
    wk,
    prob.ρ,
    prob.cp_coeffs,
    prob.k_coeffs,
    prob.dx,
    prob.dy,
    prob.z_centers,
    prob.dz,
    prob.dt;
    rtol=1e-8,
    maxiter=10,
    verbose=false,
    par="sequential"
  )
  println()
end

function check_solve_cgm(prob::SampleProblem)
  println("== @code_warntype solve_cgm! ==")
  wk = prepare_workbuffers(prob)
  params = (; max_iter=3, min_iter=1, P=2, dire_reset_every=2)
  @code_warntype IHCP_CGM.solve_cgm!(
    prob.T_init,
    prob.Y_obs,
    prob.q_init,
    wk,
    prob.dx,
    prob.dy,
    prob.z_centers,
    prob.dz,
    prob.dt,
    prob.ρ,
    prob.cp_coeffs,
    prob.k_coeffs;
    params=params,
    par="sequential"
  )
  println()
end

function main()
  prob = build_sample_problem()
  check_pbicgstab(prob)
  check_solve_dhcp(prob)
  check_solve_sensitivity(prob)
  check_solve_adjoint(prob)
  check_compute_gradient(prob)
  check_solve_cgm(prob)
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
