# Phase 1 Type Stability Report

## Scope
- Refactored `CommonSolver` to operate on a generic `T <: AbstractFloat`, introduced `Val`-based smoother dispatch, and eliminated literal `Float64` constants inside iterative kernels.
- Propagated the `T` type parameter to `DHCPSolver`, `AdjointSolver`, `SensitivitySolver`, and `CGMSolver`, ensuring hot paths (e.g. `solve_dhcp!`, `solve_adjoint_mf!`, `solve_cgm!`) accept and preserve element types instead of hard-coding `Float64`.
- Updated auxiliary routines (`tensor_dot`, `compute_step_size`, `WorkBuffers` usage) to avoid unintended Float64 promotion and keep arithmetic in the chosen precision.
- Added an automated type-stability probe script at `julia/scripts/check_type_stability.jl`.

## Key Changes
- **CommonSolver**
  - All numerical kernels now accept `AbstractArray{T,3}` and `AbstractVector{T}` inputs with `where {T <: AbstractFloat}`.
  - `smoother` keyword uses `Symbol` and dispatches through `Val{:none}` / `Val{:gs}` helpers to avoid runtime string comparisons.
  - Internal literals replaced with `zero(T)`, `one(T)`, or conversions to retain the active element type; `FloatMin` adapted via `T(FloatMin)`.
  - `rbsor!`, `resSOR`, and related helpers return values in the active precision, improving downstream stability.

- **Phase Solvers**
  - DHCP, sensitivity, and adjoint solvers expose generic signatures, allocate buffers with `zeros(T, ...)`, and pass `smoother=:gs` to the refactored `PBiCGSTAB!`.
  - CGM workflow uses promoted scalars (`T`) for optimisation parameters, keeps histories (`J_hist`) in the active precision, and avoids Float64 coercion except when formatting output.
  - `compute_gradient!` and supporting utilities respect the propagated type, enabling consistent inference through the full optimisation pipeline.

- **Utility updates**
  - `tensor_dot` now returns the promoted accumulator type instead of forcing `Float64`.
  - `compute_step_size` accepts any real arrays and promotes `eps` to the accumulator type before division.

- **Tooling**
  - New script `julia/scripts/check_type_stability.jl` builds a compact sample problem and emits `@code_warntype` output for the principal solvers (`PBiCGSTAB!`, `solve_dhcp!`, `solve_sensitivity!`, `solve_adjoint_mf!`, `compute_gradient!`, `solve_cgm!`).

## Verification
- Intended usage:
  ```bash
  julia --project=julia julia/scripts/check_type_stability.jl
  ```
- The script prints `@code_warntype` reports; no dynamic dispatch warnings are expected after the refactor.
- Automated execution was not performed in this session because the sandboxed Julia runtime cannot create its package lockfile; run the command above in a full-access environment to confirm.

## Follow-up
- Consider generalising `RHSCore` and `StoppingCriteria` to remove remaining `Float64` hard-coding if alternative precisions are required.
- After running the check script, capture the `@code_warntype` snippets as artefacts for future regressions (optional).

