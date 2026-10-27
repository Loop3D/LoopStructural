# Research: existing matrix-free interpolator work

**Feature**: `001-low-memory-interpolator` | **Date**: 2026-10-27

This note records what exists today. It is input for `spec.md` and for the
plan. Facts come from the files that this note names.

## 1. Where the work is

| Place | Branch | State |
| --- | --- | --- |
| `~/dev/loop-rs/crates/interpolation-fd` | `lachlan/drillhole` (checked out), `main` | **Most complete.** 4548 lines in `interpolator.rs`. Has all three inequality solvers. |
| `~/dev/LoopAI/packages/interpolation-fd` | `feat/interpolation/constraint-api` | **Older port** of loop-rs commit `05a45ae`. 3141 lines. ADMM only. Adds an `Interpolator` trait adapter (`adapter.rs`). |
| `~/dev/LoopStructural/packages/loop_interpolation/.../loopsolver` | `master` | Python ADMM (`admm_method.py`, `admm_solver.py`, `admm_constant_norm.py`). This is the original that the Rust code ports. |

Other worktrees were checked. The loop-rs worktrees `2d-plots` and
`interpolator-settings-page` do not change this crate in a way that matters for
this spec. The LoopAI worktree `feat+mesh+tet-isosurface-remesh` is about
meshes.

**Use loop-rs as the source.** LoopAI has the trait adapter only. Take the
adapter idea from LoopAI. Take the numerical code from loop-rs.

## 2. What loop-rs already has

Source files: `cg.rs` (177 lines), `grid.rs` (396), `stencil.rs` (798),
`interpolator.rs` (4548), `python.rs` (710, PyO3).

- **Operator.** The normal-equations operator `N = A^T W A` is never built.
  `apply_normal_operator_into(x, scratch, out)` computes `out = N x` with
  stencil sweeps. Each constraint family adds into the same `out`.
- **Constraints.** Value, gradient, normal, tangent, interface, value
  inequality, pairwise inequality. Rows are trilinear. They are stored as small
  records, not as a sparse matrix.
- **Regularisation.** Six curvature families (`dxx` to `dyz`), fused into one
  5x5x5 kernel for interior nodes (about 1.7 to 2 times faster). Centred
  second-derivative border stencils ([1, -2, 1] on layers 1 and n-2). Directional (anisotropic) regularisation.
- **Solvers.**
  - Preconditioned CG (Jacobi, or block-Jacobi along x).
  - `solve_admm`: fixed-rho ADMM. The x-update is an inner CG.
  - `solve_active_set_cg`: primal active set, with projected CG.
  - `solve_hybrid`: a fixed ADMM warm-up, then the active-set polish.
- **Design notes in the code.** CG allocates its five work vectors once.
  `Inequalities.md` (same folder) is the earlier design spec for the two
  inequality methods.
- **Tests.** Unit tests in `interpolator.rs`, edge cases in
  `tests/edge_cases.rs`, and Python parity tests in `tests_python/`.
- **Known gaps** (from its `TASKS.md`): no 2D support, no adaptive rho, no
  column scaling. The cross-derivative weights differ from the current Python
  package by a factor of 4 (`dxy`, `dxz`, `dyz`). The Python `rotation_xy` is
  never applied, so the Rust rotation has no Python reference.

## 3. Why this is not yet "low memory in a browser"

These are the gaps between the existing code and the goal. They come from
reading the code and the notes. Nobody has measured them yet.

1. **Memory is not measured.** Nothing in either repo reports bytes per
   degree of freedom (dof). By count, plain CG holds `x`, `r`, `z`, `p`, `Ap`,
   plus `rhs` and the Jacobi diagonal. That is 7 vectors of `f64`, about 56
   bytes per dof, before any scratch for the fused kernel.
2. **Precision is fixed at `f64`.** There is no `f32` path. An `f32` path
   would halve the vector memory.
3. **The active-set polish uses a dense matrix.** The working-set Gram matrix
   `Q_W Q_W^T` is factored with Cholesky (`active_set_cholesky`). Its size is
   the square of the number of active rows. With many active rows this breaks
   the low-memory goal.
4. **ADMM keeps extra vectors.** It keeps `z` and `u`, one entry per
   inequality row each. It also keeps a copy of the right-hand side.
5. **Threading uses Rayon.** In a browser, Rayon needs shared memory,
   `SharedArrayBuffer`, and cross-origin isolation. Many pages cannot meet
   these. A single-thread path must work well.
6. **No WebAssembly build exists.** The crate builds as `cdylib` for PyO3
   only. Nothing is built for `wasm32`.
7. **No progress or cancel hooks.** A browser solve must not freeze the page.
   It needs a worker, progress reports and cancel.
8. **`wasm32` memory limit.** A 32-bit module can address 4 GiB at most. In
   practice browsers refuse growth earlier. The solver must know its budget
   before it allocates.

## 4. Candidate directions (for the plan, not for the spec)

- Keep CG on the normal equations. Memory is O(dof + constraints).
- Reduce vectors: a CG variant that keeps fewer vectors, or store `r` and `z`
  in the same buffer when the preconditioner is Jacobi.
- Offer `f32` vectors with `f64` accumulation for dot products.
- Replace the dense Gram factorisation in the active-set polish with a
  matrix-free solve (for example, CG on the Gram operator) so that memory
  stays linear in the number of active rows.
- Build for `wasm32-unknown-unknown` with `wasm-bindgen`. Run in a Web Worker.
  Make Rayon optional.
- Keep the Python package as a user of the same core, through PyO3, so that
  one core is tested in two places.
