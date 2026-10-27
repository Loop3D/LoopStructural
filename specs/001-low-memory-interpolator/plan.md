# Implementation plan: Low-memory interpolator for the browser

**Feature**: `001-low-memory-interpolator` | **Date**: 2026-10-27 | **Spec**: [spec.md](spec.md)

<!-- A plan says HOW. Write in ASD-STE100 Simplified Technical English. -->

## Open points

| Id | Question | Default if the user does not answer | Tasks that wait |
|---|---|---|---|
| OP1 | Does the memory per dof of the core stay under the target of D4 (about 75 bytes per dof in `f64`)? Tasks T009, T010 and T020 measure it. | Over 95 bytes per dof after T020: stop and change the design. Between 75 and 95: ask the user. | T039 to T045, T052 (all browser tasks) |
| OP2 | Does the `.wasm` file stay under 2 MB after gzip (SC-010)? Task T044 measures it. | Review the dependencies and the size profile. | T045 |
| OP3 | The reference machine of SC-007 and the Firefox and Safari smoke tests need a person. | CI reports the time only. The release check is manual. | none (manual) |

The orchestrator starts a task with `Blocked by` only after the user
confirms the open point.

The active-set polish of `loop-rs` is **not** in this plan (see D5). The spec
does not need it (principle: no feature that the spec does not need).

## Summary

Make a new Rust core, `packages/loop_interpolation_core`. Start from the matrix-free
finite-difference code in `loop-rs` (`crates/interpolation-fd`). Change it so
that it (1) uses a known, small amount of memory, (2) builds for
`wasm32` and runs in a Web Worker, and (3) runs as a state machine that a
caller can step, watch and cancel. Ship the core in three forms from one
source: a Rust crate, a Python wheel, and a WebAssembly package.

`loop_interpolation` gets a new solver choice, `"native"`. It is a drop-in for
the existing `FiniteDifferenceInterpolator`. With this choice, the Python
class does not assemble sparse rows. It passes its stored constraints to the
core. The Python package stays pure Python. The core is an optional extra.

A web page calls the WebAssembly build directly. The page does all coordinate
transformation (FR-020). The core works in grid coordinates only.

## Technical context

- **Packages changed**
  - New: `packages/loop_interpolation_core` (Rust, PyO3, wasm-bindgen).
  - Changed: `packages/loop_interpolation` (new solver choice, new extra).
  - Changed: release and CI files (see "Project structure").
- **Source of the code**: copy `loop-rs/crates/interpolation-fd` at its current
  head on branch `lachlan/drillhole` (commit
  `bff51342da8ebff7ca3c87e07c1d4fbd31279c8b`, checked on 2026-10-27). Copy
  `src/`, `tests/edge_cases.rs` and `tests_python/`. Write the commit hash in
  the package README. The crate uses Rust edition 2024, `pyo3 0.29` and
  `numpy 0.29`. Do not use a git dependency. Do not use the `LoopAI` copy. It is
  older and it has no active-set solver. Take only the idea of a thin adapter
  from it.
- **New dependencies** (all pure Rust)
  - `wasm-bindgen` (browser glue, feature `wasm` only).
  - `pyo3` and `numpy` (feature `python` only, as in `loop-rs`).
  - `rayon` (feature `parallel` only, now off by default).
  - No `serde`, no `serde_json`, no `regex`. They add size (SC-010).
- **Python side**: no new required dependency. New optional extra
  `native = ["loop-interpolation-core"]`.
- **Build tools**: `maturin` (wheel), `wasm-pack` or `wasm-bindgen-cli`
  with `wasm-opt` (browser package), `cargo`.
- **Performance targets**: SC-001 to SC-010 of the spec. "Test strategy"
  says how to measure each one.
- **Platforms**: wheel for Windows, macOS, Linux. Browser package tested in
  Chromium on Linux (CI) first.

## Repo-fit check

This repository has no constitution. These are the rules that come from the
repository's own decisions. Each row says how the plan agrees.

| Rule | Result | How the plan agrees |
|---|---|---|
| A piece that is useful alone is a separate package (repo structure decision) | Pass | The core is useful with no Python (browser). It is its own package. |
| `loop_interpolation` stays pure Python and its behaviour does not change (FR-019) | Pass | The core is an optional extra. Import is lazy. A missing core gives an `ImportError` with an install hint when the user asks for `"native"`. |
| `release-please` paths are scoped (fix of 2026-09-27) | Pass | The new package gets its own path. The root `exclude-paths` list gets the new path. Task T001. |
| `uv sync` must work with no Rust compiler | Pass, with a change | The root `pyproject.toml` has `[tool.uv.workspace] members = ["packages/*"]`. This would make the core a member. Task T002 adds `exclude = ["packages/loop_interpolation_core"]`. |
| Keep CI small | Pass | One new workflow with path filters. The browser job runs only for this package. The big reference problem R1 runs only on manual dispatch. |
| Test first | Pass | Each phase starts with tests that fail. Task order in `tasks.md` follows this. |
| STE documentation | Pass | The package README, the user guide page and these files use STE. |
| No duplicate numerical code | Pass | After the port, one core serves Python, the browser and (later) `loop-rs`. `loop-rs` can depend on the new crate. That change is out of scope. |

## Design

### D1. What the core takes out of the `loop-rs` code

| Item | Decision | Reason |
|---|---|---|
| Grid rotation (`rotation`, `to_local`, `direction_to_local`) | **Remove** | FR-020: the caller transforms. The Python rotation is not applied today, so there is no Python behaviour to match. |
| `PyO3` bindings | Keep, behind feature `python`. Rewrite for the new API. | Needed for story 6. |
| `rayon` | Move behind feature `parallel` | Browser build runs on one thread (FR-006). |
| Dense Cholesky of the working-set Gram matrix | **Replace** (D5) | Memory grows with the square of the active rows (FR-002). |
| `main.rs` demo binary | Remove | Not needed. |
| `serde` optional feature | Remove | Not needed. SC-010. |
| Build output in the `loop-rs` crate folder (`loop_interpolation_core/*.pyd`, `__pycache__`) | Do not copy | Build output. |
| Closure-style solvers (`solve`, `solve_admm`, ...) | **Replace** with session objects (D3) | Needed for progress and cancel in a single-thread browser. |
| Cross-derivative weights | Keep the `loop-rs` meaning (used as given). The Python adapter divides `dxy`, `dxz`, `dyz` by 4 as `StructuredGrid.get_operators` does. | `loop-rs` research found a factor-4 difference. The adapter hides it. |

### D2. Types and data (Rust)

```text
Grid          origin [f64;3], counts [usize;3], step [f64;3]       (no rotation)
Problem       Grid + constraint records + RegularisationWeights + BorderWeights
              + optional directional regularisation
Constraint    fixed-size record. The `loop-rs` records already hold the 8
              (or 16) corner indices and the trilinear coefficients
              (`PointConstraint`: 8 x usize, 8 x f64, value, weight;
              `PairConstraint`: 16 + 16 + weight; the inequality records add
              lower and upper). Keep this layout. It is linear in the number
              of constraints and has no matrix. The memory estimate counts
              these bytes (up to 272 bytes for a pair record).
Precision     enum { F64, F32 }
Settings      precision, tol, max_iter, budget_bytes, ADMM settings
MemoryEstimate    { vectors, constraints, working, total }
SolveSession  state machine, see D3
SolveReport   converged, iterations, residual, max_violation, peak_bytes,
              warnings
CoreError     enum, one case for each error kind of FR-012, FR-015, FR-018
Field         the solved values and the grid; evaluate, evaluate_gradient
```

Exact fields and signatures are in D8.

The vector element type is generic over a small private trait `Scalar` (`f32`
and `f64`). The trait is in the crate. It adds no dependency. Dot products and
norms always add in `f64` (FR-013).

### D3. Solve sessions (progress, cancel, no thread)

A single-thread browser module cannot receive a "cancel" message while a long
loop runs. So the solver must give control back often. The core does this
with a session object, not with a clock or a thread.

```text
let mut s = problem.start(settings)?;     // allocates, checks the budget
loop {
    match s.step(max_iters) {              // runs at most max_iters iterations
        Step::Running(progress) => { /* report; the caller may drop `s` to cancel */ }
        Step::Done(report)      => break,
    }
}
let field = s.finish()?;                   // moves the field out, no copy
```

- `progress` has the stage, the iteration, the residual and the largest
  violation (FR-009).
- To cancel, the caller drops the session. Drop frees all memory (FR-010).
  The time to cancel is the time of one `step`. The web wrapper picks
  `max_iters` so that one step takes about 100 ms.
- `s.best_so_far()` returns the current field marked "not converged" (story 3,
  scenario 4).
- The native one-shot call `problem.solve(settings)` is a loop over `step`.
  Python uses it. One code path serves all callers.
- Each solver (CG, ADMM, active set, hybrid) is a state machine with its
  state in the session. This replaces the closures of `loop-rs`.

### D4. Memory design (FR-001, FR-002, SC-001, SC-002)

Count of full-size items, `n` = dof. Source: reading `loop-rs`
`interpolator.rs` and `cg.rs`. Nobody has measured these. Task T009 measures.

| Item (in `loop-rs`) | Bytes per dof, `f64`, desktop | Plan |
|---|---|---|
| CG: `x`, `r`, `z`, `p`, `ap` | 40 | Cut `z` (`r.z` = `sum(r^2 * inv_diag)`). Keep 4. |
| `rhs` | 8 | Build straight into `r`, then drop. Cut. |
| Jacobi diagonal `diag` | 8 | Keep as `inv_diag`. |
| `scratch`: one `dof`-length border buffer | 8 | Share with `ap` if the border pass allows (T020). Else keep. |
| `scratch`: 6 slots for each shell row (about 17% of nodes at 100^3) | about 8 | Keep. |
| `shell_row_lookup: Vec<usize>` (filled with `usize::MAX`, size `dof`) | 8 (4 on `wasm32`) | Store as `u32`, or compute the position from `(i,j,k)`. Task T010. |
| `shell_row_indices`, `shell_output_indices` | about 3 | Keep. Store as `u32`. |
| Directional weights `[Vec<f64>; 6]` | 48, only if used | Lazy. R1 does not use them. The estimate adds 48 if the problem uses them. |
| ADMM: base `rhs` copy and `rhs_admm` | 16 | Keep one (D5.1). |
| Result `x` | 0 | Moved out. |

Target after the cuts, plain CG: `x`, `r`, `p`, `ap`, `inv_diag` (40) +
border scratch (8, or 0 if shared) + shell scratch (8) + lookup (4) +
indices (3) = about 63 bytes per dof. ADMM adds one vector (8), so about
71. SC-001 allows 100 bytes per dof. The rest is for the constraint records
(0.4 MB for R1, 5.4 MB for R3), the module, and the copy of the field that
goes to JavaScript (8 bytes per dof). The margin is small. So task T009
measures the port **before** any cut, and again after the cuts of T010 and T020. If
the value is over 95 bytes per dof after the cuts, stop (OP1).

Large allocations use `Vec::try_reserve_exact`. A failure gives
`CoreError::OutOfMemory` and not an abort (FR-012, story 4, scenario 3).

`Problem::estimate_memory(settings)` is arithmetic on counts. It allocates
nothing (FR-011). The estimate adds a safety margin. SC-009 checks it.

### D5. Inequality constraints

1. **ADMM only.** Port the `loop-rs` ADMM into a session (task T022). Its
   operator `N + rho Q^T Q` is already matrix-free. It meets SC-004 at 1e-4.
   The state holds `z` and `u` (one entry per inequality row) and one
   full-size vector. Memory is linear in the dof and the rows (FR-002).
2. **Not in this plan: active-set polish and hybrid.** The `loop-rs` active
   set factors a dense Gram matrix of the active rows. Its size is the square
   of the active rows. This breaks FR-002. A matrix-free replacement would
   serve only the 1e-8 tolerance, and the spec puts that tolerance out of
   scope. So the core has no `ActiveSet` and no `Hybrid` method. A later spec
   can add them. Known risks for that later work, from the `loop-rs` notes: an
   inexact projection leaks residual into the pinned nodes; the release check
   must not run on an unconverged residual.
3. **Infeasible problems** (FR-015). ADMM reports a primal residual that does
   not fall. The session stops after a stall test (20 outer iterations
   without a fall of 1 percent) and returns `CoreError::Infeasible` with the
   20 rows that have the largest violation.

### D6. Browser package

- Feature `wasm`. Crate type `cdylib`. Build target `wasm32-unknown-unknown`.
  No threads, no `SharedArrayBuffer`.
- JavaScript API (classes made by `wasm-bindgen`):

  ```text
  new Problem(origin: Float64Array(3), counts: Uint32Array(3), step: Float64Array(3))
  problem.setRegularisation({ dxx, dyy, dzz, dxy, dxz, dyz, dx, dy, dz })
  problem.addValue(points: Float64Array(3n), values: Float64Array(n), weights: Float64Array(n))
  problem.addGradient / addNormal / addTangent / addInterface(...)
  problem.addValueBound(points, lower, upper, weights)
  problem.addPairBound(pointsA, pointsB, lower, upper, weights)
  problem.estimate(settings) -> { total, vectors, ... }
  problem.start(settings) -> Session           // consumes the problem; throws on budget error
  session.step(maxIters) -> { done, stage, iteration, residual, violation }
  session.field() -> Float32Array | Float64Array   // one copy out
  session.evaluate(points: Float64Array(3m)) -> Float64Array
  session.evaluateGradient(points) -> Float64Array(3m)
  session.free()                                   // cancel
  ```

- A small JavaScript worker wrapper (`worker.js`) owns the `Problem`. It
  receives messages, calls `step` in a loop with `setTimeout(0)` between
  steps so that a "cancel" message is handled, and posts progress. The main
  thread never loads the module (story 3, FR-008).
- Large arrays go to the worker as transferable buffers (no copy).
- Size (SC-010): `opt-level = "z"` for the browser profile only, `lto`,
  `panic = "abort"`, `wasm-opt -Oz`. Measure after phase C. If over 2 MB,
  review the dependency list.
- Peak memory in the browser is read from `WebAssembly.Memory.buffer.byteLength`
  after the solve. WebAssembly memory never shrinks, so this is the peak.

### D7. Python drop-in (`loop_interpolation`)

Changes in `packages/loop_interpolation`:

1. `_solver_strategy.resolve_solver_choice` accepts `"native"`. Today it keeps
   the list `["cg", "lsmr", "admm", None]` and logs a warning for any other
   string. Add `"native"` to the list.
2. The solver can be set at three times: `InterpolatorFactory.create_interpolator(solver=...)`
   (before setup), `InterpolatorBuilder.use_solver` (before setup if the user
   calls it first), and `solve_system(solver=...)` or
   `InterpolatorBuilder.solve(solver=...)` (after setup). So
   `FiniteDifferenceInterpolator.setup_interpolator` checks
   `self.solver == "native"` **at setup time**. If true, it stores the
   weights and keywords and does **not** call the `add_*` methods that
   assemble sparse rows. This saves the Python memory. If the solver is set
   later, the rows are already built. The result is still correct, because
   the adapter (item 4) reads the raw data with `get_value_constraints()`,
   `get_gradient_constraints()`, `get_norm_constraints()`,
   `get_tangent_constraints()`, `get_interface_constraints()`,
   `get_inequality_value_constraints()` and
   `get_inequality_pairs_constraints()`, not the rows. In that case
   `solve_system` logs one warning that says the memory saving needs
   `solver="native"` before setup. Task T037 adds the setup branch. Task T030
   adds tests for both orders. `update()` calls
   `setup_interpolator(**self._setup_kwargs)` again, so the branch must be
   safe to run twice. The setup branch must still return the diagnostics
   report (`finalize_setup_diagnostics_report`). T037 checks that the report
   does not need the rows.
3. `DiscreteInterpolator.solve_system` gets one branch for `"native"`. It
   calls a new function `_native.solve(interpolator, tol, solver_kwargs)` and
   sets `self.c` and `self.up_to_date`.
4. New module `_native.py` (the adapter, about 150 lines). It maps data to core
   calls:

   | Python data / method | Core call |
   |---|---|
   | value (`cpw`) | `add_value` |
   | gradient (`gpw`), direction only | `add_gradient` |
   | normal (`npw`), with length | `add_normal` |
   | tangent (`tpw`) | `add_tangent` |
   | interface (`ipw`), star pairing | `add_interface` |
   | value inequality | `add_value_bound` |
   | inequality pairs | `add_pair_bound` |
   | `dxx` ... `dzz`, `dxy` ... (with the factor 4) | regularisation weights |
   | `dx`, `dy`, `dz` | border weights |
   | directional regularisation | `add_directional_regularisation` |

   Names and the meaning of the gradient and normal keys are not the same in
   the two projects (see the `LoopAI` README). The parity tests (below) fix
   the mapping. Do not trust this table before the tests pass.
5. Settings that the core does not support raise `NotImplementedError` with
   the name of the setting (FR-026): a region mask, column scaling,
   `use_regularisation_weight_scale`, constant-norm iterations,
   `interface_pair_mode` other than star, 2D supports, rotated supports.
6. `solver_kwargs` for `"native"`: `tol`, `maxiter`, `precision`
   (`"f64"` or `"f32"`), `budget_bytes`, and the ADMM
   settings (`rho`, `inner_max_iter`, `inner_tol`, `admm_abs_tol`, `admm_rel_tol`).
   An unknown key raises `NotImplementedError` that names the key.
7. The node order of the core must equal the node order of `StructuredGrid`.
   A test fixes this.
8. `latest_solve_timing` gets the core's timing. The interpolator keeps the
   `SolveReport` as `latest_native_report`.

### D8. Exact Rust and Python API, and the fixture format

A sub agent must not need to make a design decision. These items are fixed.

**Rust, crate `loop_interpolation_core`.** Names are final. A field that this
list does not show does not exist.

```text
// error.rs   (Debug, Clone, PartialEq; Display and std::error::Error)
pub enum BoundKind { Value, Pair }
pub struct ViolatedConstraint { pub index: usize, pub kind: BoundKind, pub violation: f64 }
pub enum CoreError {
    InvalidInput { constraint: String, reason: String },  // constraint = "value[12]"
    GridTooSmall { axis: usize, nodes: usize },           // fewer than 3 nodes
    GridTooLarge { nodes: u64, max_nodes: u64 },          // MAX_NODES = isize::MAX as u64 / 8
    NoConstraints,
    OverBudget { estimate_bytes: u64, budget_bytes: u64, max_nodes_per_axis: usize },
    OutOfMemory { requested_bytes: u64 },
    Infeasible { worst: Vec<ViolatedConstraint> },        // at most 20, falling violation
    PointOutsideGrid { index: usize },
    Unsupported { setting: String },
    SolveRunning,                                         // used by the wasm layer
}

// settings.rs
pub enum Precision { F64, F32 }
pub struct Settings {
    pub precision: Precision,        // default F64
    pub tol: f64,                    // relative residual, default 1e-8
    pub max_iter: usize,             // default 5000 (total inner CG iterations)
    pub budget_bytes: Option<u64>,   // default None
    pub rho: f64,                    // ADMM; defaults copied from loop-rs python.rs
    pub inner_max_iter: usize,
    pub inner_tol: f64,
    pub admm_abs_tol: f64,
    pub admm_rel_tol: f64,
}

// report.rs   (Debug, Clone)
pub enum Stage { Setup, Solve, CheckConstraints }
pub struct Progress { pub stage: Stage, pub iteration: usize, pub residual: f64, pub max_violation: f64 }
pub struct SolveReport {
    pub converged: bool, pub iterations: usize, pub residual: f64,
    pub max_violation: f64, pub peak_bytes: u64, pub warnings: Vec<String>,
}
pub enum Step { Running(Progress), Done(SolveReport) }
pub struct MemoryEstimate { pub vectors: u64, pub constraints: u64, pub working: u64, pub total: u64 }

// field.rs
pub enum FieldData { F64(Vec<f64>), F32(Vec<f32>) }
pub struct Field { pub grid: StructuredGrid3D, pub data: FieldData }
impl Field {
    pub fn values_f64(&self) -> Vec<f64>;
    pub fn evaluate(&self, points: &[f64]) -> Result<Vec<f64>, CoreError>;           // 3m in, m out
    pub fn evaluate_gradient(&self, points: &[f64]) -> Result<Vec<f64>, CoreError>;  // 3m in, 3m out
}

// memory.rs
pub struct MemTracker { /* current, peak, budget */ }
impl MemTracker {
    pub fn new(budget: Option<u64>) -> Self;
    pub fn try_zeroed<T: Copy + Default>(&mut self, n: usize) -> Result<Vec<T>, CoreError>;
    pub fn release(&mut self, bytes: u64);
    pub fn current(&self) -> u64;  pub fn peak(&self) -> u64;
}
pub struct ProblemCounts { /* dof, family counts, uses_directional */ }
pub fn estimate(c: &ProblemCounts, s: &Settings) -> MemoryEstimate;
pub fn max_grid_for_budget(c: &ProblemCounts, s: &Settings, budget: u64) -> usize;

// problem.rs   (arrays are flat: points = 3n values; vectors = 3n values)
pub struct Problem { /* FiniteDifferenceInterpolator, counts, warnings */ }
impl Problem {
    pub fn new(origin: [f64; 3], counts: [usize; 3], step: [f64; 3]) -> Result<Problem, CoreError>;
    pub fn set_regularisation(&mut self, w: RegularisationWeights);
    pub fn set_border_regularisation(&mut self, w: BorderRegularisationWeights);
    pub fn add_value(&mut self, points: &[f64], values: &[f64], weights: &[f64]) -> Result<(), CoreError>;
    pub fn add_gradient(&mut self, points: &[f64], vectors: &[f64], weights: &[f64]) -> Result<(), CoreError>; // direction only
    pub fn add_normal(&mut self, points: &[f64], vectors: &[f64], weights: &[f64]) -> Result<(), CoreError>;   // with length
    pub fn add_tangent(&mut self, points: &[f64], vectors: &[f64], weights: &[f64]) -> Result<(), CoreError>;
    pub fn add_interface(&mut self, points: &[f64], ids: &[f64], weight: f64) -> Result<(), CoreError>;
    pub fn add_directional_regularisation(&mut self, directions: &[f64], weight: f64) -> Result<(), CoreError>;
    pub fn add_value_bound(&mut self, points: &[f64], lower: &[f64], upper: &[f64], weights: &[f64]) -> Result<(), CoreError>;
    pub fn add_pair_bound(&mut self, points_a: &[f64], points_b: &[f64], lower: &[f64], upper: &[f64], weights: &[f64]) -> Result<(), CoreError>;
    pub fn warnings(&self) -> &[String];
    pub fn estimate_memory(&self, s: &Settings) -> Result<MemoryEstimate, CoreError>;
    pub fn residual_f64(&self, field: &Field) -> f64;                                  // phase E
    pub fn start(self, s: Settings) -> Result<SolveSession, CoreError>;               // consumes
    pub fn solve(self, s: Settings) -> Result<(Field, SolveReport), CoreError>;       // loop over step
}

// session.rs   (the session owns the Problem, so it has no lifetime)
impl SolveSession {
    pub fn step(&mut self, max_iters: usize) -> Result<Step, CoreError>;
    pub fn best_so_far(&self) -> Result<(Field, SolveReport), CoreError>;  // converged = false until Done
    pub fn finish(self) -> Result<(Field, SolveReport), CoreError>;        // moves the field out
}
```

Notes: a bound of `-inf` or `+inf` means "no bound" and is valid. `SolveReport`
has no time field because `wasm32-unknown-unknown` has no clock; the caller
measures time. A solve that ends by `max_iter` returns `Ok` with
`converged = false`. The ported rows that drop a constraint with weight 0 or a
point outside the grid are wrapped: the wrapper checks first and gives an error
(weight) or a counted warning (outside).

**Python module `loop_interpolation_core`** (PyO3, abi3 for Python 3.9).

```text
class Problem(origin: ndarray(3), counts: ndarray(3), step: ndarray(3))
    set_regularisation(dict)        # keys dxx dyy dzz dxy dxz dyz  (core meaning, see fixture)
    set_border_regularisation(dict) # keys dx dy dz
    add_value / add_gradient / add_normal / add_tangent(points, vectors_or_values, weights)
    add_interface(points, ids, weight)
    add_directional_regularisation(directions, weight)
    add_value_bound(points, lower, upper, weights)
    add_pair_bound(points_a, points_b, lower, upper, weights)
    estimate(settings: dict) -> dict          # vectors constraints working total
    solve(settings: dict) -> (ndarray, dict)  # field, report
```

Settings keys: `precision` (`"f64"` or `"f32"`), `tol`, `max_iter`,
`budget_bytes`, `rho`, `inner_max_iter`, `inner_tol`, `admm_abs_tol`,
`admm_rel_tol`. The report keys are the `SolveReport` field names.

**Fixture file `fixtures/fixture_20.bin`** (little endian, no header other
than the magic). Bytes 0 to 3 are the ASCII text `LIF1`. Then a list of
blocks. Each block is a `u32` length `N` and `N` values of `f64`. The order is
fixed:

```text
origin(3) counts(3) step(3) reg(6: dxx dyy dzz dxy dxz dyz) border(3: dx dy dz) ridge(1)
value_points value_values value_weights
gradient_points gradient_vectors gradient_weights
normal_points normal_vectors normal_weights
tangent_points tangent_vectors tangent_weights
interface_points interface_ids interface_weight(1)
expected_field(N = counts product)
```

All weights are in core meaning (the cross-derivative weights are already
divided by 4). The generator `fixtures/make_fixture.py` writes the file and
documents the mapping from the Python data columns. The expected field is the
reference solution of the spec: the matrix solver with `solver="cg"` and
`tol=1e-10`.

## Test strategy

Order in each phase: write the tests, see them fail, then implement.

- **Unit tests (Rust, `cargo test`)**: all tests of the `loop-rs` crate move
  with the code, so the port starts from a green suite. New tests: session
  `step` gives the same answer as one-shot; dropping a session frees the
  memory; `estimate_memory` for known sizes; `try_reserve` failure path (a
  budget of 1 byte); `f32` against `f64`; no-op when no inequality rows.
- **Memory tests (Rust)**: a counting global allocator in a test binary
  measures the peak. Check `peak <= estimate <= 1.25 * peak` for sizes
  10^3 to 100^3 in CI. Sizes up to 200^3 run on manual dispatch (SC-009). The
  counting test also checks the number of full-size vectors (D4).
- **Parity tests (Python)**: in `packages/loop_interpolation/tests/
  test_native_*.py`. Skip with `pytest.importorskip` when the core is not
  installed. Start from `loop-rs/crates/interpolation-fd/tests_python/`.
  For each constraint family, build the same problem with the matrix solver
  and with `"native"` on a 30^3 grid. Compare the fields (SC-003). Include the
  inequality problems against the Python ADMM (SC-004), the node-order test,
  the unsupported-setting errors, and a test that the native path builds no
  sparse rows.
- **Cross-build fixtures**: one small problem (20^3) is stored as `fixture_20.bin` (D8)
  with its expected field. The Rust test, the Python test and the browser test
  read the same file. This proves the three builds agree (story 1, scenario
  4).
- **Browser tests (Playwright, Chromium)**: load the package in a Worker, run
  the 20^3 fixture, check the field. A second test runs a 60^3 problem and
  checks (a) the main-thread timer delay is below 100 ms (SC-006), (b) cancel
  stops in under 1 s (FR-010), (c) a budget error appears for a budget that
  is too small (story 4). Memory is `WebAssembly.Memory.buffer.byteLength`.
- **Reference problems R1 and R2 (manual dispatch)**: run in Chromium and on
  the desktop. They record peak memory, time, and field differences in the
  job summary (SC-001, 002, 003, 004, 007). The result goes in
  `specs/001-low-memory-interpolator/measurements.md`.
- **Size check**: CI prints the compressed size of the `.wasm` file and fails
  above 2 MB (SC-010).

## Phases

Each phase ends with passing tests and can merge alone.

| Phase | Content | Stories and requirements | Exit check | Tasks |
|---|---|---|---|---|
| A | Port, remove rotation, features, measure memory | FR-001, 002, 003, 005, 018, 020, 023 | Ported tests green. Memory per dof known. | T006 to T011 |
| B | Flat input, memory estimate, sessions, CG and ADMM state machines, family checks | Story 1 (core), 2, 3 (core), 4; FR-004, 009 to 012, 014, 015, 016, 024 | Memory tests. SC-004 at 1e-4. | T012 to T026 |
| D | `python` bindings, adapter, `"native"` choice, parity tests | Story 6; FR-017, 019, 021, 025, 026 | Parity tests pass. Existing tests unchanged. | T027 to T038 |
| C | `wasm` build, worker wrapper, browser tests | Story 1, 3; FR-006, 008; SC-006, SC-010 | Browser tests pass in CI. | T039 to T045 |
| E | `f32` mode | Story 5; FR-013, 022; SC-005, SC-008 | `f32` tests. | T046 to T050 |
| F | Measurement, CI, docs, final checks | SC-001, 002, 007, 008, 009, 011, 012 | `native.yml` green. `measurements.md` filled. | T051 to T056 |

Setup tasks T001 to T005 come first. The plan of the earlier draft had a phase
for the active-set polish. It is removed (D5).

Phase D can start after phase B. It does not need the browser build.

## Project structure

```text
packages/loop_interpolation_core/
  Cargo.toml                 features: parallel, python, wasm
  pyproject.toml             maturin; import name loop_interpolation_core
  README.md                  STE; source commit of loop-rs; contract
  src/
    lib.rs
    par.rs                   rayon on or off (feature parallel)
    grid.rs                  ported, no rotation
    stencil.rs               ported
    interpolator.rs          ported (not split), tables built lazily
    error.rs  settings.rs  report.rs  field.rs
    memory.rs                MemTracker, estimate, max_grid_for_budget
    session.rs               SolveSession
    cg.rs                    CG as state machine
    admm.rs                  ADMM as state machine
    problem.rs               flat input, validation, start, solve
    scalar.rs                f32 / f64 trait (phase E)
    python.rs                feature python
    wasm.rs                  feature wasm
  web/
    worker.js                step loop, messages, cancel
    index.d.ts
    README.md
  tests/
    common/mod.rs  memory.rs  memory_estimate.rs  input_validation.rs
    session_cg.rs  session_admm.rs  scale_inequality.rs  constraint_families.rs
    parity_fixture.rs  precision_f32.rs  reference_problems.rs  edge_cases.rs
  fixtures/
    make_fixture.py  fixture_20.bin
  scripts/
    build_wasm.sh  check_wasm_size.py
  tests_web/                 Playwright tests
packages/loop_interpolation/
  pyproject.toml             new extra: native
  src/loop_interpolation/
    _native.py               adapter (new)
    _solver_strategy.py      accepts "native"
    _finite_difference_interpolator.py   skip row assembly for "native"
    _discrete_interpolator.py            one branch in solve_system
  tests/test_native_*.py     parity tests (new)
.github/workflows/native.yml           new; path filters; wheel + wasm + browser
.release-please-config.json            new package path; root exclude-paths
.release-please-manifest.json          new entry
pyproject.toml (root)                  [tool.uv.workspace] exclude the core
docs/                                  user guide page for solver="native"
specs/001-low-memory-interpolator/     spec, research, plan, tasks, measurements
```

## Risks

| Risk | Effect | Mitigation |
|---|---|---|
| Real memory is above 100 bytes per dof | SC-001 fails | Phase A measures first. D4 lists the vectors to remove. |
| Node order or key meaning differs between Python and the core | Wrong field with no error | Parity tests for each family. Node-order test. |
| Python builds sparse rows before the solver is known | Memory saved only in the browser | T030 and T037 handle both call orders. A test shows that no rows are built when the solver is set first. |
| `wasm` size over 2 MB | SC-010 fails | No `serde`. Size profile. Check in CI. |
| ADMM needs too many iterations for SC-007 | Time target fails | Report the time (T051). The time target can change (assumption 2 of the spec). |
| Browser refuses memory growth | Tab abort | `try_reserve`, budget, estimate with margin. |
| Single thread is too slow for SC-007 | Time target fails | Measure. Fused kernel is already ported. `parallel` feature for desktop only. Targets can change (assumption 2). |
| Maintaining two Rust copies (`loop-rs` and here) | Drift | Record the source commit. Later, `loop-rs` depends on this crate (out of scope). |

## Complexity tracking

No rule is broken. Two items add complexity on purpose.

| Item | Why it is necessary | Simpler option, and why it is not sufficient |
|---|---|---|
| Session state machines for all solvers | A single-thread browser module must give control back to receive "cancel" and to report progress. | A closure with a callback cannot cancel in a browser with no threads and no `SharedArrayBuffer` (FR-006). |
| New package outside the uv workspace | Contributors with no Rust compiler must still run `uv sync`. | A workspace member forces Rust on all contributors. |

## Coverage of the spec

`tasks.md` has the table that maps each FR and SC to its tasks. All FR and
SC items have at least one task. Two items need a person and have no
automatic task: the Firefox and Safari smoke test (story 1, scenario 3) and
the SC-007 release check on the reference machine (OP3).

## Answers

Answers of the user during the build (2026-10-27). They record decisions that the plan did not cover.

- **OP1 (memory per dof):** measured 63.3 bytes per dof for R1 through the session (target 75). The user has not yet confirmed OP1, so the browser tasks (T039 to T045, T052) stay blocked.
- **T013, session derive:** `Clone` is removed from `SolveSession` in `session.rs`, because `Problem` now holds the interpolator.
- **T016, scaling test:** the constraint rows grow with the node count in the scaling test (same density), because a stored row is 144 bytes and the constraint part is not small for a 30^3 grid.
- **T026, parity with the matrix solver (FR-025, SC-003):** the Rust numerics follow the Python solver. The curvature rows are only at nodes that are interior on all three axes. The border rows are centred second derivatives `[1, -2, 1]` on layers 1 and n-2. The gradient strike and dip basis uses the horizontal strike of the Python solver. The one-sided first-derivative border of loop-rs is not used. Reason: a centred second derivative gives a plane no penalty, and the Python results do not change (FR-019).
- **T035, column scaling:** `apply_scaling_matrix` is `True` by default in this repository. The native adapter ignores it when it has the default value, and logs one info line that scaling is skipped. It raises `NotImplementedError` only when the user sets it explicitly to a non-default value.
