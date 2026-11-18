# Feature spec: Parallel RBF interpolator in `loop_interpolation_core`

**Feature**: `002-rbf-interpolator` | **Created**: 2026-10-28 | **Status**: Draft (OQ2 to OQ4 answered; not reviewed)

**Input**: The user wants a radial basis function (RBF) implicit-function
interpolator. It follows the `loop-rs` / `LoopAI` RBF interpolator. It goes in
the `loop_interpolation_core` package. It is written in Rust and uses `rayon`.
It must work with LoopStructural. It must allow a fast multipole method (FMM)
later.

<!--
A spec says WHAT the user needs and WHY. It does not say HOW.
Write in ASD-STE100 Simplified Technical English.
-->

## Context

- `loop_interpolation_core` today has one solver: a matrix-free finite
  difference interpolator on a regular grid (spec `001-low-memory-interpolator`).
  Spec 001 lists RBF interpolators as **out of scope**. This spec adds them.
  It does not change the finite difference interpolator.
- The RBF code exists in two places outside this repository:
  - `~/dev/loop-rs/crates/interpolation-rbf` (crate `loop_rbf_core`, about
    4 200 lines, commit `bff5134`). **Most complete.** Use it as the source.
  - `~/dev/LoopAI/packages/interpolation-rbf` (branch
    `feat/interpolation/constraint-api`, commit `c367bac`). An older copy. It
    adds `adapter.rs`, which maps the RBF engine to an `Interpolator` trait.
    Take the adapter idea from here, not the numerics.
- What `loop_rbf_core` has today (from its `TASKS.md`):
  - Hermite (dual-kriging) RBF: value, planar (gradient) and tangent
    constraints. Also "increment" constraints (`field(b) - field(a) = d`).
  - Eight kernels: cubic (`r3`), Gaussian, multiquadric, inverse multiquadric,
    Wendland C2, thin-plate spline (`r^4 log r`), Matern C4, linear.
  - Polynomial drift of degree 0, 1 or 2.
  - Anisotropy from the planar normals; greedy point reduction; regression
    smoothing; restricted range; value inequality constraints (QP active-set).
  - `evaluate` and `evaluate_gradient` at any point.
  - **Solve is a dense LU (nalgebra), single thread, O(n^3) time and O(n^2)
    memory.** Its `TASKS.md` lists this as the open performance gap (P2).
    `evaluate` is a plain loop over all sources for each point.
- More sources of ideas (inspiration only; the numerics still come from
  `loop_rbf_core`, FR-023):
  - **Surfe C++** (`~/dev/surfe`, MIT licence, Government of Canada). It is the
    reference for the constraint semantics, the kernels and the modified
    kernel. `loop_rbf_core` already cites it line by line.
  - **GeologyCore.** The user named it. It is not on this machine, so this
    spec does not describe it. The plan must read it first (see OQ8).
  - **`ferreus_rbf`** (crates.io / docs.rs, version 0.3.0, MIT). A Rust RBF
    crate for more than 1 000 000 sources. It uses an FMM (`ferreus_bbfmm`) for
    the operator and a **domain decomposition preconditioner** inside a
    **flexible GMRES (FGMRES)** solver. It claims O(N log N). It uses `rayon`
    and `faer`. It has a global-trend (anisotropy) transform, isosurface
    extraction and model save/load. From its public docs it does not list
    gradient or inequality constraints. This is the model for the iterative
    path and for the later FMM (user stories 3 and 5).
- The Python package `loop_interpolation` has `SurfeRBFInterpolator`. It wraps
  `surfepy` (C++, conda only). It is not usable on a plain `pip install`.
- The goal of this spec is a **native RBF interpolator that is portable
  (pip, no C++ dependency), uses all cores, and has a path to O(n) cost**.
  The goal is **not** to add new RBF theory.
- Words used in this spec:
  - **source**: one row of the RBF system. A value constraint gives 1 source.
    A planar constraint gives 3 sources (one for each axis). A tangent
    constraint gives 1 source.
  - **n**: the number of sources plus the number of drift terms. This is the
    size of the linear system.
  - **drift**: the polynomial part of the RBF (degree 0, 1 or 2).
  - **kernel operator**: the map `w -> M w`, where `M` is the n x n kernel
    matrix of the sources. The solver may apply it without storing `M`.
  - **direct solve**: a solve that factors a stored dense matrix.
  - **iterative solve**: a solve that only uses the kernel operator.
  - **FMM**: fast multipole method. It applies the kernel operator in about
    O(n) time instead of O(n^2).
  - **reference solution**: the field that `loop_rbf_core` (commit `bff5134`)
    gives for the same problem with its dense LU solve.
  - **evaluation grid**: a regular grid of points where the field is sampled
    from the solved RBF (for example, for an isosurface).

## User scenarios and tests

### User story 1: Solve an RBF model from Python (priority: P1)

A LoopStructural user selects the native RBF interpolator in place of Surfe.
They give contact points, orientation data and tangent data. They get a scalar
field that honours the data. They do not need conda or a C++ build.

**Why this priority**: this is the purpose. Without it nothing else has value.

**Independent test**: solve problem P1 (see Success criteria) through the
Python call. Compare with the reference solution.

**Acceptance scenarios**:

1. **Given** value, planar and tangent constraints, **when** the user solves
   with the default kernel (`r3`) and linear drift, **then** the field agrees
   with the reference solution (SC-001).
2. **Given** each of the eight kernels, **when** the user solves a small
   problem, **then** each field agrees with the reference solution (SC-001).
   A kernel that needs a `shape` and gets none gives an error that names the
   kernel (FR-014).
3. **Given** the linear kernel and a planar or tangent constraint, **when** the
   user solves, **then** the call fails with an error. It does not return a
   field (the linear kernel has no derivative).
4. **Given** the same problem with 1 thread and with N threads, **when** both
   solves end, **then** the fields agree within SC-002.
5. **Given** a field that is solved, **when** the user asks for the value and
   the gradient at many points, **then** the solver returns them. The gradient
   agrees with a central finite difference of the value (SC-003).
6. **Given** input with a NaN, an infinite value, an array of the wrong length
   or a zero-length normal, **when** the user solves, **then** the solver
   returns an error that names the constraint (FR-012).
7. **Given** a grid origin with a large coordinate value (for example
   `5 000 000`), **when** the user solves, **then** the field agrees with the
   field of the same problem moved to the origin, within SC-002. (The solver
   centres and scales the points inside, and the caller does not need to.)

### User story 2: Use the RBF interpolator in a LoopStructural model (priority: P1)

A geologist builds a `GeologicalModel` with the native RBF interpolator for one
feature. The rest of the model (faults, folds, unconformities, stratigraphic
column) works as before.

**Why this priority**: "compatible with LoopStructural" is a stated need. An
RBF that cannot be used from a model is not useful.

**Independent test**: build a small model with one conformable feature. Use the
native RBF and the FDI. Compare the horizon positions.

**Acceptance scenarios**:

1. **Given** a model feature, **when** the user selects the native RBF by name,
   **then** the feature is built with the existing constraint methods
   (`add_value_constraints`, `add_norm_constraints`, `add_gradient_constraints`,
   `add_tangent_constraints`) and the existing `evaluate_value` and
   `evaluate_gradient` calls work (FR-016).
2. **Given** the user does not install the native package, **when** the user
   uses LoopStructural, **then** all existing behaviour stays the same
   (FR-018). The existing test suite passes with no change.
3. **Given** a setting that the RBF does not support (for example, a model
   region mask, a pairwise inequality constraint or a fold constraint), **when**
   the user selects the RBF, **then** the call fails with an error that names
   the setting. It does not ignore the setting (FR-017).
4. **Given** a fault that moves the data, **when** the feature is built,
   **then** the faulted data go in as ordinary points. The RBF does not need to
   know about the fault.
5. **Given** a feature that is built, **when** the viewer or `model.evaluate_*`
   asks for values at a block of points, **then** the points are evaluated in
   parallel (FR-008).

### User story 3: Solve many sources on many cores (priority: P1)

A geologist has 20 000 to 100 000 sources. The dense solve of the loop-rs code
does not fit in memory or takes hours. The new solver uses all cores and does
not store the dense matrix.

**Why this priority**: the speed and memory limit of the dense solve is the
reason to write a new version in Rust with `rayon`.

**Independent test**: solve problem P3 with the iterative solver. Check time,
peak memory and agreement with a direct solve on a smaller version.

**Acceptance scenarios**:

1. **Given** n up to 5 000, **when** the user solves, **then** the direct solve
   is the default. Its speed on 8 threads is at least SC-005 times its speed
   on 1 thread.
2. **Given** n above the direct limit (FR-005), **when** the user solves,
   **then** the iterative solve is the default. The solver never stores the
   n x n matrix. Peak memory grows in a straight line with n (SC-006).
3. **Given** the iterative solve, **when** the solve ends, **then** the result
   reports the iteration count, the final relative residual and a converged
   flag. If the iteration limit is reached, the result says "not converged"
   (FR-013).
4. **Given** the iterative solve and the direct solve on the same problem with
   n up to 5 000, **when** both end, **then** the fields agree within SC-002.
5. **Given** a poorly conditioned problem (for example a Gaussian kernel with a
   large shape), **when** the iterative solve cannot converge, **then** the
   solver says so. It does not return a field marked as converged.
6. **Given** a solve is running, **when** the user cancels, **then** the solve
   stops in less than 1 s and frees its memory (FR-011).

### User story 4: Sample the field on a big grid fast (priority: P1)

The geologist wants the field on a 200 x 200 x 200 grid (8 million points) to
draw an isosurface. The loop-rs `evaluate` is a single-thread loop over all
sources for each point.

**Why this priority**: for most models the evaluation cost is larger than the
solve cost. It is also the simplest part to run in parallel.

**Independent test**: evaluate P2 on a 200^3 grid. Compare with the
single-thread result. Measure the speed-up.

**Acceptance scenarios**:

1. **Given** a solved field and an evaluation grid, **when** the user asks for
   the values, **then** the points are split across threads. The result does
   not depend on the thread count (SC-002).
2. **Given** a big grid, **when** the evaluation runs, **then** the solver
   works in blocks of points. Peak memory is the output plus a fixed block
   size, not the output times the number of sources.
3. **Given** the Wendland C2 kernel, **when** the user evaluates, **then** the
   solver skips the sources that are outside the support radius of the point
   (FR-009).
4. **Given** the user asks for the value and the gradient together, **when**
   the evaluation runs, **then** the solver makes one pass over the sources.

### User story 5: A path to the fast multipole method (priority: P2)

Later, a developer adds an FMM for the kernel operator. They do not rewrite the
solver, the constraint code or the Python layer.

**Why this priority**: it is a design need, not a feature for the first
release. If it is wrong now it is costly to fix later.

**Independent test**: add a second implementation of the kernel operator that
is a test double (a slow but different code path). Check that the iterative
solver and the evaluator give the same field with both.

**Acceptance scenarios**:

1. **Given** the iterative solver, **when** a second kernel operator is
   supplied, **then** the solver runs with no change (FR-006).
2. **Given** an FMM operator that is approximate, **when** the solver uses it,
   **then** the solver accepts an operator-error tolerance, and the stop
   criterion accounts for it (FR-007).
3. **Given** the same for evaluation, **when** a second evaluator is supplied,
   **then** the evaluation path uses it with no change.

### User story 6: Refuse a problem that does not fit (priority: P2)

The user gives a problem that is too big for the memory budget. The solver says
so before it starts.

**Why this priority**: the dense solve can ask for tens of GB and crash the
process. This follows the pattern of spec 001, story 4.

**Independent test**: ask for an estimate for increasing n. Compare with the
measured peak.

**Acceptance scenarios**:

1. **Given** a problem and a method, **when** the user asks for a memory
   estimate, **then** the solver returns it without allocating the large
   arrays.
2. **Given** a budget and an estimate above the budget, **when** the user
   solves with the direct method, **then** the solve does not start. The
   message gives the estimate, the budget, and the largest n that fits. If the
   iterative method would fit, the message says so.
3. **Given** an allocation fails, **when** it happens, **then** the solver
   returns an error. It does not abort.

### User story 7: Inequality and smoothing options (priority: P3)

The geologist keeps the options that `loop_rbf_core` already has: value
bounds, restricted range, regression smoothing, greedy point reduction and
anisotropy.

**Why this priority**: they exist and are tested. They matter to users of
Surfe. They are not needed for a first useful release.

**Acceptance scenarios**:

1. **Given** a value inequality constraint that the unconstrained answer
   breaks, **when** the user solves, **then** the bound holds within SC-004.
2. **Given** regression smoothing, restricted range or anisotropy, **when** the
   user solves, **then** the field agrees with the reference solution within
   SC-001.
3. **Given** greedy reduction, **when** the user solves, **then** the active
   source count is lower and the held-out points agree within SC-007.
4. **Given** an option that the iterative solve does not support (see
   Clarifications, OQ3), **when** the user selects both, **then** the call
   fails with an error that names the option (FR-017).

### Edge cases

- Fewer sources than the drift needs (for example, 2 value points with a
  linear drift). (Error that says how many are needed. The solver does not
  return a field.)
- All value points are on one plane or one line, so the drift matrix is
  singular. (Error that names the problem. It does not return NaN.)
- Two value constraints at one point with different values. (Least-squares
  compromise when `ridge_factor` or smoothing is more than zero. Otherwise an
  error. The solver adds a warning that counts the conflicts.)
- Two sources at the same point (duplicates). (The solver merges exact
  duplicates and adds a warning.)
- A constraint weight of zero or less. (Error that names the constraint.)
- A Gaussian, Wendland or Matern kernel with `shape <= 0`. (Error.)
- A planar or tangent constraint with a zero-length vector. (Error that names
  the constraint.)
- Wendland C2 with a support radius so small that no source has a neighbour.
  (Warning. The system is diagonal.)
- A very large coordinate span (for example `1e7`) or a very small span (for
  example `1e-3`). (The solver scales the coordinates inside. The kernel shape
  is given in the user's units.)
- Zero constraints. (Error: no constraints.)
- An evaluation point far outside the data. (The value follows the drift. No
  error. For compact kernels the value is the drift only.)
- Cancel during setup, before the first iteration. (Stops in less than 1 s.)
- Two solves at the same time on one object. (The second is refused.)
- The thread pool has one thread, or the `parallel` feature is off (for
  example in `wasm32`). (The solver runs on one thread and gives the same
  result.)
- n is so large that an index does not fit in `u32`. (The solver uses `usize`
  indexes. It gives an error only at the memory budget.)

## Requirements

### Functional requirements

- **FR-001**: The package must offer an RBF interpolator for the constraints
  that `loop_rbf_core` has: value, planar (full gradient vector), planar from
  strike/dip/polarity, tangent and increment.
- **FR-002**: The package must offer these kernels: cubic, Gaussian,
  multiquadric, inverse multiquadric, Wendland C2, thin-plate spline, Matern
  C4 and linear. The linear kernel must refuse constraints that need a
  derivative. The shape conventions of `loop_rbf_core` stay (they differ from
  Surfe, and the difference is documented).
- **FR-003**: The package must offer a polynomial drift of degree 0, 1 and 2.
- **FR-004**: The package must offer a **direct** solve (dense, stored matrix)
  and an **iterative** solve (kernel operator only, no stored matrix).
  The user must be able to choose, and the default must choose by n.
- **FR-005**: The default must use the direct solve up to a limit on n and the
  iterative solve above it. The limit must come from the memory estimate and
  the budget. Its default must be 5 000.
- **FR-006**: The iterative solver and the evaluator must work only through a
  **kernel operator interface** (apply the kernel matrix to a vector, apply the
  kernel matrix of the evaluation points to a vector). No other part of the
  solver may read the kernel matrix entries. A second implementation (an FMM,
  or a test double) must be possible with no change to the solver code.
- **FR-007**: The kernel operator interface must allow an approximate
  operator. It must report its error tolerance, so that the solver can set its
  stop criterion.
- **FR-008**: All work that is O(n^2) or more must run in parallel with
  `rayon` when the `parallel` feature is on: system assembly, the kernel
  operator, the direct factorisation and the evaluation. The shim in
  `par.rs` must keep the single-thread build (and `wasm32`) working.
- **FR-009**: The evaluation must work in blocks of points. For a kernel with
  compact support it must use a spatial index and skip the sources outside the
  support radius.
- **FR-010**: The result must not depend on the number of threads within SC-002.
  The solver should use a fixed split of the work (a split that depends on the
  problem size and not on the thread count) so that sums add in the same
  order.
- **FR-011**: The solver must report progress (stage, iteration, residual), and
  it must stop on request in less than 1 s. It must use the same progress and
  cancel types as the finite difference session where they fit.
- **FR-012**: Bad input (NaN, infinite, wrong array length, zero-length vector,
  bad weight, bad shape) must give an error that names the constraint. It must
  not give a wrong field.
- **FR-013**: Every result must say if the solve converged, which method ran,
  the iteration count (iterative), the final relative residual, the condition
  warning (if any), the time, and warnings.
- **FR-014**: A kernel that needs a `shape` must give an error that names the
  kernel when the `shape` is missing or not more than zero.
- **FR-015**: The package must give a memory estimate before it allocates the
  large arrays. It must accept a budget, refuse to start when the estimate is
  over, and return an error (not abort) when an allocation fails. The estimate
  must not be lower than the real peak.
- **FR-016**: The package must be usable from LoopStructural through the
  existing interpolator interface (`GeologicalInterpolator`), with the same
  constraint methods and the same `evaluate_value` and `evaluate_gradient`
  methods as the other interpolators. The user must select it by name
  (for example `interpolatortype="RBF"`), in the same way as other
  interpolators.
- **FR-017**: A setting that the RBF does not support must give an error that
  names the setting. The solver must not ignore it.
- **FR-018**: The existing behaviour of LoopStructural and of
  `loop_interpolation` must not change when the user does not select the RBF.
  The package must stay an optional dependency.
- **FR-019**: The package must keep the full precision of the points. It must
  centre and scale the coordinates inside the solver for conditioning, and
  must undo this in the evaluation.
- **FR-020**: The package must offer a Python binding (feature `python`) with
  NumPy arrays. A WebAssembly binding is **not** required for this release, but
  the core code must compile for `wasm32-unknown-unknown` without the
  `parallel` feature.
- **FR-021**: The package must keep the options of `loop_rbf_core`: value
  inequality constraints (QP active set), restricted range, regression
  smoothing, greedy point reduction and anisotropy. Each must be able to say
  that it does not work with the iterative solve (FR-017).
- **FR-022**: The package must not change the finite difference interpolator
  and its tests. The two solvers may share only the `par`, `error`, `report`
  and `memory` modules.
- **FR-023**: The numerical code must come from `loop_rbf_core` with its
  tests. The tests must move with the code, and they must pass with and
  without the `parallel` feature.

### Key entities

- **RbfProblem**: sources, kernel, drift degree, options. Built from
  constraints. It has no grid.
- **KernelOperator**: applies the kernel matrix (and the drift blocks) to a
  vector. Implementations: dense-on-the-fly (blocked, parallel), stored dense
  (direct solve), and later FMM.
- **RbfSolution**: the weights, the drift coefficients, the centring and
  scaling, and the solve report. It evaluates values and gradients.
- **Evaluator**: applies the kernel of a block of points to the weights.
  Same interface idea as KernelOperator.
- **Budget** and **Memory estimate**: as in spec 001.
- **Solve report**: see FR-013.

## Success criteria

<!-- Measurable. No implementation details. -->

Reference problem **P1**: 200 value points, 60 planar points and 40 tangent
points on a synthetic folded surface (n about 560). Reference problem **P2**:
P1 with 2 000 value points and 500 planar points (n about 3 500), evaluated on
a 200 x 200 x 200 grid. Reference problem **P3**: 30 000 value points and
10 000 planar points on the same surface (n about 60 000). Reference machine:
a laptop with 8 cores and 16 GB of RAM.

The numbers below are **proposals**. See Clarifications, OQ1.

- **SC-001**: For P1 and for each kernel, the direct solve agrees with the
  reference solution. The maximum difference is at most 1e-8 of the field
  range at 5 000 random points.
- **SC-002**: Results for different thread counts, and the iterative result
  against the direct result (P1, P2), differ by at most 1e-6 of the field range.
  For the direct solve with a fixed split (FR-010), results for different
  thread counts must be equal to 1e-12 of the field range.
- **SC-003**: The gradient from `evaluate_gradient` agrees with a central
  difference of `evaluate` to 1e-5 relative error for all kernels except
  linear.
- **SC-004**: With value inequality constraints, every bound holds with a
  violation of at most 1e-6 of the field range.
- **SC-005**: The direct solve of P2 on 8 threads is at least 4 times faster
  than on 1 thread. The evaluation of P2 on the 200^3 grid on 8 threads is at
  least 6 times faster than on 1 thread.
- **SC-006**: P3 solves with the iterative method with a peak memory of at
  most 1 GB. Peak memory against n for n = 10 000, 20 000, 40 000 fits a
  straight line with R^2 of at least 0.99.
- **SC-007**: P3 solves in at most 10 minutes on the reference machine. (The
  O(n^2) operator is the limit. FMM is the later fix.)
- **SC-008**: The evaluation of P2 on the 200^3 grid has a peak memory of at
  most 200 MB more than the output array.
- **SC-009**: The single-thread speed of the direct solve and of the
  evaluation is not more than 1.2 times slower than `loop_rbf_core` on P1 and
  P2. (No regression against the source.)
- **SC-010**: The memory estimate is never lower than the measured peak, and is
  at most 25% higher, for 8 values of n from 500 to 5 000 (direct) and for
  3 values of n from 10 000 to 40 000 (iterative).
- **SC-011**: A model of one conformable feature built with the RBF in
  LoopStructural gives horizon positions within 1% of the model extent of the
  `loop_rbf_core` result for the same data.
- **SC-012**: The existing Rust and Python test suites pass with no change.
  The ported `loop_rbf_core` tests pass with and without the `parallel`
  feature. The crate compiles for `wasm32-unknown-unknown` without `parallel`.

## Out of scope

- The fast multipole method itself. This spec only requires the interface for
  it (FR-006, FR-007, user story 5). A later spec adds the FMM.
- Other fast methods that are not the FMM (for example, hierarchical matrices,
  random Fourier features, partition of unity). The kernel operator interface
  allows them, but this release has none.
- GPU solvers.
- A WebAssembly binding and browser tests for the RBF (FR-020 only keeps the
  core compilable).
- Surfe "horizons" (stratigraphic) and "continuous property" and "vector field"
  system types, the modified kernel, and `min_stratigraphic_thickness`. They are
  open items in `loop_rbf_core`.
- Pairwise inequality constraints. `loop_rbf_core` has none, apart from
  `increment`, which is an equality. (LoopAI's adapter rejects them too.)
- A change to how LoopStructural handles faults, folds or unconformities.
- A change to the finite difference interpolator.
- Mesh and isosurface drawing.
- Saving a solved RBF to disk. (A serde feature exists in `loop_rbf_core`. A
  later spec can bring it.)
- Preconditioner research beyond what the iterative solve needs to converge on
  P3.

## Clarifications

<!-- Record each question and its answer here, with the date. -->

- 2026-10-28: The user asked for a spec for an RBF interpolator that follows
  the `loop-rs` / `LoopAI` RBF interpolator, in `loop_interpolation_core`, in
  Rust with `rayon`, compatible with LoopStructural, with a possible FMM later.
  Result of the search: `loop-rs` is the more complete source (see Context).

- 2026-10-28: The user answered: **OQ2** use a Cargo feature (`rbf`) in the
  existing crate. **OQ3** accepted: value inequality, greedy reduction and
  restricted range are direct-solve only. **OQ4** accepted: the first release
  guarantees iterative convergence only for `r3`, `tps` and Wendland C2, and
  reports "not converged" for the others. The user also named Surfe C++,
  GeologyCore and `ferreus_rbf` as sources of ideas (see Context).

### Open questions (answer before the plan)

- **OQ8: GeologyCore.** Where is it (a path or a URL)? What do you want from
  it: the constraint set, the solver, or the LoopStructural link?
- **OQ9: Iterative method.** `ferreus_rbf` uses FGMRES with a domain
  decomposition preconditioner. Proposal: the plan picks FGMRES (handles the
  indefinite drift system) and a local-neighbourhood preconditioner, and the
  later FMM spec adds the coarse level. Do you accept this?
- **OQ10: FMM source.** Proposal for the later FMM spec: study `ferreus_bbfmm`
  (MIT) and decide there whether to depend on it or to write our own. This
  spec adds no dependency on it.
- **OQ1: Numbers.** Are the SC proposals right? Main ones: n up to 60 000 in
  10 minutes (SC-007), 1 GB (SC-006), 8-core reference machine, 1e-8 against
  the reference for the direct solve (SC-001).



- **OQ5: Name in LoopStructural.** Proposal: new `InterpolatorType.RBF`
  (name `"RBF"`), separate from `SURFE`, so Surfe users are not changed. The
  alternative is to make `SURFE` use the native engine when `surfepy` is
  missing. The proposal is safer for FR-018.
- **OQ6: Shape convention.** Keep the `loop_rbf_core` convention (differs from
  Surfe for Gaussian, MQ and IMQ), or switch to the Surfe convention so
  existing Surfe settings carry over? Proposal: keep it, and convert in the
  LoopStructural adapter.
- **OQ7: Which repository has the version that you want?** `loop-rs` has the
  newer RBF code (commit `bff5134`). `LoopAI` has the older RBF code with the
  trait adapter. Proposal: port from `loop-rs`, and skip the LoopAI
  `Interpolator` trait unless you want it in this crate.
