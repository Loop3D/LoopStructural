# Feature spec: Matrix-free directional constraints in `loop_interpolation_core`

**Feature**: `004-matrix-free-directional-constraints` | **Created**: 2026-10-31 | **Status**: Draft (OQ1 to OQ5 answered; not reviewed)

**Input**: The user wants to model folds with the Rust finite difference
interpolator (FDI). This needs a gradient constraint on every grid cell and
a regularisation that is projected onto a chosen vector. The constraints must
be matrix-free, like the rest of the solver.

**Scope decision**: The core constraints are **generic directional
constraints**. They know nothing about folds. The fold builder is in
`LoopStructural` and uses them (FR-012).

<!-- Write in ASD-STE100 Simplified Technical English. WHAT and WHY, not HOW. -->

## Context

- The Python `DiscreteFoldInterpolator`
  (`LoopStructural/interpolators/_discrete_fold_interpolator.py`) adds four
  constraint families at every tetrahedron barycentre. `f` is the scalar field.
  - Fold orientation: `deformed_orientation . grad f = 0`.
  - Fold axis: `fold_axis . grad f = 0`.
  - Fold normalisation: `dgz . grad f = fold_norm` (non-zero target).
  - Fold regularisation: edge jumps of the gradient projected onto `dgz`,
    the deformed orientation and the fold axis, with three separate weights.
- `loop_interpolation_core` (spec 001) already has, on a structured grid:
  - `add_tangent_constraint` (`v . grad f = 0`, one row per call).
  - `add_norm_constraint` (all three gradient components fixed).
  - `add_directional_regularisation` (`(v . grad)^2 f = 0` with one vector per
    node, matrix-free, accumulates over calls). It is exposed in Python.
- Gaps:
  1. No constraint of the form `v . grad f = b` with `b != 0`. Fold
     normalisation needs it.
  2. Point constraints are stored as one struct per row (8 indices, 8
     coefficients, value, weight; about 150 bytes). A fold uses up to three
     rows per cell, so memory grows with the cell count and removes the
     low-memory benefit of spec 001.
  3. No Python path from a fold feature to the Rust interpolator.
- The Rust regularisation uses the squared directional second derivative. The
  Python version uses first-derivative edge jumps. The two are close but not
  identical (see FR-008).

## User stories

1. **As a modeller**, I want to build a folded foliation field on a structured
   grid with the Rust solver, so that fold models run fast and portably.
2. **As a modeller with a large grid**, I want the per-cell fold constraints to
   use little memory, so that I can use the same grid sizes as the plain FDI.
3. **As a developer**, I want one generic primitive for cell-based directional
   constraints, so that other features can use it (for example, fault frames).

## Functional requirements

- **FR-001** The core MUST provide a family of directional cell constraints.
  Core vectors and points MUST use the computational grid basis defined by
  spec 001. A family has one vector `v` per grid node (computational basis),
  an optional target
  `b` per node, a weight `w`, a quadrature rule (FR-014) and an optional cell
  mask. The family is generic: it is a constraint on `v . grad f`.
- **FR-002** In each active cell, and at each quadrature point `q`, the
  residual MUST be `w * sqrt(q_weight) * (v(q) . grad f(q) - b(q)) / |row|`.
  `v(q)` and `b(q)` are the trilinear interpolation of the node values.
  `grad f(q)` is the trilinear gradient from the 8 corner nodes. `|row|` is the
  norm of the 8 row coefficients. This is the same row normalisation as
  `add_norm_constraint`. `q_weight` is the quadrature weight, and the weights
  of one cell sum to 1.
- **FR-003** If `b` is not given, `b = 0` (homogeneous: "the gradient is
  orthogonal to `v`"). If `b` is given, the constraint fixes the gradient
  component along `v`.
- **FR-004** The core MUST NOT store one struct per row. The core MUST store
  only per-node data (`v`: 3 values; `b`: 1 value, only if given) and the
  cell mask. The 8 corner indices, the interpolated vectors and the row
  coefficients MUST be computed again when needed.
- **FR-005** Cell constraints MUST join the normal operator, the right-hand
  side and the diagonal (preconditioner) of the existing CG solve. The
  operator MUST give the same result as an explicit-row reference (see AC-001).
- **FR-006** The core MUST work with every existing solver path that uses the
  normal operator (`solve`, `solve_with`, the `SolveSession` stepping API) and
  with the memory estimate (`estimate_memory`).
- **FR-007** The core MUST allow several families in one problem, each with
  its own weight. Families that are in the same cell loop MUST share one
  gather of the 8 corner values and one scatter (the cost must not grow
  by a full pass for each family).
- **FR-008** Regularisation projected onto a vector MUST be close to the
  Python `minimise_edge_jumps` result, or the difference MUST be justified.
  - For a constant vector `v`, the local grid approximation to Python's
    edge-jump penalty is `sum_a (d_a d_v f)^2 = |H v|^2` (`H` is the
    Hessian). For a varying vector field, the derivative of `v . grad f`
    also contains terms from changes in `v`. The `|H v|^2` term below does
    not include those terms and MUST NOT be described as an exact equivalent
    to the Python penalty.
  - The existing `add_directional_regularisation` penalises
    `(v^T H v)^2`. This is a weaker, one-row term. It does not penalise
    `H v` components that are orthogonal to `v`.
  - The core MUST therefore add a second term, `|H v|^2` (three rows per node,
    built from the six second-derivative stencils with per-node `v`). This is
    the specified local approximation, including when `v` varies. Its weight
    accumulates over calls like the existing term.
  - The existing term stays. The documentation MUST state which term matches
    the constant-vector local approximation and show the AC-008 comparison
    against the Python edge-jump result for the specified varying fold field.
- **FR-009** A node with a zero `v`, or a cell outside the mask, MUST add no
  row. At a quadrature point where interpolation gives `v(q) = 0`, the row
  MUST also be skipped, including its target. This means that a zero direction
  imposes no constraint. The core MUST NOT raise an error for these cases.
- **FR-010** A vector array with the wrong length MUST be rejected with a
  `CoreError`. This differs from point constraints, which drop silently. A
  wrong length is a caller error, not a geometry case.
- **FR-011** The Python binding MUST expose the families with numpy arrays.
  It MUST accept `vectors` (N x 3, one per node), optional `targets` (N),
  `weight`, `quadrature` and an optional cell `mask`. Vector values MUST use
  the core computational grid basis from spec 001. The modelling adapter is
  responsible for sampling and projecting world directions before binding
  submission.
- **FR-012** The fold builder is in `LoopStructural`, not in the core. Its
  modelling adapter MUST sample fold geometry and masks at world positions,
  then project each direction once into the computational grid basis before
  core submission. It evaluates `deformed_orientation`, `fold_axis` and
  `dgz` at the grid nodes
  and sets one family each (orientation and axis: homogeneous;
  normalisation: `v = dgz`, `b = fold_norm`) and the regularisation terms
  with the weights of `fold_weights`. For defaults and sign conventions, the
  reference implementation MUST be
  `packages/loop_interpolation/src/loop_interpolation/_fd_fold_interpolator.py`,
  which defines the `FDFold` engine used by the golden tests. Defaults MUST be
  `fold_orientation=10.0`, `fold_axis_w=10.0`,
  `fold_regularisation=(0.1, 0.01, 0.01)`,
  `fold_normalisation=1.0`, `fold_norm=-1.0`, and
  `dgz_alignment="warn"`. `fold_normalisation` is the family weight; it MUST
  NOT scale the target. `fold_norm` is a scalar dot-product target and MUST
  remain physically unchanged under rigid rotation. It is subject to the
  reference implementation's `dgz_alignment` policy. This choice does not change the
  separate root-level `LoopStructural` implementation or the package's
  `DiscreteFoldInterpolator`. `FoldedFeatureBuilder` MUST be able to select
  the Rust implementation.
- **FR-013** The cell mask MUST replace the Python `step` and `mask_fn`
  options. A subsample mask made from a seed MUST be fixed for that seed.
- **FR-014** Two quadrature rules MUST exist:
  - `centre`: one point at the cell centre (weight 1).
  - `gauss2`: 2 x 2 x 2 Gauss-Legendre points at `0.5 +/- 1/(2 sqrt 3)` on
    each axis (weight 1/8 each). `gauss2` is exact for the integrated squared
    residual of a trilinear field with a constant `v`, because the squared
    gradient components are at most degree 2 on each axis.
  The default MUST be chosen from the AC-003 and AC-006 results.

## Non-functional requirements

- **NFR-001** For the fold configuration in FR-012 (two homogeneous families
  and one family with targets), retained family data MUST be at most 80 bytes
  per grid node: 24 bytes for each vector and 8 bytes for the target. A cell
  mask MUST use at most one byte per cell. At 1,000,000 cells on a
  100 x 100 x 100-cell grid, retained family data and masks MUST use at least
  five times less memory than the explicit-row reference: three rows per cell,
  144 bytes per row (indices, coefficients, target and weight). Compare only
  family data and masks; exclude the grid and solver workspace from both
  measurements. The `centre` rule is the baseline; `gauss2` MUST NOT add
  retained memory.
- **NFR-002** The time of one CG iteration with three families and the
  regularisation MUST be no more than 1.5 times (`centre`) and 3 times
  (`gauss2`) the time of the same grid with only the existing regularisation. Parallel execution MUST follow the existing `par.rs`
  rules.
- **NFR-003** Results MUST NOT depend on the thread count, within the
  tolerance already used for other families.
- **NFR-004** The f32 and f64 scalar paths MUST both work, as for the other
  families.

## Acceptance criteria

- **AC-001** On a small grid (at least 6x6x6) with random `v`, `b` and weights,
  the new operator `A^T A x`, the right-hand side and the diagonal agree with an
  explicit dense row matrix to a relative error of 1e-12 (f64).
- **AC-002** A homogeneous `centre` family with constant `v` gives the same
  solution as the same rows added with `add_tangent_constraint` at the cell
  centres.
- **AC-003** A deterministic synthetic fold test MUST recover the analytic
  gradient direction. Use the unit cube with a uniform 40 x 8 x 8-cell grid
  and the field `f(x,y,z) = -z + 0.1 sin(2 pi x)`. Its fold axis is the
  world-space y direction. At each grid node use the analytic vectors
  `fold_axis = (0,1,0)`, `deformed_orientation = (1,0,0.2 pi cos(2 pi x))`,
  and `dgz = (0,0,1)`, with target `fold_norm = -1`. Use all cells, the
  `gauss2` rule and the FR-012 default weights. Add unit-weight value
  constraints at every grid node on the two z-boundary planes, using the
  analytic field values, to set the field datum. Measure the axial error
  between the computed and analytic gradient directions at cell centres in
  cells that do not touch the x-boundaries. The median error MUST be below
  10 degrees and MUST be reported by the test. The adapter test MUST also
  sample the translated and rigidly rotated world geometry, directions and
  masks, convert each direction once to the computational basis, and confirm
  that the nonzero scalar dot-product target remains unchanged.
- **AC-004** The fold constraints with all three families and the regularisation
  converge with the existing CG settings. No new solver parameter is needed.
- **AC-005** The memory estimate (`estimate_memory`) includes the cell
  families and is within 10 percent of the measured peak on the AC-003 grid.
- **AC-006** NFR-001 and NFR-002 are measured and recorded in
  `measurements.md` for this feature.
- **AC-007** The Python test suite covers the binding, the wrong-length error
  (FR-010) and the zero-vector case (FR-009).
- **AC-008** On the AC-003 fold, compare each `|H v|^2` regularisation family
  with the corresponding Python `minimise_edge_jumps` family, using the same
  analytic direction field, weights, domain and value data. Use the existing
  Python edge-jump implementation on its supported mesh and evaluate both
  resulting gradient fields at the AC-003 cell centres. Report the median
  angle for each family and for the combined regularisation, and report the
  same measures for the existing `(v^T H v)^2` terms. The combined `|H v|^2`
  result SHOULD give a median angle below 5 degrees. If it does not, document
  the measured result and explain the difference, including the
  varying-vector terms that this local approximation omits. This limited
  regularisation comparison MUST NOT be presented as proof of exact
  equivalence to the tetrahedral fold result.

- **AC-009 (golden gate)** The golden suite
  (`packages/loop_interpolation/tests/golden`) MUST pass. This is a hard
  gate for merge.
  - A new engine `FDFoldRust` (Rust FDI with the generic families from FR-012)
    MUST be added to `ENGINES` and to the `fold_sparse` and `refolded_sparse`
    cases (`conftest.py`, `geomodels.py`).
  - It MUST pass the same `BOUNDS` as `FDFold` in `test_golden.py`
    (`fold_sparse`: (0.20, 15.0); `refolded_sparse`: (0.13, 12.0)), the
    relative checks (`angle_median` below 0.5 x `FD` on `fold_sparse`, below
    0.4 x `FD` on `refolded_sparse`), and `test_data_honoured`. The `BOUNDS`
    numbers MUST NOT be raised to make it pass.
  - Unaffected golden files (`golden_files/*.npz`, `*.json`) MUST stay
    unchanged. Under reviewed C-004, analytic evidence may justify correcting
    only affected old results and Python fingerprints. Numeric bounds and
    tolerances stay fixed. `--update-golden` MUST NOT regenerate old results
    because the new native engine fails. Review new baselines and any justified
    old corrections in the diff.
  - Keep unaffected Python fingerprints fixed. Under reviewed C-004,
    analytic evidence may justify correcting only an affected fingerprint.
    Review such a correction explicitly; do not regenerate it because the new
    native engine fails. For `FDFoldRust` the fingerprints are those of the
    generic families (row count, norms, weights), and they are compared with
    the explicit-row reference of AC-001.
  - Solves use the golden settings (CG, `maxiter=20000`, `rtol=1e-10`).
  - The pure Rust tests (`cargo test` in `loop_interpolation_core`, including
    `parity_fixture.rs` and `class_stencil.rs`) and
    `tests/test_native_parity.py` MUST still pass.
- **AC-010** Fix a failing golden check in the implementation (weights,
  quadrature rule, regularisation term). A reviewed C-004 correction may
  change only an affected old baseline or Python fingerprint when analytic
  evidence proves a coordinate defect. Keep unaffected results, numeric bounds
  and tolerances fixed. A failed new native engine MUST NOT trigger old-result
  regeneration. If a bound cannot be met, report the failure and reason to
  the owner before any bound is changed.

## Out of scope

- Tetrahedral and P2 meshes. The Python `DiscreteFoldInterpolator` and the P2
  fold class stay as they are.
- Numerical parity with the tetrahedral fold result. The grid is different.
- Fault frames, inequality constraints with cell families, and the ADMM path.
- Changes to the fold frame or the fold event geometry code.
- Python mesh and fold changes are out of scope for this 004 feature. Separate
  shared Python coordinate repairs are owned by source spec 005.

## Clarifications

- **2026-11-03, source 005 coordinate contract**: Core vectors and points use
  the computational grid basis from spec 001. The modelling adapter samples
  geometry and masks at world positions and projects directions once before
  core submission. Scalar dot-product targets remain unchanged under rigid
  rotation. Source 005 owns shared Python coordinate repairs; this spec has
  no dependency on its completion.
- **2026-11-03, source 005 C-004**: Analytic evidence and explicit review may
  correct only affected old golden results and Python fingerprints. Keep
  unaffected results, bounds and tolerances fixed. Do not regenerate old
  results after a failed new native engine.

## Decisions (from the open questions)

- **OQ1 (answered)** Vectors are per node and are interpolated to the
  quadrature points. This needs no per-cell or per-point storage, and it
  uses the same convention as `add_directional_regularisation`. Quadrature
  rules: FR-014.
- **OQ2 (answered)** Compute `|row|` when applying each quadrature row. Do not
  cache it per node: the interpolated vector and row norm depend on the
  quadrature point, and the memory limit in NFR-001 leaves no cache allowance.
  Measure the cost in AC-006.
- **OQ3 (answered)** The fused regularisation kernel (`build_combined_kernel`)
  needs one uniform weight for each stencil. Per-node vectors cannot use it.
  It is not applicable. The families share their own cell pass (FR-007).
- **OQ4 (answered)** The fold builder is in `LoopStructural` (FR-012). The core
  constraints are generic.
- **OQ5 (answered)** The regularisation MUST be close to Python or justified
  (FR-008, AC-008).
