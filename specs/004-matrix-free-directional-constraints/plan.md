# Implementation plan: Matrix-free directional constraints

**Feature**: `004-matrix-free-directional-constraints`  
**Date**: 2026-10-31  
**Spec**: [spec.md](spec.md)

<!-- A plan says HOW. Write in ASD-STE100 Simplified Technical English. -->

## Open points

| Id | Question | Resolution rule | Tasks that wait |
|---|---|---|---|
| OP1 | Which cell quadrature rule is the default? | Implement and measure both rules. Select the default from AC-003 and AC-006. Do not select it before those results exist. | T020 |

The accepted repository constitution exists at
`.specify/memory/constitution.md`. The 001 interpolator plan and task list are
the local format precedent.

## Summary

Add a streamed cell family to the Rust finite-difference interpolator. Each
family stores one direction at each grid node and may store one target at each
node. The solver builds a row when it applies the normal operator, right-hand
side or diagonal. It does not keep a record for each cell or quadrature point.

Add a second directional regularisation term that penalises `|H v|^2`. Keep
the current `(v^T H v)^2` term available as a separate core API. The
`FDFoldRust` adapter uses the new `|H v|^2` term for its fold regularisation
weights and compares it with the current term and the Python edge-jump result.

Expose cell families to Python. Add an `FDFoldRust` interpolator that passes
fold vectors and masks to the Rust core. Register it with the folded-feature
builder and golden-test harness. Keep existing engines and golden files
unchanged.

## Technical context

- The core stores point constraints as `PointConstraint` values in
  `packages/loop_interpolation_core/src/interpolator.rs`. Each record holds
  eight indices, eight coefficients, a target and a weight.
- `Problem` in `src/problem.rs` validates flat arrays and adds point
  constraints to the interpolator. The PyO3 API in `src/python.rs` exposes
  these methods as flat contiguous arrays.
- `apply_normal_operator_into`, `data_rhs` and `try_diagonal` in
  `src/interpolator.rs` assemble the normal system without a matrix.
- `ProblemCounts` and `estimate` in `src/memory.rs` report the solver memory.
- The Python package already has a native adapter in
  `packages/loop_interpolation/src/loop_interpolation/_native.py`.
  `FDFoldInterpolator` is in `_fd_fold_interpolator.py`.
- The package golden runner is in
  `packages/loop_interpolation/tests/golden/`. `engines.py` currently runs
  `FDFold` through the Python CG solver, and fingerprints its sparse rows.
- `LoopStructural` re-exports the package interpolators. The folded feature
  factory allow-list is in
  `LoopStructural/modelling/core/_model_feature_factory.py`.

## Repository-fit check

| Rule | Result | How the plan agrees |
|---|---|---|
| Keep the existing point and regularisation APIs | Pass | Add a separate cell family and a separate `|H v|^2` API. Do not change the existing point-row behavior or directional term. |
| Keep the normal system matrix-free | Pass | Recompute cell rows during each operator, right-hand-side and diagonal pass. |
| Keep family storage bounded by nodes and cells | Pass | Store only node vectors, optional node targets, family metadata and the cell mask. |
| Preserve Python behavior when the optional Rust core is not installed | Pass | Import the Rust core only when the Rust fold interpolator is used. |
| Keep existing golden baselines fixed | Pass | Add only `FDFoldRust` baselines for `fold_sparse` and `refolded_sparse`. |
| Test first | Pass | Each behavior test is a separate task before the implementation task that makes it pass. |
| Use the repository's build commands | Pass | Run Rust commands with `--manifest-path packages/loop_interpolation_core/Cargo.toml`; run Python tests from the package roots. |

## Design

### D1. Stored family and public API

Add a public `CellQuadrature` enum and an internal
`DirectionalCellFamily` value in `interpolator.rs`. Store families in
`FiniteDifferenceInterpolator.cell_directional_constraints`. Each family
contains:

- `vectors`: three `f64` values per node;
- optional `targets`: one `f64` value per node;
- one positive family weight;
- a quadrature enum with `Centre` and `Gauss2`;
- an optional one-byte-per-cell mask.

The family stores no row coefficients, cell corner indices or per-quadrature
values. Core vectors and points use the computational GRID basis from spec
001. The modelling adapter samples fold geometry and masks at WORLD positions,
then projects each direction once before core submission. A nonzero scalar
dot-product target remains physically unchanged under rigid rotation.
`Problem` exposes
`add_directional_cell_constraint(vectors, targets, weight, quadrature, mask)`.
The Rust API uses flat slices, consistent with the existing `Problem` methods;
`targets` and `mask` are optional slices and `quadrature` is a
`CellQuadrature`. The Python API accepts a contiguous `N x 3` vector array, an
optional `N` target array, the weight, the quadrature name and an optional
cell mask.

Check lengths, finite values and positive weights before storing the family.
Zero vectors are valid. A zero interpolated vector makes that quadrature row a
no-op, including its target.

### D2. Applying cell rows

At each active cell and quadrature point, compute the trilinear vector and
target, the gradient coefficients, and the row norm. Apply the normalized row
with scale `weight * sqrt(quadrature_weight)`. Use the same row for the normal
operator, right-hand side and Jacobi diagonal. Accumulate all families in one
cell traversal so that they share the corner gather and output scatter.

Cell writes can overlap between adjacent cells. Use eight parity colors of
the structured cells for parallel application. Cells in one color do not
share corner nodes. Keep the existing serial behavior when the parallel
feature is disabled. Follow the thresholds and iterator helpers in `par.rs`.

The direct Rust tests build an explicit row matrix for small grids. Compare
its normal operator, right-hand side and diagonal with the streamed family.

### D3. Directional regularisation

Add three residual rows per interior node for `H v`: one row for each
component of the Hessian-vector product. Implement
`FiniteDifferenceInterpolator.add_hessian_vector_regularisation` and expose
it as `Problem.add_hessian_vector_regularisation(directions, weight)` and
through PyO3. Use the six existing second-derivative stencils and the
direction at the row's node. Add the weighted normal contribution and
diagonal to the existing regularisation paths.
The fold adapter samples these directions at world positions and projects
them once into the computational grid basis before core submission.

Keep `add_directional_regularisation` and its `(v^T H v)^2` behavior. Use
`|H v|^2` for the fold adapter's three regularisation directions. Report both
forms in AC-008. Do not describe `|H v|^2` as exactly equivalent to the Python
edge-jump penalty when `v` varies.

### D4. Memory and solve paths

Add family counts to `ProblemCounts`. Estimate the retained bytes for vectors,
targets, masks and family metadata. Include them in the existing constraint
or working estimate and in the solve budget check. Do not add a per-row
allocation to `solve`, `solve_with` or `SolveSession`.

Test ordinary solve and stepped solve paths. Test `f32`, `f64`, and one-thread
versus multi-thread results. The operator implementation is shared by all
paths, so no separate solver algorithm is needed.

### D5. Python and fold integration

Extend `python.rs` with the cell-family method and explicit shape and length
errors. `FDFoldRust` setup samples geometry and masks at world positions and
projects directions once into the core computational grid basis. `_native.py`
only submits those prepared arrays; it does not sample or project again.
Create `FDFoldRustInterpolator` as a native-only finite-difference fold
adapter. Its setup samples fold geometry and masks at world positions, stores
the three cell-family descriptions and three `|H v|^2` directions, and
projects each direction once into the computational grid basis. It converts
`step` and `mask_fn` to a cell mask. If subsampling is requested, use a `seed`
parameter with default 0 to make the mask repeatable. `_native.py` submits the
prepared arrays and mask without another sampling or projection step. The
class uses the fold defaults
and `dgz_alignment` behavior in FR-012.

Register `FDFoldRust` in the package interpolator map and in
`FOLD_INTERPOLATOR_TYPES`. The factory must select it for a folded feature.
The golden runner calls the Rust CG implementation with
`maxiter=20000` and `rtol=1e-10`.

The Rust core is optional. Importing `loop_interpolation` without the core
installed must continue to work. Selecting `FDFoldRust` without the core must
give the existing native-install guidance.

### D6. Golden fingerprints and baselines

Add `FDFoldRust` only to `fold_sparse` and `refolded_sparse`. Extend the
fingerprint function to report the active row count, row norms and weights
for each cell family. Compare the family fingerprint with rows generated by
the AC-001 explicit-row helper. Keep unaffected golden files and fingerprints
fixed. Under source 005 C-004, analytic evidence and explicit review may
correct only affected old results and Python fingerprints. Keep bounds and
tolerances fixed. A failed new native engine does not permit regenerating old
results.

This task creates and reviews only:

- `fold_sparse__FDFoldRust.json` and `.npz`;
- `refolded_sparse__FDFoldRust.json` and `.npz`.

Keep the existing `FDFold` bounds and relative checks. Under C-004, an
affected old baseline or fingerprint can change only with analytic evidence
and explicit review. Such a correction belongs to its owning coordinate-fix
task. Do not automatically regenerate old results after a failed new native
engine.

### D7. Validation

- Rust unit/integration tests cover AC-001 to AC-005, quadrature rules, masks,
  and `f32`/`f64`. Python binding tests cover AC-007.
- AC-003 uses the exact domain, field, grid, directions, data and metric in
  the spec.
- AC-006 measures memory and one-iteration time for both rules and records
  the results in `measurements.md`.
- AC-008 runs the Python edge-jump reference and both Rust regularisation
  forms on the specified fold. Evaluate gradients at the same AC-003 cell
  centres. Check translated and rotated adapter inputs, world masks and the
  unchanged nonzero scalar target.
- AC-009 runs the golden fold cases and the existing Rust and native parity
  tests. Review any C-004 correction against analytic evidence before changing
  an affected old baseline or fingerprint.

## Test commands

Run Rust checks from the repository root:

```text
cargo test --manifest-path packages/loop_interpolation_core/Cargo.toml
cargo fmt --manifest-path packages/loop_interpolation_core/Cargo.toml -- --check
cargo clippy --manifest-path packages/loop_interpolation_core/Cargo.toml --all-targets -- -D warnings
```

Run Python tests from `packages/loop_interpolation`:

```text
python -m pytest tests/test_directional_cell_binding.py tests/test_native_directional_cell.py tests/test_fold_directional_regularisation.py tests/test_native_parity.py -q
python -m pytest tests/golden/test_golden.py -q -k "fold_sparse or refolded_sparse"
```

Run model-builder tests from the repository root with the project test
configuration. Before merge, run the constitution gates `pytest` and
`ruff check LoopStructural packages`. Review C-004 changes only against
analytic evidence. Do not regenerate old golden results after a failed new
native engine.

## Requirement coverage

| Requirement | Planned work |
|---|---|
| FR-001 | D1, D2; T001-T007; T012-T015 (computational basis and adapter conversion) |
| FR-002 to FR-007 | D1, D2; T001-T007 |
| FR-008 | D3; T002, T009, T014 |
| FR-009, FR-010 | D1; T003, T007 |
| FR-011 | D1, D5; T010-T013 (binding arrays use computational basis; T015 owns world sampling and one direction projection) |
| FR-012 | D3, D5; T012-T017 (world geometry and mask sampling, single projection, fold defaults) |
| FR-013, FR-014 | D1, D2; T002, T006, T020 |
| NFR-001, NFR-002 | D4, D7; T004, T008, T019 |
| NFR-003, NFR-004 | D4; T005-T006 |
| AC-001, AC-002 | D2; T002, T006 |
| AC-003, AC-004 | D4, D5, D7; T005-T006, T014-T015 (world-space translated/rotated oracle, mask and unchanged scalar target) |
| AC-005, AC-006 | D4, D7; T004, T008, T019 |
| AC-007 | D1, D5; T003, T010-T011 |
| AC-008 | D3, D7; T002, T009, T014 |
| AC-009, AC-010 | D6, D7; T015, T018, T021 (reviewed C-004 corrections only for affected results; T018 captures new baselines only) |

## Handover

Source spec 005 owns the shared Python coordinate repairs. Coordinate repair
tasks in 005 must coordinate shared factory, folded-feature builder,
interpolation-adapter and cache files with any still-open tasks here. This
spec does not depend on completion of 005.
