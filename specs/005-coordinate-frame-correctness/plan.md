# Implementation plan: Coordinate and fault-frame correctness

**Feature**: `005-coordinate-frame-correctness` | **Date**: 2026-11-04 | **Spec**: [spec.md](spec.md)

<!-- Written in ASD-STE100 Simplified Technical English.
The constitution (.specify/memory/constitution.md) has priority. -->

## Open points

- **OP1** (closed 2026-11-06, see "Answers"): This plan gives the exact design for all
  stages, S1 to S6. The spec-planner limit is 10 code files, and the maintainer selected a
  staged plan on 2026-11-04, so each stage has its own design section: S1 (shared
  transform, reference-constraint boundary and fault frame), S2 (fold, P2, directional
  regularisation and constant-norm cells), S3 (saved frames, copies and transform
  changes), S4 (display, data export, fault-centre geometry and intrusion
  classification), S5 (resize paths, cell centres, non-mutating diagnostics, residual and
  exclusion report) and S6 (auxiliary point transformation and its entry points, R2
  displacement validation, backend and capability matrix, golden review, final gates).
  S6 has its exact design since 2026-11-06 ("Stage S6 design"). Next step: the
  spec-tasks skill makes the S6 tasks. S1 to S4 are merged; S5 is in pull request #402.
- **OP2**: Run time of the mandatory R2 and R5 matrices. R2 needs 4 anchors x
  3 rotations x 2 centre modes x 2 engines = 48 solves (FDI with 20000
  requested cells). R5 needs 3 buffers x 3 rotations x 2 centre modes x 2
  engines = 36 support builds. SC-018 forbids skips of mandatory cells. Owner:
  the maintainer. Decision needed: run the full matrix in the default
  `pytest` run, or put it behind a pytest marker that a named CI job runs.
  Default until an answer: the full matrix runs in the default `pytest` run.
  The S1 test task measures the time and records it in `verification.md`. Blocks:
  nothing in S1 design; it can change the test layout only.

- **OP3**: The constant-norm variants solve with 20 fixed outer iterations.
  This plan does not know if the unrotated R1 control reaches the SC-002
  limits (`1e-8*F`). Owner: the S2 constant-norm test task (it measures
  first), then the maintainer. If the unrotated control fails, the limit is
  a solver-accuracy question, not a coordinate question, and the maintainer
  decides a separate tolerance or an iteration setting. Blocks: only the
  constant-norm test task. Default until an answer: the SC-002 limits.

- **OP4** (closed 2026-11-05, see "Answers"): The policy for a transform change after an interpolator holds data
  (FR-020). S3 uses "reject at use": the interpolator raises `StaleFrameError`
  at the next set, evaluate or save, and changes no state. It cannot stop
  `BoundingBox.set_local_transform` itself, because the box is shared by all
  interpolators of a model and does not know them. Owner: the maintainer.
  Decision needed: accept this policy, or require a hard rejection inside the
  bounding box (an observer on the box) or an automatic rebuild (it must also
  resize the support, and an unstructured mesh cannot be resized). Default
  until an answer: reject at use. Blocks: nothing in the S3 design; a
  different answer changes the lifecycle tasks only.
- **OP5** (closed 2026-11-05, see "Answers"): SC-017 asks for "1 actionable failure" when a reader from before S3
  loads a new-format file. That reader cannot be changed, because the rollback
  target is the S2 merge. It fails with `TypeError: ... unexpected keyword
  argument 'format_version'` and makes no object. Owner: the maintainer.
  Decision needed: accept this error, which names the key, as the documented
  rollback behaviour. Default until an answer: accept it. Blocks: the S3
  rollback check wording only.
- **OP7** (closed 2026-11-05, spec C-017, see "Answers"): The P2 value-only fit of R1
  (6 x 6 x 6 cells, ridge and gradient regularisation, `cg` with the SC-003
  settings) misses the SC-002 limits in every frame, also the unrotated control at
  `ANCHORS[0]`. The value error is `7.8e-8` there (limit `1e-8*F`, F = 7.5, so `7.5e-8`) and `2e-7`
  to `7.6e-7` at 30 degrees. The gradient error is `3.9e-7`
  against the limit `1e-8` (the rule is `1e-8*max(1, |gradient|)`). The figure
  `4e-9`, which earlier text gave for the gradient limit, was wrong. Also the
  value error is `7.8e-8` against the limit `7.5e-8` in the unrotated control, so
  it fails by 4 percent (an earlier figure of `1.8e-8` was a ratio to the scale of
  the query values). At 30 degrees it is `2e-7` to `7.6e-7`. FDI and PLI stay
  below `3.2e-9` (value) and `1.7e-9` (gradient). The cause is
  solver and regularisation accuracy, not a frame defect. FDI and PLI meet SC-002.
  Decision: SC-002 for P2 applies to the exact nodal field; the fit accuracy goes
  to a separate issue.
- **OP6**: The feature `copy()` methods under `LoopStructural/modelling/features`
  were not read. If one of them calls `interpolator.copy()`, the frame now
  stays on the result. Owner: the S3 copy test task (it runs `git grep -n
  'interpolator.*copy('` first). Default until an answer: no change in S3.
  Blocks: nothing.

- **OP8** (closed for the calls on 2026-11-05, see "Answers"): The four calls
  `model.rescale(...)` in `exporters.py` (lines 92, 366, 431 and 480) discard their
  result, so they change nothing. S4 removes them, because the deprecated method
  would warn on every export. The argument `real_coords=False` of the exporters is
  ignored; it cannot give model coordinates. Owner: the maintainer. Decision
  needed: document or deprecate `real_coords`. Default until an answer: leave it.
  Blocks: nothing.
- **OP11** (closed 2026-11-05, spec C-016, see "Answers"): `GeologicalModel.scale` is
  the inverse of `rescale` (world to local) and has the same "TODO move to the
  bounding box" note. Decision: S4 deprecates `scale` together with `rescale`.
- **OP9**: `LoopStructural/utils/regions.py` was not read (10-file limit). The
  abutting regions receive the world fault centre as `point`. The S1 contract
  says that region callbacks take world points. Owner: the S4 test task for the
  fault centre (it reads `regions.py` first). Default until an answer: the
  regions treat `point` and `vector` as world values. If they convert, the task
  reports it. Blocks: nothing.
- **OP10**: In `GeologicalFeature.get_data`, the `value_map` branch selects rows
  with `value_constraints[value_constraints == v, :3]`. This compares the whole
  array, not the value column. It is a pre-existing defect, not a coordinate
  defect. Owner: the maintainer. Decision needed: file an issue, as for the
  C-005 and C-006 defects. Default until an answer: S4 keeps the expression and
  tests only the branch without `value_map`. Blocks: nothing.
- **OP12**: The spec (FR-015) says that cell-centre sampling can place a centre
  beyond its cell. By reading, `BoundingBox.cell_centres` is correct for a
  rotated box. The spacing equals `node_step_vector`, and the four callers
  convert with `reproject`. The planner did not run the code. Owner: the S5
  cell-centre test task. Decision needed: if all 60 cells of the new tests pass,
  record "not reproduced" and keep the tests as guards. Default until an answer:
  no code change in `_bounding_box.py`. If a cell fails, use the closed form of
  "S5 decisions" item 7. Blocks: nothing.
- **OP13** (closed 2026-11-05, maintainer answer Q3, see "Answers"): The misfit
  methods must not change state and must not rebuild a stale model. Rule: on a
  stale builder, fault builder or interpolator, the methods write one
  `UserWarning` and return NaN ("S5 decisions" item 11).
- **OP14** (closed 2026-11-05, maintainer answer Q1, see "Answers"): The 2D
  `maximum` keeps its legacy value `origin + nsteps * step`. It counts nodes as
  cells. S5 keeps it for FR-028. The four existing assertions pass with no change.
  S5 only derives `nsteps_cells`, `n_cell_x` and `n_cell_y` from `nsteps`.
- **OP15** (closed 2026-11-05, maintainer answer Q2, see "Answers"): The old misfit
  methods keep their names and public shape. New names give the 1-D result. The
  gradient value changes meaning. `COMPAT.md` and `ROADMAP.md` record it.
- **OP16**: The old `evaluate_gradient_misfit` for two families (gradient and
  normal) of different length. Before S5 it raised `ValueError` on current NumPy.
  Equal lengths gave shape `(2, N)`. Owner: the maintainer. Decision needed: keep
  `(2, N)` for equal lengths and return one row `(1, N1 + N2)` for different
  lengths, or always return one row. Default until an answer: the first choice,
  because it keeps the shape that worked. `COMPAT.md` records it. Blocks:
  nothing.
- **OP17**: The setter fix and the extent from `set_nelements` can change numeric
  tests, including `packages/loop_interpolation/tests/golden` (the parity oracle
  for the Rust path, with its own `--update-golden` option). Owner: the S5 survey
  task. Decision needed: whether to regenerate such a file. Default until an
  answer: record each change in "S5 golden review (C-004)"; regenerate a file only
  when the cause is the intended extent change, and name the file; leave the
  gates red and report if the cause is not clear. Blocks: the merge of S5 when a
  parity file changes.

- **OP18** (closed 2026-11-06, see "Answers"): The spec (R2) named "the same existing
  test-side displacement evaluator" with trishear parameters: apex angle 50 degrees,
  propagation/slip ratio 1.5, concentration 1, 20 integration steps, tip coordinate 0,
  propagation sign -1, and 16 quadrature points for the scalar comparison. The repository
  has no such evaluator (`git grep -n -i trishear` on 2026-11-06 finds the spec only), and
  the spec scope bans a trishear algorithm in the library. Maintainer answer: amend the
  spec. S6 uses the library displacement paths on both frames, with the same algorithm:
  `FaultSegment.apply_to_points` with 20 steps (vector) and
  `FaultSegment.evaluate_displacement` (scalar), slip 3 and the default `BaseFault`
  profile. The trishear parameters and the 16 quadrature points are not used. The spec has
  the new clarification C-019 and the changed R2 fixture text. The SC-004 numbers do not
  change (1332 points, 0.03).
- **OP19** (closed 2026-11-06, see "Answers"): S5 left one decision open
  (`run-state.json`, pull request #402). The physical residual methods
  `_value_residuals_fresh` and `_gradient_cosines_fresh` of `GeologicalFeature` evaluated
  through the feature, which applies the faults a second time to rows that are already
  restored. For a faulted feature the residual compared two frames (0.080 measured, 0.0004
  to 0.0018 true). Maintainer answer: accept the fix. The two methods evaluate the
  interpolator directly at the stored positions, with no fault step. This reverses S5
  decision 8 for these two private methods. The task text of `COMPAT.md` and the entry in
  `ROADMAP.md` say so. Files: `_geological_feature.py`, `test_diagnostics_backends.py`,
  `test_diagnostics_frames.py`, `COMPAT.md`, `ROADMAP.md` ("S6 decisions" item 6).
- **OP20** (closed 2026-11-06 as "fix", see "Answers"): `EuclideanTransformation` has other
  defects that are not round-trip defects: `fit` fails for `dimensions=3`; the fitted angle
  has the wrong sign for a main axis with a negative slope; the 3D rotation matrix reverses
  z; `fit` without scikit-learn logs an error and returns `None`. First answer (default):
  do not fix them. Maintainer answer: fix all four in S6 ("S6 decisions" item 3). The
  fixes change numbers for 3D users. `ROADMAP.md` gets a decision entry (behaviour change)
  and `COMPAT.md` gets a migration notice. No caller in the repository uses the class in
  three dimensions (fact 20).
- **OP21**: The native extension (`loop_interpolation_core`) and Surfe (`surfepy`) are not
  installed in the venv, so the two extension cells of S6 never ran when this plan was
  written. Their tolerances are not measured. Owner: the maintainer. Decision needed: (a)
  whether a CI job installs the two dependencies so that the optional cells run before a
  release; (b) what to do if a cell misses SC-002. Default until an answer: (a) no CI
  change; the absence is recorded as a named exclusion. (b) The limits are SC-002
  (`1e-8*F` and `1e-8*max(1,|gradient|)`) with the native solver tolerance `1e-12`; a miss
  is reported as a fit-accuracy finding (as OP7) and the limit is not relaxed. Blocks:
  nothing; SC-018 passes with named exclusions.
- **OP22**: Scope of the diagnostics tests. Default until an answer: S6 adds P2 cells for
  the state and block tests (no analytic bound, because of issue #399) and does not add
  a residual report test for fold features; the report code does not depend on the fold
  data path. Owner: the maintainer. Blocks: nothing.
- **OP23**: Other copies of public helpers exist: `LoopStructural/utils/maths.py` and
  `loop_common/math/_maths.py` (the same function names), and the legacy modules
  `LoopStructural/interpolators/_constant_norm.py` and `_discrete_fold_interpolator.py`
  (not imported by `LoopStructural.interpolators`, which re-exports `loop_interpolation`).
  FR-025 lists the "audited" utilities, and the survey names only the `_transformation.py`
  pair. Owner: the maintainer. Decision needed: audit these copies in a later spec.
  Default until an answer: S6 does not change or test them. Blocks: nothing.
- **OP24**: The class `EuclideanTransformation` rotates about z only. The 3D `fit` finds
  the main axis from the x and y columns (PCA with 2 components) and takes the mean of all
  three columns as the translation. The z column is not part of the axis. Alternative: a
  PCA of all three columns, and the angle from the xy part of the first vector. This gives
  another angle when the z spread is large. Owner: the maintainer. Decision needed: accept
  the xy axis. Default until an answer: the xy axis ("S6 decisions" item 3a). Blocks:
  nothing.

## Answers

- 2026-11-04 (wave 3, user decision): the T006 resize repair changes the
  golden file `packages/loop_common/tests/data/structured_support_reference.npz`.
  27 values change: `after_origin.nsteps` and `after_maximum.nsteps` and
  `after_maximum.maximum`, for 3 support classes and 3 grids. The old values
  did not enclose the requested bounds (FR-014 defect). The file is added to
  the `Files` of T006. T006 regenerates it with
  `python packages/loop_common/tests/test_structured_support_reference.py --update`.
  T015 records the change under `## Golden review (C-004)`.
- 2026-11-04 (spec review, user decision): `FaultSegment.fault_ellipsoid`
  called `model.rescale` on `builder.fault_centre`, which S1 makes a world
  point. `LoopStructural/modelling/features/fault/_fault_segment.py` is added
  to the `Files` of T013. The ellipsoid uses the world centre with no
  `rescale`. A regression test goes in `test_fault_frame_analytic.py` (T011).
- 2026-11-04 (spec review, user decision): fix two minor record findings in
  S1. The golden file keeps only the 27 FR-014 changes; the 33 other keys
  get their old values back. `verification.md` gets the correct FR ranges in
  its "Pre-fix failures" titles.

- 2026-11-04 (OP2, maintainer decision, spec C-007): the R2 and R5 matrix
  stays in the default `pytest` run with no marker. Measured S1 time: about
  5 minutes for the four new test files. OP2 is closed.

- 2026-11-05 (S3 plan, OP4, maintainer decision): the policy for a transform
  change after an interpolator holds data is "reject at use". The interpolator
  raises `StaleFrameError` at the next set, evaluate or save and changes no
  state. S3 builds no rejection inside `BoundingBox` and no automatic rebuild.
  OP4 is closed.
- 2026-11-05 (S3 wave 15, maintainer decision, OP7): the P2 cell of
  `test_rebuild_in_new_frame_preserves_world_evaluation` (T037) uses the exact
  nodal control (`exact=True`). The spec allows it ("Exact evaluation control").
  No tolerance changes. The P2 fit accuracy is open point OP7. `fdi`, `pli`,
  `cn_fdi` and `cn_pli` keep the fit (`exact=False`).
- 2026-11-05 (C-018, maintainer decision, later the same day): the maintainer will
  release a major version but not yet. The legacy-load change of FR-019 stays
  backward compatible on the 1.x line. A legacy payload with no choice loads with
  the identity frame and emits one `FutureWarning`. The code carries a note to
  remove the fallback at the next major release. This replaces the S3 answer
  "no warning-only identity fallback" for 1.x and makes the C-015 entry below
  moot for the loader (no break ships in 1.x). See "S3 amendment".
- 2026-11-05 (spec gap interview, C-015, maintainer decision): the legacy-load
  break (FR-019) ships in a minor release, with a note. It does not ship as a major
  release. `ROADMAP.md` records an explicit exception to the rule that 1.x stays
  backward compatible, and `COMPAT.md` and the release notes state the break. This
  replaces the earlier text "must be a major version" (S3 answer below). The
  constitution asks for a "major-version note"; the note is the statement of this
  decision. If the maintainer reads the constitution as a demand for a major
  release, the maintainer amends it by pull request. The text of `ROADMAP.md` and
  `COMPAT.md` that stage S3 added (pull request #397) says "must be a major
  version" and needs the same change.
- 2026-11-05 (spec gap interview, C-016, maintainer decision, OP11): S4 deprecates
  `GeologicalModel.scale` together with `rescale`. Its warning names
  `model.bounding_box.project(points)`. Nothing in the library calls `scale`.
- 2026-11-05 (spec gap interview, C-017, maintainer decision, OP7): SC-002 for P2
  applies to the exact nodal field of R1. The P2 fit accuracy is a solver property.
  Issue #399 records it. The backend table of the spec says so.
- 2026-11-05 (S4 plan, maintainer request): deprecate `GeologicalModel.rescale`.
  The method keeps working and emits a `DeprecationWarning` that says it will be
  removed. The four discarded calls in `exporters.py` go (OP8, calls). `scale`
  follows in the same release (C-016).
- 2026-11-05 (S3 end review, maintainer decision): a legacy payload needs both
  `local_origin` and `local_rotation`, or `declare_identity=True`. A partial
  input raises `InterpolatorFormatError` (`legacy_no_transform`) and names the
  missing argument. The earlier default (zero origin, identity rotation) is
  removed from the legacy branch. A v2 payload still accepts a partial frame
  that equals the saved transform.
- 2026-11-05 (S3 plan, OP5, maintainer decision): the `TypeError` that names
  `format_version`, raised by a reader from before S3, is the accepted
  rollback behaviour for SC-017. S3 ships no reader patch. OP5 is closed.
- 2026-11-05 (S3 plan, maintainer decision): the legacy-load break (FR-019)
  gets a ROADMAP.md decision and a `COMPAT.md` row with a major-version note.
  No warning-only identity fallback exists. The release that contains S3 must
  be a major version, or the maintainer records another decision. (Replaced on
  2026-11-05 by C-015: a minor release with a note.)

- 2026-11-04 (S2 wave 9, user decision): T022 found an existing defect.
  `ConstantNormFDIInterpolator` always fails with `IndexError` in
  `ConstantNormInterpolator.add_constant_norm` (`_constant_norm.py`),
  because `StructuredGrid.element_size` is one scalar and the code indexes
  it per element. The defect does not depend on the frame. The statement
  "S2 adds tests only" for the constant-norm variants is changed: new task
  T032 (wave 10, `[S]`) changes only `_constant_norm.py` and uses the
  scalar for each selected element. T022 keeps the `cn_fdi` cells; they
  fail before T032 and pass after it. T030 and T031 also need T032. The
  PLI control meets SC-002 (value error 1.8e-12, gradient error 1.3e-10).

- 2026-11-05 (S5 plan): OP1 is closed for S5. After an independent review of
  the S5 design, the maintainer answered five questions on 2026-11-05. The S5
  design was revised to match.
- 2026-11-05 (Q1, 2D `maximum`, maintainer): do not change the value of the 2D
  `maximum`. Derive `nsteps_cells`, `n_cell_x` and `n_cell_y` from `nsteps`, and
  keep `maximum = origin + nsteps * step` (legacy; it counts nodes as cells). This
  is kept for FR-028. The four existing assertions in
  `packages/loop_common/tests/test_2d_discrete_support.py` (lines 14 and 22) and
  `tests/unit/interpolator/test_2d_discrete_support.py` (lines 16 and 24) stay
  unchanged and must pass. OP14 is closed.
- 2026-11-05 (Q2, misfit result, maintainer): keep the old public names and the
  old public shape of `evaluate_value_misfit` and `evaluate_gradient_misfit`. The
  new non-mutating, world-correct 1-D result has new names
  (`evaluate_gradient_cosines`, `evaluate_value_residuals`). The old methods call
  the same corrected code. The gradient value changes meaning (a cosine). This is
  a behaviour change, recorded in `COMPAT.md` (migration notice) and in a
  `ROADMAP.md` decision entry. OP15 is closed.
- 2026-11-05 (Q3, stale model, maintainer): the misfit methods change no state
  and do not rebuild a stale model. S5 chose one rule: warn once with
  `UserWarning` and return NaN. OP13 is closed.
- 2026-11-05 (Q4, `set_nelements`, maintainer): accept the exact extent after
  `set_nelements`. S5 adds a numeric-test survey, a new section "S5 golden review
  (C-004)", and a regeneration of the golden data files after review (C-004,
  SC-015). The earlier sentence that the 3D golden test passes with no change is
  withdrawn. `set_nelements` clamps `nsteps` to at least 2; the setter rejects
  non-whole values; the 2D `nsteps` assignment coerces to an int array.
- 2026-11-05 (Q5, value misfit, maintainer): the public value misfit keeps its
  normalisation by the field range. The physical residual uses an unnormalised
  private helper, so its unit "field" is true.
- 2026-11-05 (S5 build start, orchestrator, five gaps of the task texts, settled
  by reading the code): (1) `GeologicalFeature` and its base class have no
  `bounding_box` attribute. The misfit code uses `self.interpolator.bounding_box`
  (an attribute of `GeologicalInterpolator`, `None` when no frame is set), as the
  old code does. `self.model.bounding_box` is not used: the model can be `None`.
  (2) A stale model for T081: for a stale interpolator, assign
  `interpolator.up_to_date = False` after the solve (a plain attribute). For a
  stale builder, assign `builder._up_to_date = False` after the solve. Do not use
  `builder.set_not_up_to_date`: the `GeologicalFeatureBuilder` override also sets
  the interpolator stale. The test asserts the flag value before the call.
  (3) `n_elements` exists on `StructuredGrid2DGeometry` (`nsteps_cells[0] *
  nsteps_cells[1]`, `_structured_grid_2d.py:45`) and on the 3D geometry (line
  166). (4) No existing test checks the keys of `ConstraintDiagnosticsReport.to_dict`.
  T071 writes the pre-S5 key set as a literal from the code: report keys
  `interpolator_type`, `total_rows`, `families`, `region_coverage`,
  `outside_model_points`; family keys `active`, `row_count`, `dropped_rows`,
  `effective_weight_mean`, `effective_weight_min`, `effective_weight_max`,
  `source_point_count`, `outside_model_point_count`; region keys
  `total_support_nodes`, `active_region_nodes`, `inactive_region_nodes`,
  `active_fraction`; report field order `interpolator_type`, `families`,
  `region_coverage`, `outside_model_points`. (5) The gradient family adds 2
  rows for each inside point for FDI (strike and dip rows, `name="gradient"`).
  PLI is different: `P1Interpolator.add_gradient_constraints` is `pass`
  (`_p1interpolator.py:70`), so PLI assembles no `gradient` rows, and the builder
  reports `dropped_rows = 2 * source_points` for the `gradient` family of PLI.
  So `verification_report` gives an `unattributed` record for PLI gradient data.
  That record is true (PLI ignores gradient data) and it is not a defect of S5.
  T081 checks the no-false-unattributed case with FDI only, and for PLI it records
  the measured `unattributed` count (expected `2 * source_points`) and the
  `outside_model` count. The worker reports both numbers.
- 2026-11-06 (S6 plan): OP1 is closed. The S6 design is exact ("Stage S6 design"). The
  first version of this entry said that no maintainer answer existed for OP18 to OP23. The
  entries below replace that statement. The S1 to S5 closed answers above stand. P2 fit
  accuracy is issue #399 and outside spec 005; S6 does not fix it.
- 2026-11-06 (S6 plan review, maintainer, OP18): amend the spec. S6 uses the library
  displacement paths (`FaultSegment.apply_to_points` and `evaluate_displacement`) on both
  frames, with the same algorithm. `spec.md` has the new clarification C-019. It replaces
  the "existing test-side evaluator and trishear parameters" text of R2 with "library
  displacement paths, same algorithm for both frames". The R2 fixture text changed to
  match. The SC-004 numbers do not change. OP18 is closed.
- 2026-11-06 (S6 plan review, maintainer, OP19): accept. `_value_residuals_fresh` and
  `_gradient_cosines_fresh` evaluate the interpolator directly at the stored positions,
  with no fault step. This reverses S5 decision 8 for those two private methods. The task
  text of `COMPAT.md` and the entry in `ROADMAP.md` say so. OP19 is closed.
- 2026-11-06 (S6 plan review, maintainer, native and Surfe stub cells): the cells
  `native-boundary` and `surfe-boundary` are mandatory and are listed apart. The spec rows 6
  and 7 show their status from the real extension cell only (excluded in this environment).
  The table prints the stub cells as separate lines. `main` of `backend_matrix` says so, and
  `test_every_spec_row_has_a_cell` does not count a stub cell for a row.
- 2026-11-06 (S6 plan review, maintainer, duplicate class): one class lives in
  `loop_common.math._transformation` and `loop_common.math` exports it.
  `LoopStructural/utils/_transformation.py` re-exports it. The `ROADMAP.md` entry records
  that this reverses item 2c-5 ("closed as deferred"). `COMPAT.md` says that the module now
  uses the `loop_common` logger, that a pickle made after S6 stores the path
  `loop_common.math._transformation` (an older release cannot load it), and that an old
  pickle still loads on new code (checked by reading: the old module path stays importable
  and exposes the class).
- 2026-11-06 (S6 plan review, maintainer, Surfe module name): `test_p0_surfe_nans.py` uses
  `importorskip("surfe")` and the adapter imports `surfepy`. A task reads the package data
  of `surfe` on PyPI and changes the name to `surfepy` only if the installed module is
  only `surfepy`. The file goes into "S6 files" and into the exclusion ledger.
- 2026-11-06 (S6 plan review, maintainer, matrix run): `backend_matrix.main(argv)` reads the
  junit files of the three gate runs (root, `packages/loop_common/tests`,
  `packages/loop_interpolation/tests`) and merges the outcomes. The key is `classname` plus
  `name`. A module-level skip is matched by file only, and for an optional cell with an
  absent dependency it counts as "skipped".
- 2026-11-06 (S6 plan review, maintainer, 1-D input): `transform` and `inverse_transform`
  raise `ValueError` for a 1-D input.
- 2026-11-06 (S6 plan review, maintainer, OP20): fix the other defects of
  `EuclideanTransformation` in S6: the 3D `fit`, the sign of the fitted angle, the reversed
  z axis and the silent need for scikit-learn. The planner chose a lazy import with an
  `ImportError` message, and no new dependency. `ROADMAP.md` gets a decision entry and
  `COMPAT.md` gets a migration notice. OP20 is closed as "fix". OP24 records the choice of
  the xy axis for the 3D fit.
- 2026-11-06 (S6 build start, orchestrator, gaps of the task texts, settled from the plan and the spec): (1) A `Cell` has `dependency: str | tuple[str, ...] | None`. A cell with several keys is "dependency absent" when any key is absent, and the exclusion reason names every absent module. `display-viewer` uses `("viewer", "matplotlib")`. This follows decision 9 (three display dependencies) and fact 17 (two checks in one test). (2) `SPEC_ROWS` holds the 7 first-column names of the spec table "Required backend coverage", copied exactly (see T094). (3) The MATRIX has 15 cells (the plan table lists 15; earlier text said 14): 11 mandatory (including the 2 stub cells) and 4 optional. (4) The root `test_transformation.py` passes only after T102 (T100 gives the `loop_common` class; T102 gives the re-export). (5) The final matrix run uses the same `PYTHONPATH` as the pytest runs. `pyevtk` is not on the path, so `export-pyevtk` is excluded and named. `loopstructuralvisualisation` is on the path, so `display-viewer` passes; a pass with an installed dependency is correct. `verification.md` lists which optional cells passed and which were excluded, with the exact dependency name. (6) A large task (T105, T107) that fails on size is split by the orchestrator, and the final report says so.
- 2026-11-06 (S6 plan review, maintainer, R2 checks): for the 37 fault-plane points the test
  reports the largest difference and asserts finite output only. The orthogonality residual
  rows are reported with a measured bound. Only the value, gradient, normal and curvature
  rows must vanish at 1e-9 on the analytic field. The plan lists their constraint names
  (fact 18).

## Summary

Observations, fault geometry, masks and public queries stay in world
coordinates. Local assembly and supports use the computational frame of the
target interpolator (`interpolator.bounding_box`). Each boundary converts
points once and directions once (rotation only).

S1 repairs the shared boundary and the fault frame:

1. A rigid-rotation check and a corner-enclosure helper in common geometry.
2. Support resize that keeps the node and cell counts consistent.
3. Public world/local conversion methods on `GeologicalInterpolator`.
4. Reference sampling (orthogonality, equality, feature inequality, restored
   support expansion) at world locations, with directions converted into the
   target basis.
5. A fault factory that does not project an explicit centre, keeps the
   requested buffer and forced geometry, infers a 3D centre, and adds one
   extent-direction row at the centre.

S2 repairs the fold paths. S3 makes saved interpolators, copies and transform
changes keep the frame ("Stage S3 design"). S4 removes two extra conversions: the
display and export of observation data, and the intrusion classification. It also
deprecates `GeologicalModel.rescale` ("Stage S4 design"). S5 repairs the resize paths and the misfit diagnostics and adds the residual and
exclusion report ("Stage S5 design"). S6 covers the auxiliary utilities and the
validation matrix (see "Stages and rollback").

## Technical context

- **Packages changed**: `LoopStructural` (modelling builders and the feature
  factory), `packages/loop_common` (geometry), `packages/loop_interpolation`
  (base and discrete interpolators). S1 does not change
  `packages/loop_interpolation_core`.
- **Languages and versions**: Python only. Ruff target py39. New modules use
  `from __future__ import annotations` for `X | None` hints.
- **New dependencies**: none.
- **Public API and import changes**: additive only in S1. New public names:
  `loop_common.geometry.validate_rigid_rotation`,
  `BoundingBox.local_enclosure`,
  `GeologicalInterpolator.world_to_local_points`,
  `GeologicalInterpolator.local_to_world_points`,
  `GeologicalInterpolator.world_to_local_vectors`,
  `GeologicalInterpolator.local_to_world_vectors`,
  `GeologicalFeatureBuilder.excluded_rows`. The signatures of
  `GeologicalModel.create_and_add_fault`, `FaultBuilder` and
  `set_interpolation_geometry` do not change. Behaviour changes with known
  coordinate meaning (FR-028): an explicit `fault_center` is used as a world
  point; `fault_buffer` and `force_mesh_geometry` keep their effect; an
  inferred centre keeps its elevation; `BoundingBox.set_local_transform` and
  `rotation_xy` reject a non-rigid matrix.
- **Saved-data format changes**: none in S1 and S2. S3 adds `format_version: 2`, `constraint_space`, `transform` and `bounding_box` to the interpolator payload (FR-016; see "S3 saved format").
- **Coordinate contract**:

  | Boundary | Input space | Output space | Function |
  |---|---|---|---|
  | Public evaluation | world points | world values, world directions | `GeologicalInterpolator.evaluate_value`, `evaluate_gradient` (unchanged) |
  | Data to interpolator | world constraint arrays | local stored arrays | `GeologicalInterpolator.set_*_constraints` (unchanged) |
  | Support samples to reference feature | local support points | world points | `GeologicalInterpolator.local_to_world_points` |
  | Reference direction to assembly | world direction | target-local direction | `GeologicalInterpolator.world_to_local_vectors` |
  | Region and mask callbacks | world points | boolean | caller callbacks (documented contract) |
  | `DiscreteInterpolator.add_gradient_orthogonal_constraints` | local points, local unit vectors | rows | internal local contract, now documented |
  | `DiscreteInterpolator.add_inequality_feature` feature callback | world points | values | documented feature-facing contract |
  | Support geometry request | world box | local enclosure | `BoundingBox.local_enclosure` |
  | Fault restoration for support expansion | world points | world points | `FaultSegment.apply_to_points` (unchanged) |

  The geological-frame values (fault coordinates 0, 1, 2) are scalars. No
  translation applies to them.
- **Performance target**: none specified. OP2 records the test run time.
- **Platforms and backends**: S1 covers Python FDI and PLI/P1. P2, fold,
  constant-norm, native and Surfe cells belong to S2 and S6. The native
  adapter rejects rotated supports (`_native.py:170`); S1 does not change it.

## Constitution check

| Project rule | Result | How the plan agrees |
|---|---|---|
| Work stays within the existing library and package responsibilities | Pass | Changes stay in `loop_common` geometry, `loop_interpolation` base classes and `LoopStructural` builders. One new private module `_rigid.py` in `loop_common.geometry`. No new package. |
| Development commands come from .dev-tools/project.json and the constitution | Pass | Gates use `pytest` and `ruff check LoopStructural packages`, and `black .` for format. |
| Both pytest and ruff check LoopStructural packages must pass before merge | Pass | Each stage row in "Stages and rollback" names both gates. |
| Python changes use numpy-style docstrings and type hints; Ruff targets py39 | Pass | Each new function has a numpy docstring and hints; `from __future__ import annotations` keeps py39 syntax. |
| Black uses line length 100 and preserves string quoting | Pass | Run `black .`; no string normalisation. |
| ROADMAP.md owns architecture stages and decisions; update it when either changes | Pass | S1 adds no ROADMAP stage. S3 records the FR-019 legacy-load decision in ROADMAP.md. The final stage adds a "005 complete" line. |
| Stages can be tested and reversed independently | Pass | Six stages, each with an acceptance check and a rollback target. |
| The public GeologicalModel API stays compatible; a break needs a ROADMAP.md decision and a major-version note | Pass | No signature or return-shape change. Corrected behaviour of known-meaning inputs gets a `COMPAT.md` "Migration notices" entry. S3 handles the FR-019 exception. |
| The QGIS plugin stays in its separate repository | Pass | No plugin file changes. |
| Feature documents stay under specs/NNN-name/ | Pass | `plan.md` and `verification.md` are in `specs/005-coordinate-frame-correctness/`. |
| Releases use release-please with per-package paths | Pass | Conventional commits per package path; no release config change. |
| A constitution amendment requires a pull request and maintainer acceptance | Not applicable | No amendment. |

## Design

S1 to S6 have the exact design below.

### Types and data

1. **`FaultBuilder` state** (`LoopStructural/modelling/features/builders/_fault_builder.py`).
   New private fields, set in `__init__`:
   - `self._fault_buffer: float = 0.2` — last requested support buffer, as a
     fraction of the longest box side.
   - `self._force_mesh_geometry: bool = False` — last forced-geometry request.

   `create_data_from_geometry` writes both fields from its arguments before
   any other work. `_on_geometry_change` reads them when its arguments are
   `None`.

2. **`GeologicalFeatureBuilder.excluded_rows`**
   (`LoopStructural/modelling/features/builders/_geological_feature_builder.py`).
   Type `dict[str, dict[str, int]]`. Key 1: constraint name
   (`f"{feature.name}_orthogonal"` or `f"{feature.name}_equality"`). Key 2:
   reason. Value: number of rows not installed. Reasons:
   - `"invalid_reference_value"` — the reference value is not finite.
   - `"invalid_reference_direction"` — the reference gradient is not finite
     or has zero length.
   - `"evaluation_error"` — the reference evaluation raised an exception; the
     count is the number of selected rows.

   The builder resets the dictionary at the start of
   `install_gradient_constraint` (orthogonal keys) and
   `install_equality_constraints` (equality keys). S5 connects it to the
   diagnostics report.

3. **Generated fault rows** (data frame columns `X, Y, Z, feature_name, nx,
   ny, nz, val, coord, w`, world coordinates). S1 adds one row when
   `major_axis` is not `None`:

   | X, Y, Z | nx, ny, nz | val | coord | w |
   |---|---|---|---|---|
   | fault centre (world) | `strike_unit / major_axis` | NaN | 2 | `w` |

   `strike_unit = cross(n, s_perp) / |cross(n, s_perp)|`, where `n` is the
   unit normal and `s_perp` is the unit slip direction after the normal
   component is removed. For R2, `n=(1,0,0)`, `s=(0,0,-1)`, so
   `strike_unit=(0,1,0)` and the row is `(0, 1/20, 0)`. This is the only
   automatic direction row (R4, SC-013). Existing tip rows (`val=+-0.5`,
   `coord=2`) do not change.

### Functions and public API

All names below were checked with `git grep -n -w <name>` on 2026-11-04 in
`LoopStructural`, `packages` and `tests`: 0 matches for each new name.

**`packages/loop_common/src/loop_common/geometry/_rigid.py`** (new)

```python
def validate_rigid_rotation(
    matrix: np.ndarray, dimensions: int, atol: float = 1e-10
) -> np.ndarray:
    """Return a float copy of ``matrix`` if it is a proper rigid rotation.

    Raises RigidRotationError (a ValueError subclass in this module) when:
    shape != (dimensions, dimensions); any entry is not finite;
    max|R^T R - I| > atol; |det(R) - 1| > atol.
    The message names the failed condition and the action:
    "supply a proper rotation (orthonormal, det +1); scale and shear are
    not supported".
    """
```

`RigidRotationError(ValueError)` is defined in the same module. It is a
`ValueError` subclass, so existing `except ValueError` callers still work.
Export both names from `loop_common.geometry.__init__`. The module imports
only `numpy`, so it causes no import cycle.

**`BoundingBox`** (`packages/loop_common/src/loop_common/geometry/_bounding_box.py`)

- `set_local_transform(local_origin=None, rotation_matrix=None)`: call
  `validate_rigid_rotation(rotation_matrix, self.dimensions)` before any field
  changes. Also reject a non-finite `local_origin` with `LoopValueError`.
  On error the old transform stays (SC-012).
- New method:

  ```python
  def local_enclosure(
      self,
      origin: np.ndarray | None = None,
      maximum: np.ndarray | None = None,
  ) -> tuple[np.ndarray, np.ndarray]:
      """Return (local_min, local_max) that enclose the projection of all
      2**dimensions corners of the world box [origin, maximum].
      Defaults: self.origin and self.maximum. Raises LoopValueError if
      either bound is None or not finite, or if origin > maximum on any axis.
      """
  ```

  Implementation: build the corners with `itertools.product` over the two
  bounds per axis, call `self.project(corners)`, return `min` and `max` on
  axis 0.

**`StructuredGrid3DGeometry`** (`packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`)

- `rotation_xy` setter: after the existing float/int conversion, call
  `validate_rigid_rotation(rotation_xy, 3)`. `None` keeps the current value.
- New private method used by both setters:

  ```python
  def _cells_for_length(self, length: np.ndarray) -> np.ndarray:
      """Number of cells per axis that cover ``length`` with the current step.

      cells = ceil(length / step_vector - 1e-9). An axis with 0 cells gets 2
      cells (the existing degenerate rule of 3 nodes). Raises ValueError
      (existing message "Cannot resize the grid.") if any cell count < 0.
      """
  ```

- `origin` setter: `self._nsteps = self._cells_for_length(self.maximum - origin) + 1`.
  Compute `self.maximum` before the origin changes.
- `maximum` setter:
  `self._nsteps = self._cells_for_length(maximum - self.origin) + 1`.

  This repairs the 1-cell loss (node count was set to the cell count).
  Reapplying unchanged bounds keeps origin, maximum, step and counts
  (SC-011, first part).

**`GeologicalInterpolator`** (`packages/loop_interpolation/src/loop_interpolation/_geological_interpolator.py`)

New public methods. Each is a documented wrapper of the existing private
method. The private methods stay for compatibility.

```python
def world_to_local_points(self, points: np.ndarray) -> np.ndarray: ...
def local_to_world_points(self, points: np.ndarray) -> np.ndarray: ...
def world_to_local_vectors(self, vectors: np.ndarray) -> np.ndarray: ...
def local_to_world_vectors(self, vectors: np.ndarray) -> np.ndarray: ...
```

Input shape `(N, d)` or `(d,)`; output has the same shape. Vectors use
rotation only. With `bounding_box is None` each method returns its input
(existing behaviour).

**`DiscreteInterpolator`** (`packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`)

- `add_gradient_orthogonal_constraints` (abstract): docstring states the
  internal contract — `points` are target-local, `vectors` are target-local
  unit directions. No signature change.
- `add_inequality_feature(feature, lower=True, mask=None)`: call
  `feature(self.local_to_world_points(self.support.nodes))`. Docstring:
  `feature` receives world points and returns one value for each node. No
  signature change.

**`GeologicalFeatureBuilder`** (`LoopStructural/modelling/features/builders/_geological_feature_builder.py`)

- `install_gradient_constraint()`:
  1. `local_pts` = barycentres or quadrature points (existing code).
  2. `world_pts = self.interpolator.local_to_world_points(local_pts)`.
  3. `vector = feature.evaluate_gradient(world_pts, ignore_regions=True)`.
  4. `valid = all(isfinite(vector), axis=1) & (norm > 0)`.
     Record `(~valid).sum()` under `"invalid_reference_direction"`.
  5. Normalise `vector[valid]`, then
     `local_vec = self.interpolator.world_to_local_vectors(vector[valid])`.
  6. Call `add_gradient_orthogonal_constraints(local_pts[valid], local_vec,
     w=weight[valid] if weight is an array else weight, b=B, name=...)`.
  7. If `valid.sum() == 0`, install no rows and log a warning with the count.
- `install_equality_constraints()`:
  1. `world_nodes = self.interpolator.local_to_world_points(support.nodes)`.
  2. `selected = np.asarray(region(world_nodes), dtype=bool)`;
     `idc = np.flatnonzero(selected)`.
  3. `val = feature.evaluate_value(world_nodes[idc])`; `valid = isfinite(val)`.
  4. Record `(~valid).sum()` under `"invalid_reference_value"`.
  5. `self.interpolator.add_equality_constraints(idc[valid],
     val[valid] * scalefactor, name=f"{feature.name}_equality")`.
  6. The existing `except` clause stays. It also records `idc.size` (or the
     node count if `idc` is not set) under `"evaluation_error"`.
- `set_interpolation_geometry(origin, maximum, rotation=None)`: replace the
  two-point projection with
  `origin, maximum = self.interpolator.bounding_box.local_enclosure(origin, maximum)`
  when `bounding_box` is not `None`. The NaN checks and the `dof < 100` loop
  do not change.
- `check_interpolation_geometry(data, buffer=0.3)`:
  1. `pts = self.model.bounding_box.with_buffer(buffer).regular_grid(local=False)`
     (world).
  2. For each fault: `pts = f.apply_to_points(pts)` (world in, world out).
  3. `local = self.interpolator.world_to_local_points(pts)`.
  4. `origin = np.minimum(support.origin.copy(), local.min(axis=0))`;
     `maximum = np.maximum(support.maximum.copy(), local.max(axis=0))`.
  5. Assign `support.origin = origin`, then `support.maximum = maximum`.
     The copies stop the in-place change of the support arrays.

**`FaultBuilder`** (`LoopStructural/modelling/features/builders/_fault_builder.py`)

- `create_data_from_geometry(...)`: no signature change.
  - Store `self._fault_buffer = fault_buffer` and
    `self._force_mesh_geometry = force_mesh_geometry` first.
  - Inferred centre (`fault_center is None`): use the rows with `coord == 0`
    and `val == 0`, columns `X, Y, Z`; `fault_center = np.nanmean(pts, axis=0)`.
    If a component is not finite, set it to `0.0` and log a warning that names
    the axis. This keeps the elevation (FR-010).
  - Slip direction: after normalisation, `s_perp = s - (s @ n) * n`. If
    `|s_perp| < 1e-12`, raise `ValueError("fault slip vector is parallel to
    the fault normal; supply a slip direction in the fault plane")`. Use
    `s_perp / |s_perp|` for `self.fault_slip_vector`, for the strike vector,
    and for the depth rows.
  - Add the extent-direction row (see "Types and data", item 3).
  - Early-return branch (no fault trace): call
    `self.set_mesh_geometry(fault_buffer, None)` before `return`. This keeps
    the support that the factory previously set.
- `_on_geometry_change(force_mesh_geometry=None, fault_buffer=None, ...)`:
  the default of `force_mesh_geometry` changes from `False` to `None`. `None`
  means "use the stored value". This is a private method.

**`create_and_add_fault`** (`LoopStructural/modelling/core/_model_feature_factory.py`)

- Delete the call `fault_center = model.scale(fault_center, inplace=False)`.
  An explicit centre stays a world point (FR-010, U1.1).
- Delete the block
  `if "force_mesh_geometry" not in kwargs: fault_frame_builder.set_mesh_geometry(kwargs.get("fault_buffer", 0.2), 0)`.
  `create_data_from_geometry` already applies the named `fault_buffer` and
  `force_mesh_geometry` (FR-012, U1.4, U1.5).

### Messages and saved data

None in S1. S3 owns the saved-format change (FR-016 to FR-019).

### Behaviour

Call sequence for a fault through the public API after S1:

1. `GeologicalModel.create_and_add_fault(...)` calls the factory with world
   data and an optional world `fault_center`.
2. The factory makes `FaultBuilder` with `model.bounding_box` and calls
   `create_data_from_geometry` with the unchanged centre.
3. The builder makes the world rows, sets `origin`/`maximum` (forced or from
   rows), and calls `set_mesh_geometry(self._fault_buffer, None)`.
4. `set_mesh_geometry` calls `set_interpolation_geometry` on the three
   coordinate builders. Each builder projects the world box with
   `local_enclosure` and resizes its support with the repaired setters.
5. `build` sends world rows to `set_*_constraints`. The interpolator converts
   them once into local arrays (existing code).
6. `install_gradient_constraint` samples the reference coordinate at world
   points, converts the directions into the target basis, and installs rows
   at local points.
7. Evaluation and displacement take world points and return world values and
   directions (existing code).

Errors: an invalid rotation raises `RigidRotationError` before any state
change. A slip parallel to the normal raises `ValueError` before the data
frame is added to the builder.

Cache invalidation: S1 changes no cache. S3 owns transform changes after
population (FR-020).

### Stage S2 design

S2 makes the fold families, the directional regularisation, the fold masks,
the fold sign check, P2 reference rows and the constant-norm variants follow
the S1 boundary contract. Base: branch `feat/005-coordinate-frame-correctness`
(S1, pull request #392). If S1 changes before it merges, run this planner
again on the merged result.

#### S2 facts from the code (read on 2026-11-04)

- Both fold providers, `LoopStructural.modelling.features.fold.FoldEvent`
  and `loop_interpolation.FoldEvent`, give the points to fold-frame
  features. Those features take world points. So
  `get_deformed_orientation(points)` takes world points and returns world
  unit directions `(deformed_orientation, fold_axis, dgz)`, each `(N, 3)`.
  Rows with invalid geometry contain NaN.
- `FDFoldInterpolator.add_fold_constraints` calls the provider with
  `self.support.nodes` (local) and uses the result as local vectors.
  `mask_fn` gets local nodes.
- `DiscreteFoldInterpolator.add_fold_constraints` calls the provider with
  `self.support.barycentre` (local) and multiplies the result with local
  element gradients. `mask_fn` gets local barycentres. Its regularisation
  provider calls the fold provider with the sample points that the base
  method gives (local).
- `P2FoldInterpolator.add_fold_constraints` calls the provider with local
  quadrature points from `tet_quadrature_points`. `mask_fn` gets the same
  local points.
- `resolve_fold_norm_target` (`_fold_norm_alignment.py`) gets the stored
  normal constraints, which are local (`get_norm_constraints()`), and calls
  the provider with their local points. It then compares the world `dgz`
  with the local normals.
- `DiscreteInterpolator.add_directional_regularisation` calls a callable
  `direction` with `get_regularisation_sample_points()` (local: FDI nodes,
  P1 shared-face centres) and uses array or callable results as local
  vectors. The public setup argument `directional_regularisation` reaches
  this method. `_native.py` (`_add_regularisation`, lines 465 to
  485) repeats the same logic for the native core.
- `P2Interpolator.minimise_grad_steepness` calls the public
  `steepness_wtfunc` with all local barycentres, but it applies the weights
  only to the masked interior elements. The two lengths differ when the
  mask removes elements.
- `ConstantNormInterpolator.add_constant_norm` uses only local gradients of
  the current solution and a scalar norm target. A rigid rotation keeps a
  gradient norm. So it needs no code change. S2 adds tests only.
- `P2Interpolator.add_gradient_orthogonal_constraints` already follows the
  S1 local row contract, and `tet_quadrature_points` uses only the 4 corner
  nodes of a P2 element. So the S1 builder path works for P2. S2 adds tests
  only.
- `LoopStructural/interpolators/_discrete_fold_interpolator.py` is not
  imported. `LoopStructural.interpolators` re-exports the
  `loop_interpolation` class. S2 does not change this unused module.
- `_quadrature.py`, `_regularisation.py`, `_finite_difference_interpolator.py`,
  `_p1interpolator.py`, `_constant_norm.py` and `_folded_feature_builder.py`
  have no world/local crossing. S2 does not change them.

#### S2 contracts

| Boundary | Input space | Output space | Owner |
|---|---|---|---|
| Fold provider `get_deformed_orientation(points)` | world points | world unit directions | fold provider (unchanged; now documented) |
| Fold `mask_fn(points)` | world points | boolean, one for each point | caller callback (documented) |
| `DirectionalRegularisation.direction` as a callable | world points | world directions, one for each point | caller (documented) |
| `DirectionalRegularisation.direction` as an array | none | world directions, one for each sample point | caller (documented) |
| `steepness_wtfunc(points)` (P2) | world points | weights, one for each weighted element | caller (documented) |
| `resolve_fold_norm_target(normal_constraints=...)` | world points and world normals | float | `_fold_norm_alignment.py` |
| Fold rows and directional terms in the assembly | local points, local unit directions | rows | interpolator (internal) |

Fold directions are not normalised again after the rotation, because a
rigid rotation keeps the length. NaN rows stay NaN.

#### S2 functions

- `packages/loop_interpolation/src/loop_interpolation/_fold_setup.py`, two
  new functions (numpy docstrings, type hints, `from __future__ import
  annotations`):

  ```python
  def sample_fold_geometry(
      interpolator, fold, local_points: np.ndarray
  ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
      """Return (world_points, deformed_world, axis_world, dgz_world)."""
      world_points = interpolator.local_to_world_points(np.asarray(local_points, dtype=float))
      deformed, axis, dgz = fold.get_deformed_orientation(world_points)
      return (
          world_points,
          np.asarray(deformed, dtype=float),
          np.asarray(axis, dtype=float),
          np.asarray(dgz, dtype=float),
      )


  def world_norm_constraints(interpolator) -> np.ndarray:
      """Return the stored normal constraints as world xyz|nx ny nz|w rows."""
  ```

  `world_norm_constraints`: set
  `n = np.asarray(interpolator.get_norm_constraints(), dtype=float)`. If
  `n.ndim != 2`, or `n.shape[0] == 0`, or `n.shape[1] < 6`, return `n`.
  Else return a copy where columns 0:3 are
  `interpolator.local_to_world_points(n[:, 0:3])` and columns 3:6 are
  `interpolator.local_to_world_vectors(n[:, 3:6])`. The other columns do
  not change.
- `_fd_fold_interpolator.py`, `FDFoldInterpolator.add_fold_constraints`
  (no signature change):
  1. `world_nodes, deformed_w, axis_w, dgz_w = sample_fold_geometry(self, self.fold, self.support.nodes)`.
  2. `deformed = self.world_to_local_vectors(deformed_w)`. Do the same for
     `axis` and `dgz`.
  3. `weight[mask_fn(world_nodes)] = 0.0`.
  4. The orientation, axis and normalisation rows use the local `deformed`,
     `axis` and `dgz` at the local `active_pos`. The other code does not
     change.
  5. `resolve_fold_norm_target(..., normal_constraints=world_norm_constraints(self), ...)`.
  6. The directional regularisation gets the masked WORLD arrays
     (`_masked_direction(dgz_w)`, `_masked_direction(deformed_w)`,
     `_masked_direction(axis_w)`). The base method converts them (see the
     contracts).
  7. Docstring: `mask_fn` receives world node positions.
- `_discrete_fold_interpolator.py`, `DiscreteFoldInterpolator.add_fold_constraints`:
  1. `world_bc, deformed_w, axis_w, dgz_w = sample_fold_geometry(self, self.fold, self.support.barycentre)`.
     Convert the three directions to local with `world_to_local_vectors`.
     Use the local vectors in the three `einsum` rows.
  2. `weight[mask_fn(world_bc)] = 0`.
  3. `normal_constraints=world_norm_constraints(self)`.
  4. The `_masked_fold_direction` provider now receives WORLD points from
     the base method. It calls `self.fold.get_deformed_orientation(points)`
     with those points, masks with `mask_fn(points)` and returns world
     vectors. The base method converts them.
  5. Docstring: `mask_fn` receives world barycentres.
- `_p2_fold_interpolator.py`, `P2FoldInterpolator.add_fold_constraints`:
  after `tet_quadrature_points`, call
  `world_points, deformed_w, axis_w, dgz_w = sample_fold_geometry(self, self.fold, points)`.
  Convert the three directions to local. Call `mask_fn(world_points)`. Use
  `normal_constraints=world_norm_constraints(self)`. `dN` and `inside`
  stay on the local points. Docstring: `mask_fn` receives world sample
  points.
- `_fold_norm_alignment.py`, `resolve_fold_norm_target`: no code change.
  Docstring of `normal_constraints`: "Normal constraints in world
  coordinates, ``xyz|nx ny nz|w``. Use ``world_norm_constraints``."
- `_discrete_interpolator.py`, `DiscreteInterpolator.add_directional_regularisation`:
  1. Keep `sample_points` (local). Add
     `world_samples = self.local_to_world_points(sample_points)`.
  2. Callable: `vectors = term.direction(world_samples)`. Array:
     `vectors = term.direction`.
  3. Keep the shape check against `sample_points.shape`.
  4. `vectors = self.world_to_local_vectors(vectors)`, then
     `self._add_directional_regularisation(term.weight, vectors, name=term.name)`.
  5. New numpy docstring: a callable direction receives world sample points
     and returns world directions. An array gives world directions, one row
     for each sample point.
- `_native.py`, `_add_regularisation` (`_native.py:454`; the directional part is lines 465 to 485): the
  same 4 steps, with `interpolator.local_to_world_points` and
  `interpolator.world_to_local_vectors`. The rotated-support rejection
  (`_native.py:170`) does not change.
- `_p2interpolator.py`, `P2Interpolator.minimise_grad_steepness`: when
  `wtfunc` is callable, set
  `wt = wtfunc(self.local_to_world_points(self.support.barycentre[elements[mask]])) * self.support.element_size[elements[mask]]`.
  Docstring: `wtfunc` receives the world barycentres of the weighted
  elements.

#### S2 test data

`tests/fixtures/coordinate_frames.py` gets these helpers, with no tests:

- `R1_GRADIENT = np.array([0.2, -0.3, 0.4])`;
  `r1_value(points, anchor)` returns `2 + (points - anchor) @ R1_GRADIENT`;
  `r1_gradient(points, anchor)` returns `R1_GRADIENT` tiled to `(N, 3)`;
  `r1_observations(anchor) -> pd.DataFrame` returns the 27 R1 fit points
  with columns X, Y, Z, val, feature_name ("r1") and w (1).
- `class AnalyticFold` with `__init__(self, deformed_fn, axis_fn, dgz_fn)`.
  `get_deformed_orientation(points)` returns the three functions at
  `points` and appends a copy of `points` to `self.seen_points` (a list).
  Class methods: `AnalyticFold.constant(deformed, axis, dgz)` (three world
  vectors) and `AnalyticFold.r4(anchor)`, which uses `r4_frame`: deformed
  is the unit extent gradient, axis is `(0,0,1)` and dgz is the unit
  surface gradient.

#### S2 tests

- `tests/unit/modelling/test_fold_constraint_frames.py` (new). Families:
  `FDFoldInterpolator` on `StructuredGrid`, `DiscreteFoldInterpolator` on
  `TetMesh`, and `P2FoldInterpolator` on `P2TetMesh` with the sample rules
  `centroid` and `gauss4`. Box `O` to `O+(8,10,12)`, 5x6x7 cells, 4 anchors
  x 3 rotations. The support is the local enclosure of the box, and
  `interpolator.bounding_box = model.bounding_box`. Use the same pattern as
  `_builder` in `test_reference_constraint_frames.py`. Tests:
  - `test_fold_geometry_sampled_at_world_points`: with `AnalyticFold.r4`,
    each array in `fold.seen_points` equals `local_to_world_points` of the
    local sample points of the family within `Tp` (FDI nodes, PLI
    barycentres, P2 quadrature points; the norm check sees the world normal
    points).
  - `test_fold_rows_affine_residual`: field `f = y'`. Compute the exact
    local coefficients as in the S1 affine test (no large-coordinate
    cancellation). Use `AnalyticFold.constant((1,0,0), (0,0,1), (0,1,0))`,
    `fold_norm=1.0` and `dgz_alignment="none"`. The rows "fold
    orientation" and "fold axis" (target 0) and "fold normalisation"
    (target 1) have `|A @ c - b| <= 1e-10` under the normalisation of each
    assembler (SC-008).
  - `test_fold_rows_follow_curved_directions`: with `AnalyticFold.r4`. For
    FDI, spy on `add_gradient_orthogonal_constraints` and check that the
    vectors equal `world_to_local_vectors(analytic(world point))` for each
    component within `1e-10` (SC-013). For PLI and P2, check the local
    directions that `sample_fold_geometry` and `world_to_local_vectors`
    give at the sample points of the family, and check that the
    orientation row of each element is `vol * (v . eg)` (PLI) or
    `v . dN` (P2) for that local `v` within `1e-10`. If the code uses the
    centre direction at every sample, this check fails.
  - `test_fold_mask_receives_world_points`: a recording `mask_fn` that
    excludes `x' > 4` receives world points equal to the sample points
    within `Tp`. Excluded rows have zero weight, and their count equals the
    analytic count.
  - `test_fold_norm_alignment_uses_world_normals`: one stored normal
    constraint with world normal `(0,-1,0)` and `dgz = (0,1,0)`. With
    `dgz_alignment="correct"`, the target is `-fold_norm` for every
    rotation and anchor. The unrotated control and the rotated cells agree.
  - `test_directional_regularisation_world_contract`: for FDI, PLI and P2,
    use the public `directional_regularisation=[{"direction": f, "weight": 1}]`
    with a recording callable `f`. `f` receives the world sample points
    within `Tp`. A spy on `_add_directional_regularisation` receives
    `world_to_local_vectors(f(world))` within `1e-10`. Repeat with an array
    of world directions.
  - `test_p2_steepness_wtfunc_receives_world_barycentres`: a recording
    `steepness_wtfunc` receives the world barycentres of the weighted
    elements, one for each row.
- `tests/unit/modelling/test_reference_constraint_frames.py` (S1 file,
  extend): add the engine `"P2"` (support `P2TetMesh`, interpolator
  `P2Interpolator`) to each R3 test that takes `engine`. The expected
  equality set comes from the world P2 nodes. Do not reuse 168.
- `tests/unit/modelling/test_constant_norm_frames.py` (new): the FDI and
  PLI constant-norm variants, R1, 4 anchors x 3 rotations, norm target
  `sqrt(0.29)`, the 27 R1 observations, and the SC-003 solve settings
  through `solver_kwargs`. At the two R1 query offsets: value error
  `<= 1e-8*F` and gradient component error `<= 1e-8*max(1, |gradient|)`
  (SC-002). The test task runs the unrotated anchor-0 control first. If
  the control fails SC-002, the task stops and reports it (OP3).

#### S2 constitution check

| Project rule | Result | How S2 agrees |
|---|---|---|
| Work stays within the existing packages | Pass | Changes stay in `loop_interpolation`. No new module; two new functions in `_fold_setup.py`. |
| Gates `pytest` and `ruff check LoopStructural packages` | Pass | The S2 row of "Stages and rollback" names both. Also run the package tests. |
| numpy docstrings, type hints, py39 | Pass | New functions and changed docstrings use numpy style. |
| The public API stays compatible | Pass | No signature change. The callback contracts of `mask_fn`, `directional_regularisation` and `steepness_wtfunc` change from local to world. A caller with an identity frame sees no change. `COMPAT.md` gets one "Migration notices" row. |
| Stages can be tested and reversed independently | Pass | S2 reverts to the S1 merge. |
| Feature documents stay under `specs/NNN-name/` | Pass | Only `plan.md`, `tasks.md` and `verification.md` change. |

#### S2 files

| File | Change |
|---|---|
| `packages/loop_interpolation/src/loop_interpolation/_fold_setup.py` | Add `sample_fold_geometry` and `world_norm_constraints` |
| `packages/loop_interpolation/src/loop_interpolation/_fd_fold_interpolator.py` | World sampling, local rows, world mask, world norm check, world directional arrays |
| `packages/loop_interpolation/src/loop_interpolation/_discrete_fold_interpolator.py` | The same for barycentres and the regularisation provider |
| `packages/loop_interpolation/src/loop_interpolation/_p2_fold_interpolator.py` | The same for quadrature points |
| `packages/loop_interpolation/src/loop_interpolation/_fold_norm_alignment.py` | Docstring only |
| `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py` | World contract in `add_directional_regularisation` |
| `packages/loop_interpolation/src/loop_interpolation/_native.py` | The same contract in `_add_regularisation` |
| `packages/loop_interpolation/src/loop_interpolation/_p2interpolator.py` | World `steepness_wtfunc` with the correct length |
| `tests/fixtures/coordinate_frames.py` | R1 helpers and `AnalyticFold` |
| `tests/unit/modelling/test_fold_constraint_frames.py` | New |
| `tests/unit/modelling/test_reference_constraint_frames.py` | Add the P2 engine |
| `tests/unit/modelling/test_constant_norm_frames.py` | New |
| `COMPAT.md` | S2 migration notice |
| `specs/005-coordinate-frame-correctness/verification.md` | S2 pre-fix and post-fix records, backend matrix rows |

The existing package tests with fake folds (`test_fd_fold_interpolator.py`,
`test_discrete_fold_interpolator.py`, `test_p2_fold_interpolator.py`,
`test_regularisation_api.py` and `test_native_parity.py`) must pass with no
change. With `bounding_box is None`, the conversions return their input. If
one of these tests encodes the old local callback contract under a frame
that is not the identity, record it under "Golden review" (C-004) before
any change.

### Stage S3 design

S3 makes saved interpolators, copies and transform changes keep the
coordinate frame (FR-016 to FR-020; U3.1 to U3.3 and U3.6 to U3.10). Base: the
S2 branch `feat/005-s2-fold-frames` (pull request #396). If S2 changes before
it merges, run this planner again on the merged result.

#### S3 facts from the code (read on 2026-11-05)

- `GeologicalInterpolator.to_dict` (`_geological_interpolator.py:880`) writes
  `type`, `data` (the stored arrays, which are local), `up_to_date`, `valid`
  and `support`. `DiscreteInterpolator.to_dict` (`_discrete_interpolator.py:1522`)
  adds `c` and merges `**super().to_dict()`. No key holds the transform, the
  bounding box, a format version or the coordinate space.
- `InterpolatorFactory.from_dict` (`_interpolator_factory.py:74`) pops `type`,
  `support`, `data`, `c`, `up_to_date` and `valid`. It passes every other key
  as `**d` to `create_interpolator`. It sets the arrays with
  `set_*_constraints`. The loaded interpolator has `bounding_box = None`, so
  `_project_constraint_array` does nothing and the arrays are not converted
  again. But the transform is lost, so a loaded rotated model evaluates world
  points as local points. If a loader attached the bounding box before the
  `set_*_constraints` calls, the setters would convert the local arrays a
  second time (the FR-017 defect).
- `from_dict` is also a creation call: `{"type", "boundingbox", "nelements"}`
  with no `support` dict builds an empty interpolator
  (`test_from_dict_builds_interpolator`). That use has no saved arrays and
  stays as it is.
- `BoundingBox.to_dict` and `from_dict` (`_bounding_box.py:672`, `:690`)
  already carry `local_origin` and `local_rotation`. `from_dict` validates them
  with `set_local_transform` (S1). The transform maps `local = R @ (world -
  origin)`, so `world = R.T @ local + origin`.
- `interpolator.bounding_box` is a plain attribute. The factory gives every
  interpolator of a model the same `BoundingBox` object as the model
  (`create_interpolator`, line 68). `BoundingBox.set_local_transform` changes
  that object in place. The interpolator gets no notice. After such a change
  `evaluate_value` joins the new transform to the old coefficients, and a later
  `set_*_constraints` call joins the new transform to old stored rows.
- Three `copy()` methods exist: `FiniteDifferenceInterpolator.copy`
  (`:247`), `P2Interpolator.copy` (`:106`) and `P2FoldInterpolator.copy`
  (`:64`, keeps `fold`). Each returns a new interpolator on the same support
  with `bounding_box = None`. The frame is lost. The piecewise linear
  interpolator has no `copy`. The fold feature `copy` methods under
  `LoopStructural/modelling/features` copy features, not interpolators (not
  read; OP6).
- No restored-point cache exists. `git grep -i -E 'cache|_restored'` over
  `LoopStructural/modelling`, `LoopStructural/interpolators` and
  `loop_common/supports` finds only geometry caches in
  `_3d_unstructured_tetra.py:77` and `_p2_structured_tetra.py:73`. They reset
  on a geometry change, and a transform change does not touch them. The state
  that depends on the frame is `data`, `c`, `up_to_date` and
  `latest_diagnostics_report`. Fault restoration reads the fault frame
  through `evaluate_value` and `evaluate_gradient` (the test task confirms
  this in `FaultSegment.apply_to_points` before it writes the R2 test).
- `_model_serializer.py` has no `interpolator` match. S3 does not change model
  persistence. `_support_factory.py` already rebuilds a support from
  `support.to_dict()`. S3 does not change it, and it does not change
  `_fault_segment.py` or `_discrete_interpolator.py` (the survey of S3 files
  in "Later stages" was too wide).
- `BoundingBox.fit` (`_bounding_box.py:472`) sets origin, maximum and the
  local transform. The model calls its bounding-box set-up before it makes
  builders. If a library path calls `fit` or `set_local_transform` after an
  interpolator holds data, the frame check of S3 raises there. The full
  `pytest` run finds such a path (see "S3 risks").

#### S3 decisions

1. **Policy for a transform change (FR-020): reject at use.** An interpolator
   that holds data in one frame refuses every use in another frame with
   `StaleFrameError`. The refusal changes no state. The caller makes a new
   interpolator in the new frame, or calls `clean()` and sets the constraints
   again. S3 adds no automatic rebuild: a rebuild must also resize the
   support, and an unstructured mesh cannot be resized. The bounding box is
   shared and does not know its interpolators, so S3 does not reject inside
   `BoundingBox.set_local_transform` (OP4).
2. **Format version 2.** A payload with `format_version: 2` carries
   `constraint_space`, `transform` and `bounding_box`. A payload with no
   `format_version` is legacy. A payload with a higher version is not
   supported. No code wrote version 1, so version 1 is malformed.
3. **Saved payload or creation call.** A dict is a saved payload if and only
   if `support` is a dict. Other dicts keep the creation behaviour.
4. **Legacy payload.** It has no recoverable transform (old `to_dict` never
   wrote one). The caller should pass `local_origin` and `local_rotation`, or
   `declare_identity=True`. On the 1.x line, a call with no choice loads with the
   identity frame and one `FutureWarning`; from the next major release it raises
   (amendment C-018, see "S3 amendment").
5. **Copy.** A copy shares the `BoundingBox` object of its source, as the
   factory does for the interpolators of one model. A loaded interpolator
   gets its own `BoundingBox`.

#### S3 contracts

| Boundary | Input | Output | Owner |
|---|---|---|---|
| `GeologicalInterpolator.to_dict` `data` arrays | stored local arrays | same arrays, declared `constraint_space: "local"` | `_geological_interpolator.py` |
| `transform` in the payload | `bounding_box.local_origin`, `local_rotation` | `local = R @ (world - origin)` | `frame_payload` |
| Loader of a saved payload | payload, optional caller frame | interpolator, or `InterpolatorFormatError` (nothing returned) | `InterpolatorFactory.from_dict` |
| Loaded `data` arrays | local arrays | stored as they are, with no conversion | `_restore_local_constraints` |
| Copy | interpolator | empty interpolator, same support, same `bounding_box` object | `_inherit_frame` |
| Use after a frame change | stored data, new transform | `StaleFrameError`, no state change | `_check_frame` |

#### S3 saved format

Three keys follow the existing keys of `to_dict`, in this order:

```json
{
  "type": "FINITE_DIFFERENCE",
  "support": {"...": "unchanged"},
  "c": [0.0],
  "data": {"...": "unchanged, local arrays"},
  "up_to_date": true,
  "valid": true,
  "format_version": 2,
  "constraint_space": "local",
  "transform": {
    "dimensions": 3,
    "local_origin": [500000.0, 6500000.0, 100.0],
    "local_rotation": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  },
  "bounding_box": {"origin": [0.0, 0.0, 0.0], "maximum": [8.0, 10.0, 12.0], "nsteps": [6, 7, 8]}
}
```

`bounding_box` is `BoundingBox.to_dict()` without its two `local_*` keys. It
is `null` when `interpolator.bounding_box is None`; then `transform` is the
identity. `format_version` is the first new key. A reader from before S3 passes
the unknown keys to `create_interpolator` as `**d`, so it fails with `TypeError:
... unexpected keyword argument 'format_version'` before it makes an object
(OP5).

#### S3 functions and signatures

New module `packages/loop_interpolation/src/loop_interpolation/_frame_persistence.py`:

```python
FORMAT_VERSION: int = 2
CONSTRAINT_SPACE: str = "local"
REASONS: tuple[str, ...] = ("legacy_no_transform", "unsupported_version", "malformed", "conflict")

class InterpolatorFormatError(ValueError):
    reason: str  # one of REASONS
    def __init__(self, reason: str, message: str) -> None: ...

class StaleFrameError(ValueError): ...

def frame_payload(bounding_box: BoundingBox | None, dimensions: int) -> dict: ...
def payload_kind(payload: dict, caller_frame: bool) -> str: ...  # "v2" or "legacy"
def resolve_frame(
    payload: dict,
    kind: str,
    support: object,
    dimensions: int,
    *,
    local_origin: np.ndarray | None = None,
    local_rotation: np.ndarray | None = None,
    declare_identity: bool = False,
) -> BoundingBox | None: ...
def bounding_box_from_support(
    support: object, dimensions: int, local_origin: np.ndarray, local_rotation: np.ndarray
) -> BoundingBox: ...
```

`frame_payload` returns the dict `{"format_version", "constraint_space",
"transform", "bounding_box"}` in that order. Arrays become lists.

`payload_kind(payload, caller_frame)`, in this order:

1. `version = payload.get("format_version")`.
2. `version is None`: if any of `constraint_space`, `transform` or
   `bounding_box` is present, raise `malformed` ("legacy payload has frame
   keys but no format_version"). If `caller_frame` is false, raise
   `legacy_no_transform`. Else return `"legacy"`.
3. `version` is not an `int`, is a `bool`, or is below 1: `malformed`.
4. `version > FORMAT_VERSION`: `unsupported_version`.
5. `version == 1`: `malformed` ("no reader writes format_version 1").
6. Return `"v2"`.

`resolve_frame` returns the `BoundingBox` to attach, or `None` for an identity
frame with no box. `caller_frame` is `local_origin is not None or
local_rotation is not None or declare_identity`.

1. `declare_identity` with `local_origin` or `local_rotation`: `conflict`.
2. Caller frame. For kind `"legacy"` without `declare_identity`, both
   `local_origin` and `local_rotation` are required (answered 2026-11-05, spec
   review of S3: FR-019 forbids a silent identity for a missing part). If one is
   `None`, raise `legacy_no_transform` with the standard message followed by
   ` ({name} is missing)`, where `name` is the missing argument. For kind `"v2"`
   a partial caller frame is compared with the saved transform:
   `origin = zeros(dimensions)` if `declare_identity` or `local_origin is None`,
   else `local_origin`; `rotation = eye(dimensions)` if `declare_identity` or
   `local_rotation is None`, else `local_rotation`; the defaults are only a
   comparison value for step 5. Check with `validate_rigid_rotation` and finite,
   correctly shaped origin. An error is `malformed`.
3. Kind `"legacy"`: if `declare_identity`, return `None`. Else return
   `bounding_box_from_support(support, dimensions, origin, rotation)`.
4. Kind `"v2"`: read `constraint_space` (must equal `"local"`), `transform`
   (keys `dimensions`, `local_origin`, `local_rotation`; `dimensions` must
   equal `dimensions`) and `bounding_box` (a dict with `origin`, `maximum`,
   `nsteps`, or `None`). A missing or wrong item is `malformed`. A
   `RigidRotationError` or `LoopValueError` is `malformed` (raise from the
   original). If `bounding_box` is `None` and the transform is not the
   identity, `malformed`.
5. Kind `"v2"` with a caller frame: the caller transform must equal the saved
   one (`np.allclose(rtol=0, atol=1e-12)` on origin and rotation). If not,
   `conflict`. The message names the saved and the supplied value.
6. Kind `"v2"`: return `None` if `bounding_box` is `None`, else
   `BoundingBox.from_dict({**bounding_box, "local_origin": ..., "local_rotation": ...})`.

`bounding_box_from_support` finds the local box of `support.nodes` (`min` and
`max` on axis 0). It maps the `2**dimensions` corners to world with
`world = corners @ rotation + origin`. It makes
`BoundingBox(dimensions=dimensions, origin=world_min, maximum=world_max)` and
calls `set_local_transform(local_origin=origin, rotation_matrix=rotation)`.
The extent is an envelope. Only `local_enclosure` and display use it, never
evaluation.

Messages (tests match the start of each):

- `legacy_no_transform`: `legacy interpolator payload has no transform: pass local_origin and local_rotation, or declare_identity=True`
- `unsupported_version`: `unsupported interpolator format_version {v}: this reader supports version {FORMAT_VERSION}`
- `malformed`: `malformed interpolator payload: {detail}`
- `conflict`: `conflicting transform: {detail}`. If `declare_identity` is one of the
  conflicting inputs, `detail` contains `declare_identity` (and `saved` when it conflicts
  with a saved transform). Else, for a caller `local_origin` or `local_rotation` against a
  saved transform, `detail` contains `saved`, `supplied` and the name of the item.
- `StaleFrameError`: `the local transform changed after constraints were stored; the stored constraints and coefficients use the previous frame. Call clean() and set the constraints again, or build a new interpolator`

Changes to `GeologicalInterpolator` (`_geological_interpolator.py`):

```python
def frame_signature(self) -> tuple[bytes, bytes]: ...
def _note_frame(self, rows: int) -> None: ...
def _check_frame(self) -> None: ...
def _inherit_frame(self, other: GeologicalInterpolator) -> GeologicalInterpolator: ...
def _restore_local_constraints(self, data: dict) -> None: ...
@classmethod
def from_dict(cls, data, *, local_origin=None, local_rotation=None, declare_identity=False): ...
@classmethod
def from_json(cls, json_str, *, local_origin=None, local_rotation=None, declare_identity=False): ...
@classmethod
def from_yaml(cls, yaml_str, *, local_origin=None, local_rotation=None, declare_identity=False): ...
```

- `frame_signature`: `(local_origin.tobytes(), local_rotation.tobytes())` of
  `self.bounding_box`. If `bounding_box is None`, use `np.zeros(dimensions)`
  and `np.eye(dimensions)` (both `float64`). Spec 003 caches use this value as
  their key.
- `__init__` sets `self._frame_at_data = None`. `clean()` sets it to `None`.
- `_note_frame(rows)`: if `rows > 0` and `self._frame_at_data is None`, set it
  to `self.frame_signature()`.
- `_check_frame()`: if `self._frame_at_data` is not `None` and
  `self.frame_signature() != self._frame_at_data`, raise `StaleFrameError`.
  The comparison is by value, so a restored transform makes the data usable
  again.
- Each `set_*_constraints` calls `self._check_frame()` as the first line of its
  `try` block. It calls `self._note_frame(points.shape[0])` after it stores the
  array. `StaleFrameError` is not wrapped in `ValidationError`.
- `evaluate_value`, `evaluate_gradient` and `to_dict` call `self._check_frame()`
  first. `to_dict` then ends with `payload.update(frame_payload(self.bounding_box,
  self.dimensions))`. `DiscreteInterpolator.to_dict` needs no change, because
  it merges `super().to_dict()`.
- `_inherit_frame(other)`: `other.bounding_box = self.bounding_box`; return
  `other`.
- `_restore_local_constraints(data)`: for each key in `value`, `gradient`,
  `normal`, `tangent`, `interface`, `inequality`, `inequality_pairs`, if the
  array has rows, validate it with the same coercer as the setter
  (`_coerce_value_constraint`; `_coerce_gradient_constraint`, with
  `is_normal=True` for `normal`; `_coerce_interface_constraint`;
  `_coerce_inequality_constraint`; `_coerce_inequality_pair_constraint`).
  Store `array.copy()` in `self.data[key]` with no projection. Set `n_i`
  (`value`), `n_g` (`gradient`), `n_n` (`normal`) and `n_t` (`tangent`) as the
  setters do. Then call `self._note_frame(total_rows)`. The method does not
  convert, so FR-017 holds.
- `from_dict`, `from_json` and `from_yaml` pass the three keyword arguments to
  `InterpolatorFactory.from_dict` only when the caller gave them. Then the
  existing test `test_geological_interpolator_from_dict_delegates_to_factory`
  (its fake takes one argument) passes with no change.

Change to `InterpolatorFactory.from_dict` (`_interpolator_factory.py`):

```python
@staticmethod
def from_dict(d, *, local_origin=None, local_rotation=None, declare_identity=False): ...
```

1. `d = d.copy()`. `caller_frame` as in `resolve_frame`.
2. If `d.get("support")` is not a dict: if `caller_frame`, raise `malformed`
   ("frame arguments apply only to a saved payload"). Else run the existing code
   path with no change (creation call).
3. `kind = payload_kind(d, caller_frame)`. This raises before any work.
4. Pop `format_version`, `constraint_space`, `transform` and `bounding_box`
   from `d` into one dict `frame_keys` (`resolve_frame` reads it). Pop the
   other known keys as the existing code does.
5. Build the support and the interpolator as the existing code does.
6. `interpolator.bounding_box = resolve_frame(frame_keys, kind, support,
   interpolator.dimensions, local_origin=..., local_rotation=...,
   declare_identity=...)`.
7. `interpolator._restore_local_constraints(data)` if `data` is not `None`.
8. Set `c`, `up_to_date` and `valid` as the existing code does. If `c` is set
   and `_frame_at_data` is `None`, call `interpolator._note_frame(1)`. Return.

Every raise happens before the return. The input dict stays unchanged (step 1
copies it, and `resolve_frame` reads only). A caller `BoundingBox` is never
changed.

The three `copy()` methods become:

```python
def copy(self):  # FiniteDifferenceInterpolator
    return self._inherit_frame(FiniteDifferenceInterpolator(self.support))
def copy(self):  # P2Interpolator
    return self._inherit_frame(P2Interpolator(self.support))
def copy(self):  # P2FoldInterpolator
    return self._inherit_frame(P2FoldInterpolator(self.support, fold=self.fold))
```

Their docstrings say "empty copy on the same support, with the same
coordinate frame".

`packages/loop_interpolation/src/loop_interpolation/__init__.py` exports
`InterpolatorFormatError` and `StaleFrameError` (additive). All new names were
free on 2026-11-05 (`git grep -w`: zero matches each).

#### S3 tests

`tests/fixtures/coordinate_frames.py` gets only the helpers that S2 left out. A
helper builds a populated R1 interpolator for engine in `fdi`, `pli`, `p2`,
`cn_fdi`, `cn_pli`, anchor and rotation, in two ways: the exact nodal field
(set `c` from the analytic world field) and the 27-observation fit. A second
helper removes the new keys from a payload (a legacy payload). The task checks
the names that S2 added before it adds any.

`tests/unit/modelling/test_frame_persistence.py` (new):

| Test function | Cells | Checks |
|---|---|---|
| `test_payload_declares_frame` | 3 engines (`fdi`, `pli`, `p2`) x 4 anchors x 3 rotations | version 2, `constraint_space`, `transform` equals the box, `bounding_box` extent; the first new key is `format_version` (U3.8; FR-016) |
| `test_round_trip_world_values` | the same 36 | values within `1e-12*F`, gradients within `1e-12*max(1,\|g\|)` at the 2 world query offsets (U3.1; SC-006) |
| `test_local_arrays_unchanged_after_two_cycles` | the same 36 | stored arrays equal bit for bit; the payloads of cycle 1 and 2 declare the same version, transform and space (U3.2, U3.8; FR-017; SC-006) |
| `test_json_and_yaml_carry_frame` | 3 engines, 1 rotated frame | `to_json`/`from_json` and `to_yaml`/`from_yaml` keep the transform |
| `test_legacy_payload_requires_choice` | 3 engines | `reason == "legacy_no_transform"`, message names `local_origin` and `declare_identity`; the source dict equals its deep copy (U3.6; FR-019; SC-012) |
| `test_legacy_with_supplied_transform` | 3 engines x 4 anchors x 3 rotations | world values and gradients meet SC-006 (U3.7; SC-012) |
| `test_legacy_with_identity_declaration` | 3 engines, anchor 0 rotation 0 | values meet SC-006 |
| `test_conflicting_inputs` | 4 cases | `declare_identity` with `local_origin`; v2 with another origin; v2 with another rotation; v2 rotated with `declare_identity`. `reason == "conflict"`, message names both values |
| `test_malformed_metadata` | 8 cases | missing `transform`; wrong `constraint_space`; scale-2 rotation; NaN origin; wrong shape; `bounding_box` null with a rotated transform; legacy payload with frame keys; `format_version` of type `str`. `reason == "malformed"`; 1 failure; source dict unchanged (SC-012) |
| `test_unsupported_future_version` | versions 3 and 99 | `reason == "unsupported_version"`; nothing returned (SC-017) |
| `test_creation_call_unchanged` | 2 cases | `{"type","boundingbox","nelements"}` builds an empty interpolator; the same call with `declare_identity=True` raises `malformed` |
| `test_old_reader_would_name_format_version` | 1 case | pop the keys that the old reader pops; the first key left is `format_version` (SC-017, OP5) |

`tests/unit/modelling/test_frame_copy_and_lifecycle.py` (new):

| Test function | Cells | Checks |
|---|---|---|
| `test_copy_keeps_frame` | `fdi`, `p2`, `p2_fold` x 12 frames | `copy.bounding_box is source.bounding_box`; the four conversion methods agree with the source; 0 stored rows; the same support; `fold` kept for P2 fold (U3.3; FR-018) |
| `test_stale_frame_rejected_without_state_change` | `fdi`, `pli`, `p2`, `cn_fdi`, `cn_pli` | after the SC-016 change (local origin `O+(3,-2,1)`, 30 degrees about z), `evaluate_value`, `evaluate_gradient`, `set_value_constraints` and `to_dict` each raise `StaleFrameError` once; `data`, `c`, `up_to_date`, `valid`, `n_*` and `latest_diagnostics_report` equal their values before (U3.9; FR-020; SC-016) |
| `test_restored_transform_is_usable` | `fdi` | set the old transform back; evaluation works and equals the first result |
| `test_rebuild_in_new_frame_preserves_world_evaluation` | the same 5 | the documented recipe: build R1 again in the new frame; world values and gradients meet SC-002 (SC-016 acceptance) |
| `test_fault_restoration_uses_current_frame` | R2, FDI and PLI | after a transform change, restoring points through the fault frame raises `StaleFrameError`; after a new build, the restored world points equal a fresh control within `Tp` (SC-016 cache part; no restoration cache exists) |

Pre-fix evidence (FR-031). The names `InterpolatorFormatError`,
`StaleFrameError`, `frame_signature` and the keyword arguments do not exist
before S3, so those tests fail first on the missing name. The task also records
the behavioural defects with the production path and writes them in
`verification.md`: the world error of a rotated round trip, the copy that has
`bounding_box is None`, and the stale evaluation that returns a value after a
transform change. A test that fails only on a missing name does not count as
evidence of a defect.

The existing tests that must pass with no change:
`test_geological_interpolator.py::test_interpolator_json_yaml_round_trip`,
`test_geological_interpolator.py::test_geological_interpolator_from_dict_delegates_to_factory`,
`test_p2_fold_registration.py::test_to_dict_from_dict_round_trip`,
`test_interpolator_factory.py::test_from_dict_builds_interpolator` and
`::test_from_dict_missing_type_raises`, and the S1 and S2 test files. If any
existing test loads a saved payload that has no `format_version`, record it
under "S3 golden review" (C-004) before any change.

#### S3 backend and capability cells

| Cell | Result in S3 |
|---|---|
| Python FDI, PLI/P1, P2: persistence, copy (FDI, P2 and P2 fold have `copy`), transform lifecycle | Mandatory. Tests above. PLI has no `copy`. This is a documented unsupported case that S3 does not add. |
| FDI and PLI constant-norm: transform lifecycle | Mandatory. `cn_fdi` and `cn_pli` in the stale-frame tests. |
| FDI-fold, PLI discrete fold, P2-fold | `copy` for P2 fold only. The fold object is not saved (existing behaviour). |
| Common 2D and 3D geometry | No change. S1 covers it. |
| Native and Surfe | S6. The native adapter does not call `to_dict`. |

#### S3 risks

- A library path that calls `BoundingBox.fit` or `set_local_transform` after
  an interpolator holds data now raises `StaleFrameError`. Before S3 the same
  path gave a silent stale result. The full `pytest` run is the check. A
  failure there is a real defect. Report it; do not weaken the check.
- `evaluate_value` and `evaluate_gradient` now make two small `tobytes` calls
  for each call. Spec 003 measures fault evaluation time. The S3 test task
  runs the spec 003 key-time check if the 003 worktree has one, and records
  the result. Else it records that no measure exists.
- The legacy rule is a breaking change for callers that loaded old files with
  no choice (FR-019 exception, FR-028). It needs the ROADMAP.md decision and
  the note on the version below (a minor release, C-015).

#### S3 constitution check

| Project rule | Result | How S3 agrees |
|---|---|---|
| Work stays within the existing packages | Pass | One new private module in `loop_interpolation`. No new package. |
| Gates `pytest` and `ruff check LoopStructural packages` | Pass | The S3 row of "Stages and rollback" names both. Also run the package tests. |
| numpy docstrings, type hints, py39 | Pass | New functions use numpy docstrings and `from __future__ import annotations`. |
| `ROADMAP.md` owns decisions; update it when one changes | Pass | S3 adds the decision "Legacy interpolator payloads (spec 005, FR-019)" under "Decisions" and a status-log line when S3 completes. |
| The public API stays compatible; a break needs a ROADMAP.md decision and a major-version note | Pass, with the approved exception | Signatures only gain keyword arguments with defaults. Loading a legacy payload with no choice now raises. That is the FR-019 exception (C-003). The ROADMAP.md decision and a `COMPAT.md` row with a major-version note record it. |
| Stages can be tested and reversed independently | Pass | S3 reverts to the S2 merge. |
| Feature documents stay under `specs/NNN-name/` | Pass | Only `plan.md`, `tasks.md` and `verification.md` change there. |

#### S3 files

| File | Change |
|---|---|
| `packages/loop_interpolation/src/loop_interpolation/_frame_persistence.py` | New: `FORMAT_VERSION`, `CONSTRAINT_SPACE`, `REASONS`, `InterpolatorFormatError`, `StaleFrameError`, `frame_payload`, `payload_kind`, `resolve_frame`, `bounding_box_from_support` |
| `packages/loop_interpolation/src/loop_interpolation/_geological_interpolator.py` | `frame_signature`, `_note_frame`, `_check_frame`, `_inherit_frame`, `_restore_local_constraints`; checks in the setters, `evaluate_*` and `to_dict`; `to_dict` frame keys; `from_dict`, `from_json`, `from_yaml` keywords; `clean` resets the frame |
| `packages/loop_interpolation/src/loop_interpolation/_interpolator_factory.py` | `from_dict` loader steps |
| `packages/loop_interpolation/src/loop_interpolation/_finite_difference_interpolator.py` | `copy` keeps the frame |
| `packages/loop_interpolation/src/loop_interpolation/_p2interpolator.py` | `copy` keeps the frame |
| `packages/loop_interpolation/src/loop_interpolation/_p2_fold_interpolator.py` | `copy` keeps the frame |
| `packages/loop_interpolation/src/loop_interpolation/__init__.py` | Export the two errors |
| `tests/fixtures/coordinate_frames.py` | Persistence helpers (only the missing ones) |
| `tests/unit/modelling/test_frame_persistence.py` | New |
| `tests/unit/modelling/test_frame_copy_and_lifecycle.py` | New |
| `tests/fixtures/api_surface_snapshot.json` | Add the new public names, only if the file lists `loop_interpolation` names (additive) |
| `COMPAT.md` | S3 migration notice with the note that the break ships in a minor release (C-015) |
| `ROADMAP.md` | FR-019 decision under "Decisions"; status-log line when S3 completes |
| `specs/005-coordinate-frame-correctness/verification.md` | S3 pre-fix, post-fix, golden review, backend matrix, rollback and run-time sections |

Not changed (the survey in "Later stages" listed them): `_support_factory.py`,
`_fault_segment.py`, `_discrete_interpolator.py` and `_model_serializer.py`.

The `COMPAT.md` row says: loading a legacy interpolator payload (no
`format_version`) assumed an identity frame without a message. It now needs
`local_origin` and `local_rotation`, or `declare_identity=True`, and raises
`InterpolatorFormatError` otherwise. New payloads carry `format_version: 2`.
`copy()` keeps the frame. A transform change after constraints exist raises
`StaleFrameError` at the next use. Version note: by maintainer decision
(C-015) the release that contains S3 is a minor release; `ROADMAP.md`
records the exception to the 1.x rule.

#### S3 amendment: backward-compatible legacy load (C-018, 2026-11-05)

This amendment replaces the statements of "S3 decisions" item 4, of
`payload_kind` step 2, of `resolve_frame` step 3 and of the legacy tests that say
"a legacy payload with no choice raises". Everything else in S3 stays.

- Behaviour on the 1.x line. A legacy payload (no `format_version`) with no
  caller frame (`caller_frame` false) loads. `payload_kind` returns `"legacy"`
  and does not warn; `resolve_frame` returns `None` for this case (the identity
  frame, no bounding box). `InterpolatorFactory.from_dict` calls
  `warn_legacy_identity()` once, as its last step, after the load succeeded, so a
  failed load emits no warning. This is what older releases did. The
  warning text starts with `legacy interpolator payload has no transform; assuming
  the identity frame` and names `local_origin` and `local_rotation`,
  `declare_identity=True`, and the next major release (which raises an error).
- Unchanged: a partial transform, a conflict, a malformed or an unsupported payload
  still raise `InterpolatorFormatError`. These are new keyword arguments and new
  checks, so no old call depends on them.
- One switch and one note. `_frame_persistence.py` gets a module constant
  `_STRICT_LEGACY: bool = False` with the comment `REMOVE AT THE NEXT MAJOR RELEASE
  (spec 005, C-018): delete the fallback, the constant, `warn_legacy_identity` and
  `_external_stacklevel`; a legacy payload with no choice then always raises
  InterpolatorFormatError("legacy_no_transform", ...)`. When the constant is
  true, `payload_kind` raises as it did before this amendment.
- Warning category: `FutureWarning`, because the user must see it; a wrong
  placement is silent. `warn_legacy_identity()` (in `_frame_persistence.py`) computes
  the `stacklevel` so that the warning points at the first frame outside the
  `loop_interpolation` package. This holds for `from_dict`, `from_json`,
  `from_yaml` and a direct `InterpolatorFactory.from_dict` call. A fixed
  `stacklevel` was wrong for three of the four (found in the review of the change).
- Tests. `test_legacy_payload_requires_choice` becomes
  `test_legacy_payload_without_choice_warns_and_loads_identity` (3 engines): it
  expects one `FutureWarning` with the message start above, a loaded interpolator
  with `bounding_box is None`, and an unchanged source dict. New tests
  `test_legacy_warning_points_at_caller` (4 entry points) and
  `test_legacy_warning_not_emitted_when_load_fails` (3 engines) check the
  attribution to the caller and that a failed load does not warn. For `ANCHORS[0]` with
  rotation 0 it also checks the values against the source (SC-006). New test
  `test_strict_legacy_switch_raises` (3 engines): with `_STRICT_LEGACY` set to
  true through `monkeypatch`, the same call raises `InterpolatorFormatError` with
  reason `legacy_no_transform` (the behaviour of the major release). The tests of
  the partial transform, the supplied transform and the identity declaration do not
  change.
- Documents (stage S3 branch, pull request #397): the `ROADMAP.md` decision
  "Legacy interpolator payloads (spec 005, FR-019)" says: on 1.x the loader warns
  and assumes identity; at the next major release it raises. The `COMPAT.md` S3
  row moves from the migration notices to a deprecation: old path "load of a legacy
  payload with no choice", behaviour kept with a warning, earliest removal at the
  next major release. The version note ("minor release", C-015) is removed from
  both. `verification.md` records the new behaviour and corrects the P2 numbers
  (value `7.8e-8`, gradient limit `1e-8`). The text of the pull request changes the
  same way.
- Files: `packages/loop_interpolation/src/loop_interpolation/_frame_persistence.py`,
  `tests/unit/modelling/test_frame_persistence.py`, `ROADMAP.md`, `COMPAT.md`,
  `specs/005-coordinate-frame-correctness/verification.md`.
- The rollback check of S3 stays: a legacy file loads on the S3 code with the old
  path and with a supplied transform; a file written by S3 fails on the old reader
  (OP5).

#### S3 task order (for the spec-tasks skill)

1. `_frame_persistence.py` and the tests of its pure functions (no
   dependency).
2. `GeologicalInterpolator` frame methods, the checks and `to_dict`; then the
   factory loader (it needs 1 and the frame methods).
3. The three `copy()` methods (they need only `_inherit_frame`).
4. Tests: write each test file before the code it checks. The two test files
   and the fixture helpers have separate owners.
5. `COMPAT.md`, `ROADMAP.md`, `verification.md` last. The finish task runs the
   full gates and the rollback check of "Stages and rollback".

S3 rollback check (SC-017). On a throw-away worktree at the S2 merge, save R1
for FDI, PLI and P2 with its `to_dict`. On the S3 code, load each legacy file
with `declare_identity=True` (or the supplied transform) and compare with
SC-006. On the S2 reader, load a file that S3 wrote; expect the `TypeError`
that names `format_version` and no object. Record both in `verification.md`.

### Stage S4 design

S4 makes the display, the data export, the fault-centre geometry and the
intrusion classification use the stored world locations (FR-022; U3.4 and
U3.5; "world-location comparison within `Tp`"). S4 removes two extra
conversions. It also deprecates `GeologicalModel.rescale`: the method keeps its
behaviour and warns that it will be removed (maintainer request, 2026-11-05). It
adds no new public name. Base: branch
`feat/005-s3-frame-persistence` (S3, pull request #397, stacked on #396). If S3
changes before it merges, run this planner again on the merged result.

#### S4 facts from the code (read on 2026-11-05)

- The data of a feature builder is world data. `GeologicalFeatureBuilder.get_value_constraints`,
  `get_gradient_constraints`, `get_tangent_constraints` and `get_norm_constraints`
  (`_geological_feature_builder.py:395` to `466`) read `self.data`, the DataFrame
  that the caller gave. The interpolator converts the rows once, when
  `set_*_constraints` stores them. The old projection of the model data is a
  comment (`geological_model.py:593`).
- `GeologicalFeature.get_data` (`_geological_feature.py:303` to `372`) takes
  those world rows and calls `self.model.rescale(...)` on the location columns in
  four places: gradient, norm, value with `value_map`, value without `value_map`.
  `GeologicalModel.rescale` is `bounding_box.reproject` (local to world,
  `geological_model.py:1103` to `1122`). So it adds the local origin, and it
  applies the rotation, a second time. With a local origin of `(100,200,300)` and an
  identity rotation, the point `(102,203,304)` is shown at `(202,403,604)`. This is
  the display reproduction of the spec. The vectors are taken from the world rows
  and are not converted, so they are already right.
- `get_data` has two callers: `Loop2DView.add_data` (`_2d_viewer.py`) and
  `ModelExporter.save` (`_model_exporter.py:119` and `126`; it writes each object
  with `save`). `StructuralFrame.get_data` (`_structural_frame.py:160`) delegates to
  its three features. No other `get_data` converts a location: a grep for
  `rescale` finds no other call in the feature classes, and
  `IntrusionFeature.get_data` is `pass`.
- Pre-existing, outside S4: in the `value_map` branch the mask is
  `value_constraints[value_constraints == v, :3]`. It compares the whole array,
  not the value column (OP10).
- Intrusion classification. `IntrusionFrameBuilder.set_intrusion_steps_parameters`
  (`intrusion_frame_builder.py:332` to `354`) takes the X, Y and Z columns of the
  intrusion network data and passes `self.model.rescale(xyz, inplace=False)` to
  `self.model.evaluate_model`. The network data is a copy of the user's rows
  (`set_intrusion_frame_c0_data`, line 158), so it is world data.
  `GeologicalModel.evaluate_model(xyz, *, scale=True)` takes world points; `scale`
  is kept only for the signature (`geological_model.py:1212`). So the unit that
  classifies each network point is read at a wrong location. The variable
  `intrusion_network_data_xyz` is used only there.
- Fault centre. `FaultSegment.fault_centre` returns `builder.fault_centre`, a
  world point since S1. `fault_ellipsoid` uses it directly (S1 fixed it).
  `Loop2DView.add_fault_ellipse` (`_2d_viewer.py:211` to `224`) draws the ellipse at
  `center[0], center[1]` with no conversion, which is right. The abutting regions
  receive the world centre as `point` (`_fault_segment.py:463` and `469`);
  `utils/regions.py` was not read (OP9).
- Exporters. The four calls `model.rescale(verts)`, `model.rescale(points)`,
  `model.rescale(xyz)` and `model.rescale(bbox)` in `exporters.py` (lines 92, 366,
  431 and 480) discard the result, because `inplace` is `False` by default. They
  change nothing. The exporters build world inputs (`regular_grid(local=False)` and
  the world bounding box) and `evaluate_model` projects them internally, so the
  written coordinates are world coordinates. The argument `real_coords=False` cannot
  produce model coordinates; it is ignored (pre-existing, outside S4; OP8).
  `ModelExporter.get_block_model` is right: `cell_centres()` is local, `rescale`
  makes it world, and `evaluate_model` projects it again.
- Other `reproject` calls are right. `LoopStructural/interpolators/_api.py:151`,
  `169` and `187` convert stored local interpolator data to world.
  `_geological_feature.py:258`, `266` and `286` do the same for the misfit methods.
  Those methods change the stored arrays in place (`grad /= norm`) and
  `grad[:, :3:6]` selects the position columns; this belongs to S5 (FR-023).
- Callers of `rescale` (grep on 2026-11-05, the whole repository): `exporters.py`
  (4 discarded calls), `_model_exporter.py:77` (`get_block_model`),
  `_geological_feature.py` (4 live calls and 1 comment in `get_data`),
  `intrusion_frame_builder.py:353`, and two calls in
  `examples/3_fault/plot_1_faulted_intrusion.py` (lines 71 and 100). No other
  package calls it. `GeologicalModel.scale` has no caller.
- The example has the same double conversion as `get_data`. It calls
  `model.rescale(model.data[...][["X", "Y", "Z"]].values, inplace=False)` on world
  data of the model, and draws the points with `viewer.add_points`. The prefault
  points appear shifted by the local origin. This is a third display case of FR-022.
- `rescale` and `scale` are stable public API: `@public_api(tier="stable")`
  (`geological_model.py:1104`, `1125`), listed in `API.md:122` and in
  `tests/fixtures/api_surface_snapshot.json:30`. `ROADMAP.md:810` records that they
  "stay as public, signature-stable pass-throughs" to the bounding box. The
  contract test compares only the signature and the tier, so a warning in the body
  does not change the snapshot. `public_api` wraps the function and passes every
  call through.
- The repository has a deprecation pattern. `set_stratigraphic_column`
  (`geological_model.py:649`) calls `warnings.warn(message, DeprecationWarning,
  stacklevel=2)`. The versioning policy of `ROADMAP.md` ("Versioning policy") and the
  "Deprecation lifecycle" of `COMPAT.md` require, for a deprecation on the 1.x line:
  a runtime `DeprecationWarning`, an entry in the `COMPAT.md` table, and a removal
  not before 2 minor releases after the first release that warns. pytest has no
  `filterwarnings` option in `pyproject.toml`, so a warning does not fail a test.
- `tests/unit/io/test_exporters.py` has two stand-in models with a `rescale`
  method (lines 53 and 67). They may be unused after the dead exporter calls go;
  the S4 task checks that the file passes unchanged.
- Survey correction. `_3d_viewer.py`, `intrusion_builder.py` and `_fault_segment.py`
  need no change in S4. `exporters.py` and `_model_exporter.py` change only to stop
  calling the deprecated method (see the decisions).

#### S4 decisions

1. **Contract.** Every `get_data` returns world locations and world vectors. A
   feature builder holds world data, so `get_data` copies it and converts nothing.
2. **Fix 1.** In `GeologicalFeature.get_data`, remove the four calls of
   `self.model.rescale`. Keep every other expression as it is, also the
   `value_map` mask (OP10).
3. **Fix 2.** In the intrusion frame builder, the classification evaluates the model
   at the stored world points. A new private method holds the expression, so a test
   can reach it.
4. **Deprecate `GeologicalModel.rescale`.** The method keeps its behaviour and its
   signature (local to world, `bounding_box.reproject`). It emits one
   `DeprecationWarning` for each call, with `stacklevel=2`, that names the
   replacement `model.bounding_box.reproject(points)`. The tier stays `"stable"`
   until the method is removed. `GeologicalModel.scale` (world to local) is
   deprecated in the same release, with the replacement
   `model.bounding_box.project(points)` (C-016). Removal: not
   before 2 minor releases after the first release that warns, and not before the
   known callers (the QGIS plugin) have moved. The `ROADMAP.md` decision of the
   2026-07-30 status-log entry "2c-3" (line 810: "`scale()`/`rescale()` stay as public,"
   "signature-stable pass-throughs") changes for both methods: S4 records the new decision.
5. **Callers move first.** Our own code must not call the deprecated method, or a
   user would see our warning. The four live calls of `get_data` and the call of the
   intrusion classification go with fixes 1 and 2. The four discarded calls in
   `exporters.py` are removed (they do nothing). `ModelExporter.get_block_model`
   uses `model.bounding_box.reproject(...)`. The two example calls are removed,
   because they apply a second conversion to world data.
6. **No change** to the viewers (`_3d_viewer.py`, `_2d_viewer.py`),
   `_fault_segment.py` and `intrusion_builder.py`. The argument
   `real_coords` of the exporters stays as it is (OP8).

#### S4 contracts

| Boundary | Input | Output | Owner |
|---|---|---|---|
| `GeologicalFeature.get_data` | world rows of the builder | `ValuePoints.locations`, `VectorPoints.locations`: world points; `VectorPoints.vectors`: world directions | `_geological_feature.py` |
| Viewer and exporter of data | objects from `get_data` | drawn or written at their `locations` | `Loop2DView.add_data`, `ModelExporter.save` (unchanged) |
| Intrusion classification | world X, Y, Z of the network data | stratigraphic ids from `model.evaluate_model` at those points | `IntrusionFrameBuilder._network_model_values` |
| Fault-centre geometry | `FaultSegment.fault_centre` (world) | ellipsoid or ellipse at that point | `fault_ellipsoid`, `Loop2DView.add_fault_ellipse` (unchanged) |
| `GeologicalModel.rescale` (deprecated) | local points | world points, as before; one `DeprecationWarning` | `geological_model.py`; replacement `bounding_box.reproject` |
| `GeologicalModel.scale` (deprecated) | world points | local points, as before; one `DeprecationWarning` | `geological_model.py`; replacement `bounding_box.project` |

#### S4 functions and signatures

`GeologicalFeature.get_data(self, value_map: dict | None = None)`
(`_geological_feature.py`): the four location arguments become
`gradient_constraints[:, :3]`, `norm_constraints[:, :3]`,
`value_constraints[value_constraints == v, :3]` and `value_constraints[:, :3]`.
The docstring says: "The locations and the vectors are in world coordinates, as
the builder stores them." No other line changes. The commented-out block at the
end of the method stays as it is.

New private method of `IntrusionFrameBuilder` (`intrusion_frame_builder.py`):

```python
def _network_model_values(self, data: pd.DataFrame) -> np.ndarray:
    """Return the stratigraphic id of each network point.

    Parameters
    ----------
    data : pandas.DataFrame
        Network data with the world columns ``X``, ``Y`` and ``Z``.

    Returns
    -------
    numpy.ndarray
        The id of the unit at each point. The model is evaluated at the world
        points; no conversion is applied.
    """
    xyz = data.loc[:, ["X", "Y", "Z"]].to_numpy(dtype=float)
    return self.model.evaluate_model(xyz)
```

In `set_intrusion_steps_parameters`: delete the assignment of
`intrusion_network_data_xyz` (lines 333 to 335, three lines; the variable has no
other use), and replace the statement at line 352 with

```python
intrusion_network_data.loc[:, "model_values"] = self._network_model_values(
    intrusion_network_data
)
```

The name `_network_model_values` was free on 2026-11-05 (`git grep -w`: zero
matches).

`GeologicalModel.rescale` (`geological_model.py`): the signature
`rescale(self, points: np.ndarray, *, inplace: bool = False) -> np.ndarray` and the
decorator `@public_api(tier="stable")` do not change. The body becomes:

```python
warnings.warn(
    "GeologicalModel.rescale is deprecated and will be removed in a future "
    "release; use model.bounding_box.reproject(points) instead",
    DeprecationWarning,
    stacklevel=2,
)
return self.bounding_box.reproject(points, inplace=inplace)
```

The module already imports `warnings` (line 649). The docstring gets a
`.. deprecated::` note that names `bounding_box.reproject` and says "Deprecated by
spec 005, stage S4". `GeologicalModel.scale` (`geological_model.py`) gets the same change. The signature
`scale(self, points: np.ndarray, *, inplace: bool = False) -> np.ndarray` and the
decorator do not change. The body becomes:

```python
warnings.warn(
    "GeologicalModel.scale is deprecated and will be removed in a future "
    "release; use model.bounding_box.project(points) instead",
    DeprecationWarning,
    stacklevel=2,
)
return self.bounding_box.project(np.array(points).astype(float), inplace=inplace)
```

The other callers change as follows:

- `exporters.py`: delete the four statements `model.rescale(verts)`,
  `model.rescale(points)`, `model.rescale(xyz)` and `model.rescale(bbox)`. If
  `real_coords` has then no use in a function, keep the argument (OP8).
- `_model_exporter.py`, `get_block_model`: replace
  `model.rescale(model.bounding_box.cell_centres())` by
  `model.bounding_box.reproject(model.bounding_box.cell_centres())`.
- `examples/3_fault/plot_1_faulted_intrusion.py`: in the two `viewer.add_points`
  calls, pass `model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].values`
  directly (world points) and delete the `model.rescale(...)` wrapper and its
  `inplace=False`.

#### S4 tests

`tests/fixtures/coordinate_frames.py` gets one helper:
`world_foliation_model(anchor, rotation_deg) -> tuple[GeologicalModel, GeologicalFeature, pd.DataFrame]`.
It makes `model_with_frame(anchor, anchor + (8, 10, 12), rotation_deg)`, a data
frame in world coordinates with the columns `X, Y, Z, feature_name, val, nx, ny, nz, gx, gy, gz, w`
(27 value rows from `r1_observations(anchor)`, 4 gradient rows and 4 normal rows at
`anchor + (2, 3, 4)`, `(5, 6, 7)`, `(3, 8, 5)` and `(6, 4, 9)` with the unit
direction `R1_GRADIENT / norm`), and the feature
`model.create_and_add_foliation("r1", ...)` built without a solve. It returns the
model, the feature and the data frame. The test task reads
`tests/integration/test_interpolator.py` and `tests/unit/interpolator/test_outside_box.py`
for the call pattern of `create_and_add_foliation` and copies it.

`tests/unit/modelling/test_world_location_consumers.py` (new):

| Test function | Cells | Checks |
|---|---|---|
| `test_feature_data_locations_are_world` | 4 anchors x 3 rotations = 12 | `feature.get_data()` gives a `ValuePoints` (value), and `VectorPoints` (gradient and norm). Each `locations` equals the input rows within `tp(L, W)`; each `vectors` equals the input vectors within `1e-12` (U3.4) |
| `test_spec_example_observation_is_not_shifted` | `ANCHORS[1]` x 3 rotations = 3 | the observation `(102,203,304)` is at `(102,203,304)` and not at `(202,403,604)` |
| `test_exported_data_keeps_world_locations` | 4 anchors, rotation 30 = 4 | each object of `get_data()` is saved with `save` to a CSV file (the task reads `_point.py` `save` first for the supported extension), and the `x`, `y`, `z` columns equal the input within `tp(L, W)` (U3.4, export) |
| `test_intrusion_classification_uses_world_points` | 4 anchors x 3 rotations = 12 | `IntrusionFrameBuilder._network_model_values(data)` is called on `object.__new__(IntrusionFrameBuilder)` with `.model` set to a real rotated `GeologicalModel` whose `evaluate_model` is replaced by a spy; the spy receives exactly the world X, Y, Z (`np.array_equal`) (U3.5) |
| `test_2d_viewer_fault_ellipse_centre_is_world` | 2 anchors x 2 rotations (0 and 30) = 4 | `Loop2DView.add_fault_ellipse` is called with a stand-in `self` that has `ax` (a Matplotlib `Agg` axis) and `model`; the patch centre equals `fault.fault_centre[:2]` within `tp(L, W)` (U3.5, coverage, no defect). The fault comes from `_fault(anchor, rotation, "PLI")` of `test_fault_frame_analytic.py` |
| `test_rescale_warns_and_keeps_behaviour` | 4 anchors x 3 rotations = 12 | `model.rescale(points)` raises one `DeprecationWarning` (`pytest.warns`, the message contains `bounding_box.reproject`), the warning is attributed to the test file (`stacklevel`), the result equals `model.bounding_box.reproject(points)` within `tp(L, W)`, and `inplace=True` changes the passed array like `reproject` does (SC-014, FR-028) |
| `test_scale_warns_and_keeps_behaviour` | 4 anchors x 3 rotations = 12 | `model.scale(points)` raises one `DeprecationWarning` (the message contains `bounding_box.project`), attributed to the test file, with a result equal to `model.bounding_box.project(points)` within `tp(L, W)` |
| `test_internal_callers_do_not_warn` | 1 | with the filter `warnings.filterwarnings("error", message=".*(rescale|scale) is deprecated.*")`: `feature.get_data()`, `ModelExporter.get_block_model(model)` and `IntrusionFrameBuilder._network_model_values(data)` run with no warning; the exporter functions that `tests/unit/io/test_exporters.py` already calls run with no warning (the task reads that file for the call pattern) |

Pre-fix evidence (FR-031). The first four tests fail on the current tree because
of the extra conversion: the first, second and third by a location shifted by the
local origin (`(202,403,604)` for the example), and the fourth because the spy
receives the shifted points (the method does not exist yet, so the first failure is
`AttributeError`). Record also a behavioural measure with the production path: the
shift of the example point for `ANCHORS[1]` and each rotation. The fifth test
passes before the fix; record it as "coverage, no defect".

The existing tests that must pass with no change:
`tests/unit/modelling/intrusions/test_intrusions.py`,
`tests/unit/io/test_exporters.py`, `tests/unit/test_public_api_contract.py` (the
signature and the tier of `rescale` do not change), and the S1 to S3 test files. The intrusion tests build a
model. If a translated frame changes an expected value, record the change under
"S4 golden review" (C-004) with the defect and its effect, before any change.
With an identity frame (origin 0, rotation 0) nothing changes.

#### S4 backend and capability cells

| Cell | Result in S4 |
|---|---|
| Python FDI foliation data (display and export of observations) | Mandatory. The tests above use `FDI`. |
| Fault centre in the 2D viewer | Coverage test above (PLI fault). |
| Intrusion classification | Mandatory. Spy test; the intrusion model tests stay unchanged. |
| Native, Surfe, browser | Not in S4. The viewer and exporters do not call them. |

#### S4 risks

- Intrusion models in a translated frame classified their network data at wrong
  points. The corrected classification can change the result of such a model
  (a golden result). The intrusion tests decide: record each change under C-004.
- Users who corrected the display by shifting their data now see a shifted
  display. This is the corrected behaviour (`COMPAT.md` notice).
- `GeologicalModel.rescale` keeps its meaning (local to world) and now warns. The
  S4 grep (`git grep -n "\.rescale("`) in the finish task must show only the
  commented block of `get_data` and the tests of the deprecation.
- A downstream caller (the QGIS plugin, user scripts) gets a `DeprecationWarning`.
  Python hides it outside `__main__` and test runs. The message names the
  replacement. The removal waits for the lifecycle of `COMPAT.md`.

#### S4 constitution check

| Project rule | Result | How S4 agrees |
|---|---|---|
| Work stays within the existing packages | Pass | Two source files change, in `LoopStructural`. No new module. |
| Gates `pytest` and `ruff check LoopStructural packages` | Pass | The S4 row of "Stages and rollback" names both. |
| numpy docstrings, type hints, py39 | Pass | The new method has a numpy docstring and type hints. |
| The public `GeologicalModel` API stays compatible; a break needs a ROADMAP.md decision and a major-version note | Pass | No signature, tier or behaviour breaks. `rescale` is deprecated on the 1.x line with a runtime warning, as the versioning policy and the lifecycle of `COMPAT.md` require; `scale` and `evaluate_model` do not change. `get_data` keeps its signature; its locations become correct (known meaning). Removal is a later, separate change. |
| Stages can be tested and reversed independently | Pass | S4 reverts to the S3 merge. |
| `ROADMAP.md` owns decisions; update it when one changes | Pass | The decision of line 810 ("`scale()`/`rescale()` stay as public pass-throughs") changes for `rescale`. S4 adds the decision "Deprecation of GeologicalModel.rescale and scale (spec 005, S4)" under "Decisions". The "005 complete" line stays in S6. |
| Feature documents stay under `specs/NNN-name/` | Pass | Only `plan.md`, `tasks.md` and `verification.md` change there. |

#### S4 files

| File | Change |
|---|---|
| `LoopStructural/modelling/core/geological_model.py` | `rescale` and `scale`: warn with `DeprecationWarning`; docstrings |
| `LoopStructural/modelling/features/_geological_feature.py` | `get_data`: remove four `self.model.rescale` calls; docstring |
| `LoopStructural/modelling/intrusions/intrusion_frame_builder.py` | New `_network_model_values`; use it in `set_intrusion_steps_parameters`; delete `intrusion_network_data_xyz` |
| `LoopStructural/export/exporters.py` | Delete the four discarded `model.rescale` calls |
| `LoopStructural/modelling/core/_model_exporter.py` | `get_block_model`: use `bounding_box.reproject` |
| `examples/3_fault/plot_1_faulted_intrusion.py` | Pass the world points to `viewer.add_points`; remove two `model.rescale` wrappers |
| `tests/fixtures/coordinate_frames.py` | New helper `world_foliation_model` |
| `tests/unit/modelling/test_world_location_consumers.py` | New |
| `COMPAT.md` | S4 migration notice; one row in the deprecation table for each of `GeologicalModel.rescale` and `GeologicalModel.scale` |
| `ROADMAP.md` | Decision "Deprecation of GeologicalModel.rescale and scale (spec 005, S4)" under "Decisions" |
| `API.md` | Line 122: mark `rescale` as deprecated, with the replacement |
| `specs/005-coordinate-frame-correctness/verification.md` | S4 pre-fix, post-fix, golden review, backend matrix, rollback and run-time sections |

Not changed (the survey in "Later stages" listed them): `_3d_viewer.py`,
`intrusion_builder.py` and `_fault_segment.py`.

The deprecation rows of `COMPAT.md` (one for `rescale`, one for `scale`) follow
the columns of that table. The `rescale` row: old path
`GeologicalModel.rescale(points, inplace=False)`; new path
`GeologicalModel.bounding_box.reproject(points, inplace=False)`; shim added
`2026-11-05` (the date of the run); first released in "TBD (next release containing
this warning)"; introduced by the S4 commit; earliest removal "TBD (2 minor
releases after the first release that warns)"; removed in empty; status
"Active"; owner "Core maintainers". The `scale` row has the old path
`GeologicalModel.scale(points, inplace=False)` and the new path
`GeologicalModel.bounding_box.project(points, inplace=False)`; the other columns
are the same. The `API.md` line reads: "Geometry: `scale` (deprecated, use
`bounding_box.project`), `rescale` (deprecated, use `bounding_box.reproject`),
`regular_grid`, `bounding_box` (property)".

The `COMPAT.md` row says: `GeologicalFeature.get_data()` (used by
`Loop2DView.add_data` and `ModelExporter.save`) returned locations moved by the
local origin and rotated a second time, and the intrusion classification evaluated
the model at such points. Both now use the stored world locations. With an
identity frame nothing changes.

#### S4 task order (for the spec-tasks skill)

1. The fixture helper (no dependency on code changes).
2. The test file, in three parts in one file, one task after the other: display
   and export tests; the intrusion test; the 2D viewer test.
3. The pre-fix evidence in `verification.md` (after all tests, before any fix).
4. The tests of the deprecation (written before it; they must fail first), in the
   same test file, one task after the others.
5. The fixes and the deprecation in parallel (they are in separate files):
   `get_data`; the intrusion helper and its call site; the `rescale` warning; the
   exporter files; the example. The tests of step 4 pass only when all callers
   have moved, so the finish task of the group runs them together.
6. The golden review of the intrusion tests, `COMPAT.md`, `ROADMAP.md`, `API.md`,
   the rollback check, and last the gates (the orchestrator runs the full suites
   when memory is tight).

S4 rollback check. Revert to the S3 merge. Then: (a) the import
`from loopstructuralvisualisation import Loop2DView, Loop3DView` works and
`inspect.signature(Loop2DView.add_data)` is equal before and after; (b) on the S3
code, `get_data()` shows the spec example at `(202,403,604)`, on the S4 code at
`(102,203,304)` (recorded in `verification.md`); (c) with an identity frame the
two codes give the same `get_data()` output; (d) on the S3 code `model.rescale(p)`
and `model.scale(p)` give no warning, on the S4 code each gives one, and the values
are equal.

### Stage S5 design

S5 has four aims. First, S5 makes the remaining resize paths keep the node counts
and cell counts consistent (FR-014). Second, S5 checks that cell-centre sampling
agrees with the grid cells (FR-015). Third, S5 makes the feature misfit
diagnostics leave their inputs unchanged and compare quantities in one basis
(FR-023). Fourth, S5 adds a residual and exclusion report that keeps physical
residuals and assembled residuals apart (FR-024).

Base: branch `feat/005-s4-world-locations` (S4, pull request #401, not merged). If
S4 changes before it merges, run this planner again on the merged result.

The plan was revised on 2026-11-05 after an independent review and five
maintainer answers (see "Answers"). The text below is the revised design.

#### S5 facts from the code (read on 2026-11-05)

- Resize, 3D. S1 repaired the `origin` and `maximum` setters of
  `StructuredGrid3DGeometry` with `_cells_for_length`. Two setters of the same
  class were not repaired (`_structured_grid_3d.py`, `nsteps` setter at line 71,
  `step_vector` setter at line 131). The `nsteps` setter counts nodes. It does
  `step_vector /= nsteps_new / nsteps_old`. This keeps the product
  `nsteps * step` constant. It does not keep the extent `(nsteps - 1) * step`.
  Example: the node count 11 to 21 with step 1 gives step `0.524` and a new
  extent `10.48` instead of `10`.
- The `step_vector` setter sets the node count to `ceil(nsteps / factor)`. It
  treats the node count as a cell count. So the extent also drifts. It can lose
  a cell.
- `set_nelements` (`_structured_grid_3d.py:50-57`) ends with `self.nsteps =
  nsteps`. It computes `nsteps = ceil((maximum - origin) / step)` with a cubic
  step from the requested element count. So it uses the `nsteps` setter, and the
  setter fix changes what `set_nelements` leaves behind. Callers:
  `GeologicalFeatureBuilder.build` calls `interpolator.set_nelements` (line 575
  of `_geological_feature_builder.py`) and the support method is reached through
  `DiscreteInterpolator.set_nelements` (line 104). `set_interpolation_geometry`
  runs `support.step_vector = support.step_vector * 0.9` in a loop while
  `interpolator.dof < 100` (line 525 of the same file). So both setters run in
  ordinary model builds, and numeric tests that depend on a built support can
  change.
- The golden test `packages/loop_common/tests/test_structured_support_reference.py`
  runs `set_nelements(200)`, then `nsteps + 1`, then `step_vector * 2`, then
  `origin - 1`, then `maximum + 3`, in sequence on one grid (lines 171 to 180).
  Each case records the state that the previous setters left. So the cases
  `after_set_nelements`, `after_nsteps`, `after_step_vector`, `after_origin` and
  `after_maximum` all change when the two setters change. The earlier text of
  this plan said that `after_origin` and `after_maximum` stay: that was wrong.
- `StructuredGrid2DGeometry` (`_structured_grid_2d.py`) has plain attributes
  `nsteps`, `maximum`, `nsteps_cells`, `n_cell_x` and `n_cell_y`, set once in
  `__init__` (lines 23 to 30). `maximum = origin + nsteps * step_vector` counts
  nodes as cells. `nodes` uses `origin + nsteps_cells * step_vector`. The two
  disagree by one step. The three count attributes and `maximum` go stale if
  `nsteps` is assigned later. The 3D class computes `maximum` as
  `origin + nsteps_cells * step_vector`. The 2D class `set_nelements` raises
  `NotImplementedError` (`supports/_2d_structured_grid.py:50`). Four existing
  assertions read the 2D `maximum`: `packages/loop_common/tests/test_2d_discrete_support.py`
  lines 14 and 22, and `tests/unit/interpolator/test_2d_discrete_support.py`
  lines 16 and 24. They expect `origin + nsteps * step` (the legacy value).
- `RectilinearGrid` (`supports/_3d_rectilinear_grid.py`) holds node arrays
  `_xnodes`, `_ynodes`, `_znodes`. It inherits the `origin`, `nsteps` and
  `step_vector` setters of the 3D geometry. A call changes the stored counts and
  the stored origin, but not the node arrays. So the geometry and the nodes
  disagree. `maximum` is read from the node arrays and has no setter. The
  constructors of the geometry classes write the private attributes (`_origin`,
  `_nsteps`, `_step_vector`), not the setters.
- Cell centres. `BoundingBox.cell_centres` (`_bounding_box.py:647`) builds
  `nsteps - 1` centres with `linspace(origin + h, maximum - h, nsteps - 1)` and
  `h = node_step_vector / 2`. Now `maximum - origin = (nsteps - 1) * step`. So
  the spacing is exactly `node_step_vector`. Every centre is the middle of a cell
  of `regular_grid(local=False)`. The method returns local points. All four
  callers (`_base_geological_feature.py:391`, `420`, `442` and
  `_model_exporter.py:77`) call `reproject` on the result. The existing tests
  (`test_bounding_box_structured_grid.py`, `test_scalar_field_grid.py`) use
  unrotated boxes only. By reading, the default path has no defect for a rotated
  box (OP12). `structured_grid(local_coordinates=True)` builds a different grid.
  Its cells are the enclosure of the projected corners. That grid works in the
  local frame, so `grid.cell_index_at` takes local points. Its `cell_centres` are
  the centres of those cells. They are not `box.cell_centres()`. The supports'
  `cell_centres(global_index)` (3D, 2D, rectilinear) use their own local frame and
  do not apply `rotation_xy`, as their docstrings state.
- Misfit methods (`_geological_feature.py`: `evaluate_gradient_misfit` at line 240
  and `evaluate_value_misfit` at line 272 on this branch). `get_gradient_constraints()`
  and `get_norm_constraints()` return the stored arrays by reference
  (`_geological_interpolator.py`: the definitions are at lines 633 and 652, and
  the returns are at lines 640 and 660). So `grad /= norm` (line 255) changes
  `interpolator.data["gradient"]` in place. The norm is taken over all columns
  (position, vector and weight). So the position and the weight columns change
  too. `grad[:, :3:6]` (line 260) is the slice with start 0, stop 3 and step 6.
  That is column `0` only. The vector columns are `3:6`. The stored vectors are
  in the local basis. The interpolator converts them once when it stores them.
  `evaluate_gradient` returns world vectors. So a rotated frame compares two
  bases.
- Evaluation path. `GeologicalFeature.evaluate_value` (line 117) and
  `evaluate_gradient` call `self.builder.up_to_date()` (lines 141 and 191).
  `BaseBuilder.up_to_date` (`_base_builder.py:81`) first calls itself on every
  fault builder. Then, if the builder or the interpolator is stale, it calls
  `update()`, which builds again and solves. So a misfit method that calls
  `self.evaluate_value` on a stale model rebuilds it. This also holds for
  `self.max()` and `self.min()` (`_base_geological_feature.py:280` and `294`),
  because they call `evaluate_value`. The earlier text of this plan marked this
  path as not read. It is now read. No separate cache-validity flag exists in
  `GeologicalFeature`, the builder or the interpolator (grep for `cache` on
  2026-11-05). The state is `interpolator.up_to_date`, `builder._up_to_date`
  (a protected attribute; no public property exists), `interpolator.c` and
  `interpolator.latest_diagnostics_report`.
- Result shape. With one family, `np.array([list])` gives the shape `(1, N)`.
  With two families of equal length it gives `(2, N)`. With two families of
  different length, `np.array([list, ndarray])` raises `ValueError` on current
  NumPy. The value misfit divides by `self.max() - self.min()` (line 288): zero
  range gives a division by zero. It does not change its inputs. No code in
  `LoopStructural` or `packages` calls the two methods (`git grep -n -w` on
  2026-11-05). `ROADMAP.md` line 832 names them in a history note only.
- Diagnostics report. `_diagnostics.py` holds three frozen dataclasses and
  `ConstraintDiagnosticsReport` with `to_dict` and `summary`. It has
  `from __future__ import annotations`. The two builders
  (`GeologicalInterpolator._build_constraint_diagnostics_report`, line 731, and the
  `DiscreteInterpolator` override, line 221) only read state.
  `get_constraint_diagnostics_report(refresh=False)` stores the result in
  `latest_diagnostics_report`.
- `DiscreteInterpolator.calculate_residual_for_constraints()` (line 497) gives
  `A x - b` for each assembled constraint, as a dict by constraint name. The names
  are the constraint names: FDI uses `value` and `gradient`, PLI uses `value` and
  `norm`. It does not test `self.c`. The constructor sets `self.c` to zeros
  (lines 64 to 68), so `c` is not `None` after construction. Line 1355 can
  assign another value.
- Dropped rows. The `DiscreteInterpolator` builder sets `dropped_rows` as
  `source_points - rows` for `value`, `normal`, `tangent` and `inequality_value`
  (line 265), and as `2 * source_points - rows` for `gradient` (line 267). Other
  families (interface, inequality pairs, regularisation) get `None`. So
  `dropped_rows` is in matrix rows. `outside_model_point_count` is in points. The
  base builder (used by `SurfeRBFInterpolator`, which does not override it) sets
  `dropped_rows = outside` and `row_count = point count`, so for Surfe both are in
  points.
- Builder exclusions. `GeologicalFeatureBuilder.excluded_rows`
  (`dict[str, dict[str, int]]`, S1, line 88) uses constraint names as keys:
  `<feature name>_orthogonal` (line 343) and `<feature name>_equality` (lines 387
  and 393). The reasons are `invalid_reference_value`,
  `invalid_reference_direction` and `evaluation_error`. The S1 restricted-support
  case uses the feature name `restricted`, so the keys are `restricted_equality`
  and `restricted_orthogonal`. These keys are not family names. No report reads
  them.
- Python version. `pyproject.toml` and `packages/loop_common/pyproject.toml` set
  `requires-python >=3.9`. `_structured_grid_3d.py` has no `from __future__
  import annotations`. The syntax `np.ndarray | None` fails there on Python 3.9.
- Documents. `git grep` finds no mention of `ConstraintDiagnosticsReport`,
  `summarise_residuals`, `UNKNOWN_EXPLANATION` or the misfit methods in
  `docs/source` or `API.md` (checked 2026-11-05).

#### S5 decisions

1. **Resize contract.** For every structured grid, each assignment of `origin`,
   `maximum`, `nsteps` or `step_vector` keeps the extent that the caller asked
   for, or encloses it by less than one cell. `nsteps` counts nodes in the
   geometry classes. A repeated assignment of the same value changes nothing.
   Assignments of different values can add up. The sequence A, B, A can enlarge
   the extent by less than one cell of B plus less than one cell of A. A test
   measures this (`test_step_vector_a_b_a_extent_drift`).
2. **3D `nsteps` setter.** Keep origin and extent. The new step is
   `extent / (nsteps_new - 1)`. The extent is `nsteps_cells * step_vector`, taken
   before the change. A value below 2 on any axis raises `ValueError`. A value
   that is not a whole number raises `ValueError`. The setter does not round,
   because rounding hides a caller error. A float array of whole numbers is
   accepted and stored as `int`. A scalar is broadcast, as before.
3. **`set_nelements` (Q4, answered 2026-11-05).** S5 accepts the exact extent
   after `set_nelements`: `maximum` is the same before and after the call. The old
   code changed it. `set_nelements` does `nsteps = np.maximum(np.ceil(...).astype(int), 2)`.
   The clamp keeps behaviour for valid callers (a value of 2 or more is unchanged)
   and avoids a new `ValueError` when the requested element count is very small.
   The old code gave one node (zero cells) for such a request, which no caller
   can use. The setter accepts only `nsteps >= 2`, because one node has no
   extent to keep.
4. **3D `step_vector` setter.** Keep origin and the requested extent. The cell
   count is `_cells_for_length(extent, step_vector=new_step)` and the node count
   is that plus one. The stored step is a copy of the value that the caller gave.
   The grid can end up to one cell larger than the extent, never smaller. The loop
   in `set_interpolation_geometry` ends as before, because the cell count grows
   when the step gets smaller.
5. **2D grid (Q1, answered 2026-11-05).** `nsteps` becomes a property. Its setter
   coerces the value with `np.ceil(np.asarray(v)).astype(int)`, as `__init__` does
   today. `nsteps_cells`, `n_cell_x` and `n_cell_y` become read-only properties
   derived from `nsteps`. `maximum` becomes a read-only property with the legacy
   value `origin + nsteps * step_vector`. This value counts nodes as cells. S5 does
   not change it, because FR-028 asks to keep it. The four existing assertions
   of `maximum` must pass with no change. After an `nsteps` assignment, the three
   counts and `maximum` follow the new `nsteps` (before, they stayed stale).
6. **Rectilinear grid.** The `origin`, `nsteps` and `step_vector` setters raise
   `NotImplementedError` after construction. Resize needs new node arrays, and
   the grid cannot guess them. `maximum` already has no setter.
7. **Cell centres (FR-015).** No change to `BoundingBox.cell_centres` unless a new
   test fails (OP12). S5 adds tests for rotated boxes and for the cell
   association of the returned grid. If a test fails, the fix is the closed form
   `origin + (arange(nsteps - 1) + 0.5) * node_step_vector` for each axis, in
   place of `linspace`. The task records the failure first.
8. **Misfit methods (Q2, answered 2026-11-05).** The old public methods keep their
   names and their public result shape. Two new methods return the corrected
   result as a 1-D array: `evaluate_gradient_cosines` and `evaluate_value_residuals`
   (both names are free; `git grep -n -w` on 2026-11-05 gives no match). The old
   methods and the new methods call the same private code (one code path for each
   quantity). That code works on copies, converts the stored local vectors with
   `interpolator.local_to_world_vectors`, converts the stored local positions with
   `bounding_box.reproject`, selects the columns `[:, :3]` (position) and
   `[:, 3:6]` (vector), and evaluates through the feature. The in-place change, the
   `grad[:, :3:6]` column error and the local-versus-world comparison are removed
   from the old methods.
9. **Old result shape.** `evaluate_gradient_misfit` returns one row for each
   non-empty family, as before: shape `(1, N)` for one family, and `(2, N)` for two
   families of equal length. For two families of different length it returns one
   row `(1, N1 + N2)`; before, this case raised `ValueError` (OP16). With no rows it
   returns `np.array([])` of shape `(0,)`, as before. `evaluate_value_misfit`
   returns shape `(N,)`, as before.
10. **Gradient value changes meaning.** The old value was the dot product of an
    unnormalised model gradient and a column slice. The new value is the cosine
    between the world data direction and the model gradient. This is a change of
    behaviour of a public method. S5 records it in `COMPAT.md` as a migration
    notice, and in `ROADMAP.md` as a decision entry (the constitution row for the
    public API asks for it). The new name `evaluate_gradient_cosines` gives the
    same numbers in one 1-D array.
11. **Stale model (Q3, answered 2026-11-05).** The misfit methods and
    `verification_report` must not change state and must not rebuild. First they
    ask `_stale_reason()`. It returns a text if `builder._up_to_date` is `False`, if
    any fault builder of the builder is stale, or if `interpolator.up_to_date` is
    `False`. It returns `None` otherwise (and when `builder` is `None`, only the
    interpolator counts). If the text is not `None`, the method does not evaluate.
    It writes one `warnings.warn(text, UserWarning, stacklevel=2)`. The text names
    the stale object and says that the model must be built first. The new methods
    then return `np.full(n, np.nan)`. The old methods return their usual shape
    filled with NaN. Reasons: (a) the call changes no state (SC-009); (b) no
    exception stops a loop over features; (c) NaN is the marker for an unusable row
    in this design; (d) a number from old coefficients looks valid and is not;
    (e) a `UserWarning` shows by default, a log line does not. When the model is up
    to date, the methods evaluate with `self.evaluate_value` and
    `self.evaluate_gradient`. Their call to `builder.up_to_date()` then does
    nothing, because `_stale_reason()` has checked the same state (faults
    included). The number `n` is the row count of the stored array. Reading a shape
    changes nothing.
12. **Gradient misfit is a cosine.** The pure function `direction_cosines` returns
    the cosine between the data direction and the model gradient for each row. A
    row with a zero, not finite or missing vector gives `NaN`.
13. **Value misfit (Q5, answered 2026-11-05).** The public `evaluate_value_misfit`
    keeps its old meaning: `|value - field|` divided by the field range
    `self.max() - self.min()`. It divides only when the range is finite and greater
    than zero. Otherwise it returns the unnormalised absolute difference. The
    physical residual uses a private helper `_value_residuals_fresh()` that gives the
    unnormalised `|value - field|` in field units. `evaluate_value_residuals`,
    `evaluate_value_misfit` and the physical report entry `value` all call it. So
    the unit `"field"` is true.
14. **Report.** `ConstraintDiagnosticsReport` gets three new fields with defaults,
    appended after the existing ones. `to_dict` and `summary` show them only when
    they are not empty. So a report that was built before S5 gives the same
    output. A new method `GeologicalFeature.verification_report()` fills them. It
    builds a new report. It does not store it in `latest_diagnostics_report`. Its
    private call is `_build_constraint_diagnostics_report`.
15. **Separate residuals.** `physical_residuals` holds quantities in data units
    (value: field units; direction: degrees). `assembled_residuals` holds the
    `A x - b` rows of each constraint, in the assembler's normalisation. The two
    dicts keep the plain names as keys. Assembled keys are constraint names (FDI:
    `value`, `gradient`; PLI: `value`, `norm`). Physical keys are `value`, `gradient`
    and `normal`. The two sets can share a key. S5 does not prefix the keys,
    because the assembled keys must stay equal to the names that
    `calculate_residual_for_constraints` uses. Each entry has a `kind` field
    (`"physical"` or `"assembled"`). The two dicts are separate objects and never
    merge in `to_dict`. The tests check the `kind` field and the object identity,
    not a claim that the keys differ.
16. **Unsolved model.** `verification_report` fills the two residual blocks only if
    `_stale_reason()` is `None` and `interpolator.c` is not `None` and the
    interpolator has `calculate_residual_for_constraints`. If not, both blocks stay
    empty (`{}`), one warning is written (stale case), and the exclusions are still
    reported, because they read state only. `calculate_residual_for_constraints`
    itself does not change in S5.
17. **Exclusions.** Each exclusion is an `ExclusionRecord` with a family, a reason,
    a count and a plain explanation. Sources, in this order:
    (a) The builder's `excluded_rows`. A builder record keeps its own constraint key
    as `family` (for example `restricted_orthogonal`). It is not a data family. S5
    does not subtract it from any family count.
    (b) The outside-model count of each family, reason `outside_model`, in points.
    (c) For each family whose `dropped_rows` is not `None`: the rows that
    `outside_model` does not explain, reason `unattributed`, in rows. The unit of
    `dropped_rows` is the matrix row, so the outside count is first multiplied by
    `rows_per_point` (2 for `gradient` when the interpolator is a
    `DiscreteInterpolator`, 1 otherwise). The record exists when the difference is
    positive. A family with `dropped_rows` equal to `None` (interface,
    regularisation, other) gets no `unattributed` record, because it has no basis.
    For Surfe, `dropped_rows = outside` and so the difference is zero.
    A reason with no registered explanation keeps its count and gets
    `UNKNOWN_EXPLANATION`. No row is dropped from the count.
18. **Out of S5.** Residual reports for fold cells and P2 cells (the cells that the
    spec table links to FR-024 with SC-003, SC-007 and SC-008) are not in S5. S2
    covers the fold and P2 data paths. S6 covers the backend matrix. The summary
    code does not depend on the engine, so S5 adds no cell for them.
19. **Python 3.9.** `_structured_grid_3d.py` uses `Optional[np.ndarray]` from
    `typing`. It gets no `from __future__ import annotations` (the file stays
    unchanged in style). `_diagnostics.py` already has the future import and uses
    `float | None` in annotations.
20. **Exports.** `loop_interpolation.__init__` exports `ResidualSummary`,
    `ExclusionRecord`, `direction_cosines` and `build_exclusion_records`.
    `summarise_residuals`, `UNKNOWN_EXPLANATION` and `EXCLUSION_EXPLANATIONS` are
    not exported. The feature code and the tests import them from
    `loop_interpolation._diagnostics`.

#### S5 contracts

| Boundary | Input | Output | Owner |
|---|---|---|---|
| `StructuredGrid3DGeometry.nsteps` and `.step_vector` setters, `set_nelements` | node count, step, or element count | origin and extent kept (extent enclosed by less than one cell for the step setter) | `_structured_grid_3d.py` |
| `StructuredGrid2DGeometry.nsteps`, `maximum` and cell counts | `origin`, `nsteps`, `step_vector` | `nsteps` is an int array; counts derived, never stale; `maximum` legacy value `origin + nsteps * step_vector` | `_structured_grid_2d.py` |
| `RectilinearGrid` resize setters | any value | `NotImplementedError`, no state change | `_3d_rectilinear_grid.py` |
| `BoundingBox.cell_centres` | the box | local centres, one for each cell of `regular_grid(local=False)` (unchanged) | `_bounding_box.py` |
| `GeologicalFeature.evaluate_gradient_misfit` and `evaluate_value_misfit` | stored local constraint arrays | old public shape and names; world-correct; stored arrays, weights, coefficients and both up-to-date flags unchanged; stale model gives a warning and NaN | `_geological_feature.py` |
| `GeologicalFeature.evaluate_gradient_cosines` and `evaluate_value_residuals` | the same | 1-D arrays; same rules | `_geological_feature.py` |
| `direction_cosines` | two `(N, d)` arrays in one basis | `(N,)` cosines | `_diagnostics.py` |
| `GeologicalFeature.verification_report` | the feature | a new `ConstraintDiagnosticsReport` with residuals and exclusions; no state change | `_geological_feature.py` |

#### S5 functions and signatures

`StructuredGrid3DGeometry` (`_structured_grid_3d.py`; add `from typing import
Optional` at the top):

```python
def _cells_for_length(
    self, length: np.ndarray, step_vector: Optional[np.ndarray] = None
) -> np.ndarray:
    """Number of cells per axis that cover ``length``.

    ``step_vector`` defaults to the current step. The rule of S1 is unchanged:
    ``ceil(length / step - 1e-9)``; 0 cells become 2; a negative count raises
    ``ValueError("Cannot resize the grid.")``.
    """
```

`set_nelements`: the line `nsteps = np.ceil(...).astype(int)` becomes

```python
nsteps = np.maximum(
    np.ceil((self.maximum - self.origin) / step_vector).astype(int), 2
)
```

`nsteps` setter:

```python
@nsteps.setter
def nsteps(self, nsteps):
    nsteps = np.broadcast_to(np.asarray(nsteps), self._nsteps.shape)
    if not np.all(np.isfinite(nsteps)) or np.any(nsteps != np.round(nsteps)):
        raise ValueError("Cannot resize the grid. nsteps must be whole numbers.")
    nsteps = nsteps.astype(int)
    if np.any(nsteps < 2):
        raise ValueError("Cannot resize the grid. nsteps must be at least 2.")
    extent = self.nsteps_cells * self._step_vector
    self._step_vector = extent / (nsteps - 1)
    self._nsteps = np.array(nsteps)
    self._on_geometry_change()
```

`step_vector` setter:

```python
@step_vector.setter
def step_vector(self, step_vector):
    step_vector = np.array(step_vector, dtype=float)
    extent = self.nsteps_cells * self._step_vector
    cells = self._cells_for_length(extent, step_vector=step_vector)
    self._nsteps = cells + 1
    self._step_vector = step_vector
    self._on_geometry_change()
```

`StructuredGrid2DGeometry` (`_structured_grid_2d.py`): in `__init__`, replace the
assignments of `self.nsteps`, `self.maximum`, `self.nsteps_cells`, `self.n_cell_x`
and `self.n_cell_y` (lines 23, 26, 28 to 30) with `self.nsteps = nsteps` (the setter
coerces) and keep the other lines. Add:

```python
@property
def nsteps(self) -> np.ndarray:
    return self._nsteps

@nsteps.setter
def nsteps(self, nsteps):
    self._nsteps = np.ceil(np.asarray(nsteps)).astype(int)

@property
def nsteps_cells(self) -> np.ndarray:
    return self.nsteps - 1

@property
def n_cell_x(self):
    return self.nsteps[0] - 1

@property
def n_cell_y(self):
    return self.nsteps[1] - 1

@property
def maximum(self) -> np.ndarray:
    # Legacy value (FR-028): it counts nodes as cells, so it is one step
    # beyond the last node. ``nodes`` uses ``nsteps_cells``. Do not change it.
    return self.origin + self.nsteps * self.step_vector
```

The properties `n_cell_x` and `n_cell_y` return the same NumPy integer type as
the old attributes. The task greps for readers and writers of the five names
before the edit (`git grep -n "n_cell_x\|n_cell_y\|nsteps_cells\|\.maximum ="`). A
writer outside `__init__` must move to the derived form. The `__init__` order
matters: `nsteps` is set before `origin` is read.

`RectilinearGrid` (`_3d_rectilinear_grid.py`): override the three setters.

```python
_RESIZE_MESSAGE = (
    "A rectilinear grid has fixed node arrays and cannot be resized. "
    "Build a new RectilinearGrid."
)

@property
def origin(self): return self._origin
@origin.setter
def origin(self, origin): raise NotImplementedError(_RESIZE_MESSAGE)
# the same pattern for nsteps and step_vector
```

`_diagnostics.py` (new code; existing classes and fields keep their order):

```python
@dataclass(frozen=True)
class ResidualSummary:
    name: str                 # constraint or data family
    kind: str                 # "physical" or "assembled"
    unit: str                 # "field", "degree" or "assembled"
    count: int                # rows with a finite residual
    nan_count: int            # rows with no finite residual (not in the statistics)
    rms: float | None         # None when count == 0
    max_abs: float | None
    mean_abs: float | None

@dataclass(frozen=True)
class ExclusionRecord:
    family: str               # a data family, or the constraint key of a builder record
    reason: str
    count: int                # points for outside_model, rows for unattributed
    explanation: str

EXCLUSION_EXPLANATIONS: dict[str, str] = {
    "invalid_reference_value": "The reference value is NaN or not finite at the sampled point.",
    "invalid_reference_direction": "The reference direction has zero length or is not finite.",
    "evaluation_error": "The reference feature raised an error when it was evaluated.",
    "outside_model": "The point lies outside the model bounding box. The count is in points.",
    "unattributed": "The assembly has fewer rows than expected and no reason was recorded. The count is in matrix rows.",
}
UNKNOWN_EXPLANATION = "No explanation is registered for this reason."

def summarise_residuals(name: str, kind: str, unit: str, residuals: np.ndarray) -> ResidualSummary: ...
def direction_cosines(model_vectors: np.ndarray, data_vectors: np.ndarray) -> np.ndarray: ...
def build_exclusion_records(
    excluded_rows: dict[str, dict[str, int]],
    families: dict[str, ConstraintFamilyDiagnostics],
    rows_per_point: dict[str, int] | None = None,
) -> tuple[ExclusionRecord, ...]: ...
```

`summarise_residuals` flattens the array, counts non-finite values in `nan_count`
and computes `rms`, `max_abs`, `mean_abs` from the finite values. It does not
change its input. `direction_cosines` raises `ValueError` if the two shapes
differ. For each row it divides the dot product by the product of the two norms.
A row whose norm is zero or not finite gives `NaN`. Both inputs stay unchanged.

`build_exclusion_records` returns the records sorted by `(family, reason)`. The
rules are:

- One record for each reason in `excluded_rows`. The family is the key of
  `excluded_rows`. The count is not subtracted from any family.
- One `outside_model` record for each family with `outside_model_point_count > 0`.
- One `unattributed` record for each family with `dropped_rows` not `None` and
  `dropped_rows - outside_model_point_count * k > 0`. Here `k` is
  `rows_per_point.get(family, 1)` (`rows_per_point` `None` means an empty dict).
  The count of the record is the difference.
- A count of zero makes no record.

The report class gets three fields after `outside_model_points`:

```python
physical_residuals: dict[str, ResidualSummary] = field(default_factory=dict)
assembled_residuals: dict[str, ResidualSummary] = field(default_factory=dict)
exclusions: tuple[ExclusionRecord, ...] = ()
```

`to_dict` adds the keys `"physical_residuals"`, `"assembled_residuals"` (a dict of
the fields of each summary) and `"exclusions"` (a list of dicts) only when the
field is not empty. `summary()` adds one line for each non-empty block, after the
existing lines. The properties `total_rows` and `active_families` do not change.

`GeologicalFeature` (`_geological_feature.py`; the module has the future
import):

```python
def evaluate_gradient_misfit(self) -> np.ndarray:
    """Cosine between each stored direction and the model gradient (1 = equal).

    One row for each non-empty family (gradient, then normal): shape (1, N) for
    one family, (2, N) for two of equal length, (1, N1 + N2) for two of
    different length. Shape (0,) with no rows. The model gradient and the data
    direction are compared in world coordinates. The method changes no stored
    array, weight, coefficient or flag. A stale model gives a UserWarning and
    NaN. The meaning of the value changed in S5 (see COMPAT.md).
    """

def evaluate_value_misfit(self) -> np.ndarray:
    """|value - field| at each stored value row, relative to the field range.

    The result is divided by the field range only when the range is finite and
    greater than zero. Otherwise it is the absolute difference in field units.
    """

def evaluate_gradient_cosines(self) -> np.ndarray:
    """1-D cosines: gradient rows first, then normal rows. Shape (N,)."""

def evaluate_value_residuals(self) -> np.ndarray:
    """1-D |value - field| in field units, not normalised. Shape (N,)."""

def verification_report(self) -> ConstraintDiagnosticsReport:
    """Report with physical residuals, assembled residuals and exclusions."""

def _stale_reason(self) -> Optional[str]: ...
def _gradient_cosines_fresh(self) -> list: ...   # one array for each non-empty family
def _value_residuals_fresh(self) -> np.ndarray: ...
def _physical_residual_summaries(self) -> dict[str, ResidualSummary]: ...
```

`_stale_reason` checks `self.builder._up_to_date`, then the same attribute of
every builder in `self.builder.faults`, then `self.interpolator.up_to_date`. It
returns a text for the first stale object, or `None`.

Method body of `_gradient_cosines_fresh`, in order: (1) for each of the arrays
`get_gradient_constraints()` and `get_norm_constraints()` that has rows, `rows =
np.array(array, dtype=float, copy=True)`; (2) `positions = rows[:, :3]`, and
`reproject` them if `bounding_box` is not `None`; (3) `data = rows[:, 3:6]`, and
apply `self.interpolator.local_to_world_vectors(data)` when `bounding_box` is
not `None`; (4) `model = self.evaluate_gradient(positions)`; (5) append
`direction_cosines(model, data)`. The method returns the list (one array for each
non-empty family). `evaluate_gradient_cosines` calls `_stale_reason`; on a stale
model it warns and returns `np.full(n, np.nan)`; otherwise it returns
`np.concatenate` of the list (`np.zeros(0)` if the list is empty). The old
`evaluate_gradient_misfit` calls the same code and then builds the shape of
decision 9. The call `self.builder.up_to_date()` is removed from all of them.

`_value_residuals_fresh`: copy `get_value_constraints()`, reproject the positions,
return `np.abs(rows[:, 3] - self.evaluate_value(positions))`.
`evaluate_value_misfit` returns the stale NaN array, or divides that result by
`self.max() - self.min()` when the range is finite and positive.

`verification_report` body: `base = self.interpolator._build_constraint_diagnostics_report()`
(a new object; it does not set `latest_diagnostics_report`). If
`_stale_reason()` is not `None`, write the warning and skip the two residual
blocks. If `interpolator.c` is `None`, or the interpolator has no
`calculate_residual_for_constraints`, skip the assembled block. Otherwise
`assembled` = for each `(name, residual)` in
`calculate_residual_for_constraints().items()`, `summarise_residuals(name,
"assembled", "assembled", residual)`. Then `physical =
self._physical_residual_summaries()`; `excluded = getattr(self.builder,
"excluded_rows", {})`; `rows_per_point = {"gradient": 2}` if the interpolator is
a `DiscreteInterpolator`, else `None`. Return `dataclasses.replace(base,
physical_residuals=physical, assembled_residuals=assembled,
exclusions=build_exclusion_records(excluded, base.families, rows_per_point))`.
The physical summaries are `value` (from `_value_residuals_fresh`, unit
`"field"`), and `gradient` and `normal` (`degrees(arccos(clip(cosine, -1, 1)))`
of the list from `_gradient_cosines_fresh`, unit `"degree"`). So there is one
code path for each quantity.

The names `ResidualSummary`, `ExclusionRecord`, `EXCLUSION_EXPLANATIONS`,
`UNKNOWN_EXPLANATION`, `summarise_residuals`, `direction_cosines`,
`build_exclusion_records`, `verification_report`, `evaluate_gradient_cosines`,
`evaluate_value_residuals`, `_stale_reason`, `_gradient_cosines_fresh`,
`_value_residuals_fresh` and `_physical_residual_summaries` are free (`git grep
-n -w` on 2026-11-05; `_build_constraint_diagnostics_report` exists and is used).

#### S5 tests

`packages/loop_common/tests/test_resize_paths.py` (new). Grids: `G1` origin
`(0,0,0)`, `nsteps (11,6,4)`, step `(1,2,3)`; `G2` origin `(100,-20,5)`,
`nsteps (5,5,5)`, step `(0.5,0.5,0.5)`; `G3` origin `(-1e4,2e3,-30)`, `nsteps
(21,3,8)`, step `(10,100,7.5)`. The 2D grids are `H1` origin `(0,0)`, `nsteps
(5,5)`, step `(1,1)`, and `H2` origin `(10,-5)`, `nsteps (8,4)`, step `(0.5,2)`.

| Test function | Cells | Checks |
|---|---|---|
| `test_reapplying_same_value_10_times_changes_nothing_3d` | 3 grids x 4 setters (origin, maximum, nsteps, step_vector) = 12 | 10 repeats of each assignment leave `origin`, `maximum`, `step_vector` and `nsteps` equal to the start: counts exactly, lengths within `1e-12 * max(1, L)` (SC-011) |
| `test_nsteps_setter_keeps_extent` | 3 | `nsteps` doubled then halved: `maximum` equal to the start within `1e-12 * L`; node count equals the request |
| `test_step_vector_setter_encloses_extent` | 3 | step x 1.5: new `maximum` is at least the old one and less than one new step above it |
| `test_step_vector_a_b_a_extent_drift` | 3 | step A, then B = 1.5 A, then A again: the final `maximum` minus the start is at least 0 and less than `A + B` on each axis; record the measured value |
| `test_nsteps_below_two_is_rejected_without_change` | 1 | `ValueError`; state unchanged |
| `test_nsteps_not_whole_is_rejected_without_change` | 1 | `nsteps = (5.5, 5, 5)` and `NaN` raise `ValueError`; state unchanged; `(5.0, 5.0, 5.0)` is accepted and stored as `int` |
| `test_set_nelements_keeps_extent` | 3 | after `set_nelements(200)` `maximum` equals the start within `1e-12 * L`; the returned value equals `n_elements`; `nsteps` equals `max(ceil(extent / h), 2)` with `h = (volume / 200) ** (1 / 3)` computed in the test |
| `test_set_nelements_small_request_is_clamped` | 1 | `G1.set_nelements(1)` raises nothing; every `nsteps` is at least 2; the extent is kept |
| `test_2d_counts_are_derived` | 2 grids x 3 checks = 6 | after `nsteps` assignment, `nsteps_cells`, `n_cell_x`, `n_cell_y` and `n_elements` agree with the new `nsteps`; `nsteps` is an int array when the input is a float list; `maximum == origin + nsteps * step_vector` (legacy value, one step beyond `nodes.max(axis=0)`) |
| `test_rectilinear_resize_is_rejected` | 3 setters | `NotImplementedError`; nodes, origin and counts unchanged |

`packages/loop_common/tests/test_cell_centre_sampling.py` (new). Boxes: `unit`
and `offset_anisotropic` from `test_bounding_box_structured_grid.py`. Frames:
rotations 0, 30 and 90 degrees about z; local origin `(0,0,0)` and `(7,-3,11)`.

| Test function | Cells | Checks |
|---|---|---|
| `test_each_cell_has_exactly_one_centre` | 12 | `grid.cell_index_at(box.project(world_centres))` is `arange(n_cells)` with `grid = box.structured_grid(local_coordinates=True)`: no `-1`, no repeat, `np.bincount` all `1`. This is the main evidence for FR-015 and SC-011 |
| `test_local_coordinates_grid_associates_own_centres` | 12 | every `box.project(world_centres)` falls in a cell of `box.structured_grid(local_coordinates=True)`: `cell_index_at(...) != -1` for all points. The world centres are `reproject(box.cell_centres(order="F"))` |
| `test_centres_are_cell_midpoints_in_world` | 12 | `reproject(box.cell_centres(order="F"))` equals the mean of the 8 corner nodes of each cell of `box.structured_grid()`, within `Tp`. Not separate evidence for SC-011 (it restates the construction of `cell_centres`) |
| `test_centres_match_exported_cells` | 12 | world centres equal `grid.cell_centres` within `Tp`. Not separate evidence for SC-011 |
| `test_centres_inside_box` | 12 | each world centre is inside `[origin, maximum]` within `Tp`. Not separate evidence for SC-011 |

The rotation is set with `box.set_local_transform(local_origin=..., local_rotation=...)`.
The task reads the S1 test `test_rigid_transform_contract.py` for the call and for
`Tp`. Rotation 0 with origin 0 is the control. Expected result by reading: all
pass on the current tree ("coverage, no defect", OP12). Only the first two tests
count as evidence for SC-011.

`packages/loop_interpolation/tests/test_residual_exclusion_report.py` (new). It
tests the pure functions and the report class, with no model.

| Test function | Cells | Checks |
|---|---|---|
| `test_direction_cosines_pure_function` | 6 | equal, opposite and orthogonal vectors; a zero row and a NaN row give `NaN`; shape mismatch raises `ValueError`; both inputs unchanged; scale invariance |
| `test_report_without_verification_is_unchanged` | 1 | `to_dict().keys()` equals the pre-S5 key set and the constructor still takes the old fields in order (rollback check) |
| `test_builder_records_keep_their_own_key` | 1 | `excluded_rows = {"restricted_orthogonal": {"invalid_reference_direction": 3}}` gives a record with `family == "restricted_orthogonal"` and count 3; no data-family count changes; no `unattributed` record appears |
| `test_unknown_reason_keeps_its_count` | 1 | a reason `"x"` in `excluded_rows` gives a record with `UNKNOWN_EXPLANATION` and the same count |
| `test_unattributed_rows_are_reported` | 1 | a `value` family with `dropped_rows` 5 and 2 outside points gives `unattributed` 3 |
| `test_gradient_outside_points_give_no_false_unattributed` | 1 | a `gradient` family with 10 source points, 2 outside points and `dropped_rows` 4 (`2 * 10 - 16` rows) and `rows_per_point={"gradient": 2}` gives no `unattributed` record; with `dropped_rows` 6 it gives 2 |
| `test_dropped_rows_none_gives_no_unattributed` | 1 | an `interface` family with `dropped_rows=None` and no outside points gives no record |
| `test_surfe_style_dropped_equals_outside` | 1 | `rows_per_point=None`, `dropped_rows` equal to the outside count: no `unattributed` record |

`tests/unit/modelling/test_diagnostics_frames.py` (new). It uses
`world_foliation_model` (S4) and solves with the pattern of
`test_reference_constraint_frames.py`.

| Test function | Cells | Checks |
|---|---|---|
| `test_misfit_leaves_state_unchanged` | 4 anchors x 3 rotations x 2 engines (FDI, PLI) = 24 | after a solve, snapshot (deep copies) of every `interpolator.data[...]` array, `c`, `interpolator.up_to_date`, `builder._up_to_date`, the weights of `constraints`, and `id(latest_diagnostics_report)`; call the two old misfit methods; every snapshot is equal (`np.array_equal`) (SC-009). It uses no new name |
| `test_new_misfit_methods_leave_state_unchanged` | 24 | the same snapshot; call `evaluate_gradient_cosines`, `evaluate_value_residuals` and `verification_report()` |
| `test_gradient_cosine_matches_analytic_direction_exact_field` | 24 | `monkeypatch` replaces `feature.evaluate_gradient` with a stub that returns the analytic world direction of the R1 foliation at the given points. R1 gradient and normal rows: `abs(1 - cos) <= 1e-10` for both the old method (one row) and `evaluate_gradient_cosines`. This tests the basis and column logic, not the solver (SC-009) |
| `test_gradient_cosine_matches_analytic_direction_solved_model` | 24 | the real solved model: `abs(1 - cos) <= tol_solved`. The value `tol_solved` is measured on the solved models. The task records the maximum and sets the bound with a margin of 2. It is not `1e-10` |
| `test_value_residuals_are_absolute_difference` | 4 anchors x 3 rotations = 12 | `evaluate_value_residuals()` equals `abs(value - evaluate_value(world_points))` at the 27 value rows; no NaN; unit field |
| `test_value_misfit_is_normalised_by_field_range` | 12 | `evaluate_value_misfit()` equals `evaluate_value_residuals() / (feature.max() - feature.min())` within `1e-12` |
| `test_value_misfit_with_zero_range_does_not_divide` | 2 | the model has value rows with one value (a near-constant field). The solver gives noise near `1e-12`, so the test forces the range with `monkeypatch`: `feature.max` and `feature.min` return the same number (cell 1) or `NaN` (cell 2). No `RuntimeWarning` (`warnings.simplefilter("error")`); the result equals the unnormalised absolute difference |
| `test_misfit_with_no_rows_is_empty` | 2 | no gradient rows: `shape == (0,)`; no value rows: `shape == (0,)` |
| `test_old_misfit_shapes_are_kept` | 4 | one family: `(1, N)`; two families of equal length: `(2, N)`; two of different length: `(1, N1 + N2)` (OP16); the values equal `evaluate_gradient_cosines()` |
| `test_stale_misfit_is_not_rebuilt` | 2 engines x 2 kinds (interpolator stale, builder stale) = 4 | uses the old methods only. Snapshot `interpolator.up_to_date`, `builder._up_to_date`, `c` and the data arrays. After a call, every snapshot is equal. On the current tree the call rebuilds the model, so `interpolator.up_to_date` changes from `False` to `True`: the test fails there (FR-031) |
| `test_stale_misfit_warns_once_and_returns_nan` | 4 | old and new methods each write exactly one `UserWarning` that names the stale object; the new methods return `np.full(n, np.nan)`; the old methods return the old shape filled with NaN |
| `test_stale_verification_report_does_not_rebuild` | 4 | `verification_report()` on a stale model: the state snapshot is equal; the two residual blocks are empty; the exclusions are present; one warning |
| `test_report_blocks_are_separate` | 2 engines = 2 | the two dicts are different objects; every physical entry has `kind == "physical"` and every assembled entry has `kind == "assembled"`; the physical `value` entry has unit `"field"` and its `rms` equals the RMS of `evaluate_value_residuals()` within `1e-12`; the assembled entries equal `summarise_residuals` of `calculate_residual_for_constraints()`; `to_dict` has two separate blocks |
| `test_unsolved_model_gives_no_residual_blocks` | 2 | a model with `interpolator.c = None`: the two blocks are empty and the exclusions are present; no exception |
| `test_exclusions_are_counted_and_explained` | 2 engines = 2 | the restricted-support reference case of S1 (`test_reference_constraint_frames.py`): the builder records appear under their constraint keys (`restricted_equality`, `restricted_orthogonal`) with the analytic count; the explanation is not empty and not `UNKNOWN_EXPLANATION` |
| `test_gradient_outside_points_give_no_false_unattributed_model` | 2 engines = 2 | a feature with gradient data and 3 points outside the model box: the report has an `outside_model` record of 3 and no `unattributed` record for `gradient` |

Pre-fix evidence (FR-031). The state test with the old methods, the exact-field
cosine test and the stale-state test fail on the current tree. In the first, a
different array is stored in `data["gradient"]` after the call. In the second,
the cosine uses one column and local vectors. In the third, the call rebuilds the
model. The value-residual test records "not reproduced" if it passes on the
current tree (the value path does not change the stored positions). The resize
tests for `nsteps`, `step_vector`, `set_nelements` and the 2D counts fail on the
current tree (record the observed extent drift and the stale counts). The 2D
`maximum` is unchanged, so no 2D `maximum` test fails. The cell-centre tests are
coverage (OP12). The tests of the new names fail because the names do not exist.
The pre-fix part uses no new name.

Existing tests that must pass with no change:
`packages/loop_common/tests/test_2d_discrete_support.py` and
`tests/unit/interpolator/test_2d_discrete_support.py` (the four `maximum`
assertions; S5 does not edit these two files), `test_rectilinear_grid.py`,
`test_bounding_box_structured_grid.py`, `tests/unit/modelling/test_scalar_field_grid.py`,
the constraint diagnostics tests (`git grep -ln "get_constraint_diagnostics_report"
tests packages`), and the S1 to S4 test files.

Tests that change by design (Q4). The 3D golden test
`packages/loop_common/tests/test_structured_support_reference.py` and its data
file `packages/loop_common/tests/data/structured_support_reference.npz` change,
because the setters now keep the extent. The earlier sentence that this golden
test must pass with no change is withdrawn. The 2D golden test
`test_structured_grid_2d_reference.py` and `data/structured_grid_2d_reference.npz`
should not change (S5 keeps the 2D `maximum`). The task runs it, and if a key
differs, it records the key. The golden update follows this order: (1) run the
test and list each failing key; (2) review each key against the extent rule of
decision 1; (3) record each key, with the old and the new value, in the section
"S5 golden review (C-004)" of `verification.md` (SC-015: 0 tolerance changes); (4)
only then run `python packages/loop_common/tests/test_structured_support_reference.py
--update` (and the 2D script with `--update` only if a key differs); (5) run the
test again.

Numeric-test survey (Q4). The tests that build a model run `set_nelements` and
the `step_vector` loop. The task runs `pytest tests/unit/modelling` and `pytest
packages/loop_interpolation/tests` on the fixed code, and `packages/loop_common/tests`.
For each test that fails or changes a number, it records the test, the old
value, the new value and the cause in "S5 golden review (C-004)". The directory
`packages/loop_interpolation/tests/golden` has its own option `--update-golden`
and is the parity oracle for the Rust path. If a file there changes, the task
records it and does not regenerate it before the review (OP17).

#### S5 backend and capability cells

| Cell | Result in S5 |
|---|---|
| FDI and PLI misfit, residual and exclusion report | Mandatory (the cells above). |
| P2 | Not in S5. The summary code does not depend on the engine. |
| Residual reports for fold and P2 cells (FR-024 with SC-003, SC-007, SC-008) | Not in S5. S2 covers the fold and P2 data paths, and S6 covers the matrix. |
| Native and Surfe | Not in S5. The misfit methods read the Python interpolator state only. Backend matrix: S6. |
| Resize of structured, 2D and rectilinear grids | Mandatory (resize tests above). Unstructured meshes cannot resize (OP4, S3). |

#### S5 risks

- The old misfit methods keep their public shape. The gradient value changes
  meaning (decision 10). A caller outside the repository (the QGIS plugin, user
  scripts) can depend on the old number. `COMPAT.md` and `ROADMAP.md` record it.
  No caller in this repository uses the methods.
- `evaluate_gradient_misfit` used to call `builder.up_to_date()`, so it rebuilt a
  stale model. Now it warns and returns NaN. A caller that relied on the rebuild
  gets NaN and one warning, and must build the model first (OP13, closed).
- The extent change from `set_nelements` is a new risk. The old `nsteps` setter
  kept `nsteps * step`. The new one keeps the extent. A model build runs
  `set_nelements` and the `step_vector` loop, so grid steps, node coordinates and
  solved values can change slightly in tests and in user models. The golden review
  and the numeric-test survey record each change. The risk grows if an
  interpolation golden file changes (OP17).
- Repeated `step_vector` assignments can add up to more than one cell of
  enlargement (decision 1). A test measures the bound.
- The 2D `maximum` keeps the legacy value, so it still counts nodes as cells. A
  reader that needs the true last node must use `nodes`. The value now follows
  later `nsteps` assignments (it was stale before).
- `RectilinearGrid` now raises on resize. A caller that assigned `origin` or
  `step_vector` got a wrong grid before. Now it gets a clear error.
- `evaluate_gradient` and `evaluate_value` of a feature call `builder.up_to_date()`.
  The misfit methods call them only after `_stale_reason()` returned `None`. The
  stale tests check this for FDI and PLI.
- The first S5 planner run read 14 code files, by targeted line ranges. The limit
  of the spec-planner skill is 10 code files. The fix run read more. The design
  uses only those reads. The risk is that an unread part of those files has other
  readers. Each task greps for readers and writers before it edits. S5 changes 10
  code files at most (see "S5 files"). No task group changes more than 10.
- `build_exclusion_records` uses `rows_per_point` of 2 for the gradient family of
  a `DiscreteInterpolator`. If an interpolator adds a different number of rows for
  a gradient point, the `unattributed` count is wrong for that engine. The
  gradient end-to-end test covers FDI and PLI only.

#### S5 constitution check

| Project rule | Result | How S5 agrees |
|---|---|---|
| Work stays within the existing packages | Pass | Three `loop_common` files, two `loop_interpolation` files (`_diagnostics.py` and `__init__.py`) and one `LoopStructural` file change. No new module. |
| Gates `pytest` and `ruff check LoopStructural packages` | Pass | The S5 row of "Stages and rollback" names both. |
| numpy docstrings, type hints, py39 | Pass | New functions and dataclasses have numpy docstrings and hints. `_diagnostics.py` has `from __future__ import annotations`. `_structured_grid_3d.py` uses `Optional`. |
| The public `GeologicalModel` API stays compatible; a break needs a ROADMAP.md decision and a major-version note | Pass | No `GeologicalModel` signature changes. The report only gets new defaulted fields. The misfit methods keep their names and public shape. The gradient value changes meaning: this is recorded in `COMPAT.md` (migration notice) and in a `ROADMAP.md` decision entry (Q2). The 2D `maximum` does not change. |
| Stages can be tested and reversed independently | Pass | S5 reverts to the S4 merge. |
| `ROADMAP.md` owns decisions; update it when one changes | Pass | S5 adds one decision entry for the misfit methods: "Misfit diagnostics of GeologicalFeature (spec 005, S5)". No other decision changes. |
| Feature documents stay under `specs/NNN-name/` | Pass | Only `plan.md`, `tasks.md` and `verification.md` change there. |

#### S5 files

At most 10 code files (.py) change in S5: 6 source files and 4 new test files. No
task group changes more than 10 code files. Data files and documents are not
counted.

| File | Change |
|---|---|
| `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py` | `_cells_for_length` gets `step_vector`; new `nsteps` and `step_vector` setters; `set_nelements` clamp; `from typing import Optional` |
| `packages/loop_common/src/loop_common/geometry/_structured_grid_2d.py` | `nsteps` property with int coercion; `nsteps_cells`, `n_cell_x`, `n_cell_y` and `maximum` become properties (`maximum` keeps the legacy value, with a comment) |
| `packages/loop_common/src/loop_common/supports/_3d_rectilinear_grid.py` | `origin`, `nsteps`, `step_vector` setters raise `NotImplementedError` |
| `packages/loop_interpolation/src/loop_interpolation/_diagnostics.py` | `ResidualSummary`, `ExclusionRecord`, explanation table, three functions, three report fields, `to_dict` and `summary` |
| `packages/loop_interpolation/src/loop_interpolation/__init__.py` | Export four new names (`ResidualSummary`, `ExclusionRecord`, `direction_cosines`, `build_exclusion_records`). `summarise_residuals` and `UNKNOWN_EXPLANATION` are imported from `loop_interpolation._diagnostics` |
| `LoopStructural/modelling/features/_geological_feature.py` | Rewrite both misfit methods on shared private code; add `evaluate_gradient_cosines`, `evaluate_value_residuals`, `verification_report`, `_stale_reason` and the private helpers |
| `packages/loop_common/tests/test_resize_paths.py` | New |
| `packages/loop_common/tests/test_cell_centre_sampling.py` | New |
| `packages/loop_interpolation/tests/test_residual_exclusion_report.py` | New |
| `tests/unit/modelling/test_diagnostics_frames.py` | New |
| `packages/loop_common/tests/test_structured_support_reference.py` | Not edited. It must pass after the regeneration of its data file |
| `packages/loop_common/tests/data/structured_support_reference.npz` | Regenerated with `--update` after review (C-004, SC-015) |
| `packages/loop_common/tests/data/structured_grid_2d_reference.npz` | Regenerated only if the 2D reference test changes (not expected) |
| `packages/loop_common/tests/test_2d_discrete_support.py` and `tests/unit/interpolator/test_2d_discrete_support.py` | Not edited. They must pass with no change (Q1) |
| `COMPAT.md` | S5 notice: misfit gradient value (migration notice), stale model gives NaN and a warning, extent kept by the setters and `set_nelements`, rectilinear resize error, new report fields |
| `ROADMAP.md` | Decision "Misfit diagnostics of GeologicalFeature (spec 005, S5)" under "Decisions" |
| `specs/005-coordinate-frame-correctness/verification.md` | S5 pre-fix, post-fix, golden review (C-004), numeric-test survey, backend matrix, rollback and run-time sections |

Not changed (the survey in "Later stages" listed them): `_bounding_box.py`
(`cell_centres`; conditional, OP12), `_structured_grid.py`,
`supports/_3d_structured_grid.py`, `_base_geological_feature.py`,
`_model_exporter.py`, `_discrete_interpolator.py` and `API.md`. `docs/source` and
`API.md` do not mention the changed names (`git grep`, 2026-11-05), so they need
no edit.

The `COMPAT.md` notice says: `GeologicalFeature.evaluate_gradient_misfit()`
changed the stored constraint arrays in place, compared a local vector with a
world vector, and read one column. It now keeps its name and its public shape. It
changes nothing. It returns the cosine between the world data direction and the
model gradient for each row. The old value was a dot product with an unnormalised
model gradient. The new value differs. Use `evaluate_gradient_cosines()` for a
1-D array. The method no longer rebuilds a stale model: it writes a `UserWarning`
and returns NaN. `evaluate_value_misfit()` keeps its normalisation by the field
range. It no longer divides by zero. Assigning `nsteps` or `step_vector` on a 3D
grid, and `set_nelements`, keep the extent. A rectilinear grid refuses to resize.
`ConstraintDiagnosticsReport` has three new optional fields. With an identity
frame and no resize, the value misfit is the same as before.

The `ROADMAP.md` entry gives the decision date `2026-11-05` (spec 005, S5, Q2).
It says that the misfit methods keep their names and public shape, that the
gradient value is now a cosine in world coordinates, that a stale model gives a
warning and NaN, and that the new names `evaluate_gradient_cosines` and
`evaluate_value_residuals` give the 1-D result. It lists the compatibility notice
in `COMPAT.md`.

#### S5 task order (for the spec-tasks skill)

No two tasks in one wave write `verification.md`. Two tasks never write the same
file in one wave. Each task group has at most 10 code files.

1. The two `loop_common` test files (resize, cell centres), in two tasks that can
   run in parallel. They run on the current tree.
2. The pre-fix evidence for step 1 in `verification.md` (one task, own wave).
3. The 3D setter fix with the `set_nelements` clamp, the 2D properties and the
   rectilinear guard (three tasks, separate files, parallel). The 2D task greps
   the readers and writers of the five 2D names first (`git grep -n
   "n_cell_x\|n_cell_y\|nsteps_cells\|\.maximum ="`).
4. The numeric-test survey and golden review (one task, after step 3). It runs
   `tests/unit/modelling`, `packages/loop_interpolation/tests` and
   `packages/loop_common/tests`, records each change in "S5 golden review
   (C-004)", and regenerates the golden `.npz` file after review. All later tasks
   that build models depend on this step. It is the only `verification.md` writer
   of its wave.
5. The post-fix evidence for the resize and cell-centre tests (one task, own wave,
   after step 4).
6. `_diagnostics.py` (dataclasses, functions, report fields), with the pure-function
   tests in `test_residual_exclusion_report.py`. Then the `__init__` export. This
   step does not depend on steps 3 to 5 and can run in parallel with them, but its
   evidence goes to `verification.md` in a later wave.
7. The feature tests in `test_diagnostics_frames.py` (state, cosine, value, empty,
   shapes, stale, residual, exclusion), written before the misfit change. This task
   owns that file alone. The pre-fix evidence goes to `verification.md` in its own
   wave. It depends on step 4, because the solved-model tolerance and the values
   depend on the built support.
8. The misfit rewrite and `verification_report` (one file; one task). It depends on
   steps 6 and 7.
9. `COMPAT.md`, `ROADMAP.md`, the rollback check, and last the gates (the
   orchestrator runs the full suites when memory is tight). `verification.md` is
   written in the last wave only by this step.

S5 rollback check. Revert to the S4 merge. Then run these checks.

- (a) The field names and order of `ConstraintFamilyDiagnostics`,
  `RegionCoverageDiagnostics` and the first five fields of
  `ConstraintDiagnosticsReport` are equal before and after (`dataclasses.fields`).
- (b) For a report built with no verification, `to_dict()` is equal before and
  after.
- (c) `from loop_interpolation import *` on the S4 code gives a subset of the S5
  names.
- (d) The signatures of `GeologicalFeature.evaluate_gradient_misfit` and
  `evaluate_value_misfit` are equal. The result shape for one family is equal.
- (e) With an unchanged grid, the 3D `origin`, `maximum`, `nsteps` and
  `step_vector` of a structured support are equal on both trees.
- (f) The four 2D `maximum` assertions pass on both trees.

The result of each check goes to `verification.md`.

### Stage S6 design

S6 is the last stage of spec 005. It has five aims. First, S6 makes the auxiliary
point transformation correct for ordinary point clouds. It also keeps one class behind
both import paths (FR-025, SC-010). Second, S6 validates the R2 fault frame against
analytic references: values, gradients, constraint residuals and displacement
(FR-026, FR-027, SC-004). Third, S6 repairs the physical residuals of faulted
features, which S5 left open (FR-023, FR-026). Fourth, S6 builds the backend and
capability matrix, with native and Surfe (FR-030, SC-018). Fifth, S6 reviews the
golden results of all stages, closes the open findings of S5, runs the final gates
and writes the "005 complete" line (FR-029, SC-014, SC-015).

Base: branch `feat/005-s5-diagnostics-and-resize` (S5, pull request #402, open). The
planning worktree is `LoopStructural-005-s6` on the branch `docs/005-s6-plan`. If S5
changes before it merges, run this planner again on the merged result.

Review status. An independent review found problems in the first S6 design. The
maintainer answered ten questions on 2026-11-06. This text includes the fixes and the
answers (see "Answers"). OP18, OP19 and OP20 are closed.

The spec has no `[NEEDS CLARIFICATION]` mark (C-001 to C-019). The planner found two
places where the spec text and the code disagreed. First, R2 named an "existing
test-side displacement evaluator" with trishear parameters. The repository has none
(fact 9). The maintainer amended the spec: clarification C-019 replaces that text with
the library displacement paths. Second, FR-025 speaks of "duplicate public entry
points" that keep "compatible import paths". The `loop_common` copy cannot be imported
(fact 2). The maintainer decided to keep one class (decision 1).

Two terms have one meaning in all S6 text.

- **Non-vacuity check**: a check that fails when a test would pass for a wrong reason.
  For example, a comparison of two empty sets passes for any code. A non-vacuity check
  asserts that the data are not empty, or that the effect under test is large enough.
- **Measure first**: the test task runs the independent computation and records the
  maximum. The task sets the bound at 2 times that maximum (the S5 rule). If the bound
  would be above the cap in the table, the task stops and reports. The task never
  raises a cap to make a test pass.

#### S6 facts from the code (read on 2026-11-06)

1. Two copies of `EuclideanTransformation` exist. They are
   `packages/loop_common/src/loop_common/math/_transformation.py` (178 lines) and
   `LoopStructural/utils/_transformation.py` (176 lines). A `diff` shows only four
   differences: the import line of the logger, `from __future__ import annotations`,
   one type hint and the quote style. Only `LoopStructural.utils` exports the class
   (`LoopStructural/utils/__init__.py:15`). Callers (`git grep -n EuclideanTransformation`,
   2026-11-06, all file types except `specs`): `docs/source/index.rst:88` imports the
   class and does not use it. `tests/unit/utils/test_transformation.py` is the only
   test. It uses `dimensions=3` in two tests that expect a `ValueError` for too few
   columns (lines 63 and 70). No caller in the repository uses the class in three
   dimensions with real data.
2. The `loop_common` copy cannot be imported. Line 5 is `from . import getLogger`.
   `loop_common.math.__init__` does not define `getLogger`. The import raises
   `ImportError: cannot import name 'getLogger' from 'loop_common.math'` (run with the
   worktree sources on `PYTHONPATH`, 2026-11-06). `loop_common.math.__init__` does not
   import the module. So the duplicate has no working entry point today.
   `loop_common.logging.get_logger(name, level="INFO", ...)` exists
   (`loop_common/logging/logger.py:36`). `loop_common/__init__.py` imports `geometry`,
   `io`, `logging`, `math` and `supports` in that order. `geometry/_bounding_box.py`
   imports `loop_common.math` early. `loop_common.logging` imports no `loop_common.math`
   module, so an import of `get_logger` in `math/_transformation.py` has no cycle.
3. `inverse_transform` is broken in both copies. It slices `points[: self.dimensions]`.
   That slices rows, not columns. It then stacks `points[self.dimensions :]` rows below
   the result. The `einsum` contracts the last axis of the input with the rotation
   block. So the result is correct only for the shape `(1, dimensions)`. It raises
   `ValueError` when the column count differs from `dimensions`. It returns only
   `min(N, dimensions)` rows when the column count equals `dimensions` and `N > 1`. The
   existing test `test_inverse_transform_is_broken_for_normal_point_clouds_bug`
   (line 93) expects the `ValueError`. FR-029 says that a test must not keep a defect.
   S6 replaces this test.
4. `transform` copies the input with `np.copy(points)`. An integer array keeps its
   dtype, so the rotated values are truncated. A 1-D array gives `IndexError`
   (`points.shape[1]`). `fit` imports `scikit-learn` first. `scikit-learn` is a
   dependency of the root package (`pyproject.toml` lines 42 and 61). It is not a
   dependency of `loop_common` (`packages/loop_common/pyproject.toml` line 21 lists
   `numpy`, `pandas`, `pydantic`, `scipy`, `pyvista` and `pyyaml`). Without it, `fit`
   logs an error and returns `None`. It sets no state in that case.
5. Other defects of the class. The reading is confirmed by a short script outside the
   repository (scikit-learn 1.9.0, 2026-11-06). The script ran the old formulas, not the
   class.
   - `_rotation` has `-1` at index `[2, 2]`. For `dimensions=3` the matrix has
     determinant -1. It reverses z. The old output for the point `(2, 2, 5)` with
     angle 90 degrees and translation `(1, 2, 3)` is `(6e-17, 1, -2)`. The correct
     value is `(0, 1, 2)`.
   - `fit` calls `np.dot(coeffs[0, :], [1, 0])`. For `dimensions=3` the first vector
     has 3 entries and the second has 2. By reading, this raises `ValueError`. It was
     not run.
   - `fit` takes `-arccos(c0x)`. For a main axis with a negative `y` component this
     has the wrong sign. The sign rule of scikit-learn makes the largest absolute
     component of the axis positive. So the error shows for lines with a direction
     between -45 and 0 degrees. Measured: a line at -30 degrees gives the angle -30
     degrees (correct: +30) and the variance along y after the transform is 6.51. A
     line at -10 degrees has the variance 1.01. Lines at +30, +120 and -60 degrees
     have a variance of about 1e-31 (correct by chance of the sign rule).
   - `fit` needs scikit-learn and fails silently without it (fact 4).
6. Environment (`importlib.util.find_spec`, read only, 2026-11-06, venv
   `C:\Users\lachl\dev\LoopStructural\.venv`). Not installed: `loop_interpolation_core`
   (the native extension), `surfepy` (the module that `_surfe_wrapper.py:7` imports),
   `pyevtk`, `geoh5py`, `omf`, `trimesh`, `loopstructuralvisualisation`. Installed:
   `sklearn`, `scipy`, `pyvista`, `matplotlib`, `dill`, `skimage`. The venv has no `pip`
   module (`python -m pip` fails, 2026-11-06). The editable install of `loop_common` in
   this venv points to the main worktree
   (`C:\Users\lachl\dev\LoopStructural\packages\loop_common\src`). A test run in the S6
   worktree must put the worktree sources first on `PYTHONPATH`, as the S5 runs did. The
   pip package name for the module `surfepy` is not verified (fact 21).
7. Native adapter (`_native.py`, 595 lines). `check_supported` begins at line 116
   (`_native.py`) and the adapter calls it at line 533, before the core is imported. It
   rejects a rotated support (`support.rotation_xy` not the identity), a region mask,
   extra sparse rows and a class that is not `FiniteDifferenceInterpolator`. A rotated
   computational frame is not rejected: the R1 fixture builds an unrotated local support
   over the enclosure of the projected corners. The adapter reads the stored local
   arrays (`get_*_constraints`) and calls these methods of `core.Problem`: `add_value`,
   `add_gradient`, `add_normal`, `add_tangent`, `add_interface`, `add_value_bound`,
   `add_pair_bound`, `set_regularisation`, `set_border_regularisation`,
   `add_directional_regularisation`, `set_ridge_factor`, `set_start` and `solve`.
   `solve` imports `loop_interpolation_core` inside the function. Four native test
   files in `packages/loop_interpolation/tests` skip at module level with
   `pytest.importorskip("loop_interpolation_core")`: `test_native_parity.py`,
   `test_native_inequality.py`, `test_native_orthogonal_and_limits.py` and
   `test_native_warm_start_ridge.py`. The file
   `tests/unit/modelling/test_native_fault_frame.py` does the same. The file
   `test_native_no_rows.py` calls `importorskip` inside 8 test functions (lines 49 to
   166), so it has one skipped case for each function. The file `test_native_errors.py`
   has no module-level skip: its tests run without the extension, except one function
   (line 250, `test_native_accepts_normal_data_rows_after_native_setup`) that calls
   `importorskip`. `tests/unit/modelling/test_builder_solver.py` has one native test
   with `importorskip` (line 28). No test checks the rejection of a rotated support.
   The existing native test of an R3 operation is
   `test_orthogonal_constraints_give_the_same_field_as_cg` in
   `test_native_orthogonal_and_limits.py`. It checks gradient orthogonal rows and does
   not use a rotated frame. No existing native test covers the R3 world and local
   conversion.
8. Surfe adapter (`_surfe_wrapper.py`, 211 lines). The module imports `surfepy` at the
   top (line 7). `loop_interpolation.__init__` (line 92) catches the failure and exports
   a dummy class. The adapter has no support and no bounding box of its own. It gets
   local rows from the `get_*_constraints` methods and implements
   `_evaluate_value_local` and `_evaluate_gradient_local`. The public methods of
   `GeologicalInterpolator` convert world points to local points once. The tests in
   `test_surfe_rbf_interpolator.py` skip with `HAS_SURFE`. The file
   `test_p0_surfe_nans.py` skips with `pytest.importorskip("surfe", minversion=None)`
   (line 6). That is the pip name, not the module name that the adapter imports. If the
   module `surfepy` is installed and the module `surfe` is not, the file skips
   although the tests could run (fact 21).
9. Displacement code. `FaultSegment.apply_to_points` (`_fault_segment.py:309`) moves each
   masked point `steps` times along the normalised gradient of coordinate 1, with the
   step `d / steps`. `d` is `faultfunction(gx, gy, gz) * displacement`.
   `FaultSegment.evaluate_displacement` (line 277) is the scalar path.
   `FaultDisplacement.__call__` returns `scale * gx(gx) * gy(gy) * gz(gz)` with
   `scale=0.5`. The default `BaseFault` profile has a cubic profile for the hanging wall
   and the footwall, `Ones` for `gy` and `smooth_peak` for `gz`. Search results
   (`git grep -n -i`, 2026-11-06): the word "trishear" occurred in the spec only (R2, before C-019), "apex" occurs
   nowhere in code, and "propagation" occurs in one intrusion docstring only. The
   repository has no trishear evaluator, no propagation ratio and no 16-point quadrature
   for displacement. The maintainer accepted the library paths (C-019).
10. A point outside the support. `GeologicalFeature.evaluate_value` gives NaN there. In
    `apply_to_points` a start point with NaN `gx` gets `d = 0` and does not move. A point
    that is inside at the start and leaves the support during the steps gets a NaN
    position (a NaN gradient is not removed).
11. Zero-gradient protection. `apply_to_points` (lines 405 to 408) sets `g_mag` and
    divides only rows with `g_mag > 0`. This is the protection that FR-027 means.
    `FaultSegment.evaluate_gradient` begins at line 247. It divides by the norm with no
    guard (line 272). `ConstantNormInterpolator` skips elements with a zero gradient norm
    (`loop_interpolation/_constant_norm.py`, `valid = v_t_norm > 0`, line 73). The native
    adapter drops zero-norm rows (`drop_zero_norm`).
12. The S1 test `test_planar_fault_frame_matches_analytic` already checks scalar values
    (1e-6), gradient magnitudes (1e-4) and gradient directions (0.01 degree) on the 6859
    probes for FDI and PLI, 4 anchors, 3 rotations and 2 centre modes. It builds the
    fault with `displacement=0`, so no probe moves (S3 minor finding T038). No test checks
    constraint residuals of the fault frame or any displacement.
13. S5 carry-over (the major open finding). `_value_residuals_fresh` (line 362) and
    `_gradient_cosines_fresh` (line 333) of `_geological_feature.py` call
    `self.evaluate_value` and `self.evaluate_gradient`. Both apply `_apply_faults`
    (`_base_geological_feature.py:252`), and the gradient path builds a tetrahedron of
    faulted points. `GeologicalFeature.evaluate_gradient` also calls
    `self.builder.up_to_date()`, which can rebuild a stale builder. The stored positions
    of a faulted feature are already restored by its faults. So the faults act a second
    time. S5 measured 0.080 against the true residual 0.0004 to 0.0018.
    `GeologicalInterpolator.evaluate_value` and `evaluate_gradient`
    (`_geological_interpolator.py:823` and `833`) take world points, change no state and
    never rebuild. For a feature without faults, `GeologicalFeature.evaluate_gradient`
    ends in `self.interpolator.evaluate_gradient` (line 247).
14. `GeologicalFeature._first_stale` (line 252) now follows the faults of each builder
    recursively. The S5 minor finding "checks one level of faults" is therefore closed
    in the code (commit `3f20a014`). Its text in `run-state.json` is out of date.
15. Golden and baseline files. `git diff --name-status ebc7c699d5 HEAD -- . ':!specs' ':!*.py' ':!*.md'`
    (2026-11-06) lists two files: `packages/loop_common/tests/data/structured_support_reference.npz`
    (S1: 27 values; S5: 78 keys) and `packages/map2loop/pyproject.toml`. The second file
    comes from commit `a0e82096` (lint, pull request #395). It is not a spec 005 change.
    The baseline files `structured_grid_2d_reference.npz`, `tests/fixtures/api_surface_snapshot.json`,
    `tests/unit/interpolator/*.txt` and everything under
    `packages/loop_interpolation/tests/golden/golden_files` are unchanged.
16. The BOM. `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py` starts
    with the bytes `ef bb bf` on this branch. It has no BOM in the S4 merge `62c420e8`.
    S5 added it by accident. A scan of all tracked `.py`, `.toml` and `.md` files
    (2026-11-06) finds one more BOM file:
    `tests/unit/modelling/test_fault_frame_analytic.py` (created in S1, commit `03a6e58b`).
17. Existing skip sites in the files of earlier stages: `test_world_location_consumers.py`
    has `pytest.importorskip` inside two functions (viewer and `matplotlib` in one, `pyevtk` in the other).
    The other S1 to S5 test files have no skip. `tests/conftest.py` and the files in
    `tests/fixtures` have no skip call (`git grep`, 2026-11-06).
18. Names of the assembled constraints. `add_constraints_to_least_squares`
    (`_discrete_interpolator.py:340`) stores each block in `interpolator.constraints[name]`
    with the keys `matrix`, `b` and `w`. It normalises the rows first. If the name exists,
    it adds a numeric suffix. The names in use:
    - FDI (`_finite_difference_interpolator.py`): `value`, `gradient`, `norm_0`, `norm_1`,
      `norm_2`, `interface_<id>`, the curvature terms `dxx`, `dyy`, `dzz`, `dxy`, `dxz`,
      `dyz` (names from `support.get_operators`), the border terms `dx_lower`,
      `dx_upper`, `dy_lower`, `dy_upper`, `dz_lower`, `dz_upper`, and directional terms.
    - PLI (`_p1interpolator.py`): `value`, `norm`, `edge jump` (curvature). PLI
      `add_gradient_constraints` assembles no row.
    - Orthogonality rows from the geometry builder
      (`_geological_feature_builder.py`, line 331): `<reference feature name>_orthogonal`.
      Their direction comes from the reference feature, which is a solved field.
    `interpolator.constraints` is kept after the solve: the call `self._post_solve()` is
    commented out (`_discrete_interpolator.py:1451`).
19. Junit names. The test run has `junit_family=xunit2` by default (pytest 9.1.1, checked
    in `_pytest/junitxml.py`). That family has no `file` attribute on a `testcase`.
    `classname` is the node id path with dots, plus the class names. `name` is the function
    name with its parameters. The root `pyproject.toml` sets
    `addopts = "--import-mode=importlib"` and `testpaths = ["tests"]`. The import mode
    does not change the node id, so it does not change `classname`. A short run
    (2026-11-06, a scratch folder, `--import-mode=importlib`) gave
    `classname="sub.pkg.test_ok.TestC" name="test_b[1]"` for a class test. A module-level
    `importorskip` gave one skipped case with `classname=""` and `name="sub.pkg.test_mod"`.
    A run from the repository root finds the root `pyproject.toml` as the ini file,
    because `packages/loop_common/pyproject.toml` and
    `packages/loop_interpolation/pyproject.toml` have no `[tool.pytest.ini_options]`. So
    the rootdir is the repository root for all three gate runs.
20. Check of the callers for three dimensions (`git grep`, 2026-11-06): no caller in the
    repository uses `fit`, `transform` or `inverse_transform` in three dimensions with
    real data. The fixes of the class change numbers for outside users only.
21. The pip package of `surfepy`. The adapter imports `surfepy`
    (`_surfe_wrapper.py:7`). `test_p0_surfe_nans.py` imports `surfe`. Neither the
    package list in `pyproject.toml` nor the venv names the pip package. The planner could
    not check PyPI here. A task step checks it (see "S6 task order", task 4).

#### S6 decisions

1. **One class, two import paths (maintainer answer 4).** `EuclideanTransformation`
   lives in `loop_common.math._transformation`. `loop_common.math` exports it.
   `LoopStructural/utils/_transformation.py` re-exports it. `LoopStructural.utils`
   exports the same object. The old path `LoopStructural.utils._transformation` stays
   importable. This is the "compatible import paths" of FR-025. This reverses item 2c-5
   of `ROADMAP.md` ("closed as deferred"). The S6 entry of `ROADMAP.md` says so in one
   sentence. The class keeps all attribute names, so a pickle made before S6 loads on new
   code: the old module path stays importable and exposes the class (checked by reading
   `LoopStructural/utils/__init__.py:15` and the new re-export; a test builds such a
   pickle). A pickle made after S6 stores the path `loop_common.math._transformation`.
   An older release cannot load it. `COMPAT.md` says both facts. The module uses the
   `loop_common` logger (`loop_common.logging.get_logger`). `COMPAT.md` says that too.
2. **Contract of `transform` and `inverse_transform`.** The input is an array-like of
   shape `(N, d)` with `d >= dimensions` and `N >= 0`. Another rank raises
   `ValueError`, including a 1-D array (maintainer answer 7). Fewer columns than
   `dimensions` raises `ValueError` (same text as today). Both methods return a new float
   array. They never change the input. The columns after `dimensions` are copied bit for
   bit. One point is the shape `(1, d)`. A shape `(d,)` is not a documented shape; it
   raises `ValueError` (it raised `IndexError`). `inverse_transform(transform(P))`
   equals `P` within `Tp` for the transformed columns and exactly for the others.
3. **The other defects of the class are fixed (OP20, closed as "fix", maintainer answer
   8).** S6 fixes four defects.
   - (a) The 3D `fit`. `fit` finds the main axis from the x and y columns only. The
     rotation of the class turns the xy plane about z, so only the xy part of the cloud
     can align with x. The translation is the mean of the first `dimensions` columns.
     The z column keeps its meaning (OP24).
   - (b) The sign of the fitted angle. The angle is `-arctan2(c_y, c_x)` for the main
     axis `(c_x, c_y)`. The old code used `-arccos(c_x)`. For an axis with a positive
     `c_y` both give the same number. For a negative `c_y` the old sign was wrong.
   - (c) The reversed z axis. `_rotation` has `1` at index `[2, 2]`. The 3D matrix is
     now a proper rotation (determinant 1). The 2D case does not use that entry.
   - (d) The silent need for scikit-learn. `fit` imports scikit-learn only when
     `fit_rotation` is true. If the import fails, `fit` raises `ImportError` with a
     message that names the package and the install command. It changes no state. No
     dependency is added to `packages/loop_common/pyproject.toml`. Reason: only `fit` with
     rotation needs the package, and the root package already depends on it. A
     NumPy-only replacement is not chosen. It would change the sign choice of the axis,
     and the current numbers of 2D users would change.
   The physical effect: a 3D user got a mirrored z axis (`z' = -(z - tz)`). A 3D user also
   got a wrong rotation angle for lines with a direction between -45 and 0 degrees, and
   `fit` failed for them. After S6, `z' = z - tz`. `ROADMAP.md` gets a decision entry for
   this behaviour change and `COMPAT.md` gets a migration notice. The evidence for FR-029
   is in decision 10.
4. **Displacement uses library code only (OP18, closed by C-019).** The vector comparison
   uses `FaultSegment.apply_to_points` with `steps = 20` on the solved frame and on an
   analytic frame. The scalar comparison uses `FaultSegment.evaluate_displacement`. Both
   frames use the same algorithm. The slip is 3 and the profile is the default
   `BaseFault` profile. The analytic frame is a `FaultSegment` that holds three
   `AnalyticReference` features (a class of `tests/fixtures/coordinate_frames.py`) with
   the R2 values and gradients. No new displacement law is added to the library (spec
   scope: out). The SC-004 numbers do not change (1332 points, 0.03).
5. **Error sources are assessed apart (FR-027).** The tests keep four populations and
   four checks apart. (a) Frame error: the S1 test and the new residual tests. (b)
   Integration error: two checks on analytic frames. The first compares the library
   result of the R2 analytic frame with a closed form written by hand,
   `(0, 0, -s * 1.5 * f(|x'|/20))` with `f(t) = 1 - 3 t^2 + 2 t^3`. The second uses a
   curved slip field (a circle) whose explicit-Euler result has a closed form. It checks
   that the library result for 20 steps equals that closed form, and that 400 steps come
   close to the exact arc. (c) Fault-plane samples: the 37 points with `x' = 0`. The test
   reports the largest difference and asserts finite output only. (d) Outside-support
   samples: points that the support does not contain. The comparison bound of SC-004
   applies to the 1332 smooth points only.
6. **Physical residuals use the interpolator (OP19, accepted by the maintainer).**
   `_value_residuals_fresh` and `_gradient_cosines_fresh` evaluate
   `self.interpolator.evaluate_value` and `self.interpolator.evaluate_gradient` at the
   world positions of the stored rows, with no fault step. They no longer call the
   feature methods. No fault acts a second time. For a feature with no faults and no
   regions the numbers do not change. This reverses S5 decision 8 for these two private
   methods only. The task text for `COMPAT.md` and the entry in `ROADMAP.md` both say so.
   The public names, the shapes and the stale rule (S5 decisions 9 to 11) do not change.
   A region mask of the feature no longer hides a row: before, a row outside a region
   gave NaN; now it gives a number. That is correct for a diagnostic of stored data. The
   public `evaluate_gradient` of a faulted feature keeps the tetrahedron path and its
   shape error (S3 note). S6 does not fix it. `COMPAT.md` lists it as a known issue.
7. **Backend matrix is a data table with a checker.** `tests/fixtures/backend_matrix.py`
   holds the cells, the optional dependencies and pure functions (`classify_cell`,
   `outcomes_from_junit`, `rows_without_real_cell`, `render_table`). A mandatory cell is
   passed or failed. A skipped or missing mandatory cell is a violation. An optional cell
   needs a named dependency. It is "excluded" only when the dependency is absent and the
   cell was skipped. It is a violation when it passed without the dependency, and when it
   was skipped with the dependency installed. The command
   `python -m tests.fixtures.backend_matrix <junit>...` reads one or more junit files,
   prints the table for `verification.md` and exits with 1 on a violation.
8. **Native and Surfe have two kinds of cell (maintainer answer 3).** The boundary cells
   (`native-boundary`, `surfe-boundary`) use a recording stub of the extension
   (`loop_interpolation_core` or `surfepy`) that the test puts in `sys.modules`. They are
   mandatory, they need no dependency, and they test the Python side of the boundary:
   local rows, grid, and world results. They are listed apart. The spec rows 6 and 7
   ("run when available; record absence explicitly") show their status from the real
   extension cell only. In this environment both rows show "excluded". The table prints
   the stub cells as separate lines. A stub cell is never reported as a native or Surfe
   pass. `main` of `backend_matrix` states this in its output. The test
   `test_every_spec_row_has_a_cell` does not count a stub cell for a row. The extension
   cells run the real extension. They are optional. Each one calls
   `require_optional(key)`, which skips with the text "optional dependency '<module>' is
   not installed".
9. **Names of the optional dependencies.** Native: module `loop_interpolation_core`
   (extra `loop-interpolation[native]`). Surfe: module `surfepy`. The pip package name is
   unverified (fact 21). Display: module `loopstructuralvisualisation`, `matplotlib`,
   `pyevtk`. A task checks the pip name of Surfe and decides if
   `test_p0_surfe_nans.py` changes its `importorskip` name to `surfepy` (maintainer
   answer 5).
10. **Golden review.** S6 expects no change to a golden file. The review has four parts:
    a ledger of all baseline files changed since `ebc7c699d5` with their evidence, a
    check that every other baseline file is unchanged, a tolerance audit of the diff,
    and a ledger of behaviour changes that have analytic evidence but no golden file
    (FR-029 evidence). The last ledger lists: the faulted-feature residuals (decision 6;
    evidence: the independent value in `test_diagnostics_backends.py`) and the four
    `EuclideanTransformation` fixes (decision 3; evidence: the hand-computed numbers in
    `test_euclidean_transformation.py` and `test_transformation.py`). It states that
    `git grep` finds no other caller, so 0 baseline files are affected. If an S6 test run
    changes a golden number, the S5 order applies: list, review, record, then `--update`.
11. **Final gates (SC-014).** The orchestrator runs them one process at a time. Each
    pytest run writes its own junit file in a scratch folder outside the repository
    (`<junit_dir>`):
    ```
    ruff check LoopStructural packages
    pytest --junitxml=<junit_dir>/root.xml
    pytest packages/loop_common/tests --junitxml=<junit_dir>/common.xml
    pytest packages/loop_interpolation/tests --junitxml=<junit_dir>/interpolation.xml
    pytest tests/unit/test_public_api_contract.py tests/unit/test_stable_api_surface.py
    python -m tests.fixtures.backend_matrix <junit_dir>/root.xml <junit_dir>/common.xml <junit_dir>/interpolation.xml
    ```
    The earlier merge gate ran the two package folders in one process. S6 splits them so
    that each run has one junit file. The two API tests also run inside the root run. All
    runs start in the repository root, so the rootdir and the `classname` paths agree
    (fact 19). A coverage ledger maps each production file changed since `ebc7c699d5` to a
    test.
12. **Disposition of the S5 open findings (run-state.json, 12 entries and 1 major).**
    Entries 1 and 10 are the same BOM finding, so 11 distinct minor findings exist.
    S6 closes 9 of them and the major one. 2 stay. See "S6 findings disposition".
13. **Documents.** `COMPAT.md` gets the S6 notices and the move of one S5 sentence.
    `ROADMAP.md` gets two decision entries (the faulted-feature residuals and the
    `EuclideanTransformation` fixes) and the line "Spec 005 complete". `specs/README.md`
    gets the S6 status.
14. **Second BOM file.** S6 removes the BOM from
    `tests/unit/modelling/test_fault_frame_analytic.py` as well. It is a 3-byte change
    and it adds that S1 file to "S6 files". The text of the file does not change. Reason:
    the BOM has the same cause and the same fix as finding 1.

#### S6 findings disposition

| Finding (S5 open list, `run-state.json`) | S6 result | How |
|---|---|---|
| 1 and 10. BOM in `_structured_grid_3d.py` | Closed | Remove the three bytes. Nothing else changes in the file. Check: first bytes are not `ef bb bf`; `git diff -w --stat` against `62c420e8` shows no change other than S5. The S1 file `test_fault_frame_analytic.py` has the same BOM and gets the same fix (decision 14) |
| 2. Format-only hunks in `_structured_grid_2d.py` | Stays | Cosmetic. A revert is more churn than the noise. Recorded in `verification.md` |
| 3. `set_nelements` passes a cell count to the node setter (`nelements=1000` gives 729 cells) | Stays | Old behaviour, not a regression. Recorded in `verification.md` and `COMPAT.md` |
| 4. T079 sentence in `verification.md` | Already closed | The sentence is at `verification.md` line 910. S6 records the check |
| 5. Exact-field stub ignores its points | Closed | The stub records the points and the test compares them with the world positions of the stored rows |
| 6. One `TOL_SOLVED` for both engines | Closed | One constant for each engine, 2 times the measured maximum. Each constant is at most the old value 9e-4 (a tightening or no change) |
| 7. PLI gradient test prints and does not assert | Closed | Assert `unattributed == 8` for the measured case |
| 8. `id(latest_diagnostics_report)` compare | Closed | Hold the object and compare with `is` |
| 9. `_stale_reason` checks one level of faults | Already closed in code | See fact 14. S6 adds no test |
| 11. `COMPAT.md` report row has the misfit sentence | Closed | Move the sentence to the `evaluate_value_misfit` row |
| 12. `verification.md` run-time count sentence (150 against 147 + 1) | Closed | Recompute from the two run records and correct the sentence |
| Major. Physical residuals of a faulted feature | Closed | OP19, decision 6 |

Earlier-stage minor findings that S6 also closes: S2 T020 ("no test of the native adapter
world contract") by the native boundary cells; S3 T038 ("the R2 fixture displaces none of
its probes") by the displacement fixture, whose slip is 3. Findings that stay: the S1 note
that `_on_geometry_change` adds anisotropy rows (it is in `_fault_builder.py`, outside the
S6 scope), and the S3 note on a shape error in the tetrahedron gradient path of
`_geological_feature.py` for a faulted feature. The OP19 fix avoids that path in the
diagnostics. The path itself does not change.

#### S6 contracts

| Boundary | Input | Output | Owner |
|---|---|---|---|
| `EuclideanTransformation.transform` and `.inverse_transform` | `(N, d)` array-like, `d >= dimensions` | new float array, same shape; transformed columns inverse to each other within `Tp`; other columns equal bit for bit; input unchanged; other shapes (also 1-D) raise `ValueError` | `loop_common/math/_transformation.py` |
| `EuclideanTransformation.fit` | `(N, d)` array-like, `d >= dimensions` | the object; translation is the mean of the first `dimensions` columns; the angle aligns the xy main axis with x; for 3D the z column keeps its sign; `ImportError` (no state change) if rotation is requested and scikit-learn is missing | `loop_common/math/_transformation.py` |
| `loop_common.math.EuclideanTransformation`, `LoopStructural.utils.EuclideanTransformation`, `LoopStructural.utils._transformation.EuclideanTransformation` | import | the same class object | `loop_common`, `LoopStructural.utils` |
| `GeologicalFeature.evaluate_gradient_cosines`, `evaluate_value_residuals`, the two misfit methods, `verification_report` | the feature | as in S5, with no fault step in the physical residuals | `_geological_feature.py` |
| `FaultSegment.apply_to_points`, `evaluate_displacement` | world points | unchanged. Zero gradient: no move, no warning. | `_fault_segment.py` (not changed) |
| `tests/fixtures/backend_matrix.py` | cells, one or more junit files | one status text for each cell; stub cells listed apart; violations named | tests |

#### S6 functions and signatures

`packages/loop_common/src/loop_common/math/_transformation.py`. The file keeps
`from __future__ import annotations`. `__init__`, `rotation`, `inverse_rotation`,
`fit_transform`, `__call__` and `_repr_html_` do not change. The logger import changes.
`fit`, `_rotation`, `transform` and `inverse_transform` change. The helper
`_check_points` and the function `_import_decomposition` are new.

```python
from loop_common.logging import get_logger          # replaces: from . import getLogger

logger = get_logger(__name__)


def _import_decomposition():
    """Return ``sklearn.decomposition``.

    Raises
    ------
    ImportError
        If scikit-learn is not installed. The message names the package and the
        install command ("pip install scikit-learn").
    """


class EuclideanTransformation:
    def _check_points(self, points: np.ndarray) -> np.ndarray:
        """Return a float copy of ``points`` with shape (N, d), d >= dimensions.

        Raises ValueError if ``points`` is not 2D ("points must be a 2D array with
        shape (N, d)") or has fewer than ``dimensions`` columns (message:
        "Points must have at least {n} dimensions").
        """
        array = np.array(points, dtype=float)
        if array.ndim != 2:
            raise ValueError("points must be a 2D array with shape (N, d)")
        if array.shape[1] < self.dimensions:
            raise ValueError(f"Points must have at least {self.dimensions} dimensions")
        return array

    def fit(self, points: np.ndarray) -> "EuclideanTransformation":
        """Fit the translation and the angle.

        The translation is the mean of the first ``dimensions`` columns. If
        ``fit_rotation`` is true, the angle is ``-arctan2(c_y, c_x)`` for the main axis
        ``(c_x, c_y)`` of the x and y columns (PCA with 2 components). The rotation
        turns the xy plane about z, so the z column is not part of the axis. The
        method changes no attribute if it raises.
        """
        array = self._check_points(points)
        d = self.dimensions
        translation = np.mean(array[:, :d], axis=0)
        angle = 0
        if self.fit_rotation:
            decomposition = _import_decomposition()
            centred = array[:, :2] - translation[None, :2]
            axis = decomposition.PCA(n_components=2).fit(centred).components_[0]
            angle = -np.arctan2(axis[1], axis[0])
        self.translation = translation
        self.angle = angle
        logger.debug("EuclideanTransformation fitted: translation=%s angle=%s", translation, angle)
        return self

    def _rotation(self, angle):
        """Return the 3x3 rotation about z; entry [2, 2] is +1 (it was -1)."""
        return np.array(
            [
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1],
            ]
        )

    def transform(self, points: np.ndarray) -> np.ndarray:
        array = self._check_points(points)
        d = self.dimensions
        centred = array[:, :d] - self.translation[None, :]
        array[:, :d] = np.einsum("ik,jk->ij", centred, self.rotation[:d, :d])
        return array

    def inverse_transform(self, points: np.ndarray) -> np.ndarray:
        array = self._check_points(points)
        d = self.dimensions
        array[:, :d] = (
            np.einsum("ik,jk->ij", array[:, :d], self.inverse_rotation[:d, :d])
            + self.translation[None, :]
        )
        return array
```

The planner checked the formulas by hand (see "S6 tests"). The 2D and 3D matrices use
`centred @ R.T`, as the old code did. So a point `p` gives `R (p - t)` and the inverse
gives `R^T q + t`. For `dimensions=2` the old and the new `transform` give the same
numbers, because entry `[2, 2]` is not used.

`packages/loop_common/src/loop_common/math/__init__.py`: add the line
`from ._transformation import EuclideanTransformation` after `from ._maths import *`. The
name `EuclideanTransformation` is free in this module (`git grep -n -w` on 2026-11-06).
No import cycle arises (fact 2).

`LoopStructural/utils/_transformation.py` (whole new content):

```python
"""Old import path of ``EuclideanTransformation``.

The class lives in ``loop_common.math``. This module keeps the path
``LoopStructural.utils._transformation`` working.
"""
from loop_common.math._transformation import EuclideanTransformation

__all__ = ["EuclideanTransformation"]
```

`LoopStructural/modelling/features/_geological_feature.py` (two changes, nothing else).

```python
# in _gradient_cosines_fresh, replace:  model = self.evaluate_gradient(positions)
model = self.interpolator.evaluate_gradient(positions)

# in _value_residuals_fresh, replace:   self.evaluate_value(positions)
return np.abs(rows[:, 3] - self.interpolator.evaluate_value(positions))
```

Both docstrings get this sentence: "The stored positions are already restored by the
faults of the feature. The method evaluates the interpolator, not the feature, so no
fault acts a second time." The `positions` stay `bounding_box.reproject(...)` of the
stored rows when the interpolator has a frame, as in S5. `_stale_reason`,
`_physical_residual_summaries`, `verification_report` and the public methods do not
change.

`tests/fixtures/coordinate_frames.py` (changed; new names, all free by `git grep -n -w`
on 2026-11-06):

```python
R2_SLIP = 3.0
R2_STEPS = 20
R2_DISPLACEMENT_TOL = 0.03          # SC-004: 1 percent of the slip
CURVED_RADIUS = 5.0                 # start radius of the curved slip control
CURVED_GX = 0.2                     # constant coordinate-0 value of the curved control

def r2_displacement_section(anchor: np.ndarray) -> np.ndarray:
    """Return the 1369 section points y'=0, x' and z' in {-9, -8.5, ..., 9} (world)."""

def r2_section_masks(anchor: np.ndarray, section: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (smooth, plane) boolean masks. smooth is |x'| >= 0.5 (1332), plane is the rest (37)."""

def r2_analytic_fault(anchor: np.ndarray, steps: int = R2_STEPS,
                      displacement: float = R2_SLIP, gradient_override=None) -> FaultSegment:
    """Return a FaultSegment of three AnalyticReference features with the R2 fields.

    ``gradient_override(coordinate, points, gradient)`` may replace the gradient of
    one coordinate. It gets the current points of each call, not the start points.
    The default profile is BaseFault. ``anchor`` shifts the reference only.
    """

def r2_closed_form_displacement(section: np.ndarray, anchor: np.ndarray,
                                slip: float = R2_SLIP) -> np.ndarray:
    """Return the displacement vectors (N, 3) of the R2 analytic frame, by hand.

    For x' > 0 the vector is (0, 0, -0.5 * slip * f(t)); for x' < 0 it is
    (0, 0, +0.5 * slip * f(t)); for x' = 0 it is zero. Here t = |x'| / 20 and
    f(t) = 1 - 3 t^2 + 2 t^3. The code does not call the library profile.
    """

def r2_curved_slip_fault(anchor: np.ndarray, steps: int) -> FaultSegment:
    """Return a FaultSegment whose slip direction turns about a vertical axis.

    Coordinate 0 has the constant value CURVED_GX. Coordinate 2 has the constant
    value 0. Coordinate 1 has the gradient (-(y-cy), (x-cx), 0) / (20 r) about the
    axis through ``anchor`` (cx, cy). The displacement is R2_SLIP and the profile is
    BaseFault. By hand: d = 0.5 * 3 * f(0.2) = 1.344 (f(0.2) = 0.896).
    """

def euler_circle_closed_form(radius: float, d: float, steps: int) -> tuple[float, float]:
    """Return (final radius, turn angle) of the explicit Euler scheme on a circle.

    Each step has the length s = d / steps along the tangent. The radius obeys
    r_k^2 = radius^2 + k s^2. The angle gains arctan(s / r_k) in step k. The function
    sums the angle with a plain loop of ``steps`` terms and does not call the library.
    """

def build_r2_fault(anchor: np.ndarray, rotation_deg: float, interpolatortype: str,
                   explicit_center: bool = True, displacement: float = R2_SLIP) -> FaultSegment:
    """Build the R2 fault as test_fault_frame_analytic._fault does, with displacement 3.

    Copy of ``_fault`` (test_fault_frame_analytic.py lines 21 to 51). The only
    changes are ``displacement`` in place of 0 and ``fault.steps = R2_STEPS``.
    The S1 test file keeps its code.
    """

def faulted_foliation_model(anchor: np.ndarray, rotation_deg: float, engine: str):
    """Return (model, fault, strati) of the S5 faulted model, shifted to ``anchor``.

    Same data as ``_faulted_built_model`` in test_diagnostics_frames.py (unit cube):
    3 fault rows at (0.5, 0.5, 0.5) (normal (1, 0, 0) with value 0, slip direction
    (0, 0, 1), extent direction (0, 1, 0)); 8 strati value rows at z = 0.3 (value 0)
    and z = 0.7 (value 1) for x and y in {0.2, 0.8}; 1 strati normal row at
    (0.3, 0.3, 0.5) with the vector (0, 0, 1). The fault displacement is 0.1.
    Every location is shifted by ``anchor``. The model uses
    ``model_with_frame(anchor, anchor + 1, rotation_deg)``. Both features use
    ``interpolatortype=engine`` ("FDI" or "PLI") and ``nelements=1e3``. The model is
    built (``model.update()``) before it is returned.
    """

def residuals_of_field(interpolator, nodal_values: np.ndarray) -> dict[str, np.ndarray]:
    """Return ``matrix @ nodal_values - b`` for each entry of ``interpolator.constraints``."""
```

`tests/fixtures/backend_matrix.py` (new; names free). Task 2 of the order writes the
helper part. A later task writes the `MATRIX` list:

```python
Outcome = Literal["passed", "failed", "skipped", "absent"]

@dataclass(frozen=True)
class OptionalDependency:
    key: str            # "native", "surfe", "viewer", "matplotlib", "pyevtk"
    module: str         # top-level module name for find_spec
    distribution: str   # pip name or extra, for the message
    capability: str

OPTIONAL_DEPENDENCIES: dict[str, OptionalDependency]

@dataclass(frozen=True)
class OptionalStatus:
    key: str
    module: str
    available: bool
    reason: str         # "" when available, else "optional dependency '<module>' is not installed (<distribution>)"

def optional_status(key: str, find_spec: Callable[[str], object | None] = importlib.util.find_spec) -> OptionalStatus
def require_optional(key: str) -> None          # pytest.skip(status.reason) when not available

@dataclass(frozen=True)
class Cell:
    cell_id: str
    spec_row: str                       # one of SPEC_ROWS, or "extra"
    mandatory: bool
    stub: bool                          # True: boundary cell with a recording stub, never a real result
    refs: tuple[tuple[str, str, str | None], ...]   # (path with "/", test function, id selector or None)
    min_cases: int                      # least number of junit cases that the refs must match
    dependency: str | None              # key of OPTIONAL_DEPENDENCIES; None for a mandatory cell

SPEC_ROWS: tuple[str, ...]              # the 7 rows of "Required backend coverage"
MATRIX: tuple[Cell, ...]

def classify_cell(cell: Cell, outcome: Outcome, status: OptionalStatus | None) -> str
def outcomes_from_junit(paths: Sequence[str]) -> dict[str, Outcome]    # cell_id -> outcome
def rows_without_real_cell(matrix: Sequence[Cell]) -> list[str]        # rows with no cell where stub is False
def render_table(matrix: Sequence[Cell], statuses: dict[str, str]) -> str
def main(argv: list[str] | None = None) -> int     # argv: junit files; prints the table; 1 on a violation
```

`classify_cell` returns: `"passed"`, `"failed"`, `"excluded: <reason>"` or
`"VIOLATION: <text>"`. Rules: a mandatory cell that is skipped or absent is a
violation. An optional cell with an absent dependency is "excluded" when it is skipped,
a violation when it passed or is absent. An optional cell with an installed dependency is
a violation when it is skipped or absent. A failed cell is "failed" in all cases.

`outcomes_from_junit` reads every file in `paths` and merges the cases. The key of a
case is the pair (`classname`, `name`). The junit family is `xunit2`, which has no `file`
attribute (fact 19). If one key occurs in more than one file, the worse outcome wins
(failed, then skipped, then passed). A cell ref `(path, function, selector)` has the
module `path` without `.py` and with dots in place of `/`. A case matches the ref when
`classname` equals the module or starts with the module and a dot (a test class), the
part of `name` before `[` equals `function`, and the selector is `None` or occurs in
`name`. A module-level skip is one skipped case with an empty `classname` and a `name`
that equals the module. It is matched by file only: it matches every ref of that file
and counts as "skipped". For an optional cell with an absent dependency this gives
"excluded". It never gives "absent". For a mandatory cell it is a violation. The outcome
of a cell: one failed case makes it failed; else one skipped case makes it skipped; else
it needs at least `min_cases` passed cases to be passed; otherwise it is absent. The
`--import-mode=importlib` option of the root `pyproject.toml` does not change this rule.
The node id gives the `classname` and the node id does not depend on the import mode
(fact 19).

#### S6 tests

All cells are parametrised test ids. A tolerance is given with its reason. "Non-vacuity
check" and "measure first" are defined at the top of this stage design.

`tests/unit/utils/test_transformation.py` (changed; the 9 other existing tests stay).
`Tp(L, W)` is `tp` of the S1 fixture with `L = max(extent of the cloud, 1.0)` and
`W = max |coordinate|`. Reason: for one point the extent is 0 and `Tp` would only keep the
rounding part. Anchors are the four `ANCHORS`. The cloud is `anchor + rng.uniform(-50, 50)`
(seed 0). Hand numbers are in the "Checks" column. The anchor `(500000, 6500000, 100)` has
`W = 6.5e6`. A float64 number of that size has the spacing 9.3e-10. So a bound of 1e-12 is
not possible for any comparison that passes a world position through a local position.
Every bound below uses `Tp` or a value derived from it.

| Test function | Cells | Checks |
|---|---|---|
| `test_round_trip_explicit_parameters` | dimensions {2, 3} x n {1, 4, 50} x 4 anchors x columns {3, 5} x angle {30 deg, -75 deg} = 96 | transform with a set angle and translation `anchor + (1, 2, 3)`; forward result differs from the input by more than `0.1` in some transformed column; for the x and y columns the forward result equals the closed form `x' = c(x - tx) - s(y - ty)`, `y' = s(x - tx) + c(y - ty)` computed in the test, within `1e-12 * L` (the centred values are small and the formulas agree to a few `eps`); for `dimensions=3` the z column equals `z - tz` within the same bound (this fails before S6: z was reversed); `inverse_transform(forward)` equals the input within `Tp` in the transformed columns; every other column is `np.array_equal` to the input after forward and after inverse; the input array is equal to its copy; the output does not share memory with the input |
| `test_round_trip_after_fit` | 4 anchors x n {1, 4, 50} x `fit_rotation` {True, False} = 24 | `dimensions=2`. `fit` on a 50-point cloud along a line at 30 degrees, then forward and inverse of the first n points: same checks as above; the fitted translation equals `np.mean` of the cloud within `Tp` |
| `test_integer_input_is_not_truncated` | 1 | an `int64` array of 4 points: the forward result equals the closed form computed with floats within `Tp` (before S6 the result is truncated) |
| `test_one_dimensional_input_raises_value_error` | 3 (transform, inverse_transform, fit) | shape `(3,)` raises `ValueError` with the text "2D array"; before S6: `IndexError` |
| `test_too_few_columns_raises_value_error` | 2 (transform, inverse) with `dimensions=3` and 2 columns | `ValueError` with "at least 3 dimensions" |
| `test_entry_points_are_one_class` | 1 | `LoopStructural.utils.EuclideanTransformation`, `LoopStructural.utils._transformation.EuclideanTransformation`, `loop_common.math.EuclideanTransformation` and `loop_common.math._transformation.EuclideanTransformation` are the same object. Fails before S6 (the last two cannot be imported) |
| `test_pickle_with_the_old_module_path_loads` | 1 | the test sets the class attribute `__module__` to `"LoopStructural.utils._transformation"` with `monkeypatch`, pickles a fitted object, restores the attribute, and loads the bytes. The result is an instance of the class with equal `translation`, `angle`, `dimensions` and `fit_rotation`. This shows that a pickle from before S6 loads |
| `test_new_pickle_names_the_loop_common_path` | 1 | the bytes of a pickle made now contain `loop_common.math._transformation`. This documents the one-way compatibility |
| `test_fit_aligns_a_line_with_negative_slope` | 2 (line at -30 and -10 degrees through (5, 5)) | `dimensions=2`, 50 points as in `_line_points`. By hand: the axis is `(cos 30, -sin 30)`; the rotation by +30 degrees maps it to `(1, 0)`. So `t.angle` equals `+pi/6 = 0.5235987756` for -30 degrees and `+pi/18 = 0.1745329252` for -10 degrees, within `1e-9` modulo `pi` (scikit-learn may return the opposite axis for the same line, which turns the line onto the negative x axis); the variance of the y column after `fit_transform` is below `1e-20`. Before S6 (measured, scikit-learn 1.9.0): the angle is `-pi/6` and the y variance is 6.51 (for -10 degrees: 1.01) |
| `test_fit_three_dimensions` | 2 (line at +30 and -30 degrees) | `dimensions=3`, 50 points of `_line_points` (z = 0.1 i). By hand: the mean of z is `0.1 * 49 / 2 = 2.45`; the means of x and y are 5. So the translation equals `(5, 5, 2.45)` within `1e-9`; the angle is `-pi/6` for +30 degrees and `+pi/6` for -30 degrees within `1e-9` modulo `pi`; the transformed z equals `z - 2.45` within `1e-12`; the y variance is below `1e-20`; the inverse returns the input within `Tp`. Before S6: `fit` raises `ValueError` (the shapes (3,) and (2,) do not match; by reading) |

The test `test_inverse_transform_is_broken_for_normal_point_clouds_bug` is deleted (fact 3,
FR-029). `verification.md` records the deletion under the golden review.
The tolerance reason for `Tp`: the round trip is two rotations and one translation of
float64 values of size `W`. The error is about `W * eps`. `Tp` has `32 * eps * W`. So the
margin at `W = 6.5e6` is more than 10 times.

`packages/loop_common/tests/test_euclidean_transformation.py` (new). It needs no
scikit-learn, so it is a mandatory cell of `loop_common` alone. It does not import
`tests.fixtures`, because the package tests run alone. It holds a local helper
`_tp(length, width)` with the same formula as `tests/fixtures/coordinate_frames.tp`
(`1e-10 * length + 32 * eps * max(width, length)`).

| Test function | Cells | Checks |
|---|---|---|
| `test_module_imports` | 1 | `importlib.import_module("loop_common.math._transformation")` works (before S6: `ImportError`) |
| `test_export_from_loop_common_math` | 1 | `loop_common.math.EuclideanTransformation` is the class of the module above |
| `test_round_trip_without_fit` | dimensions {2, 3} x n {1, 4, 50} on the anchor `(500000, 6500000, 100)` = 6 | same round-trip checks as the root test, with fixed angle and translation, and the bound `_tp` |
| `test_two_dimensional_rotation_of_known_angle` | 1 | `dimensions=2`, angle `pi/6`, translation zero, the point `(2, 0, 7)`. By hand: `(2 cos 30, 2 sin 30, 7) = (1.7320508076, 1.0, 7)`. The result equals this within `1e-12`; the third column is exactly 7 |
| `test_three_dimensional_rotation_keeps_the_z_direction` | 2 points | `dimensions=3`, angle `pi/2`, translation `(1, 2, 3)`. By hand: the point `(2, 2, 5)` is centred to `(1, 0, 2)` and gives `(0, 1, 2)`; the point `(1, 3, 3)` is centred to `(0, 1, 0)` and gives `(-1, 0, 0)`. Both equal the result within `1e-12`. The determinant of `rotation[:3, :3]` equals 1 within `1e-12`. `inverse_transform` returns the inputs. Before S6: the first point gives z = -2 and the determinant is -1 |
| `test_fit_without_scikit_learn_raises_import_error` | 1 | `monkeypatch.setitem(sys.modules, "sklearn", None)`. `fit` with `fit_rotation=True` raises `ImportError` whose text contains "scikit-learn". `translation` and `angle` are equal before and after. Before S6: `fit` logs and returns `None` |
| `test_fit_without_rotation_needs_no_scikit_learn` | 1 | the same block of `sklearn`. `fit_rotation=False` on `[[0, 0, 5], [2, 4, 6]]`. By hand: the translation is `(1, 2)` and the angle is 0. Before S6: `fit` returns `None` and sets nothing |

`tests/unit/modelling/test_fault_frame_validation.py` (new; depends on the fixture task).
A module-scoped fixture `solved_r2` with `params=CELLS` builds one solved fault for each
cell with `build_r2_fault`. A second module-scoped fixture `solved_r2_inferred` has the
2 inferred-centre cells. Pytest keeps one param instance of a module-scoped fixture alive
at a time, so only one solved fault is in memory. `CELLS` is 4 anchors x 3 rotations x
2 engines (FDI, PLI) = 24, with the explicit centre. Every test below that says 24 cells
uses `solved_r2`, so the file does 26 R2 solves and no more. The section is
`r2_displacement_section(anchor)`. The analytic frame is `r2_analytic_fault(anchor)`.

| Test function | Cells | Checks |
|---|---|---|
| `test_vector_displacement_matches_analytic_frame` | 24 | `solved.apply_to_points(section)` and `analytic.apply_to_points(section)`. On the 1332 smooth points: both results finite, and the largest norm of the difference (`np.linalg.norm(diff, axis=1).max()`) is at most `0.03` (SC-004; 1 percent of slip 3). Non-vacuity check: the analytic displacement norm is at least `0.5` for at least 1200 of the 1332 smooth points (by hand: `d = 1.5 * f(t)` with `f(t) = 1 - 3t^2 + 2t^3`, `t = abs(x')/20` in `[0.025, 0.45]`, so `d >= 1.5 * 0.57475 = 0.86` on both walls; the footwall profile is the mirror) |
| `test_scalar_displacement_matches_analytic_frame` | 24 | same for `evaluate_displacement`, bound `0.03`, same non-vacuity check |
| `test_displacement_populations_are_counted_separately` | 24 | the section has 1369 points; the smooth mask has 1332 and the plane mask has 37 (counts from the generator, not from the solved result); the solved result is finite on all 1369 points, so the invalid-support count is 0; the three counts add up to 1369. For the 37 fault-plane points the test computes the largest difference between the solved and the analytic result, prints it, and asserts only that the solved output is finite (the hanging-wall flag `gx > 0` can switch there, so no bound applies). Precondition assertion: the support contains the analytic path (largest `abs(z')` is `9 + 1.5 = 10.5`) in world coordinates. If not, the test fails with a message instead of hiding NaN points in the invalid count |
| `test_vector_displacement_inferred_centre` | 2 (anchor 3, rotation 30, FDI and PLI) | as the first test, with `explicit_center=False` |
| `test_analytic_frame_displacement_matches_closed_form` | 4 anchors | `r2_analytic_fault(anchor).apply_to_points(section)` minus the section equals `r2_closed_form_displacement(section, anchor)`. The x and y components are exactly 0 (`np.array_equal`). The z component matches within `64 * eps * max(1, max abs(z))` (20 steps add the same vector; each addition rounds by at most `eps/2 * abs(z)`, so the sum is below `20 * eps/2 * abs(z)`; the factor 64 is the margin). The 37 plane points move by exactly 0. This check does not depend on the step count, because `gx`, `gz` and the direction are constant along the path. It tests the library profile against the independent hand formula |
| `test_integration_error_follows_the_euler_closed_form` | 4 anchors | the curved control `r2_curved_slip_fault(anchor, steps)` with 8 start points on a circle of radius `CURVED_RADIUS = 5` about `anchor` (angles k * 45 degrees, same z as the anchor). By hand: `d = 1.344`. For `steps = N` each step has the length `s = d/N` along the tangent. So `r_N^2 = 25 + d^2/N`, and the turn angle is the sum of `arctan(s / r_k)` with `r_k^2 = 25 + k s^2`. (a) The library result for `N = 20` equals `euler_circle_closed_form` within `Tp(10, W)`. (b) The library result for `N = 400` differs from the exact arc (radius 5, angle `d/5`) by at most `d^2 / (N * 5) + Tp(10, W) = 9.0e-4 + Tp` (the radial term is `d^2 / (2 * N * 5) = 4.5e-4` and the angle term is smaller). (c) The results for `N = 20` and `N = 400` differ by more than `5e-3` (by hand `r_20 = sqrt(25 + 1.806/20) = 5.00902` and `r_400 = sqrt(25 + 1.806/400) = 5.00045`, so the radial gap is 8.6e-3). So the test fails if the step count has no effect |
| `test_misplaced_frame_fails_the_displacement_bound` | 24 | the analytic frame shifted by 1 world unit along x' (`r2_analytic_fault(anchor + (1, 0, 0))`) against the solved frame: the largest smooth-point difference is more than `0.03`. By hand: shift 1 moves `gx` by 0.05; at `t = 0.45` the profile slope is -1.485, so the difference is about `1.5 * 1.485 * 0.05 = 0.11`. This shows that the bound can fail |
| `test_outside_support_points_have_no_displacement` | 24 | six probes `anchor +- (H + 5) e_axis` where `H` is the largest world offset of the 8 support corners on that axis. Precondition: the support does not contain them. `evaluate_value` and `evaluate_gradient` of each coordinate are NaN; `evaluate_displacement` is `0.0`; `apply_to_points` returns the points unchanged (`np.array_equal`). If `apply_to_points` raises for an empty mask, the task reports the exception and tests `evaluate_displacement` only |
| `test_zero_gradient_protection_remains_effective` | 4 anchors x variant {zero, nan} = 8 | analytic fault whose coordinate-1 gradient is `(0, 0, 0)` (or NaN) at each call for the points with current `z' > 5` (`gradient_override` gets the current points). The test defines three sets by the **start** position `z0'`. Set A: `z0' > 5`. Set B: `z0' <= 3.5`. Set C: `3.5 < z0' <= 5`. Reason: a footwall point moves up by at most 1.5, so a point of set B never reaches `z' = 5`; a point of set C can cross it during the steps. Zero variant, under `warnings.simplefilter("error")`: `apply_to_points` finishes with no warning; all outputs are finite; set A is `np.array_equal` to its start; set B is `np.array_equal` to the result without an override. Set C is only finite. NaN variant (no warning filter, because NumPy may warn about NaN for other reasons): the outputs of set A are NaN and the outputs of set B are finite and equal to the result without an override. Set C is not asserted. The zero variant fails if the `g_mag > 0` guard goes (0/0 gives a `RuntimeWarning`) |
| `test_fault_frame_assembled_residuals_vanish_on_the_analytic_field` | 24 | for each coordinate `i`: set the nodal values to the analytic field (`r2_frame` at `local_to_world_points(support.nodes)`) and compute `residuals_of_field`. Every constraint name must be in one of three groups, else the test fails and names it (non-vacuity check: no family is ignored silently). Group 1, must vanish: FDI `value`, `gradient`, `norm_0`, `norm_1`, `norm_2`, `dxx`, `dyy`, `dzz`, `dxy`, `dxz`, `dyz`, `dx_lower`, `dx_upper`, `dy_lower`, `dy_upper`, `dz_lower`, `dz_upper`; PLI `value`, `norm`, `edge jump` (fact 18; a numeric suffix is allowed). Their maximum residual is at most `TAU_ASSEMBLED`. Measure first. Cap `1e-9`. Reason for the cap: the field is affine, the stencils are exact for affine fields, the rows are normalised, and the nodal values carry the world rounding of the node positions (up to `4.7e-10 / 20 = 2.3e-11` at the anchor `(500000, 6500000, 100)`), so the residual is a few times 1e-11. SC-008 uses `1e-10` for the same kind of row. Group 2, reported only with a measured bound `TAU_ORTHOGONAL`: the names that end in `_orthogonal`. Their direction comes from a solved neighbour field, so they do not vanish on the analytic field. The residual is about the sine of the neighbour direction error (SC-003 allows 0.01 degree = 1.7e-4). Cap `1e-3`. Group 3, no entry may exist: any other name. The test writes the constraint names and row counts of the three coordinates as literals (non-vacuity check: each coordinate has at least one `value`, `gradient` or `norm` row, and at least 1 row in group 1) |
| `test_fault_frame_physical_residuals_match_analytic` | 24 | for each coordinate: `fault[i].verification_report().physical_residuals`. Coordinate 0: the `value` entry has 25 rows (the data rows at `x' = 0`) and `max <= 1e-6` (SC-003 scalar limit) and equals `abs(0 - fault[0].evaluate_value(rows))` within `max(1e-12, Tp(20, W) / 20)`. Reason: the field gradient is `1/20`, and the world position of a stored row passes through a local position, which moves it by about `Tp`. At `W = 6.5e6` this is 2.4e-9. Every direction family (`gradient`, `normal`): the largest angle is at most `0.01` degree (SC-003). The test also computes the angle by an independent formula, `degrees(arctan2(norm(cross(a, b)), dot(a, b)))` for the stored world direction `a` and the analytic direction `b`, and asserts the reported maximum equals the largest of these values within `1e-6` degree. The two angles come from different formulas, so the check is not the same arccos twice. Literal row counts as above |
| `test_residual_check_detects_misplaced_reference` | 24 | the assembled residual of coordinate 0 with the analytic field shifted by 1 unit along x' has a value-constraint maximum above `1e-2`. By hand: the shift moves the field by 0.05 at the 25 data rows. This shows that the residual check can fail |

`tests/unit/modelling/test_diagnostics_backends.py` (new; depends on the fixture task).
It holds the faulted-feature cells (OP19) and the P2 cells (OP22). Two module-scoped
parametrised fixtures build the models: `faulted_model` (the 24 cells of 4 anchors x 3
rotations x 2 engines, with `faulted_foliation_model`) and `p2_model` (the 12 cells of 4
anchors x 3 rotations, with `world_foliation_model`). The file builds 24 + 12 = 36 models.
Without the shared fixtures it would build one model for each test and cell (about 120).
No test in the file changes a model. Each test uses the methods that S5 proved to be
read-only. Test 5 checks that.

| Test function | Cells | Checks |
|---|---|---|
| `test_faulted_value_residuals_use_restored_positions` | 24 | `strati.evaluate_value_residuals()` equals the independent value `abs(val - interpolator.support.evaluate_value(stored local positions, interpolator.c))` within `max(1e-12, 3 * Tp(1, W))`. Reason: the feature path takes the stored position to world and back to local, which moves it by about `Tp`; the field gradient is below 3 (values 0 and 1 over 0.4 give 2.5); so the value moves by at most `3 * Tp`. The residual is at most `TOL_FAULTED_VALUE` (measure first with the independent function; cap `1e-2`; S5 measured 0.0004 to 0.0018 against 0.080 before the fix). Spy: a monkeypatch counts calls of `strati._apply_faults`, `strati.evaluate_value` and `strati.evaluate_gradient`; each count is 0 |
| `test_faulted_gradient_cosines_use_restored_positions` | 24 | `strati.evaluate_gradient_cosines()` equals the independent cosines (`interpolator.evaluate_gradient` at the stored world positions against the stored world vectors) within `max(1e-12, 10 * Tp(1, W))`. Reason: the position moves by about `Tp`, the cell size is 0.1, and the gradient of the field changes by at most its own size over one cell, so a cosine moves by at most `Tp / 0.1`. `max(1 - cos)` is at most `TOL_FAULTED_COS` (measure first; cap `1e-2`). Same spy |
| `test_fault_restores_the_stored_rows` | 24 | non-vacuity check: at least 4 stored rows differ from the raw world data by more than `1e-3` (the fault restores them). Write the measured number as a literal |
| `test_report_physical_value_matches_independent_rms` | 24 | `verification_report().physical_residuals["value"].rms` equals the RMS of the independent residuals within `max(1e-12, 3 * Tp(1, W))` (same reason as the first test), and `.kind == "physical"` |
| `test_p2_misfit_leaves_state_unchanged` | 12 (P2 foliation of `world_foliation_model`) | the S5 snapshot of `data`, `c`, flags, weights and the report object (compare with `is`) is equal after the two old methods, the two new methods and `verification_report()`. The P2 fit is not tested against the analytic field (issue #399) |
| `test_p2_report_blocks_are_separate` | 12 | the physical and assembled blocks are different objects with `kind` `"physical"` and `"assembled"`; the physical `value` rms equals the RMS of `evaluate_value_residuals()` within `1e-12` (both come from the same arrays with the same arithmetic, for any anchor) |

Pre-fix evidence for this file. The first, second and fourth tests fail on the current
tree (the spy counts are above 0, and the residual is about 0.08). The test for the
restored rows and the P2 tests are guards. They pass before the fix.

`tests/unit/modelling/test_diagnostics_frames.py` (changed; findings 5 to 8). The
exact-field test patches `feature.interpolator.evaluate_gradient` in place of
`feature.evaluate_gradient`, records the points, and asserts that they equal the world
positions of the stored rows within `Tp`. This test is a guard. For a feature without
faults, `GeologicalFeature.evaluate_gradient` ends in `self.interpolator.evaluate_gradient`
(fact 13). So the changed test passes before and after the fix. The old form of the test
(a patch of `feature.evaluate_gradient`) would fail after the fix, because the fixed code
no longer calls that method. So the change of the patch target must land with the fix
task or before it. A search of `tests` and `packages/*/tests` (2026-11-06) finds no other
test that patches `evaluate_gradient` or `evaluate_value`. `TOL_SOLVED` becomes
`TOL_SOLVED_FDI` and `TOL_SOLVED_PLI` (2 times the measured maximum of each engine; the S5
record has FDI 1.2e-4; each constant is at most the old 9e-4, so the change is a tightening
or no change). The PLI gradient test asserts `unattributed == 8`. The snapshot keeps the
report object and compares it with `is`.

`tests/unit/interpolator/test_backend_world_contract.py` (new; mandatory; no skip; depends
on the fixture task). Fixture `stub_core` puts a recording `loop_interpolation_core`
module (class `Problem`) in `sys.modules` with `monkeypatch.setitem`. Its `Problem.solve`
returns the exact R1 nodal field and `{"converged": True, "iterations": 1, "residual": 0.0}`.
Fixture `stub_surfe` puts a recording `surfepy` module (class `Surfe_API`) in
`sys.modules`, removes `loop_interpolation._surfe_wrapper` from `sys.modules` before the
import, and at teardown restores `sys.modules["surfepy"]`,
`sys.modules["loop_interpolation._surfe_wrapper"]` and the attribute `_surfe_wrapper` of the
package `loop_interpolation` to their earlier state (absent stays absent). The teardown
asserts the restoration.

| Test function | Cells | Checks |
|---|---|---|
| `test_native_adapter_sends_local_rows_and_grid` | 4 anchors x 3 rotations = 12 | R1 FDI interpolator (`r1_interpolator("fdi", anchor, rotation, exact=True)`) with 4 gradient rows added in world coordinates (the offsets of `world_foliation_model`). `_native.solve` runs with the stub. The `add_value` points equal `bounding_box.project(world observations)` within `Tp`; the `add_gradient` vectors equal `world_to_local_vectors(unit R1 gradient)` within `1e-12` (a direction passes through a rotation only, so the anchor does not matter); `Problem` gets `support.origin`, `support.nsteps` and `support.step_vector`. For rotation 30 and 90 the local points differ from the world points by more than 1 (non-vacuity check) |
| `test_native_solved_field_evaluates_in_world` | 12 | after `solve`, `interp.c` is the stub field and `interp.evaluate_value(world query)` equals `r1_value` within `1e-8 * F` (SC-002); `evaluate_gradient` equals the R1 gradient within `1e-8 * max(1, abs(g))`. The field is exact, so this tests the conversion only |
| `test_native_unsupported_settings_keep_rejection_and_state` | 3 (rotated support `rotation_xy`, region, `P1Interpolator`) | `NotImplementedError` whose text contains `rotation_xy`, `region` or `interpolator class`; `c`, `data` and `up_to_date` are equal before and after |
| `test_surfe_adapter_sends_local_points` | 12 | the R1 interpolator is the Surfe wrapper with `interp.bounding_box = box`. `setup_interpolator()` sends the interface points (27) equal to `project(world observations)` within `Tp`; the planar rows have local vectors |
| `test_surfe_adapter_returns_world_values_and_gradients` | 12 | the stub's `EvaluateInterpolantAtPoints(local)` returns `r1_value(interp.local_to_world_points(local))` and `EvaluateVectorInterpolantAtPoints(local)` returns `world_to_local_vectors(R1 gradient)`. Public `evaluate_value(world)` equals `r1_value` within `1e-8 * F` and `evaluate_gradient` equals the R1 gradient within `1e-8 * max(1, abs(g))` (the SC-002 limits). Reason: the stub gets a local point that was produced from a world point, and returns a value at the world point that it computes back from it. Two conversions at the size 6.5e6 move the position by about 1e-9 (two half-spacings of 4.7e-10). The R1 gradient has the norm 0.54, so the value moves by about 5e-10. `F` is the field range over the box (about 5 or more), so the limit `1e-8 * F` is 5e-8 or more. The margin is about 100. The task records the measured maximum error. A second conversion or a missing conversion makes the stub return wrong values by far more than this |

`tests/unit/interpolator/test_backend_optional_native_surfe.py` (new; optional cells;
depends on the fixture task and on the helper part of `backend_matrix.py`). Each test
starts with `require_optional("native")` or `require_optional("surfe")`.

| Test function | Cells | Checks |
|---|---|---|
| `test_native_extension_r1_fit_in_world` | 4 anchors x 3 rotations = 12 | R1 fit (`exact=False` support, 27 value rows) solved with `solver="native"`, `solver_kwargs={"tol": 1e-12}`. World value and gradient errors at the two query offsets: `1e-8 * F` and `1e-8 * max(1, abs(g))` (SC-002; OP21). It also equals the `cg` result within `1e-6 * F` |
| `test_surfe_extension_r1_values_and_gradients_in_world` | 12 | Surfe wrapper with `poly_order=1` and the default kernel; 27 interface rows. Same SC-002 limits at the two query offsets (OP21) |

Existing optional evidence that the matrix lists (not edited, except one line of
`test_p0_surfe_nans.py`): the native files in fact 7, `test_surfe_rbf_interpolator.py`,
`test_p0_surfe_nans.py`, `test_world_location_consumers.py::test_2d_viewer_fault_ellipse_centre_is_world`
(viewer, `matplotlib`) and `::test_internal_callers_do_not_warn` (`pyevtk`).

`packages/loop_interpolation/tests/test_p0_surfe_nans.py` (edited, one line, see task 4).
The line `pytest.importorskip("surfe", minversion=None)` uses the pip name. The adapter
imports `surfepy`. The task reads the package data of `surfe` on PyPI and decides. The
file changes to `pytest.importorskip("surfepy")` only if the installed module of that
package is `surfepy` and not `surfe`. If the task cannot decide, the file does not change.

`tests/unit/test_backend_matrix.py` (new; mandatory; no skip). It imports `MATRIX`,
`SPEC_ROWS`, `classify_cell`, `optional_status`, `require_optional`,
`outcomes_from_junit`, `rows_without_real_cell` and `render_table`. The helper tests are
written with the helper part of `backend_matrix.py`. The tests that need the `MATRIX` list
are written with the list, in a later wave.

| Test function | Cells | Checks |
|---|---|---|
| `test_every_spec_row_has_a_cell` | 7 | `rows_without_real_cell(MATRIX)` is empty: each of the 7 rows of `SPEC_ROWS` has at least one cell with `stub` false. Every cell names a row of `SPEC_ROWS` or the value `"extra"`. A stub cell alone does not satisfy a row: the test removes the two real cells of rows 6 and 7 from a copy of `MATRIX`, and `rows_without_real_cell` then returns exactly those two rows |
| `test_stub_cells_are_mandatory_and_listed_apart` | 2 | every cell with `stub` true is mandatory and has no dependency. `render_table` prints the rows 6 and 7 with the status of their real cells only, and prints the stub cells in a separate block with the heading "Boundary stub cells (Python side; not a native or Surfe result)". With the outcome "stub passed, extension skipped, dependency absent", the line of row 6 reads "excluded" and the stub line reads "passed" |
| `test_mandatory_cells_name_existing_tests` | one for each mandatory cell (11) | each ref is a file in the repository with a top-level function of that name (AST parse, no import) |
| `test_mandatory_cells_have_no_skip` | one for each mandatory cell (11) | in the referenced functions (decorators and body) and in the module-level statements and `pytestmark` of each referenced file, and in `tests/conftest.py`, every file `tests/fixtures/*.py` and the `conftest.py` of the folder of the ref, no use of `skip`, `skipif`, `importorskip`, `xfail` or `require_optional`. The one allowed exception is the definition of `require_optional` in `backend_matrix.py`. A temporary file with a skip makes the check fail (the test builds one in `tmp_path`) |
| `test_optional_status_names_exact_module` | 5 dependencies | with a fake `find_spec` that returns `None`, `available` is False and the reason contains the exact module name; with a fake that returns an object, `available` is True and the reason is empty |
| `test_require_optional_skips_with_named_dependency` | 2 (native, surfe) | with the fake `find_spec` that returns `None`, `require_optional` raises the pytest skip exception and its message contains `loop_interpolation_core` or `surfepy` |
| `test_classify_cell_rules` | 12 | the 12 combinations of (mandatory or optional with absent dependency or optional with installed dependency) x (passed, failed, skipped, absent) give the strings of the rule above; "passed without the dependency" and "mandatory skipped" are violations |
| `test_outcomes_from_junit` | 6 | a small XML text with passed, failed, skipped and missing cases gives the four outcomes; a cell with fewer passed cases than `min_cases` is "absent"; a class test (`classname="a.b.test_x.TestC"`) matches a ref of the file `a/b/test_x.py`; two files that hold the refs of one cell are merged (one passed case in each gives "passed"); a module-level skip (`classname=""`, `name="a.b.test_x"`) gives "skipped" for every ref of that file, and `classify_cell` of an optional cell with an absent dependency gives "excluded" for it, not "absent"; the same case for a mandatory cell is a violation |

#### S6 backend and capability cells

"Mandatory" cells have 0 skips (SC-018). The ids are the `cell_id` values of `MATRIX`.
The matrix task fills `refs`, the id selectors and `min_cases` from
`pytest --collect-only -q` and from the test names of S1 to S6. The selectors are
literals that occur in the collected ids.

| Cell id | Spec row | Level | Evidence (test files) | Dependency |
|---|---|---|---|---|
| `py-fdi` | Python FDI and PLI/P1 | Mandatory | R1 `test_frame_persistence.py`, `test_frame_copy_and_lifecycle.py`; R2 `test_fault_frame_analytic.py`, `test_fault_frame_validation.py`; R3 `test_reference_constraint_frames.py`; R5 `test_fault_buffer_bounds.py` (FDI ids) | none |
| `py-pli` | Python FDI and PLI/P1 | Mandatory | the same files (PLI ids) | none |
| `py-p2` | Python P2 | Mandatory | R1 exact field in the persistence and lifecycle files; `test_reference_constraint_frames.py` (P2 ids); `test_diagnostics_backends.py` (P2 cells) | none |
| `py-fold-fdi`, `py-fold-pli`, `py-fold-p2` | Python FDI-fold, PLI discrete-fold, P2-fold | Mandatory | `test_fold_constraint_frames.py` (family ids) | none |
| `py-cn-fdi`, `py-cn-pli` | Python constant-norm variants | Mandatory | `test_constant_norm_frames.py`, lifecycle files (`cn_fdi`, `cn_pli` ids) | none |
| `common-geometry` | Common 2D and 3D geometry | Mandatory | `test_rigid_transform_contract.py`, `test_resize_paths.py`, `test_cell_centre_sampling.py`, `test_euclidean_transformation.py`, `tests/unit/utils/test_transformation.py` | none |
| `native-boundary` | Existing native FDI or fold integration (Python side) | Mandatory, stub, listed apart (not a native pass; not a row result) | `test_backend_world_contract.py` (native tests) | none |
| `surfe-boundary` | Optional Surfe integration (Python side) | Mandatory, stub, listed apart (not a Surfe pass; not a row result) | `test_backend_world_contract.py` (Surfe tests) | none |
| `native-extension` | Existing native FDI or fold integration | Optional (real extension; gives the status of row 6) | `test_backend_optional_native_surfe.py`; the native files of fact 7; `test_native_fault_frame.py`; `test_native_errors.py` (13 tests run without the extension; the function at line 250 skips); `tests/unit/modelling/test_builder_solver.py` (one native test, line 28). R3: the existing test `test_orthogonal_constraints_give_the_same_field_as_cg` covers the R3 operation (gradient orthogonal rows) with no rotated frame. R3 in a rotated frame is "R3 not covered" for native (fact 7) | `loop_interpolation_core` (not installed in the venv, 2026-11-06) |
| `surfe-extension` | Optional Surfe integration | Optional (real extension; gives the status of row 7) | `test_backend_optional_native_surfe.py`; `test_surfe_rbf_interpolator.py`; `test_p0_surfe_nans.py` | `surfepy` (not installed) |
| `display-viewer` | Grid and display consumers (S4, not a spec row) | Optional | `test_world_location_consumers.py::test_2d_viewer_fault_ellipse_centre_is_world` | `loopstructuralvisualisation` (not installed) and `matplotlib` (installed) |
| `export-pyevtk` | Grid and display consumers (S4) | Optional | `test_world_location_consumers.py::test_internal_callers_do_not_warn` | `pyevtk` (not installed) |

`SPEC_ROWS` has the 7 rows of the spec table. The two display cells have the row value
`"extra"`, because the spec table has no display row. The planner kept them. S6 keeps
them. The S4 tests skip here. SC-018 needs a named dependency for each exclusion.

The matrix table prints the 7 spec rows. Rows 6 and 7 show the status of their real
extension cell only. The stub cells are printed in a separate block below the rows, with
the note that they test the Python side only. In this environment (2026-11-06) the
expected result is: rows 6 and 7 show "excluded" with the named dependency; the 2
display cells are excluded with a named dependency; 0 mandatory cells are skipped
(the 11 mandatory cells include the 2 stub cells). Unsupported capabilities keep their
documented rejection: the native adapter rejects a rotated support, a region and sparse
rows (`test_native_unsupported_settings_keep_rejection_and_state`). Surfe rejects nothing.
Its inequality hooks are no-ops by design and stay.

The exclusion ledger in `verification.md` also lists every other skipped test of the
full run with its reason. It uses the junit files of the three gate runs. It groups the
skips by dependency or reason.

- Native extension `loop_interpolation_core`: the native files of fact 7, the function at
  `test_native_errors.py:250`, and `tests/unit/modelling/test_builder_solver.py:28`.
- Surfe: `test_surfe_rbf_interpolator.py` (`HAS_SURFE`) and `test_p0_surfe_nans.py`.
- `geoh5py`, `omf` and `pyevtk` (module-level skips): `tests/unit/io/test_geoh5.py`,
  `tests/unit/io/test_omf.py` and `tests/unit/io/test_exporters.py`. No test names
  `trimesh` in a skip (`git grep`, 2026-11-06). The ledger task searches again and lists
  any file that it finds.
- `loopstructuralvisualisation`, `matplotlib` and `pyevtk` (function-level skips): the
  two cells of `test_world_location_consumers.py` and `test_api.py:121` (`matplotlib`).
- `pyvista` (function-level skips): `tests/unit/geometry/test__structured_grid.py`,
  `test__surface.py` and `tests/unit/interpolator/test_2d_discrete_support.py`.
- `dill` (installed here): `tests/integration/test_p2fold_model.py` and
  `tests/unit/modelling/test_geological_model.py`.
- Skips with no named dependency (unsupported capability, with the skip reason): for
  example the four `TetMesh ... evaluate_d2` skips in `test_p2_interpolator.py` and the
  skips in `packages/loop_interpolation/tests/golden/test_golden.py`.

The 47 skipped tests of `packages/loop_interpolation/tests` and the skips of the root
suite must all appear. None of them is a mandatory cell.

#### S6 risks

- The planner read more than 10 code files, by targeted line ranges, as S5 did. Files
  read: `math/_transformation.py`, `utils/_transformation.py` (as a `diff`),
  `_native.py` (lines 1 to 130, 234 to 368, 454 to 595), `_surfe_wrapper.py`,
  `_fault_function.py` (lines 1 to 20, 197 to 225, 305 to 458), `_fault_segment.py`
  (lines 1 to 140, 180 to 420), `_structural_frame.py` (lines 1 to 60),
  `_geological_feature.py` (lines 125 to 430), `_base_builder.py` (lines 81 to 103),
  `_fault_builder.py` (lines 405 to 560), `_geological_feature_builder.py` (lines 300 to
  365), `_discrete_interpolator.py` and the FDI and PLI interpolators (the constraint
  names), `_geological_interpolator.py` (lines 815 to 845), plus the fixture and test
  files named above. The risk is that an unread part has another reader. Each task greps
  for readers before it edits.
- Neither the native extension nor `surfepy` is installed. The two extension cells did
  not run. Their tolerances (SC-002) are not measured (OP21). The stub cells test only the
  Python side. The final report must say "excluded", never "passed", for these two cells.
  The matrix table lists the stub cells apart, so they cannot hide the exclusion.
- The pip name of Surfe is unverified (fact 21). The venv has no `pip` module. The check
  in task 4 needs a pip outside the venv and a network. If it cannot run, the file
  `test_p0_surfe_nans.py` does not change and the task reports "not decided". The tests
  of that file then skip even if only `surfepy` is installed. The ledger records this.
- The R2 analytic frame is a `FaultSegment` of `AnalyticReference` features. The planner
  read the constructor (`StructuralFrame.__init__` needs a name and a list of features)
  and `apply_to_points` (it calls `evaluate_value` and `evaluate_gradient(...,
  ignore_regions=True)` of the features, and `AnalyticReference` accepts both). The planner
  did not run it. The test task runs it first. If a method of `BaseFeature` that
  `apply_to_points` needs is missing on `AnalyticReference`, the task adds the smallest
  method to the fixture class, not to the library.
- `apply_to_points` of a solved fault with an empty mask may raise. The outside-support
  test has a fallback (see the table).
- Run time. `test_fault_frame_validation.py` does 26 R2 solves (24 + 2). S1 measured about
  5 minutes for 4 files with 48 R2 solves. `test_diagnostics_backends.py` builds 36 small
  models (24 faulted and 12 P2) with module-scoped fixtures. The task records the time (OP2
  stays closed: the default run includes them). The shared fixtures need that no test
  changes a model. Test 5 of the file checks the S5 snapshot rule.
- The OP19 change alters the numbers of the misfit methods for faulted features, and for
  features with a region mask (rows outside a region gave NaN; now they give numbers). The
  S5 version of the methods is not in a release (S5 is pull request #402). The pre-005
  methods gave wrong numbers for faulted features too. `COMPAT.md` records it.
- `EuclideanTransformation` changes for callers outside the repository. An integer array
  gives a float result. A 1-D array raises `ValueError`. `inverse_transform` works. The 3D
  numbers change: z is no longer reversed, and the fitted angle has the right sign for
  lines with a direction between -45 and 0 degrees. `fit` raises `ImportError` when
  scikit-learn is missing and rotation is requested. No caller in the repository is
  affected (fact 20). The old import path stays. A pickle made before S6 loads, because the
  old module exposes the class. A pickle made after S6 stores the path
  `loop_common.math._transformation`, and an older release cannot load it.
- The module-level cache of the Surfe stub touches `sys.modules`. The fixture restores the
  state and asserts it. A leak would make later Surfe tests use the stub.
- The matrix checker reads junit names from three files (root, `loop_common` tests,
  `loop_interpolation` tests). A change of a test id breaks the match. The test
  `test_mandatory_cells_name_existing_tests` finds a renamed function, and `min_cases`
  finds a lost parametrisation. If a gate run starts in another folder, the `classname`
  paths change and every cell becomes "absent". The final-gates task checks one
  `classname` of each file before it runs the matrix.
- A display test passes only when the viewer package is installed. In this environment it
  is excluded with its named dependency.
- The SC-004 bound (0.03) is loose. The negative control shows that it can fail. The
  frame error of S1 (1e-6) is much smaller, so a real regression of the frame shows first
  in the S1 and residual tests.
- The sign of the main axis from scikit-learn depends on the library version. The fit
  tests compare the angle modulo `pi` and compare the alignment of the line, so they do
  not depend on that sign. The hand numbers in the plan were measured with scikit-learn
  1.9.0.
- Stage size. The limit of 10 code files holds for each wave, not for the stage. The whole
  stage changes 17 code files: 5 source files and 12 test and fixture files (see "S6
  files"). The widest wave has 9. The file `backend_matrix.py` and the file
  `test_backend_matrix.py` appear in two tasks each (tasks 6 and 15). Those tasks are in
  different waves.

#### S6 constitution check

| Project rule | Result | How S6 agrees |
|---|---|---|
| Work stays within the existing packages | Pass | Three source files in `loop_common` and `LoopStructural.utils`, one in `LoopStructural/modelling`, one byte fix in `loop_common`. No new package or module in `src`. Test helpers are new files under `tests/fixtures`. |
| Gates `pytest` and `ruff check LoopStructural packages` | Pass | Decision 11 and the S6 row of "Stages and rollback". |
| numpy docstrings, type hints, py39 | Pass | New functions have numpy docstrings and hints. `_transformation.py` keeps `from __future__ import annotations`. `backend_matrix.py` uses `from __future__ import annotations` and `typing.Literal`. |
| The public `GeologicalModel` API stays compatible; a break needs a ROADMAP.md decision and a major-version note | Pass | No `GeologicalModel` signature changes. `EuclideanTransformation` keeps both import paths and its constructor. `inverse_transform` changes from an error to a result for ordinary arrays. A 1-D input changes from `IndexError` to `ValueError`. An integer input gives a float result. `fit` raises `ImportError` in place of a log message. The 3D numbers change (decision 3). `COMPAT.md` records these, and `ROADMAP.md` records the decision for the behaviour change. The physical residuals of faulted features change by OP19; `ROADMAP.md` records that decision too. |
| Stages can be tested and reversed independently | Pass | S6 reverts to the S5 merge. Test-only files and `backend_matrix.py` can stay or go. |
| `ROADMAP.md` owns decisions; update it when one changes | Pass | Two decision entries, a sentence on item 2c-5, and the line "Spec 005 complete". |
| Feature documents stay under `specs/NNN-name/` | Pass | `plan.md`, `tasks.md` and `verification.md` change there; `spec.md` gets C-019; `specs/README.md` carries the status. |

#### S6 files

The limit of 10 code files (.py) applies to each wave of tasks. In total S6 changes 17
code files: 5 source files and 12 test and fixture files. Each task changes at most 3
code files. Data files and documents are not counted.

| File | Change |
|---|---|
| `packages/loop_common/src/loop_common/math/_transformation.py` | Logger import; `_check_points`; `_import_decomposition`; new `fit`, `_rotation`, `transform`, `inverse_transform` |
| `packages/loop_common/src/loop_common/math/__init__.py` | Export `EuclideanTransformation` |
| `LoopStructural/utils/_transformation.py` | Re-export of the `loop_common` class |
| `LoopStructural/modelling/features/_geological_feature.py` | Two lines in `_gradient_cosines_fresh` and `_value_residuals_fresh`; two docstring sentences |
| `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py` | Remove the BOM (3 bytes) |
| `tests/unit/modelling/test_fault_frame_analytic.py` | Remove the BOM (3 bytes); no other change (decision 14) |
| `tests/fixtures/coordinate_frames.py` | New helpers (see "S6 functions") |
| `tests/fixtures/backend_matrix.py` | New (helper part in task 6, `MATRIX` list in task 15) |
| `tests/unit/utils/test_transformation.py` | New cells; delete the defect test |
| `packages/loop_common/tests/test_euclidean_transformation.py` | New |
| `tests/unit/modelling/test_fault_frame_validation.py` | New |
| `tests/unit/modelling/test_diagnostics_backends.py` | New |
| `tests/unit/modelling/test_diagnostics_frames.py` | Findings 5 to 8 and the patch target of the exact-field test |
| `tests/unit/interpolator/test_backend_world_contract.py` | New |
| `tests/unit/interpolator/test_backend_optional_native_surfe.py` | New |
| `tests/unit/test_backend_matrix.py` | New (helper tests in task 6, matrix tests in task 15) |
| `packages/loop_interpolation/tests/test_p0_surfe_nans.py` | Edited, one line (the `importorskip` module name), only if task 4 decides so |
| `COMPAT.md` | S6 notices; move one sentence of the S5 row |
| `ROADMAP.md` | Two decision entries; one sentence on item 2c-5; the line "Spec 005 complete" |
| `specs/005-coordinate-frame-correctness/verification.md` | S6 sections: pre-fix, post-fix, golden review (with the behaviour-change ledger), findings disposition, backend matrix and exclusion ledger, coverage ledger, rollback, run time, final gates, open points |
| `specs/005-coordinate-frame-correctness/spec.md` | Clarification C-019 and the R2 fixture text (done in the planning run) |
| `specs/README.md` | Status of the 005 line |

Not changed (the survey in "Later stages" listed them): `_native.py`, `_surfe_wrapper.py`
and `_fault_function.py`. The reading found no coordinate defect in them; their
contracts get tests only. `_fault_segment.py` is not changed. The existing S1 to S5 test
files other than `test_diagnostics_frames.py`, `tests/unit/utils/test_transformation.py`
and the BOM file `test_fault_frame_analytic.py` are not edited.

The `COMPAT.md` notices say:

- `EuclideanTransformation` is one class. `loop_common.math.EuclideanTransformation` is
  new. `LoopStructural.utils.EuclideanTransformation` is the same object.
  `inverse_transform` returns the original points for an `(N, d)` array; before it raised
  `ValueError` or returned too few rows. The methods return a float copy; an integer input
  was truncated before. A 1-D input raises `ValueError`; before it raised `IndexError`.
- The 3D behaviour changes. For `dimensions=3` the z axis is no longer reversed
  (`z' = z - tz`; before it was `-(z - tz)`). `fit` works for `dimensions=3`; before it
  failed. The fitted angle has the right sign for a main axis with a negative slope; before
  it was wrong for lines with a direction between -45 and 0 degrees. A 3D user must check
  every result that was computed with the old class: the numbers differ. `fit` with
  `fit_rotation=True` raises `ImportError` when scikit-learn is missing; before it logged
  an error and returned `None`. `fit` with `fit_rotation=False` needs no scikit-learn.
- `fit` writes its debug message through the `loop_common` logger
  (`loop_common.logging.get_logger`), not through the `LoopStructural` logger.
- A pickle of the class made before S6 still loads, because the old module path stays
  importable and exposes the class. A pickle made after S6 stores the path
  `loop_common.math._transformation`. A release from before S6 cannot load it.
- `GeologicalFeature.evaluate_gradient_misfit`, `evaluate_value_misfit`,
  `evaluate_gradient_cosines`, `evaluate_value_residuals` and `verification_report` of a
  faulted feature compare the stored restored rows with the interpolator at the same
  positions. Before, the faults acted a second time. The two private methods
  `_value_residuals_fresh` and `_gradient_cosines_fresh` now evaluate the interpolator
  directly. This reverses S5 decision 8 for those two methods. The numbers for a feature
  with no faults and no regions do not change. For a feature with a region mask, a row
  outside the region gave NaN before and gives a number now.
- Known issue, not fixed in S6: the public `evaluate_gradient` of a faulted feature keeps
  its tetrahedron path and the shape error that S3 noted.
- The S5 sentence "With an identity frame and no resize, the value misfit is the same as
  before" moves from the `ConstraintDiagnosticsReport` row to the `evaluate_value_misfit`
  row.

The `ROADMAP.md` entries are dated with the day of the S6 end gate.

- Decision "Residuals of faulted features (spec 005, S6)". It says that the residual
  diagnostics of a faulted feature evaluate the interpolator at the restored positions,
  and that S5 decision 8 changes for these two private methods. It says that the notice is
  in `COMPAT.md`.
- Decision "EuclideanTransformation (spec 005, S6)". It is a behaviour change. It names the
  four defects (the 3D `fit`, the sign of the fitted angle, the reversed z axis, the silent
  need for scikit-learn) and their physical effect (a mirrored z axis and a wrong rotation
  of the map view for some lines). It says that one class now lives in `loop_common` and
  that `LoopStructural.utils` re-exports it. It has one sentence: "This reverses item 2c-5
  of the stage sequence, which was closed as deferred." It says that the notice is in
  `COMPAT.md`.

The line "Spec 005 complete" is added to the status log: "Spec 005 complete: stages S1 to S6
implemented; the auxiliary transformation, the R2 displacement validation, the backend
matrix and the golden review are done. Gate result: root pytest <passed> passed, <skipped>
skipped (every skip names its dependency in
`specs/005-coordinate-frame-correctness/verification.md`)." The task replaces the two
placeholders with the numbers of the final gate.

#### S6 task order (for the spec-tasks skill)

No two tasks in one wave write the same file, and no task shares a wave with a task that it
depends on. `verification.md`, `COMPAT.md` and `ROADMAP.md` are written by tasks in
different waves. No two tasks in one wave write `verification.md`. Each task changes at
most 3 code files. Each wave changes at most 10 code files. The widest wave (wave 1) has
9. The file `backend_matrix.py` appears in task 6 (the helper) and task 15 (the `MATRIX`
list). The file `test_backend_matrix.py` appears in the same two tasks. They are in
different waves.

| Wave | Task | Files | Depends on |
|---|---|---|---|
| 1 | 1. Fixture helpers | `tests/fixtures/coordinate_frames.py` | none |
| 1 | 2. Transformation tests (they fail first) | `tests/unit/utils/test_transformation.py`, `packages/loop_common/tests/test_euclidean_transformation.py` | none |
| 1 | 3. Changes of `test_diagnostics_frames.py` (findings 5 to 8; measures `TOL_SOLVED_FDI` and `TOL_SOLVED_PLI`) | `tests/unit/modelling/test_diagnostics_frames.py` | none |
| 1 | 4. Name check of the Surfe package (see below) | `packages/loop_interpolation/tests/test_p0_surfe_nans.py` (only if the check decides so) | none |
| 1 | 5. BOM removal | `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`, `tests/unit/modelling/test_fault_frame_analytic.py` | none |
| 1 | 6. Matrix helper and its tests | `tests/fixtures/backend_matrix.py` (helper part, `MATRIX = ()`), `tests/unit/test_backend_matrix.py` (helper tests) | none |
| 2 | 7. R2 validation tests (measure `TAU_ASSEMBLED` and `TAU_ORTHOGONAL`) | `tests/unit/modelling/test_fault_frame_validation.py` | 1 |
| 2 | 8. Diagnostics tests (measure `TOL_FAULTED_VALUE` and `TOL_FAULTED_COS`) | `tests/unit/modelling/test_diagnostics_backends.py` | 1 |
| 2 | 9. Boundary tests with the stubs | `tests/unit/interpolator/test_backend_world_contract.py` | 1 |
| 2 | 10. Optional extension tests | `tests/unit/interpolator/test_backend_optional_native_surfe.py` | 1, 6 |
| 3 | 11. Pre-fix evidence in `verification.md` | `verification.md` | 2, 3, 7, 8, 9, 10 |
| 4 | 12. Transformation fix | `packages/loop_common/src/loop_common/math/_transformation.py`, `packages/loop_common/src/loop_common/math/__init__.py` | 11 |
| 4 | 13. Residual fix | `LoopStructural/modelling/features/_geological_feature.py` | 3, 8, 11 |
| 5 | 14. Re-export | `LoopStructural/utils/_transformation.py` | 12 |
| 5 | 15. `MATRIX` list and the matrix tests | `tests/fixtures/backend_matrix.py` (the list), `tests/unit/test_backend_matrix.py` (matrix tests) | 2, 3, 6, 7, 8, 9, 10, 12, 13 |
| 6 | 16. Post-fix evidence in `verification.md` | `verification.md` | 12, 13, 14, 15 |
| 7 | 17. Golden review, tolerance audit, coverage ledger, findings disposition | `verification.md` | 16 |
| 8 | 18. `COMPAT.md` | `COMPAT.md` | 17 |
| 9 | 19. Final gates, matrix run, rollback check, exclusion ledger, run time | `verification.md` | 18 |
| 10 | 20. `ROADMAP.md` decisions and the line "Spec 005 complete" | `ROADMAP.md` | 19 |

The orchestrator then updates `specs/README.md`.

Notes for the tasks.

- Task 1 writes `r2_curved_slip_fault` and `euler_circle_closed_form` among the other
  helpers. It runs `apply_to_points` on the analytic frame first and reports the result
  (see the risks).
- The tasks 7 and 8 depend on task 1 because they import its helpers. Both tasks also
  use module-scoped parametrised fixtures to build each model once.
- Tolerances: `TOL_FAULTED_VALUE`, `TOL_FAULTED_COS`, `TAU_ASSEMBLED`, `TAU_ORTHOGONAL`
  and the per-engine `TOL_SOLVED_*` are set inside the test tasks (3, 7 and 8), from
  independent computations, before the fix. If a bound would pass its cap, the task stops.
- Task 4. Step 1: run `pip index versions surfe` with a pip that is not in the venv (the
  venv has none). Step 2: if possible, run `pip download surfe --no-deps -d <scratch>`
  where `<scratch>` is a new folder under the scratchpad directory. Step 3: list the
  files of the wheel (`python -I -m zipfile -l <wheel>`) and read the top-level module
  names. Step 4, the decision: if the wheel holds the module `surfepy` and no module
  `surfe`, change the line to `pytest.importorskip("surfepy")`. If the wheel holds
  `surfe`, or both, or the check cannot run, do not change the file, and report the
  package name and the module names, or "not decided". Step 5: report the result for
  `verification.md` (task 11) and for the cell `surfe-extension` (task 15). The task does
  not install anything into the venv and does not unpack the wheel into the repository.
- Task 5 checks the first bytes of both files, removes `ef bb bf`, and checks that
  `git diff -w --stat` shows no other change.
- Task 11 expected failures before the fixes: the transformation tests (the measured count
  goes in the record; by reading, the old `inverse_transform` fails for `dimensions=2`
  with 3 or 5 columns, for `dimensions=3` with 5 columns, and for `dimensions=3` with 3
  columns and N greater than 3; it works for `dimensions=3`, 3 columns and N=1, but the
  z check of the forward step fails every `dimensions=3` cell), the defect tests of the
  fit (negative slope, 3D fit, the z direction, the missing scikit-learn), the three
  residual tests of faulted features (tests 1, 2 and 4 of `test_diagnostics_backends.py`),
  and `test_entry_points_are_one_class`. Expected guards that pass: the changed exact-field
  test of `test_diagnostics_frames.py` (it passes before and after the fix for a feature
  with no faults), the displacement and residual tests of the R2 frame, the P2 tests, the
  test of the restored rows, and the stub cells.
- Task 17 commands. The golden file list: `git diff --name-status ebc7c699d5 HEAD -- . ':!specs' ':!*.py' ':!*.md'`
  (fact 15). The tolerance audit has two parts. First, the removed lines:
  `git diff -U0 ebc7c699d5 HEAD -- ':!specs' | grep -iE '^-.*(tol|approx|abs=|<= *[0-9.]+e-)'`.
  The pattern ignores case, so it finds `TOL_SOLVED`, `TOL = ...`, `Tp` helpers and `tau_`
  names (the task adds `tp` and `tau` to the pattern). Second, the added lines with the
  same pattern, with `^+` in place of `^-`. The task pairs each removed line with its
  added line and compares the numbers. Every removed line with a tolerance is listed and
  explained. The expected count of relaxed tolerances is 0. The one unavoidable change is
  the rename of `TOL_SOLVED` into `TOL_SOLVED_FDI` and `TOL_SOLVED_PLI`. It is a tightening
  or no change, because each constant is at most 9e-4. The task names it for the SC-015
  review. The production-file list is
  `git diff --name-only ebc7c699d5 HEAD -- LoopStructural packages` without test paths,
  each file mapped to a test.
- Task 19 runs the commands of decision 11 one process at a time, then checks one
  `classname` of each junit file (it must start with `tests.`,
  `packages.loop_common.tests.` or `packages.loop_interpolation.tests.`), then runs the
  matrix, the rollback check and the exclusion ledger. It writes the results to
  `verification.md`.
- Task 20 reads the gate numbers from task 19.

S6 rollback check. Revert to the S5 merge commit, `<S5 merge commit>` (the task 19 fills
in the hash). Then run these checks. The result of each goes to `verification.md`.

- (a) `from LoopStructural.utils import EuclideanTransformation` and
  `import LoopStructural.utils._transformation` work on both trees, and the class
  `__init__` signature is equal. A pickle made on the S5 tree loads on the S6 tree.
- (b) Every public name (no leading underscore) of `loop_common.math` on the S5 tree is
  also in the S6 tree (the set of public names of `dir(loop_common.math)` of S5 is a
  subset). The S6 tree adds `EuclideanTransformation` only. The private submodule name
  `_transformation` is not compared.
- (c) For a float `(N, 3)` array, `transform` gives equal output on both trees
  (`np.array_equal`) for `dimensions=2` with a set angle and translation. For
  `dimensions=3` the z column differs by design.
- (d) `from loop_interpolation import *` gives the same names on both trees.
- (e) The signatures of the five `GeologicalFeature` diagnostic methods are equal. For a
  feature with no faults and no regions the residual numbers of the S5 test cells are equal.
- (f) All baseline files listed in fact 15 are byte-equal after S6
  (`git diff --stat <S5 merge commit> HEAD -- <files>` prints nothing).

### Stages and rollback

| Stage | Files and behaviour | Acceptance check | Rollback target and compatibility check |
|---|---|---|---|
| S1 | Shared transform boundary, reference relationships for FDI and PLI, fault frame, buffer and forced geometry, rotation check, resize stability (this design) | New S1 tests pass for FDI and PLI; `pytest`; `ruff check LoopStructural packages`; `pytest packages/loop_common/tests packages/loop_interpolation/tests` | Revert the S1 commits to master `ebc7c699`. Check: the existing `tests/unit/modelling/test__fault_builder.py` and `packages/loop_common/tests/test_bounding_box.py` pass on the reverted tree. No saved format changed. |
| S2 | Fold geometry, directional regularisation, quadrature, masks, sign alignment for FDI-fold, PLI discrete fold, P2 and P2-fold; constant-norm variants (FR-008; R3 for P2; R4) ("Stage S2 design") | New S2 tests (SC-002 for constant norm; SC-007, SC-008 and SC-013 for every mandatory fold and P2 cell); the existing fold package tests pass with no change; `pytest`; `ruff check LoopStructural packages`; package tests | Revert to the S1 merge. Check: the S1 tests still pass. No saved format changed. |
| S3 | Persistence format version, transform and coordinate-space metadata ("Stage S3 design"); no second projection at load; copies keep the frame; legacy load requires a transform or identity declaration; transform change policy and restored-point cache invalidation (FR-016 to FR-020); ROADMAP.md decision for FR-019; `COMPAT.md` migration note | SC-006, SC-012, SC-016, SC-017 | Revert to the S2 merge. Check: old-format R1 files load through their documented path; a new-format file fails with 1 actionable error and 0 state changes on the old reader (SC-017). |
| S4 | Display, export of observation geometry, fault-centre geometry and intrusion classification use stored world locations (FR-022; "Stage S4 design") | World-location comparison within `Tp` (U3.4, U3.5) | Revert to the S3 merge. Check: viewer import paths unchanged. |
| S5 | Remaining resize paths (3D `nsteps` and `step_vector`, `set_nelements`, 2D derived counts, rectilinear guard), golden review, cell-centre tests, non-mutating misfit diagnostics, residual and exclusion report (FR-014 remainder, FR-015, FR-023, FR-024; "Stage S5 design") | SC-009, SC-011; the new S5 tests; `pytest`; `ruff check LoopStructural packages`; package tests | Revert to the S4 merge. Check: diagnostics report fields unchanged for callers (`dataclasses.fields`, `to_dict` keys, exported names; "S5 rollback check"). |
| S6 | Auxiliary point transformation (one class, round trips, and the 3D, sign, z and scikit-learn fixes) and its entry points (FR-025); R2 validation with constraint residuals and displacement, and the physical-residual fix for faulted features (FR-026, FR-027); backend and capability matrix with native and Surfe (FR-030); golden review (FR-029); final gates and the "Spec 005 complete" line in ROADMAP.md ("Stage S6 design") | SC-004, SC-010, SC-014, SC-015, SC-018; the new S6 tests; `pytest`; `ruff check LoopStructural packages`; package tests; the matrix run over the merged junit files of the three gate runs with 0 mandatory skips | Revert to the S5 merge. Check: `LoopStructural.utils` and `loop_common.math` import paths unchanged ("S6 rollback check"). |

## Test strategy

Shared rule for every stage (FR-031): write each regression test first, run
it on the unfixed tree, and record the failure (test id, observed value) in
`specs/005-coordinate-frame-correctness/verification.md`. Then fix and record
the pass. A finding of the 1.7.0 investigation that does not reproduce on the
current tree is recorded as "not reproduced"; its test stays as a guard.

- **Unit and regression tests (S1)**:
  - `packages/loop_common/tests/test_rigid_transform_contract.py`:
    `validate_rigid_rotation` accepts rotations of 0, 30 and 90 degrees about
    z; rejects scale, shear, singular, reflection (det -1), NaN/inf and wrong
    shape (2D and 3D). SC-001 point round trips (`Tp`) and direction round
    trips (`1e-12`) for 1 point and arrays, for all 4 anchors. `local_enclosure`
    contains all 8 projected corners for each rotation. Invalid rotation:
    `set_local_transform` and `rotation_xy` leave the prior state unchanged
    (SC-012, rotation part). Reapply unchanged support bounds 10 times: 0
    changes to origin, maximum, step and counts (SC-011, resize part).
  - `tests/unit/modelling/test_reference_constraint_frames.py`: R3 for FDI
    and PLI, 4 anchors x 3 rotations. Equality region `x' > 4` selects 168 of
    336 nodes unrotated; for rotations the expected set comes from the world
    nodes. 0 missing or duplicate nodes, local index paired with the world
    value (SC-007). Orthogonality rows with `f=y'`, direction `(1,0,0)`,
    target 0, and direction `(0,1,0)`, target 1: residual at most `1e-10`
    (SC-008, FDI and PLI; barycentre, `"centroid"` and `"gauss4"`). Feature
    inequality bounds equal `h` at each node's world location. Restricted
    reference support: `excluded_rows` count equals the analytic count of
    samples outside the closed restricted box. Support expansion through a
    restoring fault: restored world samples project inside the support
    (FR-009, U2.7).
  - `tests/unit/modelling/test_fault_frame_analytic.py`: R2 through
    `GeologicalModel.create_and_add_fault` for FDI and PLI, explicit and
    inferred centre, 4 anchors x 3 rotations. SC-003 on the 6859 probes with
    CG `rtol=1e-10`, `atol=1e-12`, `maxiter=20000`. U1.1: generated rows use
    the centre with no projection. Inferred centre equals `O` in all 3
    components. R4: the extent-direction row equals the analytic extent
    gradient at the centre within `1e-10` per component; the analytic extent
    directions at `O+(0,-6,0)` and `O+(0,6,0)` differ by at least 15 degrees;
    a constant centre direction at those probes fails the check (SC-013).
    Non-orthogonal slip input: `fault_slip_vector` is orthogonal to the
    normal; a parallel slip raises `ValueError`.
  - `tests/unit/modelling/test_fault_buffer_bounds.py`: R5 for FDI and PLI,
    buffers 0, 0.2, 0.5, 3 rotations, explicit and inferred centre, forced
    geometry. Each local face is outward of the expected enclosure by at most
    1 cell plus `Tp`. The 3 buffers give distinct bounds. All generated rows
    are inside the support (SC-005). A property change after creation
    (`fault_builder.fault_normal_vector = ...`) keeps the buffer and the
    forced geometry (U1.5).
  - `tests/fixtures/coordinate_frames.py` (new helpers, no tests): anchors,
    `rot_z(deg) -> np.ndarray`, `tp(L, W) -> float`, R1 to R5 data frames,
    the R2 analytic frame, the R3 field `h`, and the R4 curved fields. Tests
    import only this module, not an external investigation workspace.
- **Integration and compatibility tests**: `tests/unit/test_public_api_contract.py`
  and `tests/unit/test_stable_api_surface.py` must pass unchanged. If
  `tests/fixtures/api_surface_snapshot.json` lists public methods of
  `GeologicalInterpolator` or `BoundingBox`, add the new names to it (additive).
- **Round-trip tests**: S1 covers point and direction round trips (SC-001).
  Persistence round trips (SC-006) are S3: `test_frame_persistence.py`.
- **Backend coverage**: S1 covers Python FDI and PLI/P1 cells for R2, R3 and
  R5. S2 covers P2, FDI-fold, PLI discrete fold, P2-fold and the FDI and PLI
  constant-norm variants ("Stage S2 design"). S6 covers the matrix with native and Surfe
  ("Stage S6 design").
- **Cross-language fixtures**: Not applicable in S1.
- **Grid and display checks**: S1 checks support enclosure and resize. S4 checks
  the display and export of data, the intrusion classification and the fault
  centre in the 2D viewer (`test_world_location_consumers.py`). Cell centres are
  S5.
- **Performance check**: none required. Record the S1 test time (OP2).
- **Merge gates**: `pytest` and `ruff check LoopStructural packages`. Also run
  `pytest packages/loop_common/tests packages/loop_interpolation/tests`,
  because the root `testpaths` is `tests` only.

## Project structure

### Common geometry (`packages/loop_common`)

**Responsibility**: rigid-rotation meaning, transformed box enclosure, stable
support resize.

**Files**:

- `packages/loop_common/src/loop_common/geometry/_rigid.py` - Create;
  `RigidRotationError`, `validate_rigid_rotation`.
- `packages/loop_common/src/loop_common/geometry/__init__.py` - Change; export
  `RigidRotationError` and `validate_rigid_rotation`.
- `packages/loop_common/src/loop_common/geometry/_bounding_box.py` - Change;
  validation in `set_local_transform`; new `local_enclosure`.
- `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py` -
  Change; `rotation_xy` validation; `_cells_for_length`; `origin` and
  `maximum` setters.
- `packages/loop_common/tests/test_rigid_transform_contract.py` - Create;
  SC-001, SC-011 (resize part), SC-012 (rotation part), enclosure.

### Interpolation (`packages/loop_interpolation`)

**Responsibility**: public world/local conversion; documented local and world
callback contracts.

**Files**:

- `packages/loop_interpolation/src/loop_interpolation/_geological_interpolator.py`
  - Change; four public conversion methods.
- `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`
  - Change; `add_inequality_feature` samples world nodes; contract docstring
  of `add_gradient_orthogonal_constraints`.

### Geological modelling (`LoopStructural`)

**Responsibility**: world fault geometry, reference sampling, support
requests.

**Files**:

- `LoopStructural/modelling/features/builders/_geological_feature_builder.py`
  - Change; `excluded_rows`; `install_gradient_constraint`;
  `install_equality_constraints`; `set_interpolation_geometry`;
  `check_interpolation_geometry`.
- `LoopStructural/modelling/features/builders/_fault_builder.py` - Change;
  stored buffer and forced geometry; 3D inferred centre; slip
  orthogonalisation; extent-direction row; early-return mesh geometry.
- `LoopStructural/modelling/core/_model_feature_factory.py` - Change; remove
  the centre projection and the buffer override in `create_and_add_fault`.
- `tests/fixtures/coordinate_frames.py` - Create; shared R1 to R5 fixtures.
- `tests/unit/modelling/test_reference_constraint_frames.py` - Create; R3,
  FR-004 to FR-007, FR-009.
- `tests/unit/modelling/test_fault_frame_analytic.py` - Create; R2, R4,
  FR-010, FR-011, SC-003, SC-013.
- `tests/unit/modelling/test_fault_buffer_bounds.py` - Create; R5, FR-012,
  FR-013, SC-005.
- `tests/unit/modelling/test__fault_builder.py` - Change only if an assertion
  encodes a defect fixed by S1 (for example an inferred centre with `z=0`).
  Record each change in `verification.md` under C-004.

### Documents

**Responsibility**: evidence, golden review, migration notes.

**Files**:

- `specs/005-coordinate-frame-correctness/verification.md` - Create;
  sections "Pre-fix failures", "Post-fix results", "Golden review (C-004)",
  "Backend matrix", "Rollback checks", "Test run time (OP2)".
- `COMPAT.md` - Change; one entry under "Migration notices": explicit
  `fault_center` is a world point; `fault_buffer` and `force_mesh_geometry`
  now keep their effect; an inferred centre keeps its elevation; non-rigid
  rotations are rejected; extra observations added as a workaround are no
  longer necessary.

### Later stages (file survey, closed)

The survey of 2026-11-04 listed the files of each stage. The planner runs read the
code and confirmed or changed each list. All stages are planned, so no later stage
remains.

- **S2**: planned. See "Stage S2 design" and "S2 files".
- **S3**: planned. See "Stage S3 design" and "S3 files". The survey listed
  `_support_factory.py` and `_fault_segment.py`; the code shows that S3 does
  not change them (no restoration cache exists).
- **S4**: planned. See "Stage S4 design" and "S4 files". The survey listed
  `_3d_viewer.py`, `intrusion_builder.py`, `_fault_segment.py` and
  `exporters.py`; the code shows that S4 does not change them. S4 changes
  `_geological_feature.py` and `intrusion_frame_builder.py`.
- **S5**: planned. See "Stage S5 design" and "S5 files". The survey listed
  `_bounding_box.py`, `_structured_grid.py`, `supports/_3d_structured_grid.py`
  and `_3d_rectilinear_grid.py`; the code shows that S5 changes only the
  rectilinear grid among them (plus `_structured_grid_3d.py` and
  `_structured_grid_2d.py`, which the survey did not list). `_diagnostics.py` and
  `_geological_feature.py` are confirmed.
- **S6**: planned. See "Stage S6 design" and "S6 files". The survey listed
  `_native.py`, `_surfe_wrapper.py` and `_fault_function.py`; the code shows that S6
  does not change them (it adds tests for their contracts). S6 changes both
  `_transformation.py` files, `loop_common/math/__init__.py`, `_geological_feature.py`
  (the S5 carry-over), the first bytes of `_structured_grid_3d.py` and of the S1 test file
  `test_fault_frame_analytic.py`, one line of `test_p0_surfe_nans.py` and `ROADMAP.md`.
  There is no later stage.

## Requirement coverage

| Requirement | Design and files | Verification | Gap |
|---|---|---|---|
| FR-001 | Existing public evaluation; S1 keeps it; fault path in `_model_feature_factory.py` | `test_fault_frame_analytic.py` (gradients in world after rotation) | None (S2 to S6 cover the other paths) |
| FR-002 | Factory keeps explicit centre in world; builder rows in world | `test_fault_frame_analytic.py` U1.1 | None |
| FR-003 | Public conversion methods; builder uses each once | `test_rigid_transform_contract.py`; `test_reference_constraint_frames.py` | None (the load path U3.2 is S3) |
| FR-004 | `install_gradient_constraint`, `install_equality_constraints`, `add_inequality_feature`; documented callbacks | `test_reference_constraint_frames.py` | None |
| FR-005 | `world_to_local_vectors` in `install_gradient_constraint`; S2 fold rows and directional terms | SC-007, SC-008 tests; S2 P2 cells | None |
| FR-006 | `install_equality_constraints` | 168-node check, rotated sets | None |
| FR-007 | `add_inequality_feature` | Bounds equal `h` at world nodes | None |
| FR-008 | "Stage S2 design": `sample_fold_geometry`, fold interpolators, `add_directional_regularisation`, native adapter, `world_norm_constraints` | `test_fold_constraint_frames.py` | None |
| FR-009 | `check_interpolation_geometry` | U2.7 test | None |
| FR-010 | Factory (no `model.scale`); 3D inferred centre | SC-003 explicit and inferred | None |
| FR-011 | Slip orthogonalisation; extent-direction row | SC-003, SC-013 | None |
| FR-012 | Stored `_fault_buffer`, `_force_mesh_geometry`; factory override removed | SC-005; U1.5 test | None |
| FR-013 | `local_enclosure`; `set_interpolation_geometry` | SC-005; enclosure unit test | None |
| FR-014 | `_cells_for_length`; S1 setters; S5: 3D `nsteps` and `step_vector` setters, `set_nelements` clamp, 2D derived counts (legacy `maximum` kept), rectilinear guard | `test_resize_paths.py` (resize 10 times; extent kept; A, B, A drift; `set_nelements`); golden review (C-004) | Extent after `set_nelements` changes by design (Q4); OP17 |
| FR-015 | S5: tests for rotated boxes and cell association; `cell_centres` unchanged unless a test fails | `test_cell_centre_sampling.py` (two tests are evidence; three are coverage) | OP12 (defect not seen by reading) |
| FR-016 | `frame_payload`; `to_dict` frame keys; format version 2 ("S3 saved format") | `test_frame_persistence.py::test_payload_declares_frame` | None |
| FR-017 | `_restore_local_constraints` stores local arrays with no conversion | `test_local_arrays_unchanged_after_two_cycles` | None |
| FR-018 | `_inherit_frame` in the three `copy()` methods | `test_frame_copy_and_lifecycle.py::test_copy_keeps_frame` | PLI has no `copy` (documented, not added); OP6 |
| FR-019 | `payload_kind`, `resolve_frame`, `InterpolatorFormatError`; ROADMAP.md decision; `COMPAT.md` row | `test_legacy_*`, `test_conflicting_inputs`, `test_malformed_metadata`, `test_unsupported_future_version` | None |
| FR-020 | `frame_signature`, `_check_frame`, `StaleFrameError` (reject at use) | `test_stale_frame_rejected_without_state_change`, `test_fault_restoration_uses_current_frame` | Policy choice OP4; no restoration cache exists |
| FR-021 | `validate_rigid_rotation` in `set_local_transform` and `rotation_xy`; saved metadata checked in `resolve_frame` | `test_rigid_transform_contract.py`; `test_malformed_metadata` | None |
| FR-022 | Remove the extra `rescale` in `GeologicalFeature.get_data` and in the intrusion classification (`_network_model_values`); the fault centre is a world point since S1 and its consumers are unchanged | `test_world_location_consumers.py` (display, export, classification, 2D ellipse); S1 `test_fault_ellipsoid_centre_is_world` | OP8 to OP10 |
| FR-023 | S5: both misfit methods keep their names and shape, work on copies, compare in world coordinates, never rebuild (warn and NaN when stale); new `evaluate_gradient_cosines` and `evaluate_value_residuals`; `direction_cosines` | `test_diagnostics_frames.py`; `test_residual_exclusion_report.py` | OP16; the gradient value changes meaning (`COMPAT.md`, `ROADMAP.md`) |
| FR-024 | `excluded_rows` (S1); S5: `ResidualSummary`, `ExclusionRecord`, `build_exclusion_records`, report fields, `verification_report` | `test_exclusions_are_counted_and_explained`, `test_physical_and_assembled_residuals_are_separate` | Interpolator-side drops with no recorded reason are reported as `unattributed` (in rows). Residual reports for fold and P2 cells (SC-003, SC-007, SC-008) are out of S5. S2 covers their data paths; S6 adds P2 state tests and no fold report test (OP22) |
| FR-025 | S6: `math/_transformation.py` (logger import, `_check_points`, `fit`, `_rotation`, new `transform` and `inverse_transform`), `math/__init__.py` export, `utils/_transformation.py` re-export (decisions 1 to 3) | `test_transformation.py` (96 explicit and 24 fitted round-trip cells, the fit defect tests, `test_entry_points_are_one_class`, the pickle tests), `test_euclidean_transformation.py` (hand-computed rotations, the scikit-learn tests) | OP23 (other copies not audited); OP24 (3D axis choice) |
| FR-026 | S1: values, gradient magnitudes and directions. S6: assembled and physical constraint residuals of the R2 frame, vector and scalar displacement, and the physical-residual fix for faulted features (OP19, accepted) | `test_planar_fault_frame_matches_analytic` (S1); `test_fault_frame_validation.py`; `test_diagnostics_backends.py` | None. C-019 replaced the trishear evaluator (OP18 closed); the 37 fault-plane points are reported, not bounded |
| FR-027 | S6 decision 5: integration error (closed form of the analytic frame and the Euler closed form of a curved control), fault-plane samples and outside-support samples assessed apart from the frame error; the zero-gradient guard of `apply_to_points` stays | `test_analytic_frame_displacement_matches_closed_form`, `test_integration_error_follows_the_euler_closed_form`, `test_displacement_populations_are_counted_separately`, `test_outside_support_points_have_no_displacement`, `test_zero_gradient_protection_remains_effective` | None |
| FR-028 | No signature change except new keyword arguments; `COMPAT.md` note; FR-019 exception recorded in ROADMAP.md; S4: `rescale` deprecated with its behaviour kept (`COMPAT.md`, `ROADMAP.md`, `API.md`) | API contract tests; existing round-trip tests unchanged; `test_rescale_warns_and_keeps_behaviour` | None |
| FR-029 | `verification.md` golden review: ledger of all baseline files changed in S1 to S5, S6 expects no change, tolerance audit, record of the deleted defect test, and a ledger of behaviour changes with analytic evidence but no golden file (faulted-feature residuals; the four `EuclideanTransformation` fixes) | Review record for each changed result | None |
| FR-030 | S1 and S2 Python cells; S6: `backend_matrix.py` (cells, `classify_cell`, merged junit files), native and Surfe boundary cells with a recording stub (mandatory, listed apart), optional extension cells (give the status of spec rows 6 and 7) | `test_backend_matrix.py`; matrix table in `verification.md` | Native and Surfe extension cells are excluded while their dependencies are absent (OP21). Native R3 in a rotated frame is not covered (fact 7) |
| FR-031 | Pre-fix rule in "Test strategy" | `verification.md` | None |
| FR-032 | "Stages and rollback" | Rollback checks in `verification.md` | None |
| SC-001 | `_rigid.py`, `BoundingBox` | `test_rigid_transform_contract.py` | None |
| SC-002 | R1 fit for the constant-norm variants (S2) | `test_constant_norm_frames.py` | Other R1 cells and lifecycle: S3. Native and Surfe R1: S6, optional cells (OP21). OP3 |
| SC-003 | Fault design | `test_fault_frame_analytic.py` | Run time OP2 |
| SC-004 | S6: displacement of the solved R2 frame against the analytic frame, library paths (C-019), 1332 smooth points | `test_vector_displacement_matches_analytic_frame`, `test_scalar_displacement_matches_analytic_frame`, `test_misplaced_frame_fails_the_displacement_bound` | None |
| SC-005 | Fault and enclosure design | `test_fault_buffer_bounds.py` | Run time OP2 |
| SC-006 | Format version 2; local arrays restored with no conversion | `test_round_trip_world_values`, `test_local_arrays_unchanged_after_two_cycles` | None |
| SC-007 | Reference sampling; S2 fold sampling and P2 | `test_reference_constraint_frames.py`; `test_fold_constraint_frames.py` | None |
| SC-008 | Reference sampling; S2 fold rows and P2 | `test_reference_constraint_frames.py` (FDI, PLI, P2); `test_fold_constraint_frames.py` | None |
| SC-009 | S5: misfit methods; `direction_cosines` | `test_misfit_leaves_state_unchanged`, `test_stale_misfit_is_not_rebuilt`, `test_gradient_cosine_matches_analytic_direction_exact_field`, `test_gradient_cosine_matches_analytic_direction_solved_model` (FDI and PLI; P2 not in S5) | None |
| SC-010 | S6: transformation fix and round-trip cells with 1, 4 and 50 points; the fit fixes with hand-computed numbers | `test_round_trip_explicit_parameters`, `test_round_trip_after_fit`, `test_round_trip_without_fit`, `test_fit_three_dimensions`, `test_three_dimensional_rotation_keeps_the_z_direction` | OP24 |
| SC-011 | `_cells_for_length` and the S5 setters (resize part); S5 cell-centre tests | `test_resize_paths.py`, `test_cell_centre_sampling.py` (first two tests) | OP12, OP17 |
| SC-012 | Rotation part in S1; metadata and legacy part in S3 | `test_rigid_transform_contract.py`; `test_legacy_*`, `test_malformed_metadata` | None |
| SC-013 | Extent-direction row; S2 fold R4 directions | `test_fault_frame_analytic.py`; `test_fold_constraint_frames.py` | None |
| SC-014 | Merge gates in every stage; S6 final gates (three pytest runs that write junit files, ruff, API tests) and the coverage ledger of production files | `pytest`, ruff, package tests, `test_public_api_contract.py`, `test_stable_api_surface.py`; ledger in `verification.md` | None |
| SC-015 | `verification.md` golden review (S1 to S5 ledger; S6 expects no change); tolerance audit (removed and added lines; the rename of `TOL_SOLVED` is the one change and it is a tightening or no change) | Review record; 0 relaxed tolerances | None |
| SC-016 | Reject at use; rebuild recipe in a new frame | `test_stale_frame_rejected_without_state_change`, `test_rebuild_in_new_frame_preserves_world_evaluation`, `test_fault_restoration_uses_current_frame` | Acceptance by automatic rebuild is not built (OP4) |
| SC-017 | Rollback column; format part S3 ("S3 rollback check") | `test_unsupported_future_version`, `test_old_reader_would_name_format_version`; `verification.md` | Old-reader error wording OP5 |
| SC-018 | S1 to S6 cells; `classify_cell`; exclusion ledger | Matrix table from `python -m tests.fixtures.backend_matrix <junit>...` (stub cells listed apart) | Native and Surfe extension cells excluded (OP21) |
| SC-019 | Gate checked on 2026-11-04: 003 spec (C-006, main and native worktree copies identical) and 004 spec (2026-11-03 source-005 entries) contain the amendments | Document check | None |

## Handover

- **Coordinate conversion (owner: `loop_interpolation`)**:
  `GeologicalInterpolator.world_to_local_points`, `local_to_world_points`,
  `world_to_local_vectors`, `local_to_world_vectors`. Spec 004 uses
  `world_to_local_vectors` to convert sampled world directions once into the
  computational basis before the low-level core. Spec 001 keeps its contract:
  the core gets local, unrotated grid inputs.
- **Rotation check (owner: `loop_common`)**:
  `loop_common.geometry.validate_rigid_rotation(matrix, dimensions, atol=1e-10)`
  and `RigidRotationError(ValueError)`. S3 uses them for saved metadata
  (`resolve_frame`).
- **Box enclosure (owner: `loop_common`)**:
  `BoundingBox.local_enclosure(origin=None, maximum=None) -> (local_min, local_max)`.
- **Callback contracts**: feature-facing reference callbacks and region
  callbacks receive world points. `add_gradient_orthogonal_constraints`
  receives target-local points and target-local unit vectors.
- **S2 callback contracts (owner: `loop_interpolation`)**: fold providers,
  fold `mask_fn`, `DirectionalRegularisation.direction` (callable or array)
  and the P2 `steepness_wtfunc` work in world coordinates. Helpers:
  `loop_interpolation._fold_setup.sample_fold_geometry(interpolator, fold, local_points)`
  and `world_norm_constraints(interpolator)`. Spec 004 builds directional
  terms on this contract. The native adapter follows it.
- **S3 saved frame and frame check (owner: `loop_interpolation`)**:
  `GeologicalInterpolator.to_dict` writes `format_version: 2`,
  `constraint_space: "local"`, `transform` and `bounding_box`. `from_dict`,
  `from_json` and `from_yaml` take `local_origin`, `local_rotation` and
  `declare_identity`. `InterpolatorFormatError(reason)` and `StaleFrameError`
  are exported. `GeologicalInterpolator.frame_signature() -> tuple[bytes,
  bytes]` is the key for any cache that depends on the frame; spec 003 keys
  its caches on it and clears them when it changes. `copy()` keeps the
  `bounding_box` object of its source.
- **S4 world data contract (owner: `LoopStructural`)**: every
  `get_data` of a feature returns `ValuePoints` and `VectorPoints` with world
  locations and world vectors; `Loop2DView.add_data` and `ModelExporter.save`
  use them as they are. `GeologicalModel.rescale` (local to world) and `scale` (world to local) are
  deprecated stable API. They warn (`DeprecationWarning`; use
  `bounding_box.reproject` and `bounding_box.project`). A caller must not pass
  `rescale` world points.
  `IntrusionFrameBuilder._network_model_values(data)` evaluates the model
  at the world X, Y and Z of the network data. Spec 004 and the QGIS plugin read
  display data through `get_data` and get world locations.
- **Exclusion record**: `GeologicalFeatureBuilder.excluded_rows:
  dict[str, dict[str, int]]`. S5 connects it to the diagnostics report.
- **S5 diagnostics (owner: `loop_interpolation`)**: `ResidualSummary`,
  `ExclusionRecord`, `direction_cosines(model_vectors, data_vectors)` and
  `build_exclusion_records(excluded_rows, families)` are exported.
  `ConstraintDiagnosticsReport` has the optional fields `physical_residuals`,
  `assembled_residuals` and `exclusions`; `to_dict` shows them only when they
  are not empty. `GeologicalFeature.verification_report()` fills them, makes
  the private call `_build_constraint_diagnostics_report`, and does not store the
  report. `summarise_residuals` and `UNKNOWN_EXPLANATION` are imported from
  `loop_interpolation._diagnostics` and are not exported. The feature misfit
  methods `evaluate_gradient_misfit` and `evaluate_value_misfit` keep their names
  and public shape, work on copies and never rebuild. A stale model gives one
  `UserWarning` and NaN. The gradient value is now a cosine in world coordinates.
  The new methods `evaluate_gradient_cosines` and `evaluate_value_residuals` give
  1-D arrays. Spec 004 can report its directional rows through
  `ResidualSummary(kind="assembled")`. Structured 3D grids keep origin and extent
  on every setter and on `set_nelements`. The 2D `maximum` keeps its legacy value.
  A rectilinear grid refuses to resize.
- **Fault geometry**: an explicit `fault_center` is a world point; the
  generated centre extent row is `coord=2` with `strike_unit / major_axis`.
  Spec 003 golden fault models can change only under C-006/C-004 with a
  record in `verification.md`. Spec 003 cache invalidation for transform
  changes uses `frame_signature` (S3).
- **S6 final handover (spec 005 complete)**: Spec 005 is complete after S6. The items
  below are for specs 003 and 004.
  - `EuclideanTransformation` is one class. It is exported by `loop_common.math` and
    re-exported by `LoopStructural.utils` and `LoopStructural.utils._transformation`.
    `transform` and `inverse_transform` take an `(N, d)` array, return a new float array,
    and keep the columns after `dimensions` bit for bit. The 3D `fit` works, the fitted
    angle has the right sign, the z axis is not reversed, and `fit` with rotation raises
    `ImportError` when scikit-learn is missing. Numbers from the old 3D class are wrong
    (`COMPAT.md`, `ROADMAP.md`).
  - The physical residual diagnostics of a faulted feature evaluate the interpolator at
    the stored, already restored positions. Spec 003 can change how points are restored
    (caches) if stored rows stay restored world positions. A diagnostic must never
    restore a stored row a second time.
  - Test fixtures for reuse: `tests/fixtures/coordinate_frames.py` (`r2_displacement_section`,
    `r2_analytic_fault`, `r2_closed_form_displacement`, `r2_curved_slip_fault`,
    `build_r2_fault`, `faulted_foliation_model`, `residuals_of_field`)
    and `tests/fixtures/backend_matrix.py` (`MATRIX`, `Cell`, `classify_cell`,
    `outcomes_from_junit`, `require_optional`, `OPTIONAL_DEPENDENCIES`). A new native or RBF
    engine of spec 003 or 004 adds one optional `Cell` with its dependency name, and a
    mandatory stub cell (listed apart) for its Python side. It must not become a
    prerequisite (FR-030).
  - Native adapter contract (unchanged in S6): local rows, a rotated frame is accepted,
    a rotated support is rejected. `test_backend_world_contract.py` is the model for the
    boundary test of a new adapter.
  - Golden rule: in spec 005 only `packages/loop_common/tests/data/structured_support_reference.npz`
    changed (S1 and S5; records in `verification.md`). The baselines and fingerprints of
    specs 003 and 004 are checked in their own worktrees under C-004. S6 did not touch
    them.
  - Open items passed on, each with a default: OP21, OP22, OP23 and OP24. Findings
    that stay: the format-only hunks in `_structured_grid_2d.py`, the node count of
    `set_nelements`, the S1 note on anisotropy rows in `_on_geometry_change`, the shape
    error in the public `evaluate_gradient` of a faulted feature (S3 note), and issue
    #399 (P2 fit accuracy).
- **Shared files**: `_geological_interpolator.py`, `_discrete_interpolator.py`
  and `_fault_builder.py` are also in the scope of specs 003 and 004.
  Coordinate the merge order with the active 003 worktree
  (`LoopStructural-native`).

## Complexity tracking

None.
