# Tasks: Coordinate and fault-frame correctness (stages S1 to S6)

**Feature**: `005-coordinate-frame-correctness` | **Plan**: [plan.md](plan.md)

<!--
Rules:
- Each task is small (less than one day) and names the files that it changes.
- Test tasks come BEFORE the implementation tasks that they test.
  Run the test and see it fail first.
- [P] = can be done at the same time as the other [P] tasks of the same wave.
- Written in ASD-STE100 Simplified Technical English.
- This file covers stages S1 (T001 to T017, done), S2 (T018 to T032), S3
  (T033 to T050, done), S4 (T051 to T068, done), S5 (T069 to T088, done) and S6 (T089 to T108) of
  plan.md. No S6 task is blocked: OP21 to OP24 have defaults.
-->

## Orchestrator notes

- Run this file with the spec-orchestrator skill. It has the rules for the
  worktree, the sub agents, the checks and the commits.
- Waves: 52 (S1: 1 to 7, done; S2: 8 to 13; S3: 14 to 22, done; S4: 23 to 31, done; S5: 32 to 42, done; S6: 43 to 52). Tiers: [S] = simple model, [M] = standard model.
- Do not start a task with "Blocked by".

## Phase 1: Setup

- [x] T001 [P] [Setup] [S] Create the S1 verification record skeleton.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: none
  - Wave: 1
  - Done when: the file exists with the 6 headings below, each followed by
    the line "No entries yet."
  - Context: Title `# Verification: 005-coordinate-frame-correctness`.
    Headings, in this order: `## Pre-fix failures`, `## Post-fix results`,
    `## Golden review (C-004)`, `## Backend matrix`, `## Rollback checks`,
    `## Test run time (OP2)`. Under the title, add one sentence: "This file
    records the evidence for FR-029, FR-031 and FR-032 of spec 005."
    Write in ASD-STE100 Simplified Technical English.

- [x] T002 [P] [Setup] [M] Create the shared analytic fixtures R1 to R5.
  - Files: `tests/fixtures/coordinate_frames.py`
  - Needs: none
  - Wave: 1
  - Done when: `pytest tests/fixtures/coordinate_frames.py` collects 0
    tests with no import error, and `ruff check LoopStructural packages`
    passes.
  - Context: Helper module, no tests. numpy docstrings and type hints; Black
    line length 100. Names: `ANCHORS` = arrays (0,0,0), (100,200,300),
    (-100,-200,-300), (500000,6500000,100); `ROTATIONS_DEG = (0.0, 30.0, 90.0)`;
    `rot_z(deg: float) -> np.ndarray` = [[c,-s,0],[s,c,0],[0,0,1]];
    `tp(L: float, W: float) -> float` = `1e-10*L + 32*np.finfo(float).eps*max(W, L)`;
    `model_with_frame(origin, maximum, rotation_deg) -> GeologicalModel` makes
    `GeologicalModel(origin, maximum)` then calls
    `model.bounding_box.set_local_transform(local_origin=origin, rotation_matrix=rot_z(rotation_deg))`.
    `r2_fault_data(anchor)`: DataFrame, columns X,Y,Z,val,coord,feature_name
    ("fault"),w; 25 rows at x'=0, y',z' in {-8,-4,0,4,8}, val 0, coord 0, w 1.
    `r5_fault_data(anchor)`: same, y',z' in {-4,-2,0,2,4}.
    `r2_probes(anchor)`: 6859 points `anchor + {-9..9}^3`.
    `r2_frame(points, anchor) -> (values (N,3), gradients (N,3,3))`:
    gx=x'/20, gy=-z'/20, gz=y'/20. `r3_value`/`r3_gradient`: h=x'+0.05*y'^2+z',
    grad (1,0.1*y',1). `r4_frame(points, anchor)`: R=40, X=R+x',
    r=sqrt(X^2+y'^2); surface (r-R)/20, slip -z'/20, extent R*atan2(y',X)/20,
    with the analytic gradients of spec.md R4. `class AnalyticReference` with
    `__init__(value_fn, gradient_fn, name, support_box=None)`,
    `evaluate_value(points)`, `evaluate_gradient(points, ignore_regions=True)`,
    `__call__ = evaluate_value`; NaN outside the closed `support_box`.

## Phase 2: Foundation (blocks all user stories)

- [x] T003 [P] [Found] [M] Test: rigid rotation, round trips, enclosure and resize stability.
  - Files: `packages/loop_common/tests/test_rigid_transform_contract.py`
  - Needs: none
  - Wave: 1
  - Done when: `pytest packages/loop_common/tests/test_rigid_transform_contract.py`
    fails before T004 (import error), then passes after T005 and T006.
  - Context: Self-contained (do not import `tests/fixtures`). Tests:
    `test_validate_rigid_rotation_accepts_proper_rotations` (0, 30, 90 deg
    about z; 2D and 3D); `test_validate_rigid_rotation_rejects_invalid`
    (scale, shear, singular, reflection det -1, NaN, inf, wrong shape) expects
    `RigidRotationError`, a `ValueError` subclass;
    `test_point_and_direction_round_trips` (SC-001: `BoundingBox.project` then
    `reproject` within `Tp = 1e-10*L + 32*eps*max(W,L)`; `project_vectors` then
    `reproject_vectors` within 1e-12; 1 point and arrays; anchors (0,0,0),
    (100,200,300), (-100,-200,-300), (500000,6500000,100));
    `test_local_enclosure_contains_projected_corners` (all 8 projected corners
    inside `BoundingBox.local_enclosure()` within Tp, for each rotation);
    `test_invalid_rotation_keeps_prior_state` (`set_local_transform` and
    `StructuredGrid3DGeometry.rotation_xy` keep the old matrix after the error);
    `test_reapply_unchanged_bounds_is_stable` (set `origin` and `maximum` of a
    `StructuredGrid3DGeometry` to their own values 10 times: 0 changes to
    origin, maximum, step_vector, nsteps). Import from `loop_common.geometry`.

- [x] T004 [Found] [S] Add `validate_rigid_rotation` and `RigidRotationError`.
  - Files: `packages/loop_common/src/loop_common/geometry/_rigid.py`,
    `packages/loop_common/src/loop_common/geometry/__init__.py`
  - Needs: T003
  - Wave: 2
  - Done when: `pytest packages/loop_common/tests/test_rigid_transform_contract.py -k validate_rigid_rotation`
    passes.
  - Context: New module, imports numpy only. `class RigidRotationError(ValueError)`.
    `def validate_rigid_rotation(matrix: np.ndarray, dimensions: int, atol: float = 1e-10) -> np.ndarray`:
    return `np.array(matrix, dtype=float)` (a copy). Raise
    `RigidRotationError` when shape != (dimensions, dimensions), any entry is
    not finite, `max|R.T @ R - I| > atol`, or `|det(R) - 1| > atol`. Message
    names the failed condition and ends with "supply a proper rotation
    (orthonormal, det +1); scale and shear are not supported". In
    `__init__.py` add `from ._rigid import RigidRotationError, validate_rigid_rotation`
    as the first import line. numpy docstrings, type hints, py39 (use
    `from __future__ import annotations`).

- [x] T005 [P] [Found] [M] Validate the local transform and add `BoundingBox.local_enclosure`.
  - Files: `packages/loop_common/src/loop_common/geometry/_bounding_box.py`
  - Needs: T004
  - Wave: 3
  - Done when: `pytest packages/loop_common/tests/test_rigid_transform_contract.py -k "round_trips or local_enclosure or prior_state"`
    and `pytest packages/loop_common/tests` pass.
  - Context: In `set_local_transform(local_origin=None, rotation_matrix=None)`
    call `validate_rigid_rotation(rotation_matrix, self.dimensions)` (import
    from `._rigid`) and reject a non-finite `local_origin` with
    `LoopValueError`, both before any field changes. New method
    `def local_enclosure(self, origin: np.ndarray | None = None, maximum: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]`:
    defaults `self.origin`, `self.maximum`; raise `LoopValueError` if a bound
    is None, not finite, or origin > maximum on an axis; build the
    `2**dimensions` corners with `itertools.product`; return
    `(min, max)` over axis 0 of `self.project(corners)`. numpy docstring.
    Existing tests in `packages/loop_common/tests/test_bounding_box.py` use
    proper rotations and must still pass.

- [x] T006 [P] [Found] [M] Repair support resize and validate `rotation_xy`.
  - Files: `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`,
    `packages/loop_common/tests/data/structured_support_reference.npz`
  - Needs: T004
  - Wave: 3
  - Done when: `pytest packages/loop_common/tests/test_rigid_transform_contract.py -k "prior_state or unchanged_bounds"`
    and `pytest packages/loop_common/tests` pass.
  - Context: `_nsteps` is a NODE count; `nsteps_cells = nsteps - 1`;
    `maximum = origin + nsteps_cells * step_vector`. Defect: the `origin` and
    `maximum` setters set `_nsteps = ceil(length/step)` (a cell count), so the
    support loses 1 cell. Add
    `def _cells_for_length(self, length: np.ndarray) -> np.ndarray`:
    `cells = np.ceil(length / self.step_vector - 1e-9).astype(np.int64)`;
    axes with 0 cells get 2 cells; raise `ValueError("Cannot resize the grid.")`
    if a count < 0. `origin` setter: compute `self.maximum` first, then
    `self._nsteps = self._cells_for_length(maximum - origin) + 1`.
    `maximum` setter: `self._nsteps = self._cells_for_length(maximum - self.origin) + 1`.
    Keep `_on_geometry_change()` calls. `rotation_xy` setter: after the
    float/int conversion call `validate_rigid_rotation(rotation_xy, 3)`
    (import from `._rigid`); `None` keeps the current value.
    Answer (plan.md "Answers", 2026-11-04): regenerate the golden file with
    `python packages/loop_common/tests/test_structured_support_reference.py --update`;
    only `after_origin`/`after_maximum` `nsteps` and `maximum` may change.

## Phase 3: User story 2 - Constrain features at the intended physical points (P1)

**Goal**: reference relationships sample at world points and install rows in
the target-local basis.
**Independent test**: `pytest tests/unit/modelling/test_reference_constraint_frames.py`.

- [x] T007 [P] [US2] [M] Test: R3 reference relationships for FDI and PLI.
  - Files: `tests/unit/modelling/test_reference_constraint_frames.py`
  - Needs: T002
  - Wave: 2
  - Done when: `pytest tests/unit/modelling/test_reference_constraint_frames.py`
    fails before T008 to T010, then passes after them. Record the failing
    test ids in the task report.
  - Context: Use `tests/fixtures/coordinate_frames.py` (ANCHORS,
    ROTATIONS_DEG, model_with_frame, AnalyticReference, r3_value, r3_gradient,
    tp). Box `O` to `O+(8,10,12)`; support 5x6x7 cells (336 nodes); engines
    "FDI" and "PLI". Tests: `test_world_local_conversion_methods`
    (`interpolator.world_to_local_points` etc. invert each other within Tp);
    `test_equality_region_selects_world_nodes` (region `x' > 4`: 168 of 336
    unrotated; for rotations compute the expected set from world nodes; 0
    missing or duplicate; value = h at the world node);
    `test_orthogonal_rows_affine_residual` (field f=y', direction (1,0,0)
    target 0, direction (0,1,0) target 1; rows `f"{name}_orthogonal"`;
    residual `A @ c_exact - b` at most 1e-10 for sample_rule None,
    "centroid", "gauss4");
    `test_inequality_feature_bounds_world` (bounds = h at world nodes);
    `test_restricted_reference_exclusion_counts` (box O+(2.1,2.1,2.1) to
    O+(5.9,7.9,9.9); `builder.excluded_rows` count equals the analytic count);
    `test_support_expansion_contains_restored_points` (FR-009);
    `test_set_interpolation_geometry_encloses_rotated_box`.

- [x] T008 [US2] [S] Add the public world/local conversion methods.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_geological_interpolator.py`
  - Needs: T007
  - Wave: 3
  - Done when: `pytest tests/unit/modelling/test_reference_constraint_frames.py -k conversion`
    passes.
  - Context: Add four public methods after `_project_constraint_array`:
    `world_to_local_points(self, points: np.ndarray) -> np.ndarray` returns
    `self._project_points(points)`; `local_to_world_points` returns
    `self._reproject_points(points)`; `world_to_local_vectors` returns
    `self._project_vectors(vectors)`; `local_to_world_vectors` returns
    `self._reproject_vectors(vectors)`. Keep the private methods. Docstrings
    (numpy style): input `(N, d)` or `(d,)`, same output shape; vectors use
    rotation only; with `bounding_box is None` the input is returned.

- [x] T009 [P] [US2] [S] Sample `add_inequality_feature` at world nodes.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`
  - Needs: T008
  - Wave: 4
  - Done when: `pytest tests/unit/modelling/test_reference_constraint_frames.py -k inequality`
    passes.
  - Context: In `add_inequality_feature(self, feature, lower=True, mask=None)`
    replace `value = feature(self.support.nodes)` with
    `value = feature(self.local_to_world_points(self.support.nodes))`. Add to
    its docstring: "`feature` receives world points and returns one value for
    each support node." In the docstring of the abstract
    `add_gradient_orthogonal_constraints(self, points, vectors, w=1.0)` state:
    "`points` are target-local support points; `vectors` are target-local
    unit directions." No signature changes.

- [x] T010 [P] [US2] [M] Sample reference features at world points in the feature builder.
  - Files: `LoopStructural/modelling/features/builders/_geological_feature_builder.py`
  - Needs: T005, T008
  - Wave: 4
  - Done when: `pytest tests/unit/modelling/test_reference_constraint_frames.py -k "not inequality"`
    and `pytest tests/unit/modelling` pass.
  - Context: Add `self.excluded_rows: dict[str, dict[str, int]] = {}` in
    `__init__`. Reasons: "invalid_reference_value",
    "invalid_reference_direction", "evaluation_error".
    `install_gradient_constraint`: local points (existing) -> `world =
    self.interpolator.local_to_world_points(pts)`; `vector =
    feature.evaluate_gradient(world, ignore_regions=True)`; `valid` = finite
    and norm > 0; record `(~valid).sum()`; normalise; `world_to_local_vectors`;
    pass `pts[valid]`, local vectors, `weight[valid]` if array. 0 valid -> no
    rows, warning. `install_equality_constraints`: world nodes, `idc =
    np.flatnonzero(region(world_nodes))`, `val =
    feature.evaluate_value(world_nodes[idc])`, finite mask, name
    `f"{feature.name}_equality"`; the `except` records "evaluation_error".
    `set_interpolation_geometry`: use
    `self.interpolator.bounding_box.local_enclosure(origin, maximum)`.
    `check_interpolation_geometry`: `regular_grid(local=False)`, restore in
    world, `world_to_local_points`, min/max on copies of support arrays.

## Phase 4: User story 1 - Build a fault frame at its supplied location (P1)

**Goal**: the public fault API makes the analytic planar frame at the supplied
world location, with the requested buffer.
**Independent test**: `pytest tests/unit/modelling/test_fault_frame_analytic.py tests/unit/modelling/test_fault_buffer_bounds.py`.

- [x] T011 [P] [US1] [M] Test: R2 planar frame, centre handling and R4 extent row.
  - Files: `tests/unit/modelling/test_fault_frame_analytic.py`
  - Needs: T002
  - Wave: 2
  - Done when: `pytest tests/unit/modelling/test_fault_frame_analytic.py`
    fails before T012 and T013, then passes after them.
  - Context: Build through `model.create_and_add_fault("fault", 0, data=...,
    interpolatortype=..., fault_center=O or None, fault_normal_vector=(1,0,0),
    fault_slip_vector=(0,0,-1), major_axis=20, minor_axis=10,
    intermediate_axis=20, fault_buffer=0.2, nelements=20000 (FDI) or 1000
    (PLI), fault_trace_anisotropy=0, fault_dip_anisotropy=0)`. Solve: CG rtol
    1e-10, atol 1e-12, maxiter 20000; orthogonality weight 1, step 10.
    Tests: `test_explicit_centre_rows_not_projected` (generated rows use O);
    `test_inferred_centre_keeps_elevation` (inferred centre == O, all 3
    components); `test_slip_vector_made_orthogonal` and
    `test_parallel_slip_raises` (ValueError);
    `test_extent_direction_row_at_centre` (one coord-2 row at O with
    n = (0,1/20,0); equals the R4 analytic extent gradient at O within 1e-10;
    R4 extent directions at O+(0,-6,0), O+(0,6,0) differ >= 15 degrees);
    `test_planar_fault_frame_matches_analytic` (SC-003 on `r2_probes`: scalar
    error <= 1e-6, relative gradient magnitude <= 1e-4, direction <= 0.01 deg;
    4 anchors x 3 rotations x 2 centre modes x FDI, PLI).

- [x] T012 [US1] [M] Repair centre, slip, extent row and buffer state in `FaultBuilder`.
  - Files: `LoopStructural/modelling/features/builders/_fault_builder.py`
  - Needs: T011, T014
  - Wave: 3
  - Done when: `pytest tests/unit/modelling/test_fault_frame_analytic.py -k "inferred_centre or slip or extent_direction"`
    passes.
  - Context: `__init__`: `self._fault_buffer: float = 0.2`,
    `self._force_mesh_geometry: bool = False`. `create_data_from_geometry`
    stores both from its arguments first. `_on_geometry_change`: default
    `force_mesh_geometry=None`; `None` args use the stored values. Inferred
    centre: rows `coord == 0` and `val == 0`, columns X,Y,Z, `np.nanmean`;
    a non-finite component -> 0.0 plus a warning naming the axis. Slip:
    `s_perp = s - (s @ n) * n`; if `|s_perp| < 1e-12` raise
    `ValueError("fault slip vector is parallel to the fault normal; supply a slip direction in the fault plane")`;
    use the unit `s_perp` for `fault_slip_vector`, strike and depth rows.
    `strike_unit = cross(n, s_perp)/|...|`. When `major_axis` is not None add
    one row: X,Y,Z = centre, nx,ny,nz = `strike_unit / major_axis`, val NaN,
    coord 2, w. Early return (no trace): call
    `self.set_mesh_geometry(fault_buffer, None)` before `return`.

- [x] T013 [US1] [S] Stop the centre projection and the buffer override in the fault factory.
  - Files: `LoopStructural/modelling/core/_model_feature_factory.py`,
    `LoopStructural/modelling/features/fault/_fault_segment.py`
  - Needs: T006, T010, T012
  - Wave: 5
  - Done when: `pytest tests/unit/modelling/test_fault_frame_analytic.py tests/unit/modelling/test_fault_buffer_bounds.py`
    passes.
  - Context: In `create_and_add_fault`, delete the two lines
    `if fault_center is not None and ~np.isnan(fault_center).any():` /
    `fault_center = model.scale(fault_center, inplace=False)`. Delete the two
    lines `if "force_mesh_geometry" not in kwargs:` /
    `fault_frame_builder.set_mesh_geometry(kwargs.get("fault_buffer", 0.2), 0)`.
    `create_data_from_geometry` already applies `fault_buffer` and
    `force_mesh_geometry`. No signature change. An explicit `fault_center`
    is a world point (FR-010).
    Answer (plan.md "Answers", 2026-11-04): in `FaultSegment.fault_ellipsoid`
    use `self.fault_centre[None, :]` with no `model.rescale`; the centre is a
    world point.

- [x] T014 [P] [US1] [M] Test: R5 buffer and forced-geometry bounds.
  - Files: `tests/unit/modelling/test_fault_buffer_bounds.py`
  - Needs: T002
  - Wave: 2
  - Done when: `pytest tests/unit/modelling/test_fault_buffer_bounds.py`
    fails before T012 and T013, then passes after them.
  - Context: Use `r5_fault_data`, `model_with_frame`, `rot_z`, `tp`. Box
    O-(10,10,10) to O+(10,10,10); normal (1,0,0), slip (0,0,-1); axes major 12,
    minor 4, intermediate 12; `force_mesh_geometry=True`; buffers b = 0, 0.2,
    0.5 through `create_and_add_fault`; FDI nelements 20000, PLI 1000;
    explicit and inferred centre; rotations 0, 30, 90 deg. Requested world
    box `O +/- (10+20*b)`; expected local box = enclosure of its 8 projected
    corners. Tests: `test_buffer_bounds_match_request` (each support face is
    outward of the expected bound by 0 to 1 cell (actual local spacing) plus
    Tp; all generated rows inside the support);
    `test_buffers_give_distinct_bounds`;
    `test_property_change_keeps_buffer_and_forced_geometry` (set
    `builder.fault_normal_vector = (1,0,0)` after creation; bounds unchanged).

- [x] T015 [US1] [M] Align `test__fault_builder.py` with the corrected behaviour.
  - Files: `tests/unit/modelling/test__fault_builder.py`,
    `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T001, T013
  - Wave: 6
  - Done when: `pytest tests/unit/modelling/test__fault_builder.py` passes.
  - Context: Change an assertion only if it encodes a defect that S1 fixed
    (for example an inferred centre with z = 0, a projected explicit centre,
    or a default buffer that replaced the requested one). For each change add
    one entry under `## Golden review (C-004)` in `verification.md`: test id,
    old value, new value, defect, expected physical effect. Do not relax a
    tolerance (SC-015). If no change is needed, write "No golden change in
    S1." under that heading.
    Answer (plan.md "Answers", 2026-11-04): also add one entry for the
    T006 change to `packages/loop_common/tests/data/structured_support_reference.npz`
    (27 values of `after_origin`/`after_maximum`; old values did not enclose
    the requested bounds; FR-014).

## Phase 5: Finish

- [x] T016 [P] [Finish] [S] Add the S1 migration notice to `COMPAT.md`.
  - Files: `COMPAT.md`
  - Needs: T013
  - Wave: 6
  - Done when: the "Migration notices" table in `COMPAT.md` has the new row
    with date 2026-11-04.
  - Context: Add one row at the end of the table under
    `## Migration notices (not breaking, no shim needed)`. Columns: Old |
    New | Since. Old: "`create_and_add_fault(fault_center=...)` projected the
    centre into local coordinates; `fault_buffer`/`force_mesh_geometry` were
    replaced by a default 0.2 buffer; an inferred centre had z = 0; any 3x3
    matrix was accepted as a rotation". New: "`fault_center` is a world point;
    the requested buffer and forced geometry are kept; an inferred centre keeps
    its elevation; `set_local_transform` and `rotation_xy` raise
    `RigidRotationError` for a non-rigid matrix. Extra observations added to
    work around these defects are no longer necessary (spec 005)". Since:
    2026-11-04.

- [x] T017 [Finish] [M] Record pre-fix evidence, gates, run time and rollback.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T003, T004, T005, T006, T007, T008, T009, T010, T011, T012, T013, T014, T015, T016
  - Wave: 7
  - Done when: `ruff check LoopStructural packages && pytest` and
    `pytest packages/loop_common/tests packages/loop_interpolation/tests`
    pass, and each heading of `verification.md` has an entry.
  - Context: Pre-fix evidence (FR-031): make a throw-away worktree at commit
    `ebc7c699` (`git worktree add <tmp> ebc7c699`), copy the 4 new test files
    and `tests/fixtures/coordinate_frames.py` into it, run them with
    `PYTHONPATH=<tmp>;<tmp>/packages/loop_common/src;<tmp>/packages/loop_interpolation/src`
    (the repo `.venv` is editable on the main checkout), record each failing
    test id and value under `## Pre-fix failures`, then
    `git worktree remove <tmp>`. A finding that does not fail is "not
    reproduced". Record pass results, the backend cells covered (Python FDI,
    PLI/P1), the S1 test time (OP2), and the rollback check: on the
    worktree at `ebc7c699`, `tests/unit/modelling/test__fault_builder.py` and
    `packages/loop_common/tests/test_bounding_box.py` pass.

## Phase 6: Stage S2 foundation

Stage S2 is in plan.md, "Stage S2 design". All S2 tasks come after S1
(T017). Environment rule for every S2 task: the shared venv is editable on
the main checkout. Before a test, set `PYTHONPATH` to the worktree root,
`<worktree>/packages/loop_common/src` and
`<worktree>/packages/loop_interpolation/src`.

- [x] T018 [P] [Found] [S] Add the R1 helpers and `AnalyticFold` to the shared fixtures.
  - Files: `tests/fixtures/coordinate_frames.py`
  - Needs: T017
  - Wave: 8
  - Done when: `pytest tests/fixtures/coordinate_frames.py` collects 0
    tests with no import error, and `ruff check LoopStructural packages`
    passes.
  - Context: Helper module, no tests. numpy docstrings, type hints, Black
    line length 100. Add: `R1_GRADIENT = np.array([0.2, -0.3, 0.4])`;
    `r1_value(points, anchor) -> np.ndarray` = `2 + (points - anchor) @ R1_GRADIENT`;
    `r1_gradient(points, anchor) -> np.ndarray` = `R1_GRADIENT` tiled to `(N, 3)`;
    `r1_observations(anchor) -> pd.DataFrame`: 27 rows at
    `anchor + (x', y', z')` with x' in {1,4,7}, y' in {1,5,9}, z' in {1,6,11};
    columns X, Y, Z, val (`r1_value`), feature_name "r1", w 1.
    `class AnalyticFold`: `__init__(self, deformed_fn, axis_fn, dgz_fn)`,
    `self.seen_points: list[np.ndarray] = []`;
    `get_deformed_orientation(self, points)` appends `np.array(points, copy=True)`
    and returns `(deformed_fn(points), axis_fn(points), dgz_fn(points))`,
    each `(N, 3)` float. `@classmethod constant(cls, deformed, axis, dgz)`
    tiles three world vectors. `@classmethod r4(cls, anchor)`: deformed =
    unit of `r4_frame(points, anchor)[1][2]` (extent gradient), axis =
    `(0, 0, 1)`, dgz = unit of `r4_frame(points, anchor)[1][0]`.

- [x] T019 [P] [Found] [S] Add `sample_fold_geometry` and `world_norm_constraints`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_fold_setup.py`
  - Needs: T017
  - Wave: 8
  - Done when: `python -c "from loop_interpolation._fold_setup import sample_fold_geometry, world_norm_constraints"`
    runs with the worktree `PYTHONPATH`, and `ruff check LoopStructural packages` passes.
  - Context: Add `from __future__ import annotations` and `import numpy as np`.
    `def sample_fold_geometry(interpolator, fold, local_points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]`:
    `world_points = interpolator.local_to_world_points(np.asarray(local_points, dtype=float))`;
    `deformed, axis, dgz = fold.get_deformed_orientation(world_points)`;
    return `(world_points, np.asarray(deformed, dtype=float), np.asarray(axis, dtype=float), np.asarray(dgz, dtype=float))`.
    `def world_norm_constraints(interpolator) -> np.ndarray`:
    `n = np.asarray(interpolator.get_norm_constraints(), dtype=float)`; if
    `n.ndim != 2 or n.shape[0] == 0 or n.shape[1] < 6` return `n`; else
    `out = n.copy()`, `out[:, 0:3] = interpolator.local_to_world_points(n[:, 0:3])`,
    `out[:, 3:6] = interpolator.local_to_world_vectors(n[:, 3:6])`, return
    `out`. numpy docstrings: the fold provider takes world points and returns
    world unit directions; NaN rows stay NaN. Keep `setup_with_fold_constraints`.

## Phase 7: User story 2 - fold, P2 and constant-norm cells (S2, FR-008)

**Goal**: fold geometry, directional terms, masks and the fold sign check
use world points and local directions, and the P2 and constant-norm cells
satisfy R3 and R1.
**Independent test**: `pytest tests/unit/modelling/test_fold_constraint_frames.py tests/unit/modelling/test_constant_norm_frames.py`.

- [x] T020 [US2] [M] Test: fold, directional-regularisation and P2 weight-callback frames (R3, R4).
  - Files: `tests/unit/modelling/test_fold_constraint_frames.py`
  - Needs: T018, T019
  - Wave: 9
  - Done when: `pytest tests/unit/modelling/test_fold_constraint_frames.py`
    fails before T023 to T028 (record the failing ids), then passes after them.
  - Context: Read plan.md "Stage S2 design", "S2 tests", for the 7 tests and
    their limits. Build each interpolator as `_builder` in
    `tests/unit/modelling/test_reference_constraint_frames.py` does: box `O`
    to `O+(8,10,12)` from `model_with_frame`, 5x6x7 cells, support from the
    local enclosure, `interp.bounding_box = model.bounding_box`. Families
    and pytest ids: `fdi_fold` = `FDFoldInterpolator(StructuredGrid)`,
    `pli_fold` = `DiscreteFoldInterpolator(TetMesh)`, `p2_fold_centroid` and
    `p2_fold_gauss4` = `P2FoldInterpolator(P2TetMesh)` with that
    `sample_rule`. Directional and weight tests use ids `fdi`, `pli`, `p2`.
    Set `interp.fold = AnalyticFold...` and call `interp.add_fold_constraints(...)`
    directly. Read each family's `add_fold_constraints` for its arguments.
    Exact field for residuals: `f = y'`, coefficients from
    `bounding_box.reproject_vectors(nodes)[:, 1] + (local_origin[1] - anchor[1])`.
    Rows are `interp.constraints[name]` with `"matrix"` and `"b"`. If P2
    rejects directional regularisation (for example `NotImplementedError`),
    test that documented rejection for `p2` instead. No tolerance above the
    plan limits (SC-015).

- [x] T021 [P] [US2] [M] Test: add the P2 engine to the R3 reference tests.
  - Files: `tests/unit/modelling/test_reference_constraint_frames.py`
  - Needs: T017
  - Wave: 8
  - Done when: `pytest tests/unit/modelling/test_reference_constraint_frames.py -k P2`
    passes. Plan.md says that the P2 reference path needs no code change, so
    this test passes on the S1 code; T031 records its pre-fix result at the
    S1 base.
  - Context: Add `"P2"` to each `engine` parametrisation. In `_builder`,
    for `"P2"` use `P2TetMesh` from `loop_common.supports` and
    `P2Interpolator` from `loop_interpolation` (same origin, `nsteps_cells=[5, 6, 7]`
    and step as the other engines). Compute the expected equality set from
    the world P2 nodes (`local_to_world_points(support.nodes)`); do not
    reuse 168. Keep all tolerances (SC-007 `Tp`, SC-008 `1e-10`). If a P2
    test fails, report the first error. Do not change production code.

- [x] T022 [P] [US2] [M] Test: constant-norm FDI and PLI fits under rotated frames (R1, SC-002).
  - Files: `tests/unit/modelling/test_constant_norm_frames.py`
  - Needs: T018
  - Wave: 9
  - Done when: `pytest tests/unit/modelling/test_constant_norm_frames.py -k cn_pli`
    passes, and the `cn_fdi` cells fail with the `IndexError` at
    `_constant_norm.py` `add_constant_norm` (T032 fixes it; plan.md
    "Answers"). If the unrotated anchor-0 PLI control fails SC-002, stop and
    report the errors (plan OP3). Do not relax a limit.
  - Context: Classes `ConstantNormFDIInterpolator` and
    `ConstantNormP1Interpolator` from `loop_interpolation` (read
    `packages/loop_interpolation/src/loop_interpolation/_constant_norm.py`
    for the norm-target attribute and the solve API). Box `O` to
    `O+(8,10,12)` from `model_with_frame`, 4 `ANCHORS` x `ROTATIONS_DEG`,
    216 requested cells (6x6x6) or the matching `TetMesh`. Data:
    `r1_observations(anchor)` as world value constraints. Norm target
    `np.sqrt(0.29)`. Solve with CG and
    `solver_kwargs={"rtol": 1e-10, "atol": 1e-12, "maxiter": 20_000}`.
    Query at `anchor + (2.2, 3.7, 4.3)` and `anchor + (6.1, 8.2, 9.1)`.
    Limits: value error `<= 1e-8 * F` with `F = max(1, max|r1_value|)`;
    gradient component error `<= 1e-8 * max(1, |R1_GRADIENT|)`. Ids
    `cn_fdi`, `cn_pli`. Run the unrotated anchor-0 control first.
    Measured PLI control (2026-11-04): value error 1.8e-12, gradient error
    1.3e-10. Do not mark the `cn_fdi` cells xfail or skip them.

- [x] T032 [P] [US2] [S] Accept a scalar `element_size` in `add_constant_norm`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_constant_norm.py`
  - Needs: T022
  - Wave: 10
  - Done when: `pytest tests/unit/modelling/test_constant_norm_frames.py`
    passes (all `cn_fdi` and `cn_pli` cells), and
    `pytest packages/loop_interpolation/tests -k constant_norm` passes.
  - Context: Existing defect found by T022 (plan.md "Answers",
    2026-11-04). In `ConstantNormInterpolator.add_constant_norm`, the line
    `volume = self.support.element_size[element_indices]` fails with
    `IndexError: invalid index to scalar variable` for a `StructuredGrid`,
    because `StructuredGrid.element_size` returns one scalar
    (`np.prod(step_vector)`, all cells have the same size). Replace it with:
    `element_size = np.asarray(self.support.element_size, dtype=float)`;
    `volume = element_size[element_indices] if element_size.ndim > 0 else np.full(element_indices.shape[0], float(element_size))`.
    No other change. Do not change `loop_common`. If the `cn_fdi` cells then
    fail an SC-002 limit, stop and report the errors (plan OP3); do not
    relax a limit or change a solver setting.

- [x] T023 [P] [US2] [S] Use the world contract in `add_directional_regularisation`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`
  - Needs: T020
  - Wave: 10
  - Done when: `pytest tests/unit/modelling/test_fold_constraint_frames.py -k directional_regularisation_world_contract`
    and `pytest packages/loop_interpolation/tests/test_regularisation_api.py` pass.
  - Context: In `DiscreteInterpolator.add_directional_regularisation`, after
    `sample_points = np.asarray(self.get_regularisation_sample_points(), dtype=float)`
    add `world_samples = self.local_to_world_points(sample_points)`. Callable:
    `vectors = term.direction(world_samples)`; array: `vectors = term.direction`.
    Keep the shape check against `sample_points.shape`. Then
    `vectors = self.world_to_local_vectors(vectors)` before
    `self._add_directional_regularisation(term.weight, vectors, name=term.name)`.
    Add a numpy docstring: a callable direction receives world sample points
    and returns world directions; an array gives world directions, one row
    for each sample point. No signature change.

- [x] T024 [P] [US2] [S] Use the same world contract in the native adapter.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_native.py`
  - Needs: T020
  - Wave: 10
  - Done when: `pytest packages/loop_interpolation/tests/test_native_parity.py`
    passes (or reports its documented skip when the native extension is not
    installed; say which in the report), and `ruff check LoopStructural packages` passes.
  - Context: In `_add_regularisation(problem, interpolator, setup_kwargs)`
    (`_native.py:454`), directional part (about lines 470 to 485): after
    `sample_points = np.asarray(interpolator.get_regularisation_sample_points(), dtype=float)`
    add `world_samples = interpolator.local_to_world_points(sample_points)`.
    Callable: `vectors = term.direction(world_samples)`. Keep the shape check
    against `sample_points.shape`. Then
    `vectors = interpolator.world_to_local_vectors(vectors)` before
    `problem.add_directional_regularisation(_flat(vectors), float(term.weight))`.
    Do not change the rotated-support rejection (`_native.py:170`).

- [x] T025 [US2] [M] Sample fold geometry at world nodes in `FDFoldInterpolator`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_fd_fold_interpolator.py`
  - Needs: T019, T020, T023
  - Wave: 11
  - Done when: `pytest tests/unit/modelling/test_fold_constraint_frames.py -k fdi_fold`
    and `pytest packages/loop_interpolation/tests/test_fd_fold_interpolator.py` pass.
  - Context: In `add_fold_constraints` (no signature change): import
    `sample_fold_geometry`, `world_norm_constraints` from `._fold_setup`.
    `world_nodes, deformed_w, axis_w, dgz_w = sample_fold_geometry(self, self.fold, self.support.nodes)`;
    local `deformed`, `axis`, `dgz` = `self.world_to_local_vectors(...)` of
    each. `weight[mask_fn(world_nodes)] = 0.0`. The three
    `add_gradient_orthogonal_constraints` calls use the local vectors at the
    local `active_pos`. `resolve_fold_norm_target(..., normal_constraints=world_norm_constraints(self), ...)`.
    The directional regularisation gets the masked WORLD arrays
    `_masked_direction(dgz_w)`, `_masked_direction(deformed_w)`,
    `_masked_direction(axis_w)`; the base method converts them (T023).
    Docstring: `mask_fn` receives world node positions.

- [x] T026 [US2] [M] Sample fold geometry at world barycentres in `DiscreteFoldInterpolator`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_discrete_fold_interpolator.py`
  - Needs: T019, T020, T023
  - Wave: 11
  - Done when: `pytest tests/unit/modelling/test_fold_constraint_frames.py -k pli_fold`
    and `pytest packages/loop_interpolation/tests/test_discrete_fold_interpolator.py` pass.
  - Context: In `add_fold_constraints`: import `sample_fold_geometry`,
    `world_norm_constraints` from `._fold_setup`.
    `world_bc, deformed_w, axis_w, dgz_w = sample_fold_geometry(self, self.fold, self.support.barycentre)`;
    convert each direction with `self.world_to_local_vectors` and use the
    local vectors in the three `np.einsum("ij,ijk->ik", v, eg)` rows.
    `weight[mask_fn(world_bc)] = 0`. `normal_constraints=world_norm_constraints(self)`.
    `_masked_fold_direction` provider: it now receives WORLD points from the
    base method; call `self.fold.get_deformed_orientation(points)` with them,
    mask with `mask_fn(points)`, return world vectors (the base converts).
    Docstring: `mask_fn` receives world barycentres.

- [x] T027 [P] [US2] [M] Sample fold geometry at world quadrature points in `P2FoldInterpolator`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_p2_fold_interpolator.py`
  - Needs: T019, T020
  - Wave: 10
  - Done when: `pytest tests/unit/modelling/test_fold_constraint_frames.py -k p2_fold`
    and `pytest packages/loop_interpolation/tests/test_p2_fold_interpolator.py` pass.
  - Context: In `add_fold_constraints`, after
    `points, tet_idx, q = tet_quadrature_points(...)`:
    `world_points, deformed_w, axis_w, dgz_w = sample_fold_geometry(self, self.fold, points)`
    (import from `._fold_setup`, with `world_norm_constraints`); convert the
    three directions with `self.world_to_local_vectors`; use the local
    vectors in `_add`. `dN` and `inside` stay on the local `points`.
    `sq = sq * (~np.asarray(mask_fn(world_points), dtype=bool))`.
    `normal_constraints=world_norm_constraints(self)`. Docstring: `mask_fn`
    receives world sample points.

- [x] T028 [P] [US2] [S] Give `steepness_wtfunc` the world barycentres of the weighted elements.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_p2interpolator.py`
  - Needs: T020
  - Wave: 10
  - Done when: `pytest tests/unit/modelling/test_fold_constraint_frames.py -k steepness_wtfunc`
    passes.
  - Context: In `P2Interpolator.minimise_grad_steepness(self, w=0.1, maskall=False, wtfunc=None)`,
    when `callable(wtfunc)`, replace
    `wt = wtfunc(self.support.barycentre) * self.support.element_size[mask]` with
    `wt = wtfunc(self.local_to_world_points(self.support.barycentre[elements[mask]])) * self.support.element_size[elements[mask]]`.
    Docstring of `wtfunc`: "receives the world barycentres of the weighted
    elements and returns one weight for each". No signature change.

- [x] T029 [P] [US2] [S] Document the world contract of `resolve_fold_norm_target`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_fold_norm_alignment.py`
  - Needs: T019
  - Wave: 9
  - Done when: `ruff check LoopStructural packages` passes and
    `pytest packages/loop_interpolation/tests -k fold_norm` passes.
  - Context: Docstring only, no code change. In the `Parameters` section of
    `resolve_fold_norm_target`, the `normal_constraints` entry becomes:
    "Normal constraints in world coordinates, ``xyz|nx ny nz|w`` format. Use
    ``loop_interpolation._fold_setup.world_norm_constraints`` to convert the
    stored local constraints." The `fold` entry adds: "Its
    ``get_deformed_orientation`` takes world points and returns world unit
    directions."

## Phase 8: Stage S2 finish

- [x] T030 [P] [Finish] [S] Add the S2 migration notice to `COMPAT.md`.
  - Files: `COMPAT.md`
  - Needs: T023, T024, T025, T026, T027, T028, T032
  - Wave: 12
  - Done when: the "Migration notices" table in `COMPAT.md` has the new row
    with the S2 date.
  - Context: Add one row at the end of the table under
    `## Migration notices (not breaking, no shim needed)`. Columns: Old |
    New | Since. Old: "Fold `mask_fn`, `directional_regularisation`
    direction callables and arrays, and the P2 `steepness_wtfunc` received
    or returned local (computational-frame) points and directions". New:
    "They receive world points and give world directions; the interpolator
    converts to its local frame once. With an identity frame nothing
    changes (spec 005, stage S2)". Since: the date of the task run
    (YYYY-MM-DD).

- [x] T031 [Finish] [M] Record the S2 evidence, backend cells and gates in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T020, T021, T022, T023, T024, T025, T026, T027, T028, T029, T030, T032
  - Wave: 13
  - Done when: `ruff check LoopStructural packages && pytest` and
    `pytest packages/loop_common/tests packages/loop_interpolation/tests`
    pass, and each heading of `verification.md` has an S2 entry.
  - Context: Add an "S2" subsection under each heading. Pre-fix (FR-031):
    make a throw-away worktree at the S2 start commit (the `start_commit`
    of stage S2 in `run-state.json`), copy the S2 test files and
    `tests/fixtures/coordinate_frames.py` into it, run them with its own
    `PYTHONPATH`, record the failing ids per test, then
    `git worktree remove --force <tmp>` and `git worktree prune`. Tests that
    pass there (P2 reference, constant norm) are "coverage, no defect".
    Backend matrix: add P2, FDI-fold, PLI discrete fold, P2-fold and the
    constant-norm cells; native: run or record the absent extension. Record
    the S2 test time (OP2) and the rollback check: the S1 tests pass at the
    S2 start commit. Golden review: list any changed existing test, or "No
    golden change in S2."

## Phase 9: Stage S3 foundation

- [x] T033 [P] [Found] [M] Add the S3 persistence and lifecycle helpers to the shared fixtures.
  - Files: `tests/fixtures/coordinate_frames.py`
  - Needs: T031
  - Wave: 14
  - Done when: `pytest tests/fixtures/coordinate_frames.py` collects 0 tests with no
    import error, and `ruff check LoopStructural packages` passes. A short
    `python -c` script (worktree `PYTHONPATH`) calls `r1_interpolator` for the 5
    engines at `ANCHORS[0]`, rotation 0, with `exact=True` for the first 3 and
    `exact=False` for the last 2, and exits 0.
  - Context: Helper module, no tests. numpy docstrings, type hints, Black line 100.
    Add `PERSIST_ENGINES = ("fdi", "pli", "p2")`, `LIFECYCLE_ENGINES = PERSIST_ENGINES + ("cn_fdi", "cn_pli")`,
    `QUERY_OFFSETS = np.array([[2.2, 3.7, 4.3], [6.1, 8.2, 9.1]])`.
    Add `r1_interpolator(engine, anchor, rotation_deg, exact=True) -> tuple[interpolator, pd.DataFrame]`.
    Steps: `model = model_with_frame(anchor, anchor + (8, 10, 12), rotation_deg)`; `bb = model.bounding_box`;
    corners of the world box through `bb.project`; local `lo`, `hi`; `cells = (5, 6, 7)` if `exact`
    else `(6, 6, 6)`; `step = (hi - lo) / cells`. Supports (copy the S2 test files): `fdi` =
    `FiniteDifferenceInterpolator(StructuredGrid(origin=lo, nsteps_cells=cells, step_vector=step))`; `pli` =
    `P1Interpolator(TetMesh(same args))`; `p2` = `P2Interpolator(P2TetMesh(origin=lo, nsteps=cells, step_vector=step))`;
    `cn_fdi` = `ConstantNormFDIInterpolator(StructuredGrid(...))`; `cn_pli` = `ConstantNormP1Interpolator(TetMesh(...))`
    (cn engines: `norm_length = np.sqrt(0.29)`). Then `interp.bounding_box = bb`; `obs = r1_observations(anchor)`;
    `interp.set_value_constraints(obs[["X", "Y", "Z", "val", "w"]].to_numpy())`. If `exact`: set
    `interp.c = r1_value(interp.local_to_world_points(interp.support.nodes), anchor)` and `interp.up_to_date = True`.
    Else: `interp.setup()` and `assert interp.solve_system(solver="cg", solver_kwargs={"rtol": 1e-10, "atol": 1e-12, "maxiter": 20000})`.
    Add `strip_frame_keys(payload: dict) -> dict`: a deep copy without the keys `format_version`,
    `constraint_space`, `transform` and `bounding_box` (a legacy payload).

- [x] T034 [P] [Found] [S] Add the S3 module stubs and export the two errors.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_frame_persistence.py`, `packages/loop_interpolation/src/loop_interpolation/__init__.py`
  - Needs: T031
  - Wave: 14
  - Done when: `python -c "from loop_interpolation import InterpolatorFormatError, StaleFrameError; from loop_interpolation._frame_persistence import FORMAT_VERSION, CONSTRAINT_SPACE, REASONS, frame_payload, payload_kind, resolve_frame, bounding_box_from_support"`
    runs with the worktree `PYTHONPATH`, and `ruff check LoopStructural packages` passes.
  - Context: New module. Start with `from __future__ import annotations`, numpy docstrings, type hints.
    `FORMAT_VERSION: int = 2`; `CONSTRAINT_SPACE: str = "local"`;
    `REASONS: tuple[str, ...] = ("legacy_no_transform", "unsupported_version", "malformed", "conflict")`.
    `class InterpolatorFormatError(ValueError)`: `__init__(self, reason: str, message: str) -> None` calls
    `super().__init__(message)` and sets `self.reason = reason`; raise `ValueError` if `reason not in REASONS`.
    `class StaleFrameError(ValueError)`: no body except a docstring.
    Stubs that raise `NotImplementedError("S3 stub")`: `frame_payload(bounding_box, dimensions) -> dict`;
    `payload_kind(payload: dict, caller_frame: bool) -> str`;
    `resolve_frame(payload: dict, kind: str, support: object, dimensions: int, *, local_origin=None, local_rotation=None, declare_identity: bool = False)`;
    `bounding_box_from_support(support: object, dimensions: int, local_origin, local_rotation)`.
    In `__init__.py`: add `from ._frame_persistence import InterpolatorFormatError, StaleFrameError` after
    the `._diagnostics` import, and add both names to `__all__` in alphabetical order.

## Phase 10: User story 3 - save, copy and keep one frame (S3, FR-016 to FR-020)

**Goal**: a saved interpolator keeps its transform and the space of its constraints. A copy keeps the
frame. A legacy file needs an explicit choice. A frame change never gives a stale result.
**Independent test**: `pytest tests/unit/modelling/test_frame_persistence.py tests/unit/modelling/test_frame_copy_and_lifecycle.py`.

- [x] T035 [P] [US3] [M] Test: payload frame keys, round trips and two save/load cycles (U3.1, U3.2, U3.8; SC-006).
  - Files: `tests/unit/modelling/test_frame_persistence.py`
  - Needs: T033, T034
  - Wave: 15
  - Done when: `pytest tests/unit/modelling/test_frame_persistence.py -k "payload_declares or round_trip_world or two_cycles or json_and_yaml"`
    fails before T043 (record the failing ids and the first error of each) and passes after T043.
  - Context: New file. Import `GeologicalInterpolator` and `InterpolatorFactory` from `loop_interpolation`
    and the helpers from `tests.fixtures.coordinate_frames` (T033). Cells: `engine in PERSIST_ENGINES`
    x `ANCHORS` x `ROTATIONS_DEG` = 36 per test (use `r1_interpolator(engine, anchor, rot, exact=True)`).
    `test_payload_declares_frame`: `payload = interp.to_dict()`; `payload["format_version"] == 2`;
    `payload["constraint_space"] == "local"`; `payload["transform"]` has `dimensions == 3`, `local_origin`
    and `local_rotation` equal to `interp.bounding_box` values; `payload["bounding_box"]` has `origin`,
    `maximum`, `nsteps` and no `local_*` key; `format_version` is the first of the new keys.
    `test_round_trip_world_values`: `loaded = GeologicalInterpolator.from_dict(payload)`; at `anchor + QUERY_OFFSETS`
    values differ by at most `1e-12 * F` (F = max(1, max abs of `r1_value`)) and gradients by at most
    `1e-12 * max(1, abs(R1_GRADIENT).max())`.
    `test_local_arrays_unchanged_after_two_cycles`: save, load, save, load; `loaded.data[k]` equals
    `interp.data[k]` bit for bit for each key; both payloads equal; same version, transform and space.
    `test_json_and_yaml_carry_frame`: `engine in PERSIST_ENGINES`, `ANCHORS[3]`, 30 degrees; use `to_json`/`from_json`
    and `to_yaml`/`from_yaml`; the loaded transform equals the source.

- [x] T036 [US3] [M] Test: legacy, conflicting, malformed and unsupported payloads, and the creation call (U3.6, U3.7; SC-012, SC-017).
  - Files: `tests/unit/modelling/test_frame_persistence.py`
  - Needs: T035
  - Wave: 16
  - Done when: `pytest tests/unit/modelling/test_frame_persistence.py -k "legacy or conflicting or malformed or unsupported or creation_call or old_reader"`
    fails before T043 (record the failing ids) and passes after T043. A cell that passes before T043
    (the creation call) is recorded as "coverage, no defect".
  - Context: Add to the file of T035. Import `InterpolatorFormatError` from `loop_interpolation`.
    `legacy = strip_frame_keys(payload)`. Use `copy.deepcopy` to prove that the source dict is unchanged.
    `test_legacy_payload_requires_choice` (3 engines): `pytest.raises(InterpolatorFormatError)`, `reason == "legacy_no_transform"`,
    message starts with `legacy interpolator payload has no transform` and names `local_origin` and `declare_identity`.
    `test_legacy_with_supplied_transform` (3 engines x 4 anchors x 3 rotations): load with `local_origin=bb.local_origin`,
    `local_rotation=bb.local_rotation`; values and gradients meet SC-006 (`1e-12` limits) against the source.
    `test_legacy_with_identity_declaration` (3 engines, anchor 0, rotation 0): `declare_identity=True`.
    `test_conflicting_inputs` (4 cases): `declare_identity` with `local_origin`; v2 payload with another origin; with another rotation;
    a rotated v2 payload with `declare_identity=True`. `reason == "conflict"`, message starts with `conflicting transform`.
    `test_malformed_metadata` (8 cases): missing `transform`; `constraint_space` "world"; rotation scaled by 2; NaN origin; origin of length 2;
    `bounding_box` None with a rotated transform; legacy payload that keeps `transform`; `format_version` equal to "2". `reason == "malformed"`.
    `test_unsupported_future_version` (3 and 99): `reason == "unsupported_version"`. `test_creation_call_unchanged`:
    `InterpolatorFactory.from_dict({"type": "FDI", "boundingbox": bb, "nelements": 1000})` builds an empty interpolator; the same call
    with `declare_identity=True` raises `malformed`. `test_old_reader_would_name_format_version`: pop `type`, `support`, `data`, `c`,
    `up_to_date`, `valid` from a payload copy; `next(iter(rest)) == "format_version"`.

- [x] T037 [P] [US3] [M] Test: copies keep the frame, a frame change is rejected without state change, and the rebuild recipe (U3.3, U3.9; SC-016).
  - Files: `tests/unit/modelling/test_frame_copy_and_lifecycle.py`
  - Needs: T033, T034
  - Wave: 15
  - Done when: `pytest tests/unit/modelling/test_frame_copy_and_lifecycle.py -k "copy_keeps_frame or stale_frame or restored_transform or rebuild_in_new_frame"`
    fails before T042 to T046 (record the failing ids) and passes after them. The rebuild cells may pass
    before: record them as "coverage, no defect".
  - Context: New file. First run `git grep -n 'interpolator.*copy(' -- LoopStructural/modelling/features` and put any hit in your
    report (OP6). `test_copy_keeps_frame`: `engine in ("fdi", "p2", "p2_fold")` x 12 frames. Build `fdi` and `p2` with
    `r1_interpolator(engine, anchor, rot)`. Build `p2_fold` as `P2FoldInterpolator(P2TetMesh(...), fold=AnalyticFold.constant(d, a, g))`
    (unit vectors) with the same box, support and `bounding_box` recipe as `r1_interpolator`. `c = interp.copy()`:
    `c.bounding_box is interp.bounding_box`; the four conversion methods agree with the source on a point array; the
    copy has 0 stored rows; same support; for `p2_fold`, `c.fold is interp.fold`.
    `test_stale_frame_rejected_without_state_change`: engine in `LIFECYCLE_ENGINES`, `ANCHORS[1]`, rotation 0, `exact=False`. Take a snapshot
    (deep copies of `data`, `c`, `up_to_date`, `valid`, `n_g`, `n_i`, `n_n`, `n_t`; the object `latest_diagnostics_report`). Then
    `bb.set_local_transform(local_origin=anchor + (3, -2, 1), rotation_matrix=rot_z(30.0))`. Each of `evaluate_value`,
    `evaluate_gradient`, `set_value_constraints` (27 rows) and `to_dict` raises `StaleFrameError` (message starts with `the local transform changed`)
    once. The snapshot is unchanged. `test_restored_transform_is_usable` (`fdi`): set the old transform back; `evaluate_value` equals the first result.
    `test_rebuild_in_new_frame_preserves_world_evaluation` (5 engines): the documented recipe, `r1_interpolator(engine, anchor, 30.0, exact=False)`
    (a new interpolator in the new frame, `ANCHORS[1]`); values and gradients at `anchor + QUERY_OFFSETS` meet SC-002
    (`1e-8 * F`; `1e-8 * max(1, abs(R1_GRADIENT).max())`). Exception (OP7, decided 2026-11-05): the `p2` cell uses `exact=True` because the P2 fit misses SC-002 in every frame; do not change the limit.

- [x] T038 [P] [US3] [M] Test: fault restoration after a transform change (SC-016, cache part).
  - Files: `tests/unit/modelling/test_frame_copy_and_lifecycle.py`
  - Needs: T037
  - Wave: 16
  - Done when: `pytest tests/unit/modelling/test_frame_copy_and_lifecycle.py -k fault_restoration`
    fails before T042 (record the first error) and passes after T042.
  - Context: Add to the file of T037. First read `FaultSegment.apply_to_points` in
    `LoopStructural/modelling/features/fault/_fault_segment.py` (read only) and write in your report which public
    interpolator methods it calls. `test_fault_restoration_uses_current_frame`: `engine in ("FDI", "PLI")` (the R2 fault of
    `tests/unit/modelling/test_fault_frame_analytic.py`; copy its model set-up, anchor `ANCHORS[1]`, rotation 0). Restore
    the 6859-probe subset `r2_probes(anchor)[:200]` through the fault frame. Change the local transform of `model.bounding_box`
    as in T037: restoring the points raises `StaleFrameError`. Then build the same fault in a fresh model with the new frame
    (`model_with_frame(...)`, rotation 30) and with a second fresh model (uncached control); the restored world points of the
    two agree within `tp(L, W)`. No restoration cache exists, so the test also records that fact in its docstring.

- [x] T039 [US3] [M] Record the S3 pre-fix evidence in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T036, T038
  - Wave: 17
  - Done when: `verification.md` has the heading `## S3 pre-fix failures` with a table of the failing test ids and
    their first error for the 2 new test files, and the 3 behavioural records below, each with a measured value.
  - Context: Run the 2 new test files on the current tree with the worktree `PYTHONPATH` (check that `loop_interpolation.__file__`
    comes from the worktree). Record per test function: failed / total, one representative id and value. Mark a function whose cells
    all fail only on a missing name ("missing name"). Also record 3 behavioural defects from the production path with a short script: (1) a rotated
    `fdi` round trip (`ANCHORS[3]`, 30 degrees, `exact=True`): `from_dict(to_dict())` gives world values with maximum error E (give E);
    (2) `interp.copy().bounding_box is None` for `fdi`; (3) after `bb.set_local_transform(...)` (T037 change) `evaluate_value`
    returns a value with no error and differs from the earlier value by D (give D). Add `## S3 pre-fix failures` after the S2 sections.
    Write in ASD-STE100 Simplified Technical English.

## Phase 11: Stage S3 implementation

- [x] T040 [US3] [M] Implement `frame_payload`, `payload_kind` and `bounding_box_from_support`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_frame_persistence.py`
  - Needs: T039
  - Wave: 18
  - Done when: a `python -c` script (worktree `PYTHONPATH`) checks that `frame_payload(None, 3)` has the key order
    `format_version`, `constraint_space`, `transform`, `bounding_box`, with an identity transform and `bounding_box` None;
    that `payload_kind({}, False)` raises `InterpolatorFormatError` with reason `legacy_no_transform`; and that
    `ruff check LoopStructural packages` passes.
  - Context: Read the section "S3 functions and signatures" of `specs/005-coordinate-frame-correctness/plan.md`. Replace the stubs of
    these 3 functions only. `frame_payload(bounding_box, dimensions)`: `{"format_version": FORMAT_VERSION, "constraint_space": CONSTRAINT_SPACE,
    "transform": {"dimensions", "local_origin" (list), "local_rotation" (list of lists)}, "bounding_box": dict or None}`; for a box,
    `bounding_box` = `BoundingBox.to_dict()` without `local_origin` and `local_rotation`; for `None`, identity of `dimensions`.
    `payload_kind`: the 6 numbered steps of the plan; messages exactly as in the plan. `bounding_box_from_support`: local min and max of
    `support.nodes`; the `2**dimensions` corners; `world = corners @ rotation + origin`; `BoundingBox(dimensions=dimensions, origin=world.min(0),
    maximum=world.max(0))` then `set_local_transform(local_origin=origin, rotation_matrix=rotation)`. Import `BoundingBox` from
    `loop_common.geometry`. Rule: py39 syntax with `from __future__ import annotations`; numpy docstrings.
    Also fix two docstrings from T034: the module docstring must not mention stubs; `StaleFrameError` is raised when the local transform changed after constraints were stored.

- [x] T041 [P] [US3] [M] Implement `resolve_frame`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_frame_persistence.py`
  - Needs: T040
  - Wave: 19
  - Done when: a `python -c` script (worktree `PYTHONPATH`) calls `resolve_frame` on a payload from `frame_payload` of a rotated box and
    gets a box with the same `local_origin` and `local_rotation`; a call with `declare_identity=True` and `local_origin=[1, 0, 0]` raises
    `InterpolatorFormatError` with reason `conflict`; and `ruff check LoopStructural packages` passes.
  - Context: Read the section "S3 functions and signatures" of `specs/005-coordinate-frame-correctness/plan.md`. Replace the stub of `resolve_frame` only.
    Follow its 6 numbered steps and the message list. `caller_frame` = `local_origin is not None or local_rotation is not None or declare_identity`.
    Use `validate_rigid_rotation` from `loop_common.geometry`. Wrap `RigidRotationError`, `LoopValueError`, `KeyError`, `TypeError` and
    `ValueError` from parsing as `InterpolatorFormatError("malformed", ...)` with `raise ... from exc`. Use `BoundingBox.from_dict` for the saved box.
    The function reads `payload` only; it never changes it or a caller object.
    Before it calls `bounding_box_from_support`, check that `support.nodes` has shape `(n, dimensions)` and that the origin has shape `(dimensions,)`; else raise `malformed`.
    A conflict message after `conflicting transform: ` contains `declare_identity` (and `saved` against a saved transform) if `declare_identity` is involved; else it contains `saved`, `supplied` and the item name (`local_origin` or `local_rotation`).

- [x] T042 [P] [US3] [M] Add the frame check, local restore and frame keys to `GeologicalInterpolator`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_geological_interpolator.py`
  - Needs: T040
  - Wave: 19
  - Done when: `pytest tests/unit/modelling/test_frame_copy_and_lifecycle.py -k "stale_frame or restored_transform or rebuild_in_new_frame or fault_restoration"`
    passes, and `pytest packages/loop_interpolation/tests/test_geological_interpolator.py` passes with no change, except
    `test_interpolator_json_yaml_round_trip`. That test loads a payload that `to_dict` now writes with frame keys, and the old loader
    rejects them until T043. Record it as "fails until T043" in your report.
  - Context: Read the section "S3 functions and signatures" of `specs/005-coordinate-frame-correctness/plan.md` (class `GeologicalInterpolator`). Add:
    `self._frame_at_data = None` in `__init__`; `clean()` sets it to `None`; `frame_signature()`; `_note_frame(rows)`; `_check_frame()`;
    `_inherit_frame(other)`; `_restore_local_constraints(data)`. Call `self._check_frame()` first in each of the 7 `set_*_constraints`
    (inside the `try`) and call `self._note_frame(points.shape[0])` after the array is stored; also first in `evaluate_value`, `evaluate_gradient`
    and `to_dict`. `to_dict` ends with `payload.update(frame_payload(self.bounding_box, self.dimensions))`. Add the keyword arguments
    `local_origin=None, local_rotation=None, declare_identity=False` to `from_dict`, `from_json`, `from_yaml`; pass them to
    `InterpolatorFactory.from_dict` only when given (the existing test uses a fake that takes one argument). Import `StaleFrameError` and
    `frame_payload` from `._frame_persistence`. `StaleFrameError` is not wrapped in `ValidationError`.
    Rule: numpy docstrings, type hints, py39.

## Phase 12: Stage S3 loader and copies

- [x] T043 [US3] [M] Rewrite `InterpolatorFactory.from_dict` for saved payloads.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_interpolator_factory.py`
  - Needs: T041, T042
  - Wave: 20
  - Done when: `pytest tests/unit/modelling/test_frame_persistence.py tests/unit/interpolator/test_interpolator_factory.py`
    and, in a separate call (the conftest plugins clash), `pytest packages/loop_interpolation/tests/test_p2_fold_registration.py packages/loop_interpolation/tests/test_geological_interpolator.py`
    pass, and `ruff check LoopStructural packages` passes.
  - Context: Read the section "S3 functions and signatures" of `specs/005-coordinate-frame-correctness/plan.md` (function `InterpolatorFactory.from_dict`). New
    signature: `from_dict(d, *, local_origin=None, local_rotation=None, declare_identity=False)`. A dict is a saved payload if and only if
    `d.get("support")` is a dict. Other dicts keep the existing code path with no change; frame arguments with such a dict raise
    `InterpolatorFormatError("malformed", "frame arguments apply only to a saved payload")`. For a saved payload follow the 8 steps of the plan:
    `payload_kind` first (before any work), pop `format_version`, `constraint_space`, `transform`, `bounding_box` into one dict, build the support and the
    interpolator as now, set `interpolator.bounding_box = resolve_frame(...)`, call `_restore_local_constraints(data)` (never `set_*_constraints`), then `c`,
    `up_to_date`, `valid`. Do not change the input dict. Keep the existing `rotation_xy` `TypeError` retry.
    Wrap an error of `_restore_local_constraints` (`ValueError`, `ValidationError`) as `InterpolatorFormatError("malformed", "malformed interpolator payload: ...")` with `raise ... from exc`; the half-built interpolator is dropped.

- [x] T044 [P] [US3] [S] Keep the frame in `FiniteDifferenceInterpolator.copy`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_finite_difference_interpolator.py`
  - Needs: T042
  - Wave: 20
  - Done when: `pytest tests/unit/modelling/test_frame_copy_and_lifecycle.py -k "copy_keeps_frame and fdi"` passes.
  - Context: In `FiniteDifferenceInterpolator.copy` (about line 247) replace the body with
    `return self._inherit_frame(FiniteDifferenceInterpolator(self.support))`. Change the docstring to: "Create an empty copy on the same
    support, with the same coordinate frame." and "Returns" to "A new empty interpolator that shares the support and the bounding box.".
    Change nothing else in the file.

- [x] T045 [P] [US3] [S] Keep the frame in `P2Interpolator.copy`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_p2interpolator.py`
  - Needs: T042
  - Wave: 20
  - Done when: `pytest tests/unit/modelling/test_frame_copy_and_lifecycle.py -k "copy_keeps_frame and p2 and not p2_fold"` passes.
  - Context: In `P2Interpolator.copy` (about line 106) replace the body with `return self._inherit_frame(P2Interpolator(self.support))`.
    Add this docstring: "Create an empty copy on the same support, with the same coordinate frame." Change nothing else in the file.

- [x] T046 [P] [US3] [S] Keep the frame in `P2FoldInterpolator.copy`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_p2_fold_interpolator.py`
  - Needs: T042
  - Wave: 20
  - Done when: `pytest tests/unit/modelling/test_frame_copy_and_lifecycle.py -k "copy_keeps_frame and p2_fold"` passes.
  - Context: In `P2FoldInterpolator.copy` (about line 64) replace the body with
    `return self._inherit_frame(P2FoldInterpolator(self.support, fold=self.fold))`. Add this docstring: "Create an empty copy on the same
    support, with the same fold and the same coordinate frame." Change nothing else in the file.

## Phase 13: Stage S3 finish

- [x] T047 [P] [Finish] [S] Add the S3 migration notice to `COMPAT.md`.
  - Files: `COMPAT.md`
  - Needs: T043
  - Wave: 21
  - Done when: `COMPAT.md` has one new row at the end of the "Migration notices" table, and `grep -c "spec 005, stage S3" COMPAT.md` prints 1.
  - Context: Copy the 3-column layout of the last row (the S2 row). Column 1: "Loading a legacy interpolator payload (no `format_version`) assumed an
    identity frame without a message; `copy()` returned an interpolator with no frame; a transform change after constraints existed gave a stale result".
    Column 2: "A legacy payload needs `local_origin` and `local_rotation`, or `declare_identity=True`, and raises `InterpolatorFormatError` otherwise. New payloads carry
    `format_version: 2`, `constraint_space`, `transform` and `bounding_box`. `copy()` keeps the frame. A transform change after constraints exist raises `StaleFrameError`
    at the next set, evaluate or save (spec 005, stage S3). Major-version note: the release that contains S3 must be a major version, or the maintainer records another
    decision in `ROADMAP.md`". Column 3: the date of the run in the form YYYY-MM-DD. Write in ASD-STE100 Simplified Technical English.

- [x] T048 [P] [Finish] [S] Record the FR-019 decision in `ROADMAP.md`.
  - Files: `ROADMAP.md`
  - Needs: T043
  - Wave: 21
  - Done when: `grep -c "Legacy interpolator payloads (spec 005, FR-019)" ROADMAP.md` prints 1.
  - Context: Add a new `###` subsection as the last one of `## Decisions`, directly before `## Stage sequence`. Heading: `### Legacy interpolator payloads (spec 005, FR-019)`.
    Body (ASD-STE100): the decision (2026-11-03, C-003); a legacy interpolator payload has no recoverable transform; the loader needs a caller-supplied
    transform (`local_origin`, `local_rotation`) or `declare_identity=True` and raises `InterpolatorFormatError` otherwise; no warning-only identity fallback
    exists; new payloads carry `format_version: 2`; this is the only approved break of the public loading behaviour in spec 005 (FR-028); the release that contains S3
    must be a major version, or the maintainer records another decision here; the policy for a transform change is "reject at use" (`StaleFrameError`, decided 2026-11-05).
    Keep the style of the neighbouring subsections. Do not change other text.

- [x] T049 [Finish] [M] Record the S3 rollback check in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T043, T044, T045, T046
  - Wave: 21
  - Done when: `verification.md` has the heading `## S3 rollback checks` with the results of the 2 checks below, and
    `git worktree list` shows no throw-away worktree.
  - Context: Use a throw-away worktree at the S3 start commit (the `start_commit` of stage S3 in `run-state.json`; if the key is missing, use the parent of the first S3 test commit):
    `git worktree add --detach ../LoopStructural-005-s3-rollback <commit>`. With its own `PYTHONPATH` (check `loop_interpolation.__file__`): (1) save `r1_interpolator`
    payloads for `fdi`, `pli`, `p2` (anchor 3, rotation 30, `exact=True`) with the old `to_dict`; write them to a scratch folder. On the S3 worktree, load each with the
    supplied transform and with `declare_identity=True` (anchor 0, rotation 0 only) and record the SC-006 errors. (2) On the old reader, load a payload that the S3 code wrote; record the exact
    `TypeError` text (it names `format_version`) and that no interpolator object exists (OP5). Then `git worktree remove --force <tmp>` and `git worktree prune`.

- [x] T050 [Finish] [M] Run the S3 gates and complete `verification.md` and the `ROADMAP.md` status log.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`, `ROADMAP.md`
  - Needs: T047, T048, T049
  - Wave: 22
  - Done when: `ruff check LoopStructural packages && pytest` and `pytest packages/loop_common/tests packages/loop_interpolation/tests` pass, `verification.md`
    has the S3 sections named below, and `ROADMAP.md` has one new status-log line for S3.
  - Context: Run the root `pytest` and the package tests in separate calls (the conftest plugins clash) with the worktree `PYTHONPATH`. Add to `verification.md`: `## S3 post-fix results`
    (commands, counts, times, the test-time of the 2 new files), `## S3 golden review (C-004)` (list each existing test file that changed, or "No golden change in S3."; record
    any existing test that loads a payload with no `format_version`), `## S3 backend matrix` (engines `fdi`, `pli`, `p2`, `cn_fdi`, `cn_pli`, `p2_fold` copy; PLI has no `copy`;
    native and Surfe are S6), `## S3 test run time (OP2)` and a line for OP6 (the grep result of T037). Measure the effect of the frame check: time `evaluate_value` on 1e5 points for
    `fdi` before and after the check (a loop of 1000 calls on 100 points is enough) and record it; if the 003 worktree has a key-time check, record its result, else write "no 003 measure". If a library
    path now raises `StaleFrameError` in `pytest`, report it; do not weaken the check. Add to the status log of `ROADMAP.md` a line dated today: "Spec 005 stage S3 implemented: saved interpolator frame,
    legacy-load rule, frame-keeping copies, stale-frame check." Write in ASD-STE100 Simplified Technical English.

## Phase 14: Stage S4 foundation

- [x] T051 [Found] [M] Add the S4 helper `world_foliation_model` to the shared fixtures.
  - Files: `tests/fixtures/coordinate_frames.py`
  - Needs: T050
  - Wave: 23
  - Done when: `python -c "from tests.fixtures.coordinate_frames import world_foliation_model"` prints no error, and `pytest tests/unit/modelling/test_frame_persistence.py -k payload_declares` still passes.
  - Context: Add `world_foliation_model(anchor, rotation_deg) -> tuple[GeologicalModel, GeologicalFeature, pd.DataFrame]`. Build the model with the existing
    `model_with_frame(anchor, anchor + (8, 10, 12), rotation_deg)`. The data frame is in WORLD coordinates with the columns
    `X, Y, Z, feature_name, val, nx, ny, nz, gx, gy, gz, w`: 27 value rows from `r1_observations(anchor)`; 4 gradient rows and 4 normal rows at
    `anchor + (2, 3, 4)`, `(5, 6, 7)`, `(3, 8, 5)` and `(6, 4, 9)` with the unit direction `R1_GRADIENT / norm(R1_GRADIENT)`.
    `feature_name` is `"r1"` for every row. Set `model.data` to the frame and make the feature with `model.create_and_add_foliation("r1", ...)`. Do NOT solve
    (no call that runs the interpolation). Copy the call pattern of `create_and_add_foliation` from `tests/integration/test_interpolator.py` and
    `tests/unit/interpolator/test_outside_box.py` (read both). Return the model, the feature and the data frame. Do not change existing helpers.
    Rule: input data is world data; the interpolator converts it once.

## Phase 15: User story 3 - display, export and classification use world locations (S4, FR-022)

**Goal**: `get_data`, the data export and the intrusion classification use the stored world locations, with no second conversion. `rescale` and `scale` warn.
**Independent test**: `pytest tests/unit/modelling/test_world_location_consumers.py`.

- [x] T052 [US3] [M] Test: display and export locations are world (U3.4; FR-022).
  - Files: `tests/unit/modelling/test_world_location_consumers.py`
  - Needs: T051
  - Wave: 24
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k "feature_data_locations or example_observation or exported_data"` fails before T057 (record the failing ids and the first error of each) and passes after T057.
  - Context: New file. Use `world_foliation_model` (T051), `ANCHORS`, `ROTATIONS_DEG` and `tp(L, W)` from `tests.fixtures.coordinate_frames` (read the file for the exact names).
    `test_feature_data_locations_are_world`: 4 anchors x 3 rotations = 12 cells. `feature.get_data()` gives a `ValuePoints` (value) and `VectorPoints` (gradient and norm).
    Each `locations` equals the input rows within `tp(L, W)`; each `vectors` equals the input vectors within `1e-12`.
    `test_spec_example_observation_is_not_shifted`: `ANCHORS[1]` x 3 rotations = 3 cells. The observation at `anchor + (2, 3, 4)` is at that point, and not at the point shifted by the local origin
    (for origin `(100,200,300)` and identity rotation: `(102,203,304)` not `(202,403,604)`).
    `test_exported_data_keeps_world_locations`: 4 anchors, rotation 30 = 4 cells. Save each object of `get_data()` with its `save` method to a CSV file in `tmp_path`
    (read `LoopStructural/datatypes/_point.py` `save` first for the extension it supports; read only). The `x`, `y`, `z` columns equal the input within `tp(L, W)`.
    The test file must not call `model.rescale`. Pre-fix expectation: the three tests fail by a location shifted by the local origin.

- [x] T053 [US3] [M] Test: intrusion classification evaluates the model at world points (U3.5; FR-022).
  - Files: `tests/unit/modelling/test_world_location_consumers.py`
  - Needs: T052
  - Wave: 25
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k intrusion_classification` fails before T058 with `AttributeError` (the method does not exist yet), then passes after T058.
  - Context: Append to the file of T052. `test_intrusion_classification_uses_world_points`: 4 anchors x 3 rotations = 12 cells. Make a real rotated `GeologicalModel` (use `model_with_frame`).
    Make the builder with `object.__new__(IntrusionFrameBuilder)` (import from `LoopStructural.modelling.intrusions.intrusion_frame_builder`) and set `builder.model = model`.
    Replace `model.evaluate_model` with a spy that stores its argument and returns `np.zeros(len(xyz))`. Call `builder._network_model_values(data)` on a data frame with the world
    columns `X, Y, Z` (at least 5 rows, taken from the anchor). The spy gets exactly the world X, Y, Z: `np.array_equal(spy_arg, data[["X","Y","Z"]].to_numpy(float))`.
    The signature of the method under test is `_network_model_values(self, data: pd.DataFrame) -> np.ndarray`.

- [x] T054 [US3] [M] Test: the 2D viewer draws the fault ellipse at the world centre (U3.5, coverage; no defect).
  - Files: `tests/unit/modelling/test_world_location_consumers.py`
  - Needs: T053
  - Wave: 26
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k fault_ellipse` passes before and after the fixes (record it as "coverage, no defect").
  - Context: Append to the file of T053. `test_2d_viewer_fault_ellipse_centre_is_world`: 2 anchors x 2 rotations (0 and 30) = 4 cells. Call `Loop2DView.add_fault_ellipse`
    (`packages/loopstructural_visualisation/src/loopstructuralvisualisation/_2d_viewer.py`, lines 211 to 224; read it) as an unbound function with a stand-in `self` that has
    `ax` (a Matplotlib axis from the `Agg` backend) and `model`. The patch centre equals `fault.fault_centre[:2]` within `tp(L, W)`.
    Build the fault with `_fault(anchor, rotation, "PLI")` from `tests/unit/modelling/test_fault_frame_analytic.py` (import it; read the file for the signature).
    Use `matplotlib.use("Agg")`. Do not change the viewer.

- [x] T055 [US3] [M] Record the S4 pre-fix evidence in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T054
  - Wave: 27
  - Done when: `verification.md` has the heading `## S4 pre-fix failures` with the failing test ids, the first error of each, and the measured shift; no fix commit exists yet.
  - Context: Run `pytest tests/unit/modelling/test_world_location_consumers.py` on the unfixed code (the root suite only; one call). Record: the tests that fail and their first error; the fault ellipse test
    as "coverage, no defect". Add a behavioural measure with the production path: for `ANCHORS[1]` and each of the 3 rotations, the distance between the observation
    `anchor + (2, 3, 4)` and the location in `feature.get_data()` (use `world_foliation_model`). Use the layout of the S3 pre-fix section of the same file. FR-031.
    Write in ASD-STE100 Simplified Technical English.

- [x] T056 [US3] [M] Test: `rescale` and `scale` warn, keep behaviour, and no library caller warns (U3.11; FR-033).
  - Files: `tests/unit/modelling/test_world_location_consumers.py`
  - Needs: T054
  - Wave: 27
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k "rescale_warns or scale_warns or internal_callers"` fails before T057 to T062 (record the failing ids) and passes after T062.
  - Context: Append to the file of T054. `test_rescale_warns_and_keeps_behaviour`: 4 anchors x 3 rotations = 12 cells. `with pytest.warns(DeprecationWarning, match="bounding_box.reproject")`:
    `model.rescale(points)` raises one warning; the warning is attributed to the test file (`record[0].filename == __file__`); the result equals `model.bounding_box.reproject(points)` within `tp(L, W)`;
    `inplace=True` changes the passed array like `reproject(points, inplace=True)` does. `test_scale_warns_and_keeps_behaviour`: same 12 cells with `match="bounding_box.project"` and `model.scale`
    against `model.bounding_box.project`. `test_internal_callers_do_not_warn`: 1 cell. Use `warnings.filterwarnings("error", message=".*(rescale|scale) is deprecated.*")`; run `feature.get_data()`,
    `ModelExporter.get_block_model(model)` and `IntrusionFrameBuilder._network_model_values(data)` with no warning. Also call the exporter functions that
    `tests/unit/io/test_exporters.py` already calls (read it for the call pattern) with no warning. `rescale(self, points, *, inplace=False)` and `scale(self, points, *, inplace=False)` keep their signatures.

## Phase 16: Stage S4 implementation

- [x] T057 [P] [US3] [M] Remove the extra conversion from `GeologicalFeature.get_data`.
  - Files: `LoopStructural/modelling/features/_geological_feature.py`
  - Needs: T052, T055
  - Wave: 28
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k "feature_data_locations or example_observation or exported_data"` passes.
  - Context: In `get_data(self, value_map: dict | None = None)` (lines 303 to 372) remove the four `self.model.rescale(...)` calls. The four location arguments become
    `gradient_constraints[:, :3]`, `norm_constraints[:, :3]`, `value_constraints[value_constraints == v, :3]` (keep this mask as it is; pre-existing defect OP10) and `value_constraints[:, :3]`.
    Add to the docstring: "The locations and the vectors are in world coordinates, as the builder stores them." Change no other line. The commented-out block at the end of the method stays.
    Do not touch the misfit methods (S5).

- [x] T058 [P] [US3] [S] Evaluate the intrusion network at the stored world points.
  - Files: `LoopStructural/modelling/intrusions/intrusion_frame_builder.py`
  - Needs: T053, T055
  - Wave: 28
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k intrusion_classification` passes.
  - Context: Add the private method to `IntrusionFrameBuilder`, with a numpy docstring ("The id of the unit at each point. The model is evaluated at the world points; no conversion is applied."):
    `def _network_model_values(self, data: pd.DataFrame) -> np.ndarray:` with the body `xyz = data.loc[:, ["X", "Y", "Z"]].to_numpy(dtype=float)` then `return self.model.evaluate_model(xyz)`.
    In `set_intrusion_steps_parameters`: delete the assignment of `intrusion_network_data_xyz` (lines 333 to 335; no other use), and replace the statement at line 352 (the one with
    `self.model.rescale(...)`) by `intrusion_network_data.loc[:, "model_values"] = self._network_model_values(intrusion_network_data)`. Keep `pd` and `np` imports as they are.

- [x] T059 [P] [US3] [S] Deprecate `GeologicalModel.rescale` and `GeologicalModel.scale`.
  - Files: `LoopStructural/modelling/core/geological_model.py`
  - Needs: T056
  - Wave: 28
  - Done when: `pytest tests/unit/modelling/test_world_location_consumers.py -k "rescale_warns or scale_warns"` passes, and `pytest tests/unit/test_public_api_contract.py` passes.
  - Context: Keep both signatures and the decorator `@public_api(tier="stable")`. `warnings` is already imported (line 649). New body of `rescale(self, points: np.ndarray, *, inplace: bool = False) -> np.ndarray`:
    `warnings.warn("GeologicalModel.rescale is deprecated and will be removed in a future release; use model.bounding_box.reproject(points) instead", DeprecationWarning, stacklevel=2)` then
    `return self.bounding_box.reproject(points, inplace=inplace)`. New body of `scale(...)`: the same warning with `scale` and `model.bounding_box.project(points)`, then
    `return self.bounding_box.project(np.array(points).astype(float), inplace=inplace)`. Add to each docstring a `.. deprecated::` note that names the replacement and says
    "Deprecated by spec 005, stage S4". Change no other line.

- [x] T060 [P] [US3] [S] Delete the four discarded `model.rescale` calls in the exporters.
  - Files: `LoopStructural/export/exporters.py`
  - Needs: T056
  - Wave: 28
  - Done when: `grep -c "model.rescale" LoopStructural/export/exporters.py` prints 0, and `pytest tests/unit/io/test_exporters.py` passes.
  - Context: The four statements are `model.rescale(verts)`, `model.rescale(points)`, `model.rescale(xyz)` and `model.rescale(bbox)` (near lines 92, 366, 431 and 480). Each discards its result,
    so each does nothing. Delete only these statements. Keep the argument `real_coords` in every signature (it is ignored; pre-existing). Remove a local name only if the deleted line was its last use and ruff reports it unused.

- [x] T061 [P] [US3] [S] Use the bounding box in `ModelExporter.get_block_model`.
  - Files: `LoopStructural/modelling/core/_model_exporter.py`
  - Needs: T056
  - Wave: 28
  - Done when: `grep -c "model.rescale" LoopStructural/modelling/core/_model_exporter.py` prints 0, and `pytest tests/unit/io/test_exporters.py` passes.
  - Context: At line 77 replace `model.rescale(model.bounding_box.cell_centres())` by `model.bounding_box.reproject(model.bounding_box.cell_centres())`. The result is the same (local to world).
    Change no other line.

- [x] T062 [P] [US3] [S] Fix the double conversion in the faulted intrusion example.
  - Files: `examples/3_fault/plot_1_faulted_intrusion.py`
  - Needs: T056
  - Wave: 28
  - Done when: `grep -c "rescale" examples/3_fault/plot_1_faulted_intrusion.py` prints 0, and `python -m py_compile examples/3_fault/plot_1_faulted_intrusion.py` passes.
  - Context: Two `viewer.add_points` calls (near lines 71 and 100) wrap the world data of the model in `model.rescale(..., inplace=False)`. Pass
    `model.data[model.data["feature_name"] == "strati"][["X", "Y", "Z"]].values` directly (world points) and delete the `model.rescale(` wrapper with its `inplace=False`. Keep the other arguments of each call.
    Do not run the example (it needs the viewer).

## Phase 17: Stage S4 finish

- [x] T063 [Finish] [M] Review the existing tests for golden changes and record the S4 post-fix results.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T057, T058, T059, T060, T061, T062
  - Wave: 29
  - Done when: `verification.md` has the headings `## S4 post-fix results` and `## S4 golden review (C-004)`, and `pytest tests/unit/modelling/test_world_location_consumers.py` passes (all cells).
  - Context: Run, in separate calls: the new test file; `tests/unit/modelling/intrusions/test_intrusions.py`; `tests/unit/io/test_exporters.py`; `tests/unit/test_public_api_contract.py`.
    Record commands, counts and times. Under the golden review list each existing test whose expected value must change, with the defect and its effect (C-004), or write
    "No golden change in S4." Do not change an existing test; report one that fails. Record the new cells (count) and that the fault ellipse test is "coverage, no defect". Write in ASD-STE100 Simplified Technical English.

- [x] T064 [P] [Finish] [M] Record the S4 notices in `COMPAT.md`.
  - Files: `COMPAT.md`
  - Needs: T057, T059
  - Wave: 29
  - Done when: `grep -c "GeologicalModel.rescale(points, inplace=False)" COMPAT.md` prints at least 1 and `grep -c "GeologicalModel.scale(points, inplace=False)" COMPAT.md` prints at least 1.
  - Context: Read `COMPAT.md` first, the S1 to S3 notices and the deprecation table, and copy their layout. (1) Add a migration notice: `GeologicalFeature.get_data()` (used by `Loop2DView.add_data` and
    `ModelExporter.save`) returned locations moved by the local origin and rotated a second time, and the intrusion classification evaluated the model at such points. Both now use the stored world locations.
    With an identity frame nothing changes. (2) Add two rows to the deprecation table. `rescale`: old path `GeologicalModel.rescale(points, inplace=False)`; new path
    `GeologicalModel.bounding_box.reproject(points, inplace=False)`. `scale`: old path `GeologicalModel.scale(points, inplace=False)`; new path `GeologicalModel.bounding_box.project(points, inplace=False)`.
    Both: shim added `2026-11-05`; first released "TBD (next release containing this warning)"; introduced by the S4 commit; earliest removal "TBD (2 minor releases after the first release that warns)"; removed empty; status "Active"; owner "Core maintainers".
    (3) Raise the count of "Active shims" by 2. Do not change other text.

- [x] T065 [P] [Finish] [S] Record the S4 deprecation decision in `ROADMAP.md`.
  - Files: `ROADMAP.md`
  - Needs: T059
  - Wave: 29
  - Done when: `grep -c "Deprecation of GeologicalModel.rescale and scale (spec 005, S4)" ROADMAP.md` prints 1.
  - Context: Add a new `###` subsection as the last one of `## Decisions`, directly before `## Stage sequence`. Heading: `### Deprecation of GeologicalModel.rescale and scale (spec 005, S4)`.
    Body (ASD-STE100): the decision (2026-11-05, C-013 and C-016); `rescale` (local to world) and `scale` (world to local) keep their behaviour and signature and emit a `DeprecationWarning`
    that names `bounding_box.reproject` and `bounding_box.project`; this replaces the status-log entry "2c-3" that said both stay as signature-stable pass-throughs; the tier stays "stable" until removal;
    removal not before 2 minor releases after the first release that warns, and not before known callers (the QGIS plugin) have moved; the library no longer calls them. Keep the style of the neighbouring subsections.
    Do not change other text.

- [x] T066 [P] [Finish] [S] Mark the two methods as deprecated in `API.md`.
  - Files: `API.md`
  - Needs: T059
  - Wave: 29
  - Done when: `grep -c "deprecated, use .bounding_box.reproject" API.md` prints 1.
  - Context: Line 122 lists the geometry methods of `GeologicalModel`. Change that entry to read: "Geometry: `scale` (deprecated, use `bounding_box.project`), `rescale` (deprecated, use `bounding_box.reproject`),
    `regular_grid`, `bounding_box` (property)". Keep the line format of its neighbours. Do not change other text.

- [x] T067 [Finish] [M] Record the S4 rollback check in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T063
  - Wave: 30
  - Done when: `verification.md` has the heading `## S4 rollback checks` with the results of the 4 checks below, and `git worktree list` shows no throw-away worktree.
  - Context: Use a throw-away worktree at the S4 start commit (the `start_commit` of stage S4 in `run-state.json`; if the key is missing, use the parent of the first S4 code commit):
    `git worktree add --detach ../LoopStructural-005-s4-rollback <commit>`, with its own `PYTHONPATH` (check `LoopStructural.__file__`). Checks: (a) `from loopstructuralvisualisation import Loop2DView, Loop3DView` works and
    `inspect.signature(Loop2DView.add_data)` is equal before and after; (b) on the old code `get_data()` shows the spec example at `(202,403,604)`, on the S4 code at `(102,203,304)`; (c) with an identity frame the two codes
    give the same `get_data()` output; (d) on the old code `model.rescale(p)` and `model.scale(p)` give no warning, on the S4 code each gives one, and the values are equal.
    Then `git worktree remove --force <tmp>` and `git worktree prune`.

- [x] T068 [Finish] [M] Run the S4 gates and complete `verification.md` and the `ROADMAP.md` status log.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`, `ROADMAP.md`
  - Needs: T064, T065, T066, T067
  - Wave: 31
  - Done when: `ruff check LoopStructural packages && pytest` and `pytest packages/loop_common/tests packages/loop_interpolation/tests` pass, `verification.md` has the S4 sections named below, and `ROADMAP.md` has one new status-log line for S4.
  - Context: Run the root `pytest` and the package tests in separate calls (the conftest plugins clash) with the worktree `PYTHONPATH`. Add to `verification.md`: `## S4 backend matrix` (Python FDI foliation data: mandatory; fault centre in the 2D viewer: coverage;
    intrusion classification: mandatory; native, Surfe, browser: not in S4), `## S4 test run time (OP2)`, and the result of `git grep -n "\.rescale("` (it must show only the commented block of `get_data` and the tests of the deprecation; report any other hit).
    Add the S4 gate counts to the post-fix section. Add to the status log of `ROADMAP.md` a line dated today: "Spec 005 stage S4 implemented: display, export and intrusion classification use stored world locations; rescale and scale deprecated."
    Report any library warning with the text "is deprecated" in the pytest summary. Write in ASD-STE100 Simplified Technical English.

## Phase 18: Stage S5 tests (FR-014, FR-015, FR-023, FR-024)

**Goal**: Tests for the remaining resize paths, the cell-centre sampling and the residual and exclusion report are written before the code that they test.
**Independent test**: `pytest packages/loop_common/tests/test_resize_paths.py packages/loop_common/tests/test_cell_centre_sampling.py` and `pytest packages/loop_interpolation/tests/test_residual_exclusion_report.py`.

- [x] T069 [US4] [M] Test: the resize paths keep the extent and the counts (U4.2; FR-014; SC-011).
  - Files: `packages/loop_common/tests/test_resize_paths.py`
  - Needs: T068
  - Wave: 32
  - Done when: `pytest packages/loop_common/tests/test_resize_paths.py` fails before T075, T076 and T077 (record the failing ids and the first error of each) and passes after T077.
  - Context: New file. Classes: `StructuredGrid3DGeometry` (`loop_common.geometry._structured_grid_3d`), `StructuredGrid2DGeometry` (`_structured_grid_2d`), `RectilinearGrid` (`loop_common.supports._3d_rectilinear_grid`). Read the golden test `test_structured_support_reference.py` for the constructor calls.
    3D grids: `G1` origin `(0,0,0)`, nsteps `(11,6,4)`, step `(1,2,3)`; `G2` origin `(100,-20,5)`, nsteps `(5,5,5)`, step `(0.5,0.5,0.5)`; `G3` origin `(-1e4,2e3,-30)`, nsteps `(21,3,8)`, step `(10,100,7.5)`. 2D: `H1` origin `(0,0)`, nsteps `(5,5)`, step `(1,1)`; `H2` origin `(10,-5)`, nsteps `(8,4)`, step `(0.5,2)`.
    `nsteps` counts nodes. Extent `L = nsteps_cells * step_vector`. Length tolerance `1e-12 * max(1, L)`; counts must be exactly equal.
    `test_reapplying_same_value_10_times_changes_nothing_3d`: 3 grids x 4 setters (origin, maximum, nsteps, step_vector) = 12 cells; after 10 repeats `origin`, `maximum`, `step_vector`, `nsteps` equal the start (SC-011).
    `test_nsteps_setter_keeps_extent`: 3 cells; nsteps doubled then halved; `maximum` equals the start; node count equals the request. `test_step_vector_setter_encloses_extent`: 3 cells; step x 1.5; new `maximum` >= old and < old + one new step.
    `test_step_vector_a_b_a_extent_drift`: 3 cells; step A, then B = 1.5 A, then A; final `maximum` minus start is >= 0 and < `A + B` on each axis; print the measured value. `test_nsteps_below_two_is_rejected_without_change`: 1 cell; `ValueError`; state equal.
    `test_nsteps_not_whole_is_rejected_without_change`: 1 cell; `(5.5, 5, 5)` and NaN raise `ValueError`, state equal; `(5.0, 5.0, 5.0)` is accepted and stored as `int`.
    `test_set_nelements_keeps_extent`: 3 cells; after `set_nelements(200)` `maximum` equals the start; the return value equals `n_elements`; `nsteps == max(ceil(L / h), 2)` with `h = (volume / 200) ** (1 / 3)` computed in the test. `test_set_nelements_small_request_is_clamped`: `G1.set_nelements(1)` raises nothing; every nsteps >= 2; extent kept.
    `test_2d_counts_are_derived`: 2 grids x 3 checks = 6 cells; after `nsteps` assignment (a float list), `nsteps_cells`, `n_cell_x`, `n_cell_y` and `n_elements` agree with the new `nsteps`; `nsteps` is an int array; `maximum == origin + nsteps * step_vector` (legacy value). `test_rectilinear_resize_is_rejected`: 3 setters; `NotImplementedError`; nodes, origin and counts unchanged.
    `StructuredGrid2DGeometry.n_elements` exists (`nsteps_cells[0] * nsteps_cells[1]`). Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T070 [P] [US4] [M] Test: cell-centre sampling agrees with the grid cells for rotated boxes (U4.3; FR-015; SC-011).
  - Files: `packages/loop_common/tests/test_cell_centre_sampling.py`
  - Needs: T068
  - Wave: 32
  - Done when: `pytest packages/loop_common/tests/test_cell_centre_sampling.py` passes on the current tree (60 cells, "coverage, no defect", OP12). If a cell fails, record the failing ids in the report and keep the test.
  - Context: New file. Boxes: `unit` and `offset_anisotropic` from `test_bounding_box_structured_grid.py` (read it for the construction). Frames: rotation 0, 30 and 90 degrees about z; local origin `(0,0,0)` and `(7,-3,11)` (6 frames; with 2 boxes = 12 cells for each test).
    Set the frame with `box.set_local_transform(local_origin=..., local_rotation=...)`. Read `test_rigid_transform_contract.py` for the call and for the tolerance `Tp`. Rotation 0 with origin 0 is the control.
    `test_each_cell_has_exactly_one_centre` (12 cells; the main evidence for FR-015 and SC-011): `grid = box.structured_grid(local_coordinates=True)`; `world = box.reproject(box.cell_centres(order="F"))`; `grid.cell_index_at(box.project(world))` equals `arange(n_cells)`: no `-1`, no repeat; `np.bincount` is all 1.
    `test_local_coordinates_grid_associates_own_centres` (12 cells): every `box.project(world)` falls in a cell of that grid (`cell_index_at(...) != -1`). `test_centres_are_cell_midpoints_in_world` (12): `world` equals the mean of the 8 corner nodes of each cell of `box.structured_grid()` within `Tp`.
    `test_centres_match_exported_cells` (12): world centres equal `grid.cell_centres` within `Tp`. `test_centres_inside_box` (12): each world centre is inside `[origin, maximum]` within `Tp`.
    Only the first two tests are evidence for SC-011. Say so in their docstrings. Do not change `_bounding_box.py`. Write in ASD-STE100 Simplified Technical English.

- [x] T071 [P] [US4] [M] Test: the pure residual and exclusion functions and the report class (U4.4, U4.5; FR-023, FR-024).
  - Files: `packages/loop_interpolation/tests/test_residual_exclusion_report.py`
  - Needs: T068
  - Wave: 32
  - Done when: `pytest packages/loop_interpolation/tests/test_residual_exclusion_report.py` fails before T073 with `ImportError` (the names do not exist yet), then passes after T073.
  - Context: New file. No model. Import from `loop_interpolation._diagnostics`: `ConstraintDiagnosticsReport`, `ConstraintFamilyDiagnostics`, `RegionCoverageDiagnostics` (read the file for their fields and for how to build a family), and the new `direction_cosines`, `build_exclusion_records`, `summarise_residuals`, `ExclusionRecord`, `ResidualSummary`, `UNKNOWN_EXPLANATION`.
    Signatures: `direction_cosines(model_vectors, data_vectors) -> ndarray (N,)`; `build_exclusion_records(excluded_rows: dict[str, dict[str, int]], families: dict[str, ConstraintFamilyDiagnostics], rows_per_point: dict[str, int] | None = None) -> tuple[ExclusionRecord, ...]`.
    `test_direction_cosines_pure_function` (6 cells): equal, opposite and orthogonal vectors give 1, -1 and 0; a zero row and a NaN row give NaN; shape mismatch raises `ValueError`; both inputs unchanged; scale invariance.
    `test_report_without_verification_is_unchanged` (1): `to_dict().keys()` equals the pre-S5 key set (no existing test has it: write this literal: report keys `interpolator_type`, `total_rows`, `families`, `region_coverage`, `outside_model_points`; family keys `active`, `row_count`, `dropped_rows`, `effective_weight_mean`, `effective_weight_min`, `effective_weight_max`, `source_point_count`, `outside_model_point_count`; region keys `total_support_nodes`, `active_region_nodes`, `inactive_region_nodes`, `active_fraction`; build the report with a family and a region coverage) and the old fields keep their order (`dataclasses.fields`).
    `test_builder_records_keep_their_own_key` (1): `{"restricted_orthogonal": {"invalid_reference_direction": 3}}` gives a record with `family == "restricted_orthogonal"`, count 3, no `unattributed` record, no family count changed. `test_unknown_reason_keeps_its_count` (1): reason `"x"` keeps its count and gets `UNKNOWN_EXPLANATION`.
    `test_unattributed_rows_are_reported` (1): `value` family, `dropped_rows` 5, outside points 2: `unattributed` count 3. `test_gradient_outside_points_give_no_false_unattributed` (1): `gradient` family, 10 source points, 2 outside, `dropped_rows` 4, `rows_per_point={"gradient": 2}`: no record; with `dropped_rows` 6: `unattributed` 2.
    `test_dropped_rows_none_gives_no_unattributed` (1): `interface` family, `dropped_rows=None`, no outside points: no record. `test_surfe_style_dropped_equals_outside` (1): `rows_per_point=None`, `dropped_rows` equal to the outside count: no record. Write in ASD-STE100 Simplified Technical English.

- [x] T072 [US4] [M] Record the S5 pre-fix evidence for the resize and cell-centre tests in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T069, T070
  - Wave: 33
  - Done when: `verification.md` has the heading `## S5 pre-fix failures` with the failing test ids, the first error of each, and the measured values; no S5 fix commit exists yet.
  - Context: Run the two `loop_common` test files on the unfixed code, in one call: `pytest packages/loop_common/tests/test_resize_paths.py packages/loop_common/tests/test_cell_centre_sampling.py` (PowerShell; the worktree on `PYTHONPATH`; the venv python is `C:\Users\lachl\dev\LoopStructural\.venv\Scripts\python.exe`).
    Record: the failing resize tests and the first error of each; the measured extent drift of the `nsteps` setter (example: nodes 11 to 21 with step 1 gives `maximum` 10.48 instead of 10); the extent drift of the `step_vector` setter; the stale 2D counts. Record the cell-centre tests as "coverage, no defect" if they pass, or "reproduced" with ids if not (OP12).
    The 2D `maximum` is unchanged, so no 2D `maximum` test must fail. Use the layout of `## S4 pre-fix failures` in the same file. FR-031. Write in ASD-STE100 Simplified Technical English.

## Phase 19: Stage S5 diagnostics types and resize code

- [x] T073 [P] [US4] [M] Add the residual and exclusion types and functions to `_diagnostics.py`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/_diagnostics.py`
  - Needs: T071
  - Wave: 33
  - Done when: `pytest packages/loop_interpolation/tests/test_residual_exclusion_report.py` passes, and `pytest packages/loop_interpolation/tests -k diagnostics` still passes.
  - Context: The file has `from __future__ import annotations`; keep the existing classes and field order. Add `ResidualSummary` (frozen dataclass: `name: str`, `kind: str` ("physical" or "assembled"), `unit: str` ("field", "degree" or "assembled"), `count: int` (finite rows), `nan_count: int`, `rms`, `max_abs`, `mean_abs`: `float | None`, `None` when `count == 0`).
    Add `ExclusionRecord` (frozen: `family: str`, `reason: str`, `count: int`, `explanation: str`). Add `EXCLUSION_EXPLANATIONS: dict[str, str]` with the keys `invalid_reference_value`, `invalid_reference_direction`, `evaluation_error`, `outside_model`, `unattributed` (plain one-sentence texts; `outside_model` says the count is in points, `unattributed` says it is in matrix rows) and `UNKNOWN_EXPLANATION = "No explanation is registered for this reason."`.
    `summarise_residuals(name, kind, unit, residuals) -> ResidualSummary`: flatten; count non-finite values in `nan_count`; statistics from finite values only; do not change the input.
    `direction_cosines(model_vectors, data_vectors) -> np.ndarray`: shape mismatch raises `ValueError`; per row dot / (norm * norm); a row with a zero or not finite norm gives NaN; inputs unchanged.
    `build_exclusion_records(excluded_rows, families, rows_per_point=None) -> tuple[ExclusionRecord, ...]`, sorted by `(family, reason)`: one record per reason in `excluded_rows` (family = its key; subtract nothing); one `outside_model` record per family with `outside_model_point_count > 0`; one `unattributed` record per family with `dropped_rows` not `None` and `dropped_rows - outside_model_point_count * k > 0`, with `k = rows_per_point.get(family, 1)`; the count is the difference. Count 0 makes no record. An unknown reason gets `UNKNOWN_EXPLANATION`.
    Report class: append after `outside_model_points`: `physical_residuals: dict[str, ResidualSummary] = field(default_factory=dict)`, `assembled_residuals: dict[str, ResidualSummary] = field(default_factory=dict)`, `exclusions: tuple[ExclusionRecord, ...] = ()`. `to_dict` adds `"physical_residuals"`, `"assembled_residuals"` (dicts of the fields of each summary) and `"exclusions"` (list of dicts) only when not empty. `summary()` adds one line per non-empty block after the existing lines. `total_rows` and `active_families` do not change. numpy docstrings and type hints (ruff `D`, `ANN`); py39.

- [x] T074 [P] [US4] [S] Export four new names from `loop_interpolation`.
  - Files: `packages/loop_interpolation/src/loop_interpolation/__init__.py`
  - Needs: T073
  - Wave: 34
  - Done when: `python -c "from loop_interpolation import ResidualSummary, ExclusionRecord, direction_cosines, build_exclusion_records"` prints no error.
  - Context: Read the file: it already imports from `._diagnostics` and has an `__all__`. Add `ResidualSummary`, `ExclusionRecord`, `direction_cosines` and `build_exclusion_records` to the same import and to `__all__`, in the style of the neighbours. Do NOT export `summarise_residuals`, `UNKNOWN_EXPLANATION` or `EXCLUSION_EXPLANATIONS`. Change no other line.

- [x] T075 [P] [US4] [M] Keep the extent in the 3D `nsteps` and `step_vector` setters and in `set_nelements`.
  - Files: `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`
  - Needs: T069, T072
  - Wave: 34
  - Done when: `pytest packages/loop_common/tests/test_resize_paths.py -k "3d or nsteps or step_vector or set_nelements or a_b_a"` passes.
  - Context: Add `from typing import Optional` (no `from __future__ import annotations`; py39). `_cells_for_length(self, length, step_vector: Optional[np.ndarray] = None)`: `step_vector` defaults to the current step; the S1 rule is unchanged (`ceil(length / step - 1e-9)`; 0 cells become 2; a negative count raises `ValueError("Cannot resize the grid.")`).
    `set_nelements`: `nsteps = np.maximum(np.ceil((self.maximum - self.origin) / step_vector).astype(int), 2)`. `nsteps` setter: `nsteps = np.broadcast_to(np.asarray(nsteps), self._nsteps.shape)`; if not all finite or any `nsteps != np.round(nsteps)` raise `ValueError("Cannot resize the grid. nsteps must be whole numbers.")`; `nsteps = nsteps.astype(int)`; if any `< 2` raise `ValueError("Cannot resize the grid. nsteps must be at least 2.")`;
    `extent = self.nsteps_cells * self._step_vector`; `self._step_vector = extent / (nsteps - 1)`; `self._nsteps = np.array(nsteps)`; `self._on_geometry_change()`. No state changes before a raise.
    `step_vector` setter: `step_vector = np.array(step_vector, dtype=float)`; `extent = self.nsteps_cells * self._step_vector`; `cells = self._cells_for_length(extent, step_vector=step_vector)`; `self._nsteps = cells + 1`; `self._step_vector = step_vector`; `self._on_geometry_change()`.
    Before the edit run `git grep -n "\.nsteps = \|\.step_vector = \|set_nelements"` and read each writer; report one that depends on the old product rule. Keep the `origin` and `maximum` setters of S1. Principle: change no public name.

- [x] T076 [P] [US4] [M] Derive the 2D counts from `nsteps` and keep the legacy `maximum`.
  - Files: `packages/loop_common/src/loop_common/geometry/_structured_grid_2d.py`
  - Needs: T069, T072
  - Wave: 34
  - Done when: `pytest packages/loop_common/tests/test_resize_paths.py -k 2d_counts` passes, and `pytest packages/loop_common/tests/test_2d_discrete_support.py` and `pytest tests/unit/interpolator/test_2d_discrete_support.py` still pass with no change (the four `maximum` assertions; separate calls).
  - Context: First run `git grep -n "n_cell_x\|n_cell_y\|nsteps_cells\|\.maximum ="` and read each writer outside `__init__`; a writer must move to the derived form (report it). In `__init__` (lines 23 to 30) replace the assignments of `self.nsteps`, `self.maximum`, `self.nsteps_cells`, `self.n_cell_x` and `self.n_cell_y` with `self.nsteps = nsteps` (the setter coerces); keep the other lines. `nsteps` is set before `origin` is read.
    Add: `nsteps` property returning `self._nsteps`; setter `self._nsteps = np.ceil(np.asarray(nsteps)).astype(int)`; read-only properties `nsteps_cells` (`self.nsteps - 1`), `n_cell_x` (`self.nsteps[0] - 1`), `n_cell_y` (`self.nsteps[1] - 1`); read-only `maximum` returning `self.origin + self.nsteps * self.step_vector`.
    Keep a comment on `maximum`: legacy value (FR-028); it counts nodes as cells, so it is one step beyond the last node; `nodes` uses `nsteps_cells`; do not change it. Keep the NumPy integer type of the old attributes. Use numpy docstrings and type hints. Rule: the four existing `maximum` assertions pass with no change (decision Q1).

- [x] T077 [P] [US4] [S] Make the resize setters of `RectilinearGrid` raise.
  - Files: `packages/loop_common/src/loop_common/supports/_3d_rectilinear_grid.py`
  - Needs: T069, T072
  - Wave: 34
  - Done when: `pytest packages/loop_common/tests/test_resize_paths.py -k rectilinear` passes, and `pytest packages/loop_common/tests/test_rectilinear_grid.py` still passes.
  - Context: `RectilinearGrid` inherits the `origin`, `nsteps` and `step_vector` setters of the 3D geometry; a call changes the counts but not the node arrays. Add a module constant `_RESIZE_MESSAGE = "A rectilinear grid has fixed node arrays and cannot be resized. Build a new RectilinearGrid."`.
    Override the three properties: each getter returns the private attribute (`self._origin`, `self._nsteps`, `self._step_vector`); each setter does `raise NotImplementedError(_RESIZE_MESSAGE)` and changes no state. Pattern: `@property def origin(self): return self._origin` then `@origin.setter def origin(self, origin): raise NotImplementedError(_RESIZE_MESSAGE)`. The constructors write the private attributes, not the setters, so construction is unchanged. `maximum` has no setter: leave it. Use numpy docstrings and type hints.

## Phase 20: Stage S5 golden review and feature tests

- [x] T078 [US4] [M] Survey the numeric tests, review the golden changes and regenerate the 3D golden data (C-004).
  - Files: `specs/005-coordinate-frame-correctness/verification.md`, `packages/loop_common/tests/data/structured_support_reference.npz`, `packages/loop_common/tests/data/structured_grid_2d_reference.npz`
  - Needs: T075, T076, T077
  - Wave: 35
  - Done when: `pytest packages/loop_common/tests/test_structured_support_reference.py packages/loop_common/tests/test_structured_grid_2d_reference.py` passes, and `verification.md` has the heading `## S5 golden review (C-004)` with one entry for each changed key or test (old value, new value, cause).
  - Context: Order (C-004, SC-015: 0 tolerance changes). (1) Run the 3D golden test and list each failing key. (2) Check each key against the rule: each setter keeps the extent or encloses it by less than one cell. (3) Write each key with its old and new value in the section. (4) Only then run `python packages/loop_common/tests/test_structured_support_reference.py --update`. (5) Run the test again.
    Run the 2D golden test. It must not change (the 2D `maximum` stays). Run `--update` for it only if a key differs, and record the key. Numeric survey, in separate calls: `pytest tests/unit/modelling`; `pytest packages/loop_interpolation/tests`; `pytest packages/loop_common/tests` (PowerShell; worktree on `PYTHONPATH`; `packages/loopstructural_visualisation/src` also on `PYTHONPATH`; `pyevtk` is only in a scratch folder: add it to `PYTHONPATH` or report the skip).
    For each test that fails or changes a number, record the test, the old value, the new value and the cause. Do not edit a test file. `packages/loop_interpolation/tests/golden` is the parity oracle (OP17): record each change there and do NOT regenerate it. If the cause is the intended extent change, name the file for the maintainer; if the cause is not clear, report it. Write in ASD-STE100 Simplified Technical English.

- [x] T079 [P] [US4] [M] Record the S5 post-fix evidence for the resize and cell-centre tests.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T078
  - Wave: 36
  - Done when: `verification.md` has the heading `## S5 resize and cell-centre results` and `pytest packages/loop_common/tests/test_resize_paths.py packages/loop_common/tests/test_cell_centre_sampling.py` passes.
  - Context: Run the two files in one call. Record the commands, the counts (resize: 12 + 3 + 3 + 3 + 1 + 1 + 3 + 1 + 6 + 3 = 36 cells; cell centres: 5 tests x 12 = 60 cells), the times, and the measured extent drift of the A, B, A test for each grid (the bound is `A + B`). Compare with `## S5 pre-fix failures`: for each id that failed, write that it now passes.
    Write that `test_each_cell_has_exactly_one_centre` and `test_local_coordinates_grid_associates_own_centres` are the evidence for FR-015 and SC-011, and the other three are coverage. Record OP12: "not reproduced" if all 60 cells passed before the fix, and no change to `_bounding_box.py`. Record the `ruff check LoopStructural packages` result for the three `loop_common` source files. Write in ASD-STE100 Simplified Technical English.

- [x] T080 [US4] [M] Test: misfit methods keep state, compare in world coordinates, and keep their shapes (U4.4; FR-023; SC-009).
  - Files: `tests/unit/modelling/test_diagnostics_frames.py`
  - Needs: T073, T078
  - Wave: 36
  - Done when: `pytest tests/unit/modelling/test_diagnostics_frames.py -k "leaves_state_unchanged or exact_field or value_ or empty or old_misfit_shapes"` fails before T083 (record the failing ids and the first error of each) and passes after T084.
  - Context: New file. Use `world_foliation_model`, `ANCHORS`, `ROTATIONS_DEG` from `tests.fixtures.coordinate_frames`; solve with the pattern of `tests/unit/modelling/test_reference_constraint_frames.py` (read both). Engines: FDI and PLI. Cell grid: 4 anchors x 3 rotations x 2 engines = 24 (value tests: 12 cells).
    `test_misfit_leaves_state_unchanged` (24): after a solve, take deep copies of every `interpolator.data[...]` array, `c`, `interpolator.up_to_date`, `builder._up_to_date`, the weights of `constraints` and `id(latest_diagnostics_report)`; call the two OLD methods `evaluate_gradient_misfit`, `evaluate_value_misfit`; every snapshot is equal (`np.array_equal`). It uses no new name. `test_new_misfit_methods_leave_state_unchanged` (24): the same, with `evaluate_gradient_cosines`, `evaluate_value_residuals` and `verification_report()`.
    `test_gradient_cosine_matches_analytic_direction_exact_field` (24): `monkeypatch` replaces `feature.evaluate_gradient` with a stub that returns the analytic world direction of the R1 foliation at the given points; `abs(1 - cos) <= 1e-10` for the old method (one row) and for `evaluate_gradient_cosines`. `..._solved_model` (24): the real solved model; `tol_solved` is measured: record the maximum and set the bound with a margin of 2 (it is not `1e-10`).
    `test_value_residuals_are_absolute_difference` (12): equals `abs(value - evaluate_value(world_points))` at the 27 value rows; no NaN. `test_value_misfit_is_normalised_by_field_range` (12): equals `evaluate_value_residuals() / (feature.max() - feature.min())` within `1e-12`. `test_value_misfit_with_zero_range_does_not_divide` (2): `monkeypatch` makes `feature.max` and `feature.min` return the same number (cell 1) or NaN (cell 2); `warnings.simplefilter("error")`; the result equals the unnormalised absolute difference.
    `test_misfit_with_no_rows_is_empty` (2): no gradient rows gives shape `(0,)`; no value rows gives shape `(0,)`. `test_old_misfit_shapes_are_kept` (4): one family `(1, N)`; two families of equal length `(2, N)`; two of different length `(1, N1 + N2)` (OP16, default); values equal `evaluate_gradient_cosines()`.

- [x] T081 [US4] [M] Test: a stale model is not rebuilt, and the report keeps residuals and exclusions apart (U4.4, U4.5; FR-023, FR-024; SC-009).
  - Files: `tests/unit/modelling/test_diagnostics_frames.py`
  - Needs: T080
  - Wave: 37
  - Done when: `pytest tests/unit/modelling/test_diagnostics_frames.py -k "stale or report_blocks or unsolved or exclusions or outside_points"` fails before T083 and T084 (record the failing ids; `test_stale_misfit_is_not_rebuilt` must fail on the current tree because the call rebuilds the model) and passes after T084.
  - Context: Append to the file of T080 (same fixtures; engines FDI and PLI). Stale recipe (settled, plan "Answers"): after the solve, for a stale interpolator assign `interpolator.up_to_date = False`; for a stale builder assign `builder._up_to_date = False` directly (not `builder.set_not_up_to_date`: it also stales the interpolator); assert the flag before the call. PLI note (settled): `P1Interpolator.add_gradient_constraints` is `pass`, so PLI assembles no `gradient` rows and the builder reports `dropped_rows = 2 * source_points` for it; `unattributed` for PLI gradient data is true. Run the no-false-unattributed check with FDI only; for PLI record in the test output the measured `unattributed` count and the `outside_model` count (print them) and assert `outside_model == 3`. `test_stale_misfit_is_not_rebuilt` (2 engines x 2 kinds, interpolator stale or builder stale = 4): uses the OLD methods only. Snapshot `interpolator.up_to_date`, `builder._up_to_date`, `c` and the data arrays; after the call every snapshot is equal (on the current tree `up_to_date` changes from `False` to `True`: the test fails; FR-031).
    `test_stale_misfit_warns_once_and_returns_nan` (4): the old and the new methods each write exactly one `UserWarning` that names the stale object; the new methods return `np.full(n, np.nan)`; the old methods return the old shape filled with NaN. `test_stale_verification_report_does_not_rebuild` (4): snapshot equal; both residual blocks empty; exclusions present; one warning.
    `test_report_blocks_are_separate` (2): the two dicts are different objects; physical entries have `kind == "physical"`, assembled entries `kind == "assembled"`; the physical `value` entry has unit `"field"` and its `rms` equals the RMS of `evaluate_value_residuals()` within `1e-12`; the assembled entries equal `summarise_residuals(name, "assembled", "assembled", r)` for `calculate_residual_for_constraints()`; `to_dict` has two separate blocks. Import `summarise_residuals` and `UNKNOWN_EXPLANATION` from `loop_interpolation._diagnostics`.
    `test_unsolved_model_gives_no_residual_blocks` (2): `interpolator.c = None`; both blocks `{}`; exclusions present; no exception. `test_exclusions_are_counted_and_explained` (2): the restricted-support reference case of `test_reference_constraint_frames.py`: builder records under `restricted_equality` and `restricted_orthogonal` with the analytic count; the explanation is not empty and not `UNKNOWN_EXPLANATION`.
    `test_gradient_outside_points_give_no_false_unattributed_model` (2): gradient data with 3 points outside the box: an `outside_model` record of 3 and no `unattributed` record for `gradient`.

- [x] T082 [US4] [M] Record the S5 pre-fix evidence for the diagnostics tests in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T081
  - Wave: 38
  - Done when: `verification.md` has the heading `## S5 diagnostics pre-fix failures` with the failing ids, the first error of each, and the measured values; no misfit fix commit exists yet.
  - Context: Run `pytest tests/unit/modelling/test_diagnostics_frames.py` on the code before T083 (root suite only; one call). Record: the state test with the old methods (a different array in `data["gradient"]` after the call); the exact-field cosine test (one column, local vectors); the stale-state test (the call rebuilds the model: `up_to_date` goes from `False` to `True`); the tests of the new names (they fail because the names do not exist).
    Record `test_value_residuals_are_absolute_difference` as "not reproduced" if it passes (the value path does not change the stored positions). Add one behavioural measure with the production path: for `ANCHORS[1]` and rotation 30, the largest absolute change of `interpolator.data["gradient"]` after one call of `evaluate_gradient_misfit`. Use the layout of `## S4 pre-fix failures`. FR-031. Write in ASD-STE100 Simplified Technical English.

## Phase 21: Stage S5 misfit code

- [x] T083 [US4] [M] Rewrite the misfit methods on shared non-mutating code and add the new methods.
  - Files: `LoopStructural/modelling/features/_geological_feature.py`
  - Needs: T073, T082
  - Wave: 39
  - Done when: `pytest tests/unit/modelling/test_diagnostics_frames.py -k "not report_blocks and not unsolved and not exclusions and not outside_points and not verification_report and not new_misfit_methods"` passes.
  - Context: The module has `from __future__ import annotations`. Add (import `direction_cosines` from `loop_interpolation._diagnostics`; `warnings`): `_stale_reason(self) -> Optional[str]`: text for the first stale object among `self.builder._up_to_date` is False, each builder in `self.builder.faults`, `self.interpolator.up_to_date` is False; `None` otherwise (builder `None`: only the interpolator). The text names the object and says to build the model first.
    `_gradient_cosines_fresh(self) -> list`: for each of `get_gradient_constraints()` and `get_norm_constraints()` with rows: `rows = np.array(array, dtype=float, copy=True)`; `positions = rows[:, :3]`, reproject with `self.interpolator.bounding_box.reproject` when `self.interpolator.bounding_box` is not `None` (settled: `GeologicalFeature` has no `bounding_box`; do not use `self.model.bounding_box`); `data = rows[:, 3:6]`, convert with `self.interpolator.local_to_world_vectors(data)` when that bounding box is not `None`; `model = self.evaluate_gradient(positions)`; append `direction_cosines(model, data)`.
    `_value_residuals_fresh(self) -> np.ndarray`: copy `get_value_constraints()`, reproject positions, return `np.abs(rows[:, 3] - self.evaluate_value(positions))`.
    `evaluate_gradient_cosines() -> ndarray (N,)`: on a stale model `warnings.warn(text, UserWarning, stacklevel=2)` and return `np.full(n, np.nan)` (`n` = row count of the stored arrays); else `np.concatenate` of the list (`np.zeros(0)` if empty). `evaluate_value_residuals() -> ndarray (N,)`: same rule with `_value_residuals_fresh`.
    `evaluate_gradient_misfit()` (old name, new body): same code; one row per non-empty family: `(1, N)` for one family, `(2, N)` for two of equal length, `(1, N1 + N2)` for two of different length, `(0,)` with no rows; stale gives that shape filled with NaN. `evaluate_value_misfit()`: stale gives NaN `(N,)`; else the residuals divided by `self.max() - self.min()` only when the range is finite and > 0, else unnormalised.
    Remove every `self.builder.up_to_date()` call and the in-place `grad /= norm` from these methods. Never write to a stored array, weight, coefficient or flag. numpy docstrings (state that the gradient value is a cosine in world coordinates; see COMPAT.md) and type hints. Read the whole old methods first (about lines 240 to 300).

- [x] T084 [US4] [M] Add `GeologicalFeature.verification_report` with physical and assembled residuals and exclusions.
  - Files: `LoopStructural/modelling/features/_geological_feature.py`
  - Needs: T083
  - Wave: 40
  - Done when: `pytest tests/unit/modelling/test_diagnostics_frames.py` passes (all cells), and `pytest tests/unit/modelling/test_scalar_field_grid.py` still passes.
  - Context: Add `_physical_residual_summaries(self) -> dict[str, ResidualSummary]`: `value` from `_value_residuals_fresh()` (kind `"physical"`, unit `"field"`); `gradient` and `normal` from `_gradient_cosines_fresh()`: `np.degrees(np.arccos(np.clip(cosine, -1, 1)))`, unit `"degree"`; use `summarise_residuals(name, "physical", unit, values)`. Add the public `verification_report(self) -> ConstraintDiagnosticsReport`.
    Body: `base = self.interpolator._build_constraint_diagnostics_report()` (a new object; do NOT set `latest_diagnostics_report`). `physical = {}`, `assembled = {}`. If `_stale_reason()` is not `None`: `warnings.warn(text, UserWarning, stacklevel=2)` and skip both blocks. Else `physical = self._physical_residual_summaries()`; and only if `interpolator.c` is not `None` and the interpolator has `calculate_residual_for_constraints`: for `(name, residual)` in `calculate_residual_for_constraints().items()` set `assembled[name] = summarise_residuals(name, "assembled", "assembled", residual)`.
    Exclusions are always built (they read state only): `excluded = getattr(self.builder, "excluded_rows", {})`; `rows_per_point = {"gradient": 2}` if the interpolator is a `DiscreteInterpolator`, else `None`. Return `dataclasses.replace(base, physical_residuals=physical, assembled_residuals=assembled, exclusions=build_exclusion_records(excluded, base.families, rows_per_point))`.
    An unsolved model (`c` is `None`) leaves BOTH residual blocks empty with no warning (plan decision 16; settled after the wave 37 review: guard the physical block on `interpolator.c is not None` too, because `evaluate_value` with `c = None` fails; the test `test_unsolved_model_gives_no_residual_blocks` requires both `{}`). The stale warning text must contain the word `interpolator` or `builder` (the stale object). The `gradient` rows_per_point is 2 for any `DiscreteInterpolator` (FDI gives 2 rows for each point; PLI gives 0 rows, so PLI gradient data gives an `unattributed` record, measured 8). Do not change `calculate_residual_for_constraints`. numpy docstring and type hints. Import `ConstraintDiagnosticsReport`, `ResidualSummary`, `summarise_residuals`, `build_exclusion_records` from `loop_interpolation._diagnostics`.

## Phase 22: Stage S5 finish

- [x] T085 [P] [Finish] [M] Record the S5 notices in `COMPAT.md`.
  - Files: `COMPAT.md`
  - Needs: T084
  - Wave: 41
  - Done when: `grep -c "evaluate_gradient_cosines" COMPAT.md` prints at least 1 and `grep -c "RectilinearGrid" COMPAT.md` prints at least 1.
  - Context: Read `COMPAT.md` first, the S1 to S4 notices, and copy their layout. Add migration notices. (1) `GeologicalFeature.evaluate_gradient_misfit()` changed the stored constraint arrays in place, compared a local vector with a world vector, and read one column. It keeps its name and its public shape and changes nothing. It returns the cosine between the world data direction and the model gradient for each row. The old value was a dot product with an unnormalised model gradient, so the new value differs. Use `evaluate_gradient_cosines()` for a 1-D array.
    Two families of different length now give one row `(1, N1 + N2)`; before, this raised `ValueError` (OP16). (2) Both misfit methods no longer rebuild a stale model: they write one `UserWarning` and return NaN; build the model first. (3) `evaluate_value_misfit()` keeps its normalisation by the field range and no longer divides by zero (zero or not finite range gives the absolute difference). `evaluate_value_residuals()` gives the unnormalised 1-D result.
    (4) Assigning `nsteps` or `step_vector` on a 3D structured grid, and `set_nelements`, keep the extent (a new step, not a product rule); `nsteps` must be whole and at least 2. `set_nelements` keeps `maximum`. Grid steps and solved values can change slightly. (5) A `RectilinearGrid` raises `NotImplementedError` on assignment of `origin`, `nsteps` or `step_vector`. (6) The 2D grid derives `nsteps_cells`, `n_cell_x`, `n_cell_y` from `nsteps`; `maximum` keeps its legacy value. (7) `ConstraintDiagnosticsReport` has three new optional fields (`physical_residuals`, `assembled_residuals`, `exclusions`); with an identity frame and no resize the value misfit is the same as before.
    Do not change other text. Write in ASD-STE100 Simplified Technical English.

- [x] T086 [P] [Finish] [S] Record the S5 misfit decision in `ROADMAP.md`.
  - Files: `ROADMAP.md`
  - Needs: T084
  - Wave: 41
  - Done when: `grep -c "Misfit diagnostics of GeologicalFeature (spec 005, S5)" ROADMAP.md` prints 1.
  - Context: Add a new `###` subsection as the last one of `## Decisions`, directly before `## Stage sequence`. Heading: `### Misfit diagnostics of GeologicalFeature (spec 005, S5)`. Body (ASD-STE100): decision date 2026-11-05 (spec 005, S5, Q2); the misfit methods `evaluate_gradient_misfit` and `evaluate_value_misfit` keep their names and public shape; the gradient value is now a cosine in world coordinates; a stale model gives a `UserWarning` and NaN (no rebuild);
    the new names `evaluate_gradient_cosines` and `evaluate_value_residuals` give the 1-D result; the notice is in `COMPAT.md`. Keep the style of the neighbouring subsections. Do not change other text.

- [x] T087 [Finish] [M] Record the S5 rollback check in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T074, T084
  - Wave: 41
  - Done when: `verification.md` has the heading `## S5 rollback checks` with the results of the 6 checks below, and `git worktree list` shows no throw-away worktree.
  - Context: Use a throw-away worktree at the S4 end commit (the commit before the first S5 code commit): `git worktree add --detach ../LoopStructural-005-s5-rollback <commit>`, with its own `PYTHONPATH` (check `LoopStructural.__file__` and `loop_interpolation.__file__`). Checks: (a) `dataclasses.fields` of `ConstraintFamilyDiagnostics`, `RegionCoverageDiagnostics` and the first five fields of `ConstraintDiagnosticsReport` are equal on both trees;
    (b) `to_dict()` of a report built with no verification is equal on both trees; (c) `from loop_interpolation import *` on the S4 tree gives a subset of the S5 names; (d) the signatures of `GeologicalFeature.evaluate_gradient_misfit` and `evaluate_value_misfit` are equal, and the result shape for one family is `(1, N)` on both;
    (e) with an unchanged grid, `origin`, `maximum`, `nsteps` and `step_vector` of a structured support are equal on both trees; (f) the four 2D `maximum` assertions pass on both trees. Then `git worktree remove --force <tmp>` and `git worktree prune`. Write in ASD-STE100 Simplified Technical English.

- [x] T088 [Finish] [M] Run the S5 gates and complete `verification.md` and the `ROADMAP.md` status log.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`, `ROADMAP.md`
  - Needs: T079, T085, T086, T087
  - Wave: 42
  - Done when: `ruff check LoopStructural packages && pytest` and `pytest packages/loop_common/tests packages/loop_interpolation/tests` pass, `verification.md` has the S5 sections named below, and `ROADMAP.md` has one new status-log line for S5.
  - Context: Run the root `pytest` and the package tests in separate calls (the conftest plugins clash) with the worktree `PYTHONPATH`, plus `packages/loopstructural_visualisation/src`; `pyevtk` is only in a scratch folder (add it or report the skip). Add to `verification.md`: `## S5 post-fix results` (commands, counts, times, the cells of `test_diagnostics_frames.py` and `test_residual_exclusion_report.py`, the measured `tol_solved`),
    `## S5 backend matrix` (FDI and PLI misfit, residual and exclusion report: mandatory; resize of structured, 2D and rectilinear grids: mandatory; P2, native, Surfe, residual reports for fold and P2 cells: not in S5; S2 and S6 cover them), `## S5 test run time (OP2)`, and a line for each open point: OP12 (cell-centre result), OP16 (default shape kept), OP17 (parity files changed or not, from T078).
    Add to the status log of `ROADMAP.md` a line dated today: "Spec 005 stage S5 implemented: resize paths keep the extent, cell-centre guard tests, non-mutating world-correct misfit diagnostics, residual and exclusion report." Report any failing gate and any file under `packages/loop_interpolation/tests/golden` that changed. Write in ASD-STE100 Simplified Technical English.

## Phase 23: Stage S6 tests and guards (FR-025 to FR-027, FR-030)

**Goal**: Fixture helpers, the transformation tests, the R2 validation tests, the faulted-feature residual tests, the backend boundary and optional cells, and the matrix helper are written before the code that they test. A short list of S5 minor findings is closed in the test files.
**Independent test**: `pytest tests/unit/modelling/test_fault_frame_validation.py tests/unit/modelling/test_diagnostics_backends.py tests/unit/interpolator/test_backend_world_contract.py tests/unit/test_backend_matrix.py` (root run) and `pytest packages/loop_common/tests/test_euclidean_transformation.py` (separate call).
**Environment (all S6 tasks that run tests)**: run pytest from PowerShell with the venv python `C:\Users\lachl\dev\LoopStructural\.venv\Scripts\python.exe -m pytest`. Set `PYTHONPATH` to `<worktree>;<worktree>\packages\loop_common\src;<worktree>\packages\loop_interpolation\src;<worktree>\packages\loopstructural_visualisation\src` (worktree first, because the editable install points to another worktree). `pyevtk` is an optional dependency: do not put it on the path. Run the root tests and the package tests in separate calls. Do not install anything into the venv. Plan changes in the task graph: the first S6 wave tasks `Need` T088 (the S5 end) and no S6 task has `Blocked by` (OP21 to OP24 have defaults).

- [x] T089 [P] [Found] [M] Add the S6 fixture helpers for the R2 displacement, the curved control and the faulted model (FR-026, FR-027, SC-004).
  - Files: `tests/fixtures/coordinate_frames.py`
  - Needs: T088
  - Wave: 43
  - Done when: `pytest tests/fixtures/coordinate_frames.py` collects 0 tests with no import error, `ruff check LoopStructural packages` passes, and the report gives the measured result of `r2_analytic_fault(anchor).apply_to_points(r2_displacement_section(anchor))` (shape, finite count, largest z shift).
  - Context: Edit the S1 helper module. Keep all old names. Add these names (all free; numpy docstrings, type hints, line length 100): constants `R2_SLIP = 3.0`, `R2_STEPS = 20`, `R2_DISPLACEMENT_TOL = 0.03`, `CURVED_RADIUS = 5.0`, `CURVED_GX = 0.2`.
    `r2_displacement_section(anchor) -> (1369, 3)`: world points with y' = 0 and x', z' in {-9, -8.5, ..., 9} (37 x 37). `r2_section_masks(anchor, section) -> (smooth, plane)`: smooth is `abs(x') >= 0.5` (1332), plane is the rest (37).
    `r2_analytic_fault(anchor, steps=R2_STEPS, displacement=R2_SLIP, gradient_override=None) -> FaultSegment`: a `FaultSegment` of three `AnalyticReference` features with the R2 values and gradients (gx = x'/20, gy = -z'/20, gz = y'/20; see `r2_frame`); default profile `BaseFault`. `gradient_override(coordinate, points, gradient)` may replace the gradient of one coordinate and gets the current points of each call. Run `apply_to_points` on this fault FIRST and report the result. If a `BaseFeature` method that `apply_to_points` needs is missing on `AnalyticReference`, add the smallest method to the fixture class, never to the library.
    `r2_closed_form_displacement(section, anchor, slip=R2_SLIP) -> (N, 3)`: by hand, no library profile: for x' > 0 the vector is `(0, 0, -0.5*slip*f(t))`, for x' < 0 it is `(0, 0, +0.5*slip*f(t))`, for x' = 0 zero; `t = abs(x')/20`, `f(t) = 1 - 3 t^2 + 2 t^3`.
    `r2_curved_slip_fault(anchor, steps) -> FaultSegment`: coordinate 0 constant `CURVED_GX`, coordinate 2 constant 0, coordinate 1 gradient `(-(y-cy), (x-cx), 0) / (20 r)` about the vertical axis through `anchor`; displacement `R2_SLIP`, profile `BaseFault`. By hand: `d = 0.5 * 3 * f(0.2) = 1.344` (`f(0.2) = 0.896`).
    `euler_circle_closed_form(radius, d, steps) -> (final_radius, turn_angle)`: step length `s = d/steps`; `r_k^2 = radius^2 + k s^2`; the angle gains `arctan(s/r_k)` in step k; sum with a plain loop of `steps` terms; no library call.
    `build_r2_fault(anchor, rotation_deg, interpolatortype, explicit_center=True, displacement=R2_SLIP) -> FaultSegment`: a copy of `_fault` in `tests/unit/modelling/test_fault_frame_analytic.py` (lines 21 to 51); only `displacement` (3, not 0) and `fault.steps = R2_STEPS` change. Do not edit the S1 test file.
    `faulted_foliation_model(anchor, rotation_deg, engine) -> (model, fault, strati)`: same data as `_faulted_built_model` in `test_diagnostics_frames.py` (unit cube): 3 fault rows at (0.5, 0.5, 0.5) (normal (1,0,0) value 0, slip direction (0,0,1), extent direction (0,1,0)); 8 strati value rows at z = 0.3 (value 0) and z = 0.7 (value 1) for x, y in {0.2, 0.8}; 1 strati normal row at (0.3, 0.3, 0.5) vector (0,0,1); fault displacement 0.1; every location shifted by `anchor`; `model_with_frame(anchor, anchor + 1, rotation_deg)`; both features use `interpolatortype=engine` ("FDI" or "PLI") and `nelements=1e3`; call `model.update()` before return.
    `residuals_of_field(interpolator, nodal_values) -> dict[str, ndarray]`: `matrix @ nodal_values - b` for each entry of `interpolator.constraints` (keys `matrix`, `b`, `w`). Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T090 [P] [US5] [M] Test: round trips and defect cases of `EuclideanTransformation` (FR-025, FR-029, SC-010; the tests fail first).
  - Files: `tests/unit/utils/test_transformation.py`, `packages/loop_common/tests/test_euclidean_transformation.py`
  - Needs: T088
  - Wave: 43
  - Done when: `pytest packages/loop_common/tests/test_euclidean_transformation.py` and `pytest tests/unit/utils/test_transformation.py` (two calls) fail before T100 and T102 (record the failing ids and the first error of each) and pass after T102.
  - Context: Root file: keep the 9 other existing tests; DELETE `test_inverse_transform_is_broken_for_normal_point_clouds_bug` (it keeps a defect; FR-029). `Tp(L, W)` is `tp` of `tests/fixtures/coordinate_frames.py` with `L = max(extent of the cloud, 1.0)`, `W = max abs coordinate`. Anchors: the four `ANCHORS`; cloud = `anchor + rng.uniform(-50, 50)` (seed 0). Every bound uses `Tp` or a value derived from it (at `W = 6.5e6` a float64 spacing is 9.3e-10).
    Root tests (cells): `test_round_trip_explicit_parameters` (dimensions {2,3} x n {1,4,50} x 4 anchors x columns {3,5} x angle {30, -75 deg} = 96): translation `anchor + (1,2,3)`; forward differs from input by more than 0.1 in a transformed column; x, y equal `x' = c(x-tx) - s(y-ty)`, `y' = s(x-tx) + c(y-ty)` within `1e-12 * L`; for `dimensions=3` z equals `z - tz`; `inverse_transform(forward)` equals the input within `Tp`; the other columns are `np.array_equal`; the input is unchanged; no shared memory. `test_round_trip_after_fit` (4 anchors x n {1,4,50} x `fit_rotation` {True, False} = 24; 2D; 50-point line at 30 deg; fitted translation equals `np.mean` within `Tp`). `test_integer_input_is_not_truncated` (int64 array of 4 points; float closed form within `Tp`). `test_one_dimensional_input_raises_value_error` (3: transform, inverse_transform, fit; shape `(3,)`; text "2D array"). `test_too_few_columns_raises_value_error` (2; `dimensions=3`, 2 columns; "at least 3 dimensions").
    `test_entry_points_are_one_class` (`LoopStructural.utils.EuclideanTransformation`, `LoopStructural.utils._transformation.EuclideanTransformation`, `loop_common.math.EuclideanTransformation`, `loop_common.math._transformation.EuclideanTransformation` are one object). `test_pickle_with_the_old_module_path_loads` (monkeypatch the class `__module__` to `"LoopStructural.utils._transformation"`, pickle a fitted object, restore, load; equal `translation`, `angle`, `dimensions`, `fit_rotation`). `test_new_pickle_names_the_loop_common_path` (bytes contain `loop_common.math._transformation`).
    `test_fit_aligns_a_line_with_negative_slope` (lines at -30 and -10 deg through (5,5), 50 points as `_line_points`): angle is `+pi/6` (0.5235987756) and `+pi/18` (0.1745329252) within 1e-9 modulo pi; y variance after `fit_transform` below 1e-20 (before: angle `-pi/6`, variance 6.51 and 1.01). `test_fit_three_dimensions` (`dimensions=3`, lines at +30 and -30 deg, z = 0.1 i): translation `(5, 5, 2.45)` within 1e-9; angle `-pi/6` and `+pi/6` within 1e-9 modulo pi; z' = z - 2.45 within 1e-12; y variance below 1e-20; inverse returns the input within `Tp` (before: `fit` raises `ValueError`).
    Package file (new; no `tests.fixtures` import, no scikit-learn; local `_tp(length, width) = 1e-10*length + 32*eps*max(width, length)`): `test_module_imports` (`importlib.import_module("loop_common.math._transformation")` works), `test_export_from_loop_common_math`, `test_round_trip_without_fit` (dimensions {2,3} x n {1,4,50}, anchor (500000, 6500000, 100) = 6), `test_two_dimensional_rotation_of_known_angle` (`dimensions=2`, angle pi/6, translation 0, point (2,0,7) gives (1.7320508076, 1.0, 7) within 1e-12; third column exactly 7), `test_three_dimensional_rotation_keeps_the_z_direction` (angle pi/2, translation (1,2,3): (2,2,5) gives (0,1,2); (1,3,3) gives (-1,0,0); within 1e-12; determinant of `rotation[:3,:3]` is 1; inverse returns inputs), `test_fit_without_scikit_learn_raises_import_error` (`monkeypatch.setitem(sys.modules, "sklearn", None)`; `ImportError` text contains "scikit-learn"; `translation` and `angle` unchanged), `test_fit_without_rotation_needs_no_scikit_learn` (same block; `fit_rotation=False` on `[[0,0,5],[2,4,6]]` gives translation (1,2), angle 0). Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T091 [P] [US5] [M] Change `test_diagnostics_frames.py` for the S5 minor findings 5 to 8.
  - Files: `tests/unit/modelling/test_diagnostics_frames.py`
  - Needs: T088
  - Wave: 43
  - Done when: `pytest tests/unit/modelling/test_diagnostics_frames.py` passes before and after T101 (the changed tests are guards), and the report gives the measured maximum of the solved residual for FDI and for PLI.
  - Context: Four changes. (5) The exact-field test patches `feature.interpolator.evaluate_gradient` in place of `feature.evaluate_gradient`, records the points, and asserts that they equal the world positions of the stored rows within `Tp`. For a feature with no faults `GeologicalFeature.evaluate_gradient` ends in `self.interpolator.evaluate_gradient`, so the test passes before and after T101. The old patch target would fail after T101, so this change must land first. No other test patches `evaluate_gradient` or `evaluate_value`.
    (6) `TOL_SOLVED` becomes `TOL_SOLVED_FDI` and `TOL_SOLVED_PLI`: measure first (run the cells, record the largest residual of each engine), set each constant to 2 times that maximum. Each constant must be at most the old `9e-4`; the S5 record has FDI 1.2e-4. If a bound would pass 9e-4, stop and report. Never raise a bound to make a test pass. (7) The PLI gradient test asserts `unattributed == 8` (it printed before). (8) The test that compares `id(latest_diagnostics_report)` holds the report object and compares with `is`.
    Do not change other tests. Tolerances: 0 relaxed. Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T092 [P] [Found] [M] Find the module name of the Surfe package and fix the `importorskip` name only if it is wrong.
  - Files: `packages/loop_interpolation/tests/test_p0_surfe_nans.py`
  - Needs: T088
  - Wave: 43
  - Done when: the report names the pip package, the module names found, and the decision ("changed to `surfepy`", "unchanged: module is `surfe`" or "not decided"), and `pytest packages/loop_interpolation/tests/test_p0_surfe_nans.py` runs with no error (it skips here).
  - Context: The adapter `_surfe_wrapper.py:7` imports `surfepy`. Line 6 of the test file is `pytest.importorskip("surfe", minversion=None)` (the pip name). Steps: (1) run `pip index versions surfe` with a pip that is NOT in the venv (the venv has none); (2) if it works, `pip download surfe --no-deps -d <scratch>` where `<scratch>` is a new folder under the scratchpad directory; (3) list the wheel files with `python -I -m zipfile -l <wheel>` and read the top-level module names. Decision: if the wheel holds the module `surfepy` and no module `surfe`, change the line to `pytest.importorskip("surfepy")`. If it holds `surfe`, or both, or the check cannot run (no pip or no network), do not change the file and report the names or "not decided". Do not install anything in the venv. Do not unpack the wheel in the repository. The result goes to T099 and to the cell `surfe-extension` (T103). Write in ASD-STE100 Simplified Technical English.

- [x] T093 [P] [Found] [S] Remove the byte order mark from two files (S5 minor findings 1 and 10).
  - Files: `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`, `tests/unit/modelling/test_fault_frame_analytic.py`
  - Needs: T088
  - Wave: 43
  - Done when: `python -I -c "import sys; sys.exit(any(open(p,'rb').read(3)==b'\xef\xbb\xbf' for p in sys.argv[1:]))" packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py tests/unit/modelling/test_fault_frame_analytic.py` exits 0, and `git diff --numstat` for the two files prints `1 1` for each.
  - Context: Both files start with the bytes `ef bb bf`. Remove those 3 bytes and nothing else. Do not use an editor tool that rewrites line endings: read the bytes, assert that they start with `ef bb bf`, write `data[3:]`. Check the line endings of each file stay as they were. Check that `git diff -w --stat packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py` against the S4 merge `62c420e8` shows no change other than the S5 changes. Write in ASD-STE100 Simplified Technical English.

- [x] T094 [P] [Found] [M] Write the backend-matrix helper and its tests (FR-030, SC-018; the tests fail first).
  - Files: `tests/fixtures/backend_matrix.py`, `tests/unit/test_backend_matrix.py`
  - Needs: T088
  - Wave: 43
  - Done when: `pytest tests/unit/test_backend_matrix.py` fails with an import error before `backend_matrix.py` exists (write the tests first), then passes; `MATRIX = ()` in this task.
  - Context: New helper module (`from __future__ import annotations`, `typing.Literal`, numpy docstrings, py39). Names: `Outcome = Literal["passed","failed","skipped","absent"]`; frozen dataclass `OptionalDependency(key, module, distribution, capability)`; `OPTIONAL_DEPENDENCIES` with keys `native` (module `loop_interpolation_core`, extra `loop-interpolation[native]`), `surfe` (module `surfepy`; pip name unverified, see T092), `viewer` (`loopstructuralvisualisation`), `matplotlib`, `pyevtk`; frozen `OptionalStatus(key, module, available, reason)` (reason is "" when available, else `"optional dependency '<module>' is not installed (<distribution>)"`); `optional_status(key, find_spec=importlib.util.find_spec) -> OptionalStatus`; `require_optional(key) -> None` (`pytest.skip(status.reason)` when not available); frozen `Cell(cell_id, spec_row, mandatory, stub, refs: tuple[tuple[str, str, str | None], ...], min_cases, dependency)` where `dependency: str | tuple[str, ...] | None` (a tuple holds several keys of `OPTIONAL_DEPENDENCIES`; the cell is "dependency absent" when ANY key is not available; `classify_cell(cell, outcome, status)` takes `status` as one `OptionalStatus`, a sequence of them, or None; the "excluded" reason then names the exact module of EACH absent key, joined with "; "; add one helper `cell_dependency_keys(cell) -> tuple[str, ...]`; `render_table(matrix, statuses)` takes `statuses: dict[str, OptionalStatus]`). Add to `test_classify_cell_rules` 2 cases with a tuple of two statuses (one absent: excluded when skipped, reason names only the absent module; both available and skipped: violation) (a ref is (path with "/", test function, id selector or None)); `SPEC_ROWS` (the 7 rows of "Required backend coverage" in `spec.md`, exactly these strings in this order: `"Python FDI and PLI/P1"`, `"Python P2"`, `"Python FDI-fold, PLI discrete-fold, and P2-fold"`, `"Python FDI and PLI constant-norm variants"`, `"Common 2D and 3D geometry"`, `"Existing native FDI or fold integration"`, `"Optional Surfe integration"`; the extra value `"extra"` is not in `SPEC_ROWS`); `MATRIX: tuple[Cell, ...] = ()`; `classify_cell(cell, outcome, status) -> str`; `outcomes_from_junit(paths) -> dict[str, Outcome]`; `rows_without_real_cell(matrix) -> list[str]`; `render_table(matrix, statuses) -> str`; `main(argv=None) -> int` (prints the table, returns 1 on a violation).
    `classify_cell` returns `"passed"`, `"failed"`, `"excluded: <reason>"` or `"VIOLATION: <text>"`. A failed cell is "failed" in all cases. A mandatory cell that is skipped or absent is a violation. An optional cell with an absent dependency is "excluded" when skipped, a violation when passed or absent. An optional cell with an installed dependency is a violation when skipped or absent.
    `outcomes_from_junit` reads every file and merges cases by key (`classname`, `name`); same key in two files: the worse outcome wins (failed, skipped, passed). Ref match: module = path without `.py`, "/" to "."; `classname` equals the module or starts with module + "."; the part of `name` before `[` equals the function; the selector is None or occurs in `name`. A module-level skip is one skipped case with `classname=""` and `name` equal to the module: it matches every ref of that file as "skipped" (for an optional cell with an absent dependency this gives "excluded", never "absent"; for a mandatory cell a violation). Cell outcome: one failed case makes it failed; else one skipped case makes it skipped; else at least `min_cases` passed cases make it passed; else absent.
    Tests now (helper only; the MATRIX tests come in T103): `test_optional_status_names_exact_module` (5; fake `find_spec` returning None: not available and the reason holds the exact module name; returning an object: available, reason empty), `test_require_optional_skips_with_named_dependency` (2: native, surfe; message holds `loop_interpolation_core` or `surfepy`), `test_classify_cell_rules` (12 combinations: mandatory / optional absent / optional installed x passed, failed, skipped, absent), `test_outcomes_from_junit` (6: the four outcomes from a small XML text; fewer than `min_cases` gives absent; a class test `classname="a.b.test_x.TestC"` matches a ref of `a/b/test_x.py`; two files merged; a module-level skip gives skipped for every ref of the file, "excluded" for an optional cell with an absent dependency and a violation for a mandatory cell). The command `python -m tests.fixtures.backend_matrix <junit>...` must work from the repository root. Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T095 [P] [US5] [M] Test: the R2 fault frame against the analytic frame, in displacement and in constraint residuals (FR-026, FR-027, SC-004, SC-003).
  - Files: `tests/unit/modelling/test_fault_frame_validation.py`
  - Needs: T089
  - Wave: 44
  - Done when: `pytest tests/unit/modelling/test_fault_frame_validation.py` passes on the current tree (these are guards, OP18 and OP22 closed; record the run time; about 26 solves), and the report gives the measured `TAU_ASSEMBLED` and `TAU_ORTHOGONAL`.
  - Context: New file; import helpers from `tests.fixtures.coordinate_frames` (T089). Module-scoped fixture `solved_r2` with `params=CELLS` (4 anchors x 3 rotations x 2 engines FDI/PLI = 24, explicit centre) builds one fault with `build_r2_fault`; a second module-scoped fixture `solved_r2_inferred` has 2 inferred-centre cells (anchor 3, rotation 30, FDI and PLI). Only these 26 solves; no test changes a model. Section = `r2_displacement_section(anchor)`; analytic frame = `r2_analytic_fault(anchor)`.
    Tests (24 cells unless noted): `test_vector_displacement_matches_analytic_frame` (`apply_to_points` on both; 1332 smooth points finite; largest row norm of the difference at most 0.03; non-vacuity: analytic norm at least 0.5 for at least 1200 of 1332; by hand `d = 1.5 f(t)`, `t` in [0.025, 0.45], so `d >= 0.86`); `test_scalar_displacement_matches_analytic_frame` (`evaluate_displacement`, same bound and check); `test_displacement_populations_are_counted_separately` (1369 = 1332 smooth + 37 plane + 0 invalid; counts from the generator; solved output finite on all 1369; the largest plane difference is printed, asserted finite only; precondition: the support holds the path, largest `abs(z')` is 10.5, else fail with a message); `test_vector_displacement_inferred_centre` (2 cells, as the first, `explicit_center=False`).
    `test_analytic_frame_displacement_matches_closed_form` (4 anchors): displacement minus section equals `r2_closed_form_displacement`; x and y exactly 0; z within `64 * eps * max(1, max abs z)`; 37 plane points move by exactly 0. `test_integration_error_follows_the_euler_closed_form` (4 anchors): `r2_curved_slip_fault(anchor, N)`, 8 start points on radius 5 at k*45 deg; `d = 1.344`; (a) N = 20 equals `euler_circle_closed_form` within `Tp(10, W)`; (b) N = 400 differs from the exact arc (radius 5, angle d/5) by at most `d^2/(N*5) + Tp(10, W)` (9.0e-4 + Tp); (c) the N = 20 and N = 400 results differ by more than 5e-3 (radial gap 8.6e-3). `test_misplaced_frame_fails_the_displacement_bound` (24): analytic frame at `anchor + (1,0,0)`: largest smooth difference above 0.03 (by hand about 0.11).
    `test_outside_support_points_have_no_displacement` (24): six probes `anchor +- (H+5) e_axis`, `H` = largest world offset of the 8 support corners on that axis; precondition: outside the support; value and gradient NaN; `evaluate_displacement` is 0.0; `apply_to_points` returns the points unchanged (`np.array_equal`); if `apply_to_points` raises on an empty mask, report the exception and test `evaluate_displacement` only. `test_zero_gradient_protection_remains_effective` (4 anchors x {zero, nan} = 8): `gradient_override` sets the coordinate-1 gradient to (0,0,0) or NaN for current `z' > 5`; sets by start position: A `z0' > 5`, B `z0' <= 3.5`, C in between. Zero variant under `warnings.simplefilter("error")`: no warning, all finite, A equals start, B equals the no-override result, C finite only. NaN variant: A is NaN, B finite and equal to the no-override result, C not asserted.
    `test_fault_frame_assembled_residuals_vanish_on_the_analytic_field` (24): per coordinate, `residuals_of_field` with nodal values = analytic field at `local_to_world_points(support.nodes)`. Group 1 must vanish (max at most `TAU_ASSEMBLED`; measure first, 2 x max, cap 1e-9): FDI `value`, `gradient`, `norm_0`, `norm_1`, `norm_2`, `dxx`, `dyy`, `dzz`, `dxy`, `dxz`, `dyz`, `dx_lower`, `dx_upper`, `dy_lower`, `dy_upper`, `dz_lower`, `dz_upper`; PLI `value`, `norm`, `edge jump` (a numeric suffix is allowed). Group 2, reported with `TAU_ORTHOGONAL` (measure first, cap 1e-3): names ending `_orthogonal`. Group 3: any other name fails the test and is named. Non-vacuity: each coordinate has at least 1 row in group 1; write the names and row counts as literals.
    `test_fault_frame_physical_residuals_match_analytic` (24): `fault[i].verification_report().physical_residuals`; coordinate 0 `value` has 25 rows, max at most 1e-6 and equals `abs(0 - fault[0].evaluate_value(rows))` within `max(1e-12, Tp(20, W)/20)`; every direction family largest angle at most 0.01 degree, and equal within 1e-6 degree to the independent `degrees(arctan2(norm(cross(a,b)), dot(a,b)))`. `test_residual_check_detects_misplaced_reference` (24): the field shifted 1 unit along x' gives a value maximum above 1e-2 (by hand 0.05 at 25 rows). If a measured bound would pass its cap, stop and report. Never raise a cap. Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T096 [P] [US5] [M] Test: residual diagnostics of faulted features and P2 state (FR-023, FR-026; tests 1, 2 and 4 fail first).
  - Files: `tests/unit/modelling/test_diagnostics_backends.py`
  - Needs: T089
  - Wave: 44
  - Done when: `pytest tests/unit/modelling/test_diagnostics_backends.py -k "restored_positions or independent_rms"` fails before T101 (record the failing ids and the residual, about 0.08) and passes after T101; the other tests pass now. The report gives the measured `TOL_FAULTED_VALUE` and `TOL_FAULTED_COS`.
  - Context: New file. Two module-scoped parametrised fixtures: `faulted_model` (24 cells = 4 anchors x 3 rotations x 2 engines, `faulted_foliation_model`) and `p2_model` (12 cells = 4 anchors x 3 rotations, `world_foliation_model`). 36 models in all. No test changes a model.
    Tests: `test_faulted_value_residuals_use_restored_positions` (24): `strati.evaluate_value_residuals()` equals `abs(val - interpolator.support.evaluate_value(stored local positions, interpolator.c))` within `max(1e-12, 3*Tp(1, W))`; at most `TOL_FAULTED_VALUE` (measure first with the independent function, 2 x max, cap 1e-2; S5 measured 0.0004 to 0.0018; 0.080 before the fix); a spy (monkeypatch) counts calls of `strati._apply_faults`, `strati.evaluate_value`, `strati.evaluate_gradient`: each count 0. `test_faulted_gradient_cosines_use_restored_positions` (24): `evaluate_gradient_cosines()` equals the independent cosines (`interpolator.evaluate_gradient` at the stored world positions against the stored world vectors) within `max(1e-12, 10*Tp(1, W))`; `max(1 - cos)` at most `TOL_FAULTED_COS` (measure first, cap 1e-2); same spy. `test_fault_restores_the_stored_rows` (24; guard): at least 4 stored rows differ from the raw world data by more than 1e-3; write the measured number as a literal. `test_report_physical_value_matches_independent_rms` (24): `verification_report().physical_residuals["value"].rms` equals the RMS of the independent residuals within `max(1e-12, 3*Tp(1, W))`; `.kind == "physical"`.
    `test_p2_misfit_leaves_state_unchanged` (12; guard): the S5 snapshot of `data`, `c`, flags, weights and the report object (compare with `is`) is equal after the two old misfit methods, the two new methods and `verification_report()`. `test_p2_report_blocks_are_separate` (12; guard): the physical and assembled blocks are different objects with `kind` "physical" and "assembled"; the physical `value` rms equals the RMS of `evaluate_value_residuals()` within 1e-12. The P2 fit is not tested against the analytic field (issue #399; OP22). Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T097 [P] [US5] [M] Test: the native and Surfe adapters send local data and return world results, with recording stubs (FR-030, SC-018, S2 T020).
  - Files: `tests/unit/interpolator/test_backend_world_contract.py`
  - Needs: T089
  - Wave: 44
  - Done when: `pytest tests/unit/interpolator/test_backend_world_contract.py` passes with 0 skips (guards; the adapters are not changed) and the report gives the measured maximum world error of the Surfe cells.
  - Context: New file; mandatory cells `native-boundary` and `surfe-boundary` (stub cells; they test the Python side only). No skip call in the file. Fixture `stub_core` puts a recording module `loop_interpolation_core` (class `Problem`) in `sys.modules` with `monkeypatch.setitem`; `Problem.solve` returns the exact R1 nodal field and `{"converged": True, "iterations": 1, "residual": 0.0}`. Fixture `stub_surfe` puts a recording `surfepy` module (class `Surfe_API`) in `sys.modules`, removes `loop_interpolation._surfe_wrapper` from `sys.modules` before the import, and at teardown restores `sys.modules["surfepy"]`, `sys.modules["loop_interpolation._surfe_wrapper"]` and the attribute `_surfe_wrapper` of the package `loop_interpolation` (absent stays absent) and asserts it.
    Tests: `test_native_adapter_sends_local_rows_and_grid` (12 = 4 anchors x 3 rotations): R1 FDI interpolator (`r1_interpolator("fdi", anchor, rotation, exact=True)`) plus 4 gradient rows in world coordinates; `_native.solve` runs with the stub; `add_value` points equal `bounding_box.project(world observations)` within `Tp`; `add_gradient` vectors equal `world_to_local_vectors(unit R1 gradient)` within 1e-12; `Problem` gets `support.origin`, `support.nsteps`, `support.step_vector`; for rotation 30 and 90 the local points differ from the world points by more than 1. `test_native_solved_field_evaluates_in_world` (12): after `solve`, `interp.c` is the stub field; `evaluate_value(world)` equals `r1_value` within `1e-8*F`; `evaluate_gradient` equals the R1 gradient within `1e-8*max(1, abs(g))`. `test_native_unsupported_settings_keep_rejection_and_state` (3: rotated support `rotation_xy`, region, `P1Interpolator`): `NotImplementedError` whose text holds `rotation_xy`, `region` or `interpolator class`; `c`, `data`, `up_to_date` equal before and after.
    `test_surfe_adapter_sends_local_points` (12): the Surfe wrapper with `interp.bounding_box = box`; `setup_interpolator()` sends the 27 interface points equal to `project(world observations)` within `Tp`; planar rows have local vectors. `test_surfe_adapter_returns_world_values_and_gradients` (12): the stub's `EvaluateInterpolantAtPoints(local)` returns `r1_value(interp.local_to_world_points(local))` and `EvaluateVectorInterpolantAtPoints(local)` returns `world_to_local_vectors(R1 gradient)`; public `evaluate_value(world)` equals `r1_value` within `1e-8*F` and `evaluate_gradient` equals the R1 gradient within `1e-8*max(1, abs(g))` (SC-002). Two conversions at size 6.5e6 move a position by about 1e-9, a value by about 5e-10; `F` is about 5 or more, so the margin is about 100. Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T098 [P] [US5] [M] Test: the optional native and Surfe extension cells (FR-030, SC-018, OP21).
  - Files: `tests/unit/interpolator/test_backend_optional_native_surfe.py`
  - Needs: T089, T094
  - Wave: 44
  - Done when: `pytest tests/unit/interpolator/test_backend_optional_native_surfe.py -rs` shows 24 skipped cases, each with the text "optional dependency '<module>' is not installed" (the extensions are absent here), and no error.
  - Context: New file; optional cells `native-extension` and `surfe-extension`. Each test starts with `require_optional("native")` or `require_optional("surfe")` from `tests.fixtures.backend_matrix` (T094). `test_native_extension_r1_fit_in_world` (12 = 4 anchors x 3 rotations): R1 fit (`exact=False` support, 27 value rows) solved with `solver="native"`, `solver_kwargs={"tol": 1e-12}`; world value error at most `1e-8*F` and gradient error at most `1e-8*max(1, abs(g))` at the two query offsets (SC-002); the result equals the `cg` result within `1e-6*F`. `test_surfe_extension_r1_values_and_gradients_in_world` (12): Surfe wrapper with `poly_order=1` and the default kernel; 27 interface rows; same SC-002 limits. The extensions are not installed, so the tolerances are not measured (OP21): never relax a limit; a miss is a fit-accuracy finding. Do not stub the extension in this file. Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

- [x] T099 [US5] [M] Record the S6 pre-fix evidence in `verification.md` (FR-031).
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T090, T091, T092, T093, T095, T096, T097, T098
  - Wave: 45
  - Done when: `verification.md` has the heading `## S6 pre-fix failures` with the commands, the failing ids and the first error of each, and the heading has a line for each guard test that passes before the fix.
  - Context: Run the new tests on the tree before T100 and T101 (see "Environment" in Phase 23; root and package runs in separate calls). Expected to fail: the transformation tests (the old `inverse_transform` fails for `dimensions=2` with 3 or 5 columns, for `dimensions=3` with 5 columns, and with 3 columns for N > 3; every `dimensions=3` cell fails the z check), the fit defect tests (negative slope, 3D fit, the z direction, missing scikit-learn), `test_entry_points_are_one_class`, and tests 1, 2 and 4 of `test_diagnostics_backends.py`. Expected guards that pass: the R2 displacement and residual tests, the P2 tests, the restored-rows test, the stub cells, the changed exact-field test. Record the T092 result (pip name, module names, decision). Record that `test_fault_frame_analytic.py` and `_structured_grid_3d.py` lost their BOM (T093). Record the measured tolerances: `TOL_SOLVED_FDI`, `TOL_SOLVED_PLI`, `TAU_ASSEMBLED`, `TAU_ORTHOGONAL`, `TOL_FAULTED_VALUE`, `TOL_FAULTED_COS`. Record the deleted test `test_inverse_transform_is_broken_for_normal_point_clouds_bug`. Write only a new section; do not change other text. Write in ASD-STE100 Simplified Technical English.

## Phase 24: Stage S6 code (FR-025, FR-026)

**Goal**: One working `EuclideanTransformation` behind both import paths, and residual diagnostics of faulted features that evaluate the interpolator at the restored positions.
**Independent test**: `pytest tests/unit/utils/test_transformation.py tests/unit/modelling/test_diagnostics_backends.py` (root run) and `pytest packages/loop_common/tests/test_euclidean_transformation.py`.

- [ ] T100 [P] [US5] [M] Fix `EuclideanTransformation` in `loop_common` and export it (FR-025, OP20, OP24).
  - Files: `packages/loop_common/src/loop_common/math/_transformation.py`, `packages/loop_common/src/loop_common/math/__init__.py`
  - Needs: T099
  - Wave: 46
  - Done when: `pytest packages/loop_common/tests/test_euclidean_transformation.py` passes (separate call), and `ruff check LoopStructural packages` passes. The root test `tests/unit/utils/test_transformation.py` passes after T102.
  - Context: Keep `from __future__ import annotations`. `__init__`, `rotation`, `inverse_rotation`, `fit_transform`, `__call__` and `_repr_html_` do not change. Replace `from . import getLogger` by `from loop_common.logging import get_logger` and `logger = get_logger(__name__)` (no import cycle). New `_import_decomposition() -> module`: returns `sklearn.decomposition`; raises `ImportError` whose message names "scikit-learn" and "pip install scikit-learn". New method `_check_points(self, points) -> np.ndarray`: `array = np.array(points, dtype=float)`; `ndim != 2` raises `ValueError("points must be a 2D array with shape (N, d)")`; `shape[1] < self.dimensions` raises `ValueError(f"Points must have at least {self.dimensions} dimensions")`; return the copy.
    `fit(points)`: `array = self._check_points(points)`; `d = self.dimensions`; `translation = np.mean(array[:, :d], axis=0)`; `angle = 0`; if `self.fit_rotation`: `decomposition = _import_decomposition()`, `centred = array[:, :2] - translation[None, :2]`, `axis = decomposition.PCA(n_components=2).fit(centred).components_[0]`, `angle = -np.arctan2(axis[1], axis[0])`; then set `self.translation`, `self.angle`; `logger.debug(...)`; return `self`. No attribute changes if it raises. Fewer than 2 rows (n = 1): PCA cannot fit, so `angle = 0` and the PCA is not called (settled at S6 build start; test `test_round_trip_after_fit` has n = 1 with `fit_rotation=True`). `_rotation(angle)`: 3x3 about z with entry `[2, 2] = +1` (it was -1; determinant 1).
    `transform(points)`: `array = self._check_points(points)`; `centred = array[:, :d] - self.translation[None, :]`; `array[:, :d] = np.einsum("ik,jk->ij", centred, self.rotation[:d, :d])`; return `array`. `inverse_transform(points)`: `array[:, :d] = np.einsum("ik,jk->ij", array[:, :d], self.inverse_rotation[:d, :d]) + self.translation[None, :]`; return `array`. Columns after `d` stay bit for bit. The 2D numbers do not change. `math/__init__.py`: add `from ._transformation import EuclideanTransformation` after `from ._maths import *`. No new dependency in `packages/loop_common/pyproject.toml`. numpy docstrings (numpy style) for all new items. Write in ASD-STE100 Simplified Technical English.

- [ ] T101 [P] [US5] [M] Fix the physical residual methods of faulted features to evaluate the interpolator (FR-023, FR-026, OP19).
  - Files: `LoopStructural/modelling/features/_geological_feature.py`
  - Needs: T091, T096, T099
  - Wave: 46
  - Done when: `pytest tests/unit/modelling/test_diagnostics_backends.py tests/unit/modelling/test_diagnostics_frames.py` passes (the tests that failed in T099 pass now), and `ruff check LoopStructural packages` passes.
  - Context: Change two lines and two docstrings, nothing else. In `_gradient_cosines_fresh` (about line 333) replace `model = self.evaluate_gradient(positions)` by `model = self.interpolator.evaluate_gradient(positions)`. In `_value_residuals_fresh` (about line 362) replace `self.evaluate_value(positions)` by `self.interpolator.evaluate_value(positions)`: `return np.abs(rows[:, 3] - self.interpolator.evaluate_value(positions))`. Add this sentence to both docstrings: "The stored positions are already restored by the faults of the feature. The method evaluates the interpolator, not the feature, so no fault acts a second time." The `positions` stay `bounding_box.reproject(...)` of the stored rows when the interpolator has a frame. Do not change `_stale_reason`, `_physical_residual_summaries`, `verification_report` or the public methods. The public `evaluate_gradient` of a faulted feature keeps its tetrahedron path and its shape error (known issue). Before the edit, grep `LoopStructural packages` for other readers of the two private methods. Also run `pytest tests/unit/modelling/test_scalar_field_grid.py`. Write in ASD-STE100 Simplified Technical English.

- [ ] T102 [P] [US5] [S] Re-export the `loop_common` class from the old module path (FR-025).
  - Files: `LoopStructural/utils/_transformation.py`
  - Needs: T100
  - Wave: 47
  - Done when: `pytest tests/unit/utils/test_transformation.py` passes (root call), and `ruff check LoopStructural packages` passes.
  - Context: Replace the whole content of the file by exactly:
    ```
    """Old import path of ``EuclideanTransformation``.

    The class lives in ``loop_common.math``. This module keeps the path
    ``LoopStructural.utils._transformation`` working.
    """
    from loop_common.math._transformation import EuclideanTransformation

    __all__ = ["EuclideanTransformation"]
    ```
    `LoopStructural/utils/__init__.py` line 15 already exports the name from this module; do not change that file. Rule: keep both import paths (FR-025). Write in ASD-STE100 Simplified Technical English.

- [ ] T103 [P] [US5] [M] Write the `MATRIX` list and the matrix tests (FR-030, SC-018).
  - Files: `tests/fixtures/backend_matrix.py`, `tests/unit/test_backend_matrix.py`
  - Needs: T090, T091, T094, T095, T096, T097, T098, T100, T101
  - Wave: 47
  - Done when: `pytest tests/unit/test_backend_matrix.py` passes with 0 skips.
  - Context: Fill `MATRIX` (15 cells; 11 mandatory) and add tests to `test_backend_matrix.py`. Fill `refs`, id selectors and `min_cases` from `pytest --collect-only -q` of the test files; selectors are literals that occur in the collected ids. Cells (id, spec row, level, evidence): `py-fdi` and `py-pli` (Python FDI and PLI/P1; mandatory; `test_frame_persistence.py`, `test_frame_copy_and_lifecycle.py`, `test_fault_frame_analytic.py`, `test_fault_frame_validation.py`, `test_reference_constraint_frames.py`, `test_fault_buffer_bounds.py`; FDI or PLI ids); `py-p2` (Python P2; mandatory; R1 exact field in the persistence and lifecycle files, `test_reference_constraint_frames.py` P2 ids, `test_diagnostics_backends.py` P2 cells); `py-fold-fdi`, `py-fold-pli`, `py-fold-p2` (`test_fold_constraint_frames.py` family ids); `py-cn-fdi`, `py-cn-pli` (`test_constant_norm_frames.py`, lifecycle files; ids `cn_fdi`, `cn_pli`); `common-geometry` (`test_rigid_transform_contract.py`, `test_resize_paths.py`, `test_cell_centre_sampling.py`, `test_euclidean_transformation.py`, `tests/unit/utils/test_transformation.py`).
    `native-boundary` and `surfe-boundary`: mandatory, `stub=True`, no dependency, `test_backend_world_contract.py` (native tests or Surfe tests). `native-extension`: optional, dependency `native`, `test_backend_optional_native_surfe.py`, the native files (`test_native_parity.py`, `test_native_inequality.py`, `test_native_orthogonal_and_limits.py`, `test_native_warm_start_ridge.py`, `test_native_no_rows.py`, `test_native_errors.py`, `tests/unit/modelling/test_native_fault_frame.py`, `tests/unit/modelling/test_builder_solver.py` native test); native R3 in a rotated frame is "R3 not covered". `surfe-extension`: optional, dependency `surfe`, the optional test file, `test_surfe_rbf_interpolator.py`, `test_p0_surfe_nans.py`. `display-viewer` (row "extra", optional, dependencies `viewer` and `matplotlib`; `test_world_location_consumers.py::test_2d_viewer_fault_ellipse_centre_is_world`) and `export-pyevtk` (row "extra", optional, `pyevtk`; `::test_internal_callers_do_not_warn`). The 11 mandatory cells include the 2 stub cells. A cell may hold several dependency keys (settled at S6 build start, see "Answers" of plan.md): `display-viewer` has `dependency=("viewer", "matplotlib")`. The viewer test imports `loopstructuralvisualisation._2d_viewer`; with the viewer package source on `PYTHONPATH` it passes (the run environment puts `packages/loopstructural_visualisation/src` on the path), and that is a correct "passed" for an optional cell with an installed dependency. The matrix has 15 cells (the task title said 14; the list of this task and the plan table give 15): 11 mandatory (`py-fdi`, `py-pli`, `py-p2`, `py-fold-fdi`, `py-fold-pli`, `py-fold-p2`, `py-cn-fdi`, `py-cn-pli`, `common-geometry`, `native-boundary`, `surfe-boundary`) and 4 optional (`native-extension`, `surfe-extension`, `display-viewer`, `export-pyevtk`). Tests assert these literal counts. Add a test `test_cell_with_two_dependencies` (the `display-viewer` cell is excluded when `matplotlib` is absent even if the viewer is present, and the reason names `matplotlib`).
    Tests to add: `test_every_spec_row_has_a_cell` (7: `rows_without_real_cell(MATRIX)` is empty; every cell row is in `SPEC_ROWS` or "extra"; removing the two real cells of rows 6 and 7 from a copy returns exactly those two rows), `test_stub_cells_are_mandatory_and_listed_apart` (every stub cell is mandatory with no dependency; `render_table` prints rows 6 and 7 from their real cells only and the stub cells in a separate block headed "Boundary stub cells (Python side; not a native or Surfe result)"; with "stub passed, extension skipped, dependency absent" row 6 reads "excluded" and the stub line reads "passed"), `test_mandatory_cells_name_existing_tests` (11; each ref is a file with a top-level function of that name; AST parse, no import), `test_mandatory_cells_have_no_skip` (11; no `skip`, `skipif`, `importorskip`, `xfail` or `require_optional` in the referenced functions, the module-level statements and `pytestmark` of each file, `tests/conftest.py`, every `tests/fixtures/*.py` and the `conftest.py` of the folder of the ref, except the definition of `require_optional`; a temporary file with a skip makes the check fail). Make `main` print the note that stub cells are not a native or Surfe result. Rule: tests must show a wrong answer (principle V). Write in ASD-STE100 Simplified Technical English.

## Phase 25: Stage S6 finish (FR-029, SC-014, SC-015, SC-018)

- [ ] T104 [US5] [M] Record the S6 post-fix results in `verification.md`.
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T100, T101, T102, T103
  - Wave: 48
  - Done when: `verification.md` has the heading `## S6 post-fix results` with the commands, the counts and the run time of the new S6 test files, and each failing id of T099 is shown as passing.
  - Context: Run, with the "Environment" of Phase 23: `pytest tests/unit/utils/test_transformation.py tests/unit/modelling/test_fault_frame_validation.py tests/unit/modelling/test_diagnostics_backends.py tests/unit/modelling/test_diagnostics_frames.py tests/unit/interpolator/test_backend_world_contract.py tests/unit/interpolator/test_backend_optional_native_surfe.py tests/unit/test_backend_matrix.py` (root call) and `pytest packages/loop_common/tests/test_euclidean_transformation.py` (separate call). Record the run time of each file (OP2 stays closed: the default run includes them; `test_fault_frame_validation.py` does 26 R2 solves). List the id of each T099 failure with its new result. Write only a new section; do not change other text. Write in ASD-STE100 Simplified Technical English.

- [ ] T105 [US5] [M] Do the S6 golden review, tolerance audit, coverage ledger and findings disposition in `verification.md` (FR-029, SC-014, SC-015).
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T104
  - Wave: 49
  - Done when: `verification.md` has the headings `## S6 golden review (C-004)`, `## S6 behaviour-change ledger`, `## S6 tolerance audit`, `## S6 coverage ledger` and `## S6 findings disposition`, and states 0 relaxed tolerances and 0 unaffected baseline files changed.
  - Context: Golden list: `git diff --name-status ebc7c699d5 HEAD -- . ':!specs' ':!*.py' ':!*.md'`. Expected: `packages/loop_common/tests/data/structured_support_reference.npz` (S1: 27 values; S5: 78 keys) and `packages/map2loop/pyproject.toml` (commit `a0e82096`, pull request #395, not a spec 005 change). Check that these baseline files are unchanged: `structured_grid_2d_reference.npz`, `tests/fixtures/api_surface_snapshot.json`, `tests/unit/interpolator/*.txt`, and everything under `packages/loop_interpolation/tests/golden/golden_files`. S6 expects no golden change. If a number changes, the S5 order applies: list, review, record, then `--update`.
    Behaviour-change ledger (no golden file; analytic evidence): faulted-feature residuals (evidence: the independent value in `test_diagnostics_backends.py`) and the four `EuclideanTransformation` fixes (evidence: hand-computed numbers in `test_euclidean_transformation.py` and `test_transformation.py`); state that `git grep` finds no other caller, so 0 baseline files are affected. Tolerance audit: `git diff -U0 ebc7c699d5 HEAD -- ':!specs' | grep -iE '^-.*(tol|tp|tau|approx|abs=|<= *[0-9.]+e-)'`, then the same with `^+`; pair each removed line with its added line; list and explain each; expected relaxed count 0; name the rename `TOL_SOLVED` to `TOL_SOLVED_FDI` and `TOL_SOLVED_PLI` (at most 9e-4, a tightening or no change). Record the deleted defect test.
    Coverage ledger: `git diff --name-only ebc7c699d5 HEAD -- LoopStructural packages` without test paths; map each production file to a test (SC-014). Findings disposition: 1 and 10 BOM closed (T093); 2 (format hunks in `_structured_grid_2d.py`) stays, cosmetic; 3 (`set_nelements` passes a cell count to the node setter) stays, old behaviour; 4 already closed (sentence at `verification.md` line 910; record the check); 5 to 8 closed (T091); 9 closed in code (`_first_stale` follows faults recursively, commit `3f20a014`; no test); 11 closed in T106; 12 closed here: recompute the run-time count sentence from the two run records and correct it; the major finding closed (T101). Also closed: S2 T020 (native boundary cells, T097) and S3 T038 (slip 3 in the R2 fixture, T089). Stays: the S1 note on `_on_geometry_change` anisotropy rows, and the S3 note on the tetrahedron gradient path of a faulted feature. Write in ASD-STE100 Simplified Technical English.

- [ ] T106 [Finish] [M] Record the S6 notices in `COMPAT.md` (FR-025, FR-029, OP19, OP20).
  - Files: `COMPAT.md`
  - Needs: T105
  - Wave: 50
  - Done when: `grep -c "loop_common.math.EuclideanTransformation" COMPAT.md` prints at least 1, `grep -c "_value_residuals_fresh" COMPAT.md` prints at least 1, and the sentence "With an identity frame and no resize, the value misfit is the same as before" is in the `evaluate_value_misfit` row only.
  - Context: Read `COMPAT.md` first and copy the layout of the S1 to S5 notices. Add: (1) `EuclideanTransformation` is one class. `loop_common.math.EuclideanTransformation` is new; `LoopStructural.utils.EuclideanTransformation` is the same object; `inverse_transform` returns the original points for an `(N, d)` array (before: `ValueError` or too few rows); the methods return a float copy (an integer input was truncated); a 1-D input raises `ValueError` (before: `IndexError`). (2) The 3D behaviour changes: for `dimensions=3` z is no longer reversed (`z' = z - tz`; before `-(z - tz)`); `fit` works for `dimensions=3` (before it failed); the fitted angle has the right sign for a main axis with a negative slope (before it was wrong for lines with a direction between -45 and 0 degrees); a 3D user must check every result of the old class; `fit` with `fit_rotation=True` raises `ImportError` when scikit-learn is missing (before it logged an error and returned `None`); `fit` with `fit_rotation=False` needs no scikit-learn. (3) `fit` logs through the `loop_common` logger (`loop_common.logging.get_logger`). (4) A pickle made before S6 still loads (the old module path stays importable and exposes the class). A pickle made after S6 stores `loop_common.math._transformation`; a release from before S6 cannot load it.
    (5) `evaluate_gradient_misfit`, `evaluate_value_misfit`, `evaluate_gradient_cosines`, `evaluate_value_residuals` and `verification_report` of a faulted feature compare the stored restored rows with the interpolator at the same positions (before, the faults acted a second time). `_value_residuals_fresh` and `_gradient_cosines_fresh` now evaluate the interpolator directly; this reverses S5 decision 8 for these two methods. Numbers for a feature with no faults and no regions do not change. For a feature with a region mask, a row outside the region gave NaN before and gives a number now. (6) Known issue, not fixed in S6: the public `evaluate_gradient` of a faulted feature keeps its tetrahedron path and the shape error that S3 noted. (7) Move the S5 sentence "With an identity frame and no resize, the value misfit is the same as before" from the `ConstraintDiagnosticsReport` row to the `evaluate_value_misfit` row (S5 finding 11). Also note the stays of S5 finding 3 (`set_nelements`). Do not change other text. Write in ASD-STE100 Simplified Technical English.

- [ ] T107 [Finish] [M] Run the S6 final gates, the backend matrix, the rollback check and the exclusion ledger; complete `verification.md` (FR-030, SC-014, SC-018).
  - Files: `specs/005-coordinate-frame-correctness/verification.md`
  - Needs: T106
  - Wave: 51
  - Done when: `verification.md` has the headings `## S6 final gates`, `## S6 backend matrix`, `## S6 exclusion ledger`, `## S6 rollback checks`, `## S6 test run time (OP2)` and a line for each of OP21 to OP24, and `python -m tests.fixtures.backend_matrix <junit>...` exits 0.
  - Context: Run one process at a time from the repository root, with the "Environment" of Phase 23 (no `pyevtk` on the path; the venv python). Each pytest run writes its own junit file in a new scratch folder `<junit_dir>` outside the repository: `ruff check LoopStructural packages`; `pytest --junitxml=<junit_dir>/root.xml`; `pytest packages/loop_common/tests --junitxml=<junit_dir>/common.xml`; `pytest packages/loop_interpolation/tests --junitxml=<junit_dir>/interpolation.xml`; `pytest tests/unit/test_public_api_contract.py tests/unit/test_stable_api_surface.py`; `python -m tests.fixtures.backend_matrix <junit_dir>/root.xml <junit_dir>/common.xml <junit_dir>/interpolation.xml`. First check one `classname` of each junit file: it must start with `tests.`, `packages.loop_common.tests.` or `packages.loop_interpolation.tests.`.
    Optional dependencies (settled at S6 build start): `python -m tests.fixtures.backend_matrix` must run with the SAME `PYTHONPATH` as the three pytest runs (its `optional_status` uses `find_spec` in its own process), so the statuses agree with the runs. `pyevtk` is NOT on the path: `export-pyevtk` must be "excluded: optional dependency 'pyevtk' is not installed". `loopstructuralvisualisation` IS on the path through `packages/loopstructural_visualisation/src` and `matplotlib` is installed: `display-viewer` must be "passed" if the viewer test passes (this is correct for an installed dependency); if it is skipped although the dependency is found, that is a VIOLATION: report it, do not hide it. In `verification.md` list each optional cell with its classification and the exact dependency name (`loop_interpolation_core`, `surfepy`, `loopstructuralvisualisation`, `matplotlib`, `pyevtk`), and say which optional cells passed and which were excluded. Matrix: expected here: rows 6 and 7 "excluded" with the named dependency; the 2 display cells excluded with a named dependency; 0 mandatory cells skipped; the 2 stub cells printed apart and never counted as a native or Surfe pass. If `loopstructuralvisualisation` is found on the path, the display cell must pass (report it). Exclusion ledger: every other skipped test of the three junit files, grouped by dependency or reason: native extension (the native files, `test_native_errors.py:250`, `test_builder_solver.py:28`); Surfe (`test_surfe_rbf_interpolator.py`, `test_p0_surfe_nans.py`); `geoh5py`, `omf`, `pyevtk` (`tests/unit/io/test_geoh5.py`, `test_omf.py`, `test_exporters.py`; search for `trimesh` again); viewer, `matplotlib`, `pyevtk` (two cells of `test_world_location_consumers.py`, `test_api.py:121`); `pyvista` (`test__structured_grid.py`, `test__surface.py`, `test_2d_discrete_support.py`); `dill` (`test_p2fold_model.py`, `test_geological_model.py`); skips with no named dependency (the four `TetMesh ... evaluate_d2` skips in `test_p2_interpolator.py`, `golden/test_golden.py`) with the skip reason. The 47 skipped tests of `packages/loop_interpolation/tests` must all appear.
    Rollback check: throw-away worktree `git worktree add --detach ../LoopStructural-005-s6-rollback <S5 merge commit>` (write the hash in the record), own `PYTHONPATH` (check `LoopStructural.__file__`); (a) `from LoopStructural.utils import EuclideanTransformation` and `import LoopStructural.utils._transformation` work on both trees, the `__init__` signatures are equal, and a pickle made on the S5 tree loads on the S6 tree; (b) the public names of `loop_common.math` on S5 are a subset of S6, which adds `EuclideanTransformation` only; (c) for a float `(N, 3)` array `transform` is `np.array_equal` on both trees for `dimensions=2` with set angle and translation (the z column differs by design for `dimensions=3`); (d) `from loop_interpolation import *` gives the same names; (e) the signatures of the five `GeologicalFeature` diagnostic methods are equal and the residual numbers of a feature with no faults and no regions are equal; (f) the baseline files are byte-equal (`git diff --stat <S5 merge commit> HEAD -- <files>` prints nothing). Then `git worktree remove --force` and `git worktree prune`. Record: OP21 (extension cells excluded; SC-002 limits not measured), OP22 (P2 state cells only), OP23 (other copies not audited), OP24 (xy axis kept). Write in ASD-STE100 Simplified Technical English.

- [ ] T108 [Finish] [M] Record the S6 decisions and "Spec 005 complete" in `ROADMAP.md`.
  - Files: `ROADMAP.md`
  - Needs: T107
  - Wave: 52
  - Done when: `grep -c "Residuals of faulted features (spec 005, S6)" ROADMAP.md` prints 1, `grep -c "EuclideanTransformation (spec 005, S6)" ROADMAP.md` prints 1, and `grep -c "Spec 005 complete" ROADMAP.md` prints 1.
  - Context: Read `ROADMAP.md` first and copy the style of the neighbouring subsections. Date the entries with the day of the S6 end gate. Add two `###` subsections at the end of `## Decisions`, before `## Stage sequence`. (1) `### Residuals of faulted features (spec 005, S6)`: the residual diagnostics of a faulted feature evaluate the interpolator at the restored positions; S5 decision 8 changes for the two private methods `_value_residuals_fresh` and `_gradient_cosines_fresh`; the notice is in `COMPAT.md`. (2) `### EuclideanTransformation (spec 005, S6)`: a behaviour change; name the four defects (the 3D `fit`, the sign of the fitted angle, the reversed z axis, the silent need for scikit-learn) and their physical effect (a mirrored z axis and a wrong rotation of the map view for some lines); one class now lives in `loop_common` and `LoopStructural.utils` re-exports it; include the sentence "This reverses item 2c-5 of the stage sequence, which was closed as deferred."; the notice is in `COMPAT.md`. Also add a sentence on item 2c-5 in the stage sequence if that item is listed there.
    Add to the status log: "Spec 005 complete: stages S1 to S6 implemented; the auxiliary transformation, the R2 displacement validation, the backend matrix and the golden review are done. Gate result: root pytest <passed> passed, <skipped> skipped (every skip names its dependency in `specs/005-coordinate-frame-correctness/verification.md`)." Replace the two placeholders with the numbers of the root run in the `## S6 final gates` section written by T107. Do not change other text. Write in ASD-STE100 Simplified Technical English.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001, T002, T003 | 1 | 2 |
| 2 | T004, T007, T011, T014 | 1 | 3 |
| 3 | T005, T006, T008, T012 | 1 | 3 |
| 4 | T009, T010 | 1 | 1 |
| 5 | T013 | 1 | 0 |
| 6 | T015, T016 | 1 | 1 |
| 7 | T017 | 0 | 1 |
| 8 | T018, T019, T021 | 2 | 1 |
| 9 | T020, T022, T029 | 1 | 2 |
| 10 | T032, T023, T024, T027, T028 | 4 | 1 |
| 11 | T025, T026 | 0 | 2 |
| 12 | T030 | 1 | 0 |
| 13 | T031 | 0 | 1 |
| 14 | T033, T034 | 1 | 1 |
| 15 | T035, T037 | 0 | 2 |
| 16 | T036, T038 | 0 | 2 |
| 17 | T039 | 0 | 1 |
| 18 | T040 | 0 | 1 |
| 19 | T041, T042 | 0 | 2 |
| 20 | T043, T044, T045, T046 | 3 | 1 |
| 21 | T047, T048, T049 | 2 | 1 |
| 22 | T050 | 0 | 1 |
| 23 | T051 | 0 | 1 |
| 24 | T052 | 0 | 1 |
| 25 | T053 | 0 | 1 |
| 26 | T054 | 0 | 1 |
| 27 | T055, T056 | 0 | 2 |
| 28 | T057, T058, T059, T060, T061, T062 | 5 | 1 |
| 29 | T063, T064, T065, T066 | 2 | 2 |
| 30 | T067 | 0 | 1 |
| 31 | T068 | 0 | 1 |
| 32 | T069, T070, T071 | 0 | 3 |
| 33 | T072, T073 | 0 | 2 |
| 34 | T074, T075, T076, T077 | 2 | 2 |
| 35 | T078 | 0 | 1 |
| 36 | T079, T080 | 0 | 2 |
| 37 | T081 | 0 | 1 |
| 38 | T082 | 0 | 1 |
| 39 | T083 | 0 | 1 |
| 40 | T084 | 0 | 1 |
| 41 | T085, T086, T087 | 1 | 2 |
| 42 | T088 | 0 | 1 |
| 43 | T089, T090, T091, T092, T093, T094 | 1 | 5 |
| 44 | T095, T096, T097, T098 | 0 | 4 |
| 45 | T099 | 0 | 1 |
| 46 | T100, T101 | 0 | 2 |
| 47 | T102, T103 | 1 | 1 |
| 48 | T104 | 0 | 1 |
| 49 | T105 | 0 | 1 |
| 50 | T106 | 0 | 1 |
| 51 | T107 | 0 | 1 |
| 52 | T108 | 0 | 1 |

## Critical path

T003 -> T004 -> T005 -> T010 -> T013 -> T015 -> T017 -> T018 -> T020 -> T023 -> T025 -> T030 -> T031 -> T033 -> T035 -> T036 -> T039 -> T040 -> T041 -> T043 -> T047 -> T050 -> T051 -> T052 -> T053 -> T054 -> T055 -> T057 -> T063 -> T067 -> T068 -> T069 -> T072 -> T075 -> T078 -> T080 -> T081 -> T082 -> T083 -> T084 -> T085 -> T088 -> T089 -> T095 -> T099 -> T100 -> T102 -> T104 -> T105 -> T106 -> T107 -> T108

## Dependencies

- Setup (T001, T002) and the foundation test T003 start at once.
- Foundation: T003 -> T004 -> T005, T006.
- US2: T007 -> T008 -> T009, T010 (T010 also needs T005).
- US1: T011, T014 -> T012 -> T013 (T013 also needs T006 and T010, because
  the fault frame uses the repaired reference sampling and resize).
- Finish: T015 and T016 after T013; T017 last.
- S2 starts after T017. Foundation: T018 (fixtures), T019 (fold helpers).
- S2 tests: T020 (fold frames) needs T018 and T019; T021 (P2 reference
  cells) and T022 (constant norm) are coverage tests.
- S2 code: T023 (directional contract) -> T025, T026 (FDI and PLI fold use
  it); T024 (native), T027 (P2 fold), T028 (P2 weights) after T020; T029
  (docstring) after T019. T032 (scalar `element_size` fix, found by T022)
  after T022.
- Finish: T030 after the code tasks; T031 last.
- S3 starts after T031 (S2 is its base). Foundation: T033 (fixtures) and T034
  (module stubs and exports) in wave 14.
- S3 tests: T035 (payload and round trips) then T036 (loader failures), both in
  `test_frame_persistence.py`; T037 (copies, stale frame) then T038 (fault
  restoration), both in `test_frame_copy_and_lifecycle.py`. T039 records the
  pre-fix evidence after all four.
- S3 code: T040 -> T041 (`_frame_persistence.py`); T042 (base class) needs T040;
  T043 (factory loader) needs T041 and T042; T044, T045, T046 (copies) need T042.
- Finish: T047 (`COMPAT.md`), T048 (`ROADMAP.md`) and T049 (rollback check) after
  the code; T050 (gates, `verification.md`, ROADMAP status) last.
- S4 starts after T050. Foundation: T051 (fixture helper). Tests, one file in a chain:
  T052 (display and export), T053 (intrusion), T054 (2D viewer); T055 (pre-fix
  evidence) and T056 (deprecation tests) after T054.
- S4 code, wave 28, separate files: T057 (`get_data`), T058 (intrusion helper), T059
  (`rescale` and `scale` warn), T060 and T061 (exporters), T062 (example).
- Finish: T063 (golden review) after the code; T064 to T066 (`COMPAT.md`,
  `ROADMAP.md`, `API.md`); T067 (rollback check); T068 (gates) last.
- S5 starts after T068 (S4 is its base). Tests, wave 32: T069 (resize), T070 (cell centres), T071 (residual and exclusion functions). T072 (pre-fix evidence) after T069 and T070.
- S5 code: T073 (`_diagnostics.py`) -> T074 (export); T075 (3D setters), T076 (2D counts) and T077 (rectilinear guard) after T072. T078 (numeric survey, golden review, regenerate golden data) after the three; T079 (post-fix evidence) after T078.
- S5 feature tests, one file in a chain: T080 -> T081 (after T073 and T078); T082 (pre-fix evidence). Misfit code: T083 -> T084 (one file; after T073 and T082).
- Finish: T085 (`COMPAT.md`), T086 (`ROADMAP.md`) and T087 (rollback check) after T084; T088 (gates) last.
- S6 starts after T088 (S5 is its base); each wave-43 task Needs T088. Tests and guards, wave 43: T089 (fixture helpers), T090 (transformation tests, fail first), T091 (`test_diagnostics_frames.py` findings), T092 (Surfe module name), T093 (BOM removal), T094 (matrix helper and tests). Wave 44: T095 (R2 validation) , T096 (faulted-feature residuals), T097 (boundary stubs) need T089; T098 (optional cells) needs T089 and T094. T099 (pre-fix evidence) after all of them except T089 and T094 (they are reached through T095 to T098).
- S6 code: T100 (`loop_common` transformation) and T101 (residual fix) after T099; T102 (re-export) after T100; T103 (`MATRIX` list and matrix tests) after the tests and the two fixes. Finish: T104 (post-fix evidence), T105 (golden review and ledgers), T106 (`COMPAT.md`), T107 (final gates, matrix, rollback, exclusion ledger), T108 (`ROADMAP.md`, "Spec 005 complete") in a chain, one file in each wave.
- No open point blocks a task. OP12, OP16 and OP17 have defaults (OP17 can block the S5 merge if a parity file changes: T078 reports it). OP18 to OP20 are closed; OP21 to OP24 have defaults and block nothing. The only change to the plan graph: the six first S6 tasks Need T088, so that they start after S5.
