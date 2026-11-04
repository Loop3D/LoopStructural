# Tasks: Coordinate and fault-frame correctness (stages S1 and S2)

**Feature**: `005-coordinate-frame-correctness` | **Plan**: [plan.md](plan.md)

<!--
Rules:
- Each task is small (less than one day) and names the files that it changes.
- Test tasks come BEFORE the implementation tasks that they test.
  Run the test and see it fail first.
- [P] = can be done at the same time as the other [P] tasks of the same wave.
- Written in ASD-STE100 Simplified Technical English.
- This file covers stages S1 (T001 to T017, done) and S2 (T018 to T032) of
  plan.md. Stages S3 to S6 wait for OP1 (a later spec-planner run for each
  stage, then a new spec-tasks run).
-->

## Orchestrator notes

- Run this file with the spec-orchestrator skill. It has the rules for the
  worktree, the sub agents, the checks and the commits.
- Waves: 13 (S1: 1 to 7, done; S2: 8 to 13). Tiers: [S] = simple model, [M] = standard model.
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

- [ ] T030 [P] [Finish] [S] Add the S2 migration notice to `COMPAT.md`.
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

- [ ] T031 [Finish] [M] Record the S2 evidence, backend cells and gates in `verification.md`.
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
| 10 | T023, T024, T027, T028, T032 | 4 | 1 |
| 11 | T025, T026 | 0 | 2 |
| 12 | T030 | 1 | 0 |
| 13 | T031 | 0 | 1 |

## Critical path

T003 -> T004 -> T005 -> T010 -> T013 -> T015 -> T017 -> T018 -> T020 -> T023 -> T025 -> T030 -> T031

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
- Stages S3 to S6 have no tasks yet (OP1 in plan.md).
