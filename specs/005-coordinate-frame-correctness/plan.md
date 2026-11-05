# Implementation plan: Coordinate and fault-frame correctness

**Feature**: `005-coordinate-frame-correctness` | **Date**: 2026-11-04 | **Spec**: [spec.md](spec.md)

<!-- Written in ASD-STE100 Simplified Technical English.
The constitution (.specify/memory/constitution.md) has priority. -->

## Open points

- **OP1**: This plan gives the exact design only for stages S1, S2 and S3 (S1:
  shared transform, reference-constraint boundary and fault frame; S2: fold, P2,
  directional-regularisation and constant-norm cells; S3: saved frames, copies
  and transform changes). The spec-planner limit is 10 code files, and the
  maintainer selected a staged plan on 2026-11-04. S2 has its exact design
  since 2026-11-04 ("Stage S2 design"). S3 has its exact design since
  2026-11-05 ("Stage S3 design"). Stages S4 to S6 have a stage map, a file
  list from a grep survey and a requirement allocation, but no exact
  signatures. Owner: the spec-planner skill, one run for each stage. Blocks:
  the tasks of S4 to S6. The spec-tasks skill makes the S3 tasks next (S1 is
  merged; S2 is in pull request #396).
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
- **OP7**: The P2 value-only fit of R1 (6 x 6 x 6 cells, ridge and gradient
  regularisation, `cg` with the SC-003 settings) misses the SC-002 limits in every
  frame, also the unrotated control at `ANCHORS[0]`: value error `1.8e-8` (limit
  `1e-8*F`, F = 7.5) and gradient error `3.9e-7` (limit `4e-9`). The cause is solver
  and regularisation accuracy, not a frame defect. FDI and PLI meet SC-002. Owner:
  the maintainer or the solver owner. Decision needed: improve the P2 fit set-up
  (production solver code, a separate task), or record P2 fitting as outside
  SC-002. Default until an answer: the S3 P2 cells use the exact nodal control.
  Blocks: nothing in S3.
- **OP6**: The feature `copy()` methods under `LoopStructural/modelling/features`
  were not read. If one of them calls `interpolator.copy()`, the frame now
  stays on the result. Owner: the S3 copy test task (it runs `git grep -n
  'interpolator.*copy('` first). Default until an answer: no change in S3.
  Blocks: nothing.

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
  be a major version, or the maintainer records another decision.

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
changes keep the frame ("Stage S3 design"). S4 to S6 repair display
consumers, grid and diagnostics, and the auxiliary utilities and validation
matrix (see "Stages and rollback").

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

S1, S2 and S3 have the exact design below. S4 to S6 have scope only (OP1).

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
   wrote one). The caller must pass `local_origin` and `local_rotation`, or
   `declare_identity=True`. The loader never guesses.
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
  the major-version note below.

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
| `COMPAT.md` | S3 migration notice with the major-version note |
| `ROADMAP.md` | FR-019 decision under "Decisions"; status-log line when S3 completes |
| `specs/005-coordinate-frame-correctness/verification.md` | S3 pre-fix, post-fix, golden review, backend matrix, rollback and run-time sections |

Not changed (the survey in "Later stages" listed them): `_support_factory.py`,
`_fault_segment.py`, `_discrete_interpolator.py` and `_model_serializer.py`.

The `COMPAT.md` row says: loading a legacy interpolator payload (no
`format_version`) assumed an identity frame without a message. It now needs
`local_origin` and `local_rotation`, or `declare_identity=True`, and raises
`InterpolatorFormatError` otherwise. New payloads carry `format_version: 2`.
`copy()` keeps the frame. A transform change after constraints exist raises
`StaleFrameError` at the next use. Major-version note: the release that
contains S3 must be a major version, or the maintainer records another
decision in `ROADMAP.md`.

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

### Stages and rollback

| Stage | Files and behaviour | Acceptance check | Rollback target and compatibility check |
|---|---|---|---|
| S1 | Shared transform boundary, reference relationships for FDI and PLI, fault frame, buffer and forced geometry, rotation check, resize stability (this design) | New S1 tests pass for FDI and PLI; `pytest`; `ruff check LoopStructural packages`; `pytest packages/loop_common/tests packages/loop_interpolation/tests` | Revert the S1 commits to master `ebc7c699`. Check: the existing `tests/unit/modelling/test__fault_builder.py` and `packages/loop_common/tests/test_bounding_box.py` pass on the reverted tree. No saved format changed. |
| S2 | Fold geometry, directional regularisation, quadrature, masks, sign alignment for FDI-fold, PLI discrete fold, P2 and P2-fold; constant-norm variants (FR-008; R3 for P2; R4) ("Stage S2 design") | New S2 tests (SC-002 for constant norm; SC-007, SC-008 and SC-013 for every mandatory fold and P2 cell); the existing fold package tests pass with no change; `pytest`; `ruff check LoopStructural packages`; package tests | Revert to the S1 merge. Check: the S1 tests still pass. No saved format changed. |
| S3 | Persistence format version, transform and coordinate-space metadata ("Stage S3 design"); no second projection at load; copies keep the frame; legacy load requires a transform or identity declaration; transform change policy and restored-point cache invalidation (FR-016 to FR-020); ROADMAP.md decision for FR-019; `COMPAT.md` migration note | SC-006, SC-012, SC-016, SC-017 | Revert to the S2 merge. Check: old-format R1 files load through their documented path; a new-format file fails with 1 actionable error and 0 state changes on the old reader (SC-017). |
| S4 | Display, export of observation geometry, fault-centre geometry and intrusion classification use stored world locations (FR-022) | World-location comparison within `Tp` (U3.4, U3.5) | Revert to the S3 merge. Check: viewer import paths unchanged. |
| S5 | Cell-centre sampling, remaining resize paths, non-mutating diagnostics, residual and exclusion report (FR-014 remainder, FR-015, FR-023, FR-024) | SC-009, SC-011 | Revert to the S4 merge. Check: diagnostics report fields unchanged for callers. |
| S6 | Auxiliary point transformations and duplicate entry points (FR-025); displacement validation R2 (FR-026, FR-027); backend matrix incl. native and Surfe (FR-030); golden review (FR-029); 005 completion line in ROADMAP.md | SC-004, SC-010, SC-014, SC-015, SC-018 | Revert to the S5 merge. Check: `LoopStructural.utils` and `loop_common.math` import paths unchanged. |

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
  constant-norm variants ("Stage S2 design"). S6 covers native and Surfe
  (OP1).
- **Cross-language fixtures**: Not applicable in S1.
- **Grid and display checks**: S1 checks support enclosure and resize. Cell
  centres and display are S4 and S5.
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

### Later stages (file survey only, OP1)

The next planner runs read and confirm these files. The lists come from a
grep survey on 2026-11-04 and are not final.

- **S2**: planned. See "Stage S2 design" and "S2 files".
- **S3**: planned. See "Stage S3 design" and "S3 files". The survey listed
  `_support_factory.py` and `_fault_segment.py`; the code shows that S3 does
  not change them (no restoration cache exists).
- **S4**: `packages/loopstructural_visualisation/src/loopstructuralvisualisation/_3d_viewer.py`,
  `LoopStructural/modelling/intrusions/intrusion_builder.py`,
  `LoopStructural/modelling/intrusions/intrusion_frame_builder.py`,
  `LoopStructural/modelling/features/fault/_fault_segment.py`,
  `LoopStructural/export/exporters.py`.
- **S5**: `packages/loop_common/src/loop_common/geometry/_bounding_box.py`
  (`cell_centres`), `_structured_grid.py`, `_structured_grid_2d.py`,
  `packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`,
  `_3d_rectilinear_grid.py`,
  `packages/loop_interpolation/src/loop_interpolation/_diagnostics.py`.
- **S6**: `packages/loop_common/src/loop_common/math/_transformation.py`,
  `LoopStructural/utils/_transformation.py`,
  `packages/loop_interpolation/src/loop_interpolation/_native.py`,
  `_surfe_wrapper.py`, `LoopStructural/modelling/features/fault/_fault_function.py`,
  `ROADMAP.md`.

## Requirement coverage

| Requirement | Design and files | Verification | Gap |
|---|---|---|---|
| FR-001 | Existing public evaluation; S1 keeps it; fault path in `_model_feature_factory.py` | `test_fault_frame_analytic.py` (gradients in world after rotation) | None for S1 paths; other paths OP1 |
| FR-002 | Factory keeps explicit centre in world; builder rows in world | `test_fault_frame_analytic.py` U1.1 | None |
| FR-003 | Public conversion methods; builder uses each once | `test_rigid_transform_contract.py`; `test_reference_constraint_frames.py` | Load path (U3.2) is S3, OP1 |
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
| FR-014 | `_cells_for_length`; setters | Resize-10-times test | Other supports S5, OP1 |
| FR-015 | S5 | S5 | OP1 |
| FR-016 | `frame_payload`; `to_dict` frame keys; format version 2 ("S3 saved format") | `test_frame_persistence.py::test_payload_declares_frame` | None |
| FR-017 | `_restore_local_constraints` stores local arrays with no conversion | `test_local_arrays_unchanged_after_two_cycles` | None |
| FR-018 | `_inherit_frame` in the three `copy()` methods | `test_frame_copy_and_lifecycle.py::test_copy_keeps_frame` | PLI has no `copy` (documented, not added); OP6 |
| FR-019 | `payload_kind`, `resolve_frame`, `InterpolatorFormatError`; ROADMAP.md decision; `COMPAT.md` row | `test_legacy_*`, `test_conflicting_inputs`, `test_malformed_metadata`, `test_unsupported_future_version` | None |
| FR-020 | `frame_signature`, `_check_frame`, `StaleFrameError` (reject at use) | `test_stale_frame_rejected_without_state_change`, `test_fault_restoration_uses_current_frame` | Policy choice OP4; no restoration cache exists |
| FR-021 | `validate_rigid_rotation` in `set_local_transform` and `rotation_xy`; saved metadata checked in `resolve_frame` | `test_rigid_transform_contract.py`; `test_malformed_metadata` | None |
| FR-022 | S4 | S4 | OP1 |
| FR-023 | S5 | S5 | OP1 |
| FR-024 | `excluded_rows` (S1 part) | Restricted-support exclusion count | Report integration S5, OP1 |
| FR-025 | S6 | S6 | OP1 |
| FR-026 | S1: values, gradients, residuals (SC-003, SC-008); displacement S6 | `test_fault_frame_analytic.py` | Displacement OP1 |
| FR-027 | S6 | S6 | OP1 |
| FR-028 | No signature change except new keyword arguments; `COMPAT.md` note; FR-019 exception recorded in ROADMAP.md | API contract tests; existing round-trip tests unchanged | None |
| FR-029 | `verification.md` golden review section | Review record per changed result | None |
| FR-030 | S1 FDI and PLI cells; S2 P2, fold and constant-norm cells | S1 and S2 tests | Native and Surfe cells S6, OP1 |
| FR-031 | Pre-fix rule in "Test strategy" | `verification.md` | None |
| FR-032 | "Stages and rollback" | Rollback checks in `verification.md` | None |
| SC-001 | `_rigid.py`, `BoundingBox` | `test_rigid_transform_contract.py` | None |
| SC-002 | R1 fit for the constant-norm variants (S2) | `test_constant_norm_frames.py` | Other R1 cells and lifecycle S3, OP1; OP3 |
| SC-003 | Fault design | `test_fault_frame_analytic.py` | Run time OP2 |
| SC-004 | S6 | S6 | OP1 |
| SC-005 | Fault and enclosure design | `test_fault_buffer_bounds.py` | Run time OP2 |
| SC-006 | Format version 2; local arrays restored with no conversion | `test_round_trip_world_values`, `test_local_arrays_unchanged_after_two_cycles` | None |
| SC-007 | Reference sampling; S2 fold sampling and P2 | `test_reference_constraint_frames.py`; `test_fold_constraint_frames.py` | None |
| SC-008 | Reference sampling; S2 fold rows and P2 | `test_reference_constraint_frames.py` (FDI, PLI, P2); `test_fold_constraint_frames.py` | None |
| SC-009 | S5 | S5 | OP1 |
| SC-010 | S6 | S6 | OP1 |
| SC-011 | `_cells_for_length` (resize part) | Resize test | Cell centres S5, OP1 |
| SC-012 | Rotation part in S1; metadata and legacy part in S3 | `test_rigid_transform_contract.py`; `test_legacy_*`, `test_malformed_metadata` | None |
| SC-013 | Extent-direction row; S2 fold R4 directions | `test_fault_frame_analytic.py`; `test_fold_constraint_frames.py` | None |
| SC-014 | Merge gates in every stage | `pytest`, ruff, package tests | None |
| SC-015 | `verification.md` golden review | Review record; 0 tolerance changes | None |
| SC-016 | Reject at use; rebuild recipe in a new frame | `test_stale_frame_rejected_without_state_change`, `test_rebuild_in_new_frame_preserves_world_evaluation`, `test_fault_restoration_uses_current_frame` | Acceptance by automatic rebuild is not built (OP4) |
| SC-017 | Rollback column; format part S3 ("S3 rollback check") | `test_unsupported_future_version`, `test_old_reader_would_name_format_version`; `verification.md` | Old-reader error wording OP5 |
| SC-018 | S1 and S2 cells; matrix S6 | `verification.md` backend matrix | Native and Surfe S6, OP1 |
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
- **Exclusion record**: `GeologicalFeatureBuilder.excluded_rows:
  dict[str, dict[str, int]]`. S5 connects it to the diagnostics report.
- **Fault geometry**: an explicit `fault_center` is a world point; the
  generated centre extent row is `coord=2` with `strike_unit / major_axis`.
  Spec 003 golden fault models can change only under C-006/C-004 with a
  record in `verification.md`. Spec 003 cache invalidation for transform
  changes uses `frame_signature` (S3).
- **Shared files**: `_geological_interpolator.py`, `_discrete_interpolator.py`
  and `_fault_builder.py` are also in the scope of specs 003 and 004.
  Coordinate the merge order with the active 003 worktree
  (`LoopStructural-native`).

## Complexity tracking

None.
