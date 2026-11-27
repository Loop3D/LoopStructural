# Implementation plan: native fault application

**Feature**: `007-native-fault-kernel` | **Date**: 2026-11-27 | **Spec**: [spec.md](spec.md)

<!-- A plan says HOW. Written in ASD-STE100 Simplified Technical English. -->

**Code base of this plan.** The plan was written on the spec 003 branch at
`a558a1bd`. Spec 003 merged to master as #391 (`6874364c`) on 2026-11-27. The
OP1 check found that each planned file is the same on master, except
`benchmarks/fault_evaluation/run.py` (a new medium-case repeat tier; `main()`
and its options do not change). So the plan applies to master `6874364c`. The build branch is
`feat/007-native-fault-kernel` (worktree
`/home/lgrose/dev/LoopStructural-wt/feat-007-native-fault-kernel`), which adds
the 94-fault benchmark in `907acf42`.

## Open points

All open points are answered (see "Answers").

- **OP1 Base of the plan.** Spec 003 is not merged. When it merges, compare
  the private names of "Names that this plan uses" with master. If a name
  changed, update this plan before a task of stage S2 starts. Default: stage
  S1 can start on the spec 003 branch, because it changes only
  `_fault_block.py` and two callers.
- **OP2 The 94-fault model.** The model is not in the repository. SC-001 and
  SC-004 need it. Default: the maintainer gives a saved model file. The
  benchmark case `network94` loads it from the path in the environment
  variable `LOOPSTRUCTURAL_NETWORK94`. If the variable is not set, the case is
  skipped and the stage report says that SC-001 or SC-004 was not measured.
- **OP3 Cause of the 40 ms.** The spec says that most of the 40 ms fixed cost
  of a cut fault frame is the restore order. The code shows a second cost: a
  cut frame uses `_evaluate_frame_legacy`, which restores the same points 4
  times for each step (gx, gy, gz and the gradient of gy). Default: task S1-a
  measures first. D1a (order cache) is always done (FR-001). D1b (one shared
  restoration for the members of a cut frame) is done only if SC-001 fails
  after D1a.
- **OP4 Regions on a fault frame.** An abutting fault adds a
  `PositiveRegion` or `NegativeRegion` to member 0 of the frame. Default: the
  native path supports these two region types when their parent is a fault
  of the same chain. A fault frame with any other region (for example an
  unconformity) is not supported, so the feature uses the Python path
  (C-005).
- **OP5 Gradient of a faulted feature (FR-004).** The tetrahedron method
  checks the fault side of each corner with `fillnan='nearest'` (a k-d tree
  of SciPy). Default: the loop of the method (the straddle check, at most 10
  refinements, the NaN rule) stays in Python. The native path does each
  numerical part of the loop: the frame values of the side check, the
  restoration of the corners and the field values. A full native loop is not
  part of this spec.
- **OP6 Data kept between calls (FR-009, FR-010).** The native data of a
  feature is kept on the feature until the next call that finds a changed
  state. The next call never uses old data, but the memory of the old state
  is released only at that call. Default: accept this. The alternative is a
  weak reference to each array, which NumPy arrays do not allow on all types.
- **OP7 Copy of the points.** The native call copies the input points to Rust
  memory so that it can release the GIL (24 bytes for each point). Default:
  accept the copy. Task S2-h measures SC-007 at 10 million points. If SC-007
  fails, read the points in place under the GIL-free call (unsafe code), in a
  new task.
- **OP8 In-place change of node values.** The state token identifies the node
  values of a grid by the array object, not by its content (see D2c). A
  change in place (`interpolator.c[:] = ...`) is not found. Default: accept,
  because the solvers make a new array. If the maintainer knows code that
  changes `c` in place, the token must add a digest of `c`.
- **OP9 The proof-of-concept crate.** The crate is outside the repository
  (`/home/lgrose/dev/LoopStructural-wt/_runs003/rust_poc/loop-snapshot`, see
  "Answers"). Default: task S2-a reads it (read
  only) to name the cause of its difference (C-007). The new kernel is
  written from the Python code, not copied from the crate.

## Answers

Answered by the maintainer on 2026-11-27.

- **OP1**: done 2026-11-27 after spec 003 merged (#391, `6874364c`). No
  name of "Names that this plan uses" changed. Base of the plan: master.
- **OP2**: answered 2026-11-27 by `specs/007-native-fault-kernel/HANDOFF.md`
  (commit `907acf42` on `feat/007-native-fault-kernel`). The model is
  `benchmarks/fault_evaluation/many_faults.py` (deterministic, seed 100; the
  build takes about 7 minutes; not in the standard suite):
  `python -m benchmarks.fault_evaluation.many_faults --output RESULT.json [--points 500 2000 5000] [--save-reference REF.npz | --compare REF.npz]`.
  `--compare` exits with 1 if a model id or group value differs bit for bit
  (SC-001). The Python-path reference from master is in
  `/home/lgrose/dev/LoopStructural-wt/_runs007/` (`many_faults_python.json`,
  `many_faults_python_ref.npz`): cutting fault mean 40.0 ms at 500 points,
  `evaluate_model` 2.09 s at 5000 points. This plan adds no benchmark case.
- **OP3**: measure, then fix (the default). D1a always; D1b only if SC-001
  fails after D1a.
- **OP4**: support abutting regions **and** unconformity regions on a fault
  frame (not only the default). D2a `region.rs`, D2b and D2c are changed.
- **OP5**: the whole tetrahedron gradient is native (not the default),
  including the side check with the nearest-point fill and the 10
  refinements. D3 is changed; new file `gradient.rs` and `nearest.rs`.
- **OP6**: accept (the default).
- **OP7**: copy the points and measure SC-007 (the default).
- **OP8**: identity of the array (the default).
- **OP9**: read the proof-of-concept crate only to find the cause (the
  default). Its path (from `HANDOFF.md`) is
  `/home/lgrose/dev/LoopStructural-wt/_runs003/rust_poc/` (`README.md`,
  `loop-snapshot/`, `pyext/`, `export_snapshot.py`), not the kiridashi
  repository.

## Summary

The work has four stages. Each stage can merge alone.

- **S1 (story 1, Python only).** Compute the restore order of a feature once
  for each state of its faults, and keep it on the feature.
- **S2 (stories 2, 4 and 5).** A Rust fault kernel in the existing package
  `loop_interpolation_core` restores points through the faults of a supported
  feature. LoopStructural uses it automatically. A setting, an environment
  variable, a kernel version check and a thread count control it.
- **S3 (story 3).** The same kernel evaluates the field values and grid
  gradients of a supported feature, and the whole tetrahedron gradient of a
  faulted feature (OP5).
- **S4 (story 4, packages).** The existing wheel workflow builds and smoke
  tests the fault kernel on all platforms of FR-014.

The native path lives in `loop_interpolation_core` because that package
already has PyO3, rayon, an abi3-py39 build, a wheel workflow for the five
platforms of FR-014 (`.github/workflows/core-wheels.yml`) and the extra
`loopstructural[native]`. A new package would repeat all of this.

## Technical context

- Python 3.9 or later. Rust edition 2024, PyO3 0.29, numpy 0.29, rayon 1.10
  (all already in `packages/loop_interpolation_core/Cargo.toml`). No new
  dependency.
- The kernel is an optional import. Without it, LoopStructural behaves as
  before (FR-006c).
- The Python path stays the reference. The native path agrees within
  1e-8 × the size of the model box (spec 003 FR-004, C-007).
- Present restore path (spec 003 branch):
  `EvaluationContext.restore` and `BaseFeature._apply_faults` call
  `restore_through` (`LoopStructural/utils/_fault_block.py`). It calls
  `with_cutters` and `cutters_first` at each call, and then
  `FaultSegment.apply_to_points` for each fault. A cut fault frame uses
  `FaultSegment._evaluate_frame_legacy`, which calls `evaluate_value` of each
  member, which calls `restore_through` again.
- Inside `restore_through`, each fault that cuts fault k is in the list
  `done` before fault k. So the members of a cut frame are evaluated at
  points that are already restored, and they use the plain grid gradient, not
  the tetrahedron method. The native restoration needs no tetrahedron method.

## Constitution check

| Rule | Result |
|---|---|
| Purpose and scope: a Python library, staged rebuild | Pass. The kernel is an optional part of the existing package `loop_interpolation_core`. |
| Development commands: `pytest`, `black .`, `ruff check LoopStructural packages` | Pass. Each task runs them. Rust tasks also run `cargo fmt --check`, `cargo clippy --all-targets -- -D warnings` and `cargo test --release` in `packages/loop_interpolation_core` (`.dev-tools/project.json`). |
| Merge gate: `pytest` and Ruff pass before merge | Pass. FR-013 adds a test run with the native path on and one with it off. |
| Code conventions: numpy docstrings, type hints, py39 | Pass. New Python modules use `from __future__ import annotations` for `X | None`. |
| `ROADMAP.md` is the source of truth; update it at each stage | Pass. FR-015: one task in each stage updates `ROADMAP.md`. |
| Work ships in stages that can be tested and reverted alone | Pass. S1 to S4 each merge alone. Setting `LoopStructuralConfig.native_faults = False` turns S2 and S3 off at run time. |
| The public `GeologicalModel` API is a hard boundary | Pass. No signature of `GeologicalModel` or of a feature changes (FR-007). `LoopStructuralConfig` gets two new fields with defaults; this adds to the API and breaks nothing. `tests/fixtures/api_surface_snapshot.json` does not list `LoopStructuralConfig` (checked 2026-11-27), so it does not change. |
| The QGIS plugin stays in a separate repository | Pass. No plugin change. |
| Specs live under `specs/NNN-name/` | Pass. |
| Releases use release-please for each package | Pass. The kernel ships in the next release of `loop-interpolation-core`. No new package. |

## Design

### D1. Restore order once for each state (stage S1, FR-001, FR-002)

**D1a. Order cache.** New function in `LoopStructural/utils/_fault_block.py`:

```python
def restore_order(faults: Iterable, *, reverse: bool = False, owner: object | None = None) -> list:
    """Return the faults and their chain of cutting faults in restore order."""
```

- Without `owner`, it returns `cutters_first(with_cutters(faults))`, or for
  `reverse=True` `cutters_first(with_cutters(faults)[::-1])[::-1]`, as
  `restore_through` does today.
- With `owner`, it keeps the result in `owner._restore_order`, a tuple
  `(key, chain, ordered)`:
  - `chain = with_cutters(faults)` (a list; it keeps a reference to each
    fault, so an `id` in `key` cannot be used again by a new object while
    the entry exists);
  - `key = (bool(reverse), tuple((id(f), tuple(id(c) for c in _cutters(f))) for f in chain))`.
- At each call it computes `chain` and `key` (linear in the chain) and
  returns the kept `ordered` list if `key` is equal. Otherwise it computes the
  order (quadratic) and replaces the entry. So a change of the cutting faults
  of any fault of the chain (FR-002) changes the key. A change of the build
  state that changes no cutter does not change the order, so the order stays
  correct.
- `owner._restore_order` is set with one assignment, so two threads that set
  it at the same time keep a correct entry (FR-011a).

`restore_through` gets a new keyword `owner: object | None = None` and calls
`restore_order(faults, reverse=reverse, owner=owner)` in place of its two
lines of order code. No other change.

Callers that pass `owner`:

- `BaseFeature._apply_faults` (`_base_geological_feature.py`): `owner=self`.
- `EvaluationContext.restore` (`_evaluation.py`), both calls of
  `restore_through`: `owner=feature`.

The builder calls (`_geological_feature_builder.py`) do not change. They run
once at build time.

**D1b. One restoration for the members of a cut frame (only if OP3 needs
it).** In `FaultSegment._evaluate_frame_legacy`: if the three members have
the same fault list (compare with `is` for each item) and the list is not
empty, restore the points once with
`restore_through(points, member0.faults, block=block, restored_through=restored_through, owner=member0)`,
and then call each member with the restored points and
`restored_through=tuple(restored_through) + tuple(restore_order(member0.faults, owner=member0))`.
The masks of the members use the original points, as today. The result must
be bit-identical (SC-001). If a test shows a difference, D1b is dropped and
the stage report says so.

### D2. Native fault kernel (stage S2)

#### D2a. Rust module

New folder `packages/loop_interpolation_core/src/fault/`:

| File | Content |
|---|---|
| `mod.rs` | `pub struct FaultChain { faults: Vec<Fault>, block: BlockSettings }`, `pub struct Fault { frame: [GridField; 3], cutters: Vec<usize>, regions: Vec<Region>, profile: Profile, displacement: f64, steps: usize }`, `pub struct BlockSettings { side_epsilon: f64, extension: Extension, exact: f64, limit: f64, step: f64 }`, `pub enum Extension { Smooth, Mirror }`, `pub enum Side { Free, Absent, Hanging, Foot }` (Python 2, 0, +1, -1). |
| `grid_field.rs` | `pub struct GridField { origin: [f64; 3], nsteps: [usize; 3], step_vector: [f64; 3], local_origin: Option<[f64; 3]>, local_rotation: Option<[[f64; 3]; 3]>, values: Vec<f64> }` with `fn value(&self, p: [f64; 3]) -> f64` and `fn gradient(&self, p: [f64; 3]) -> [f64; 3]`. Port of the grid methods `_lookup_points`, `trilinear`, `_values_from_coefficients` and `_evaluate_gradient_from_lookup` (`packages/loop_common/src/loop_common/supports/_3d_structured_grid.py`, `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`), and of the interpolator methods `_project_points` and `_reproject_vectors` (the world-to-local frame of the bounding box). NaN for a NaN coordinate and outside the grid, as the Python lookup. |
| `profile.rs` | `pub enum Profile1 { Ones, Zeros, SmoothPeak, Cubic { w: [f64; 4], lim: [f64; 2], min: f64, max: f64, value_points: Vec<(f64, f64)> }, Composite { positive: Box<Profile1>, negative: Box<Profile1> } }`, `pub struct Profile { scale: f64, gx: Profile1, gy: Profile1, gz: Profile1 }`. Port of `smooth_peak`, `CubicFunction.__call__` and `_limit_value`, `Composite`, `Ones`, `Zeros`, `FaultDisplacement.__call__`, `_saturate`, `_taylor`, `continue_profile` and `block_displacement` (`_fault_function.py`). |
| `region.rs` | `pub enum Region { Sign { positive: bool, fault: usize, vector: Option<[f64; 3]>, point: Option<[f64; 3]> }, Unconformity { field: GridField, faults: Vec<usize>, value: f64, sign: bool } }`. `Sign`: port of `BaseSignRegion._calculate_value_and_distance`, `PositiveRegion.__call__` and `NegativeRegion.__call__` (`LoopStructural/utils/regions.py`), with the value of frame coordinate 0 of the parent fault from this chain (OP4). `Unconformity`: port of `UnconformityFeature.evaluate` (`LoopStructural/modelling/features/_unconformity_feature.py`): restore the point through `faults` (indices in this chain, in restore order; the faults in `restored` are skipped, as `restored_through`), read `field`, then `value_at <= value` if `sign`, else `value_at >= value`; NaN gives false. |
| `nearest.rs` | `pub struct KdTree` over a set of 3D points, `fn nearest(&self, p: [f64; 3]) -> usize`. Exact Euclidean nearest point; for equal distances the lowest index. Port of the `fillnan='nearest'` rule of `GeologicalFeature.evaluate_value` (OP5). |
| `gradient.rs` | `fn feature_gradient(chain, field, points, order, sides, restored, jumping, scale) -> (Vec<[f64; 3]>, usize)`: port of the faulted branch of `GeologicalFeature.evaluate_gradient` (OP5): a regular tetrahedron for each point (`regular_tetraherdron_for_points`, `packages/loop_common/src/loop_common/math/_maths.py`); the side check with frame coordinate 0 of each fault in `jumping`, with NaN corner values filled from the nearest known corner of the same batch (`nearest.rs`); at most 10 refinements, each with half the scale, only for the tetrahedra that still cross; the corners restored with `restore`; the field read at the restored corners; `gradient_from_tetrahedron` (same file). A point that still crosses a fault gets NaN. The second result is the number of points refined at least once. |
| `restore.rs` | `fn apply(&self, fault: usize, p: [f64; 3], reverse: bool, sides: &[Side], restored: &mut Vec<bool>) -> [f64; 3]`: port of `FaultSegment.apply_to_points` for one point, including `force_side`, `_displacement_fraction`, the steps, the NaN rules (#413) and the step direction of #425 (a NaN gradient of a cut frame uses the grid gradient of member 1 at the point). `fn restore(&self, order: &[usize], p: [f64; 3], reverse: bool, sides: &[Side], restored: &[bool]) -> [f64; 3]`: port of the loop of `restore_through`. `fn frame_value(&self, fault: usize, coordinate: usize, p, sides, restored) -> f64`: one frame coordinate with the region mask, as `FaultSegment._evaluate_frame` gives it. |
| `pool.rs` | `pub fn pool(threads: usize) -> Arc<rayon::ThreadPool>`: one pool for each thread count, kept in a `OnceLock<Mutex<HashMap<usize, Arc<ThreadPool>>>>`. `threads == 0` means all cores. |

Each point is independent, so the batch runs as a parallel loop over chunks
of 4096 points on `pool(threads)`. The kernel is immutable after it is made,
so it is `Send + Sync` (FR-011a).

The kernel is written from the Python code, not from the proof-of-concept
crate (OP9).

#### D2b. Python binding (`src/python.rs`, `src/lib.rs`)

`src/lib.rs` adds `pub mod fault;`. `src/python.rs` adds to the module
`loop_interpolation_core`:

```python
FAULT_KERNEL_VERSION: int = 1

class FaultChain:
    def __init__(self, spec: dict) -> None: ...
    def restore(self, points: np.ndarray, order: np.ndarray, reverse: bool,
                sides: np.ndarray, restored: np.ndarray, threads: int) -> np.ndarray: ...
    def frame_value(self, fault: int, coordinate: int, points: np.ndarray,
                    sides: np.ndarray, restored: np.ndarray, threads: int) -> np.ndarray: ...
    def feature_gradient(self, field: GridField, points: np.ndarray, order: np.ndarray,
                         sides: np.ndarray, restored: np.ndarray, jumping: np.ndarray,
                         scale: float, threads: int) -> tuple[np.ndarray, int]: ...

class GridField:
    def __init__(self, spec: dict) -> None: ...
    def values(self, points: np.ndarray, threads: int) -> np.ndarray: ...
    def gradients(self, points: np.ndarray, threads: int) -> np.ndarray: ...
```

- `points`: float64, shape (N, 3). `order`: int64 indices into `spec["faults"]`
  in restore order. `sides`: int8, one for each fault of the spec (2 = the
  block does not name the fault). `restored`: bool, one for each fault of the
  spec. Returns a new float64 array, (N, 3) for `restore`, (N,) for
  `frame_value` and `values`, (N, 3) for `gradients`.
- Each call copies `points` to a Rust vector, releases the GIL
  (`Python::detach`; `allow_threads` before PyO3 0.26) for the work, and
  returns a new array (FR-011). The model data is copied only in `__init__`.
- `feature_gradient`: `jumping` is bool, one for each fault of the spec
  (true: the side check uses this fault). `scale` is the first
  `element_scale_parameter`. Returns the (N, 3) gradients and the number of
  refined points (D3).
- An error raises `ValueError` (bad spec or shape) or `RuntimeError`
  (internal error) with a message. No partial result is returned.
- An empty `points` array returns an empty array.

**Spec of a `FaultChain` (Python dict; shown as JSON).** Arrays are float64
NumPy arrays; lists in JSON below stand for them.

```json
{
  "version": 1,
  "block": {"side_epsilon": 1e-9, "extension": "smooth", "exact": 0.25, "limit": 0.5, "step": 1e-3},
  "faults": [
    {
      "name": "f1",
      "displacement": 100.0,
      "steps": 10,
      "frame": [
        {"origin": [0, 0, 0], "nsteps": [20, 20, 20], "step_vector": [50, 50, 50],
         "local_origin": [0, 0, 0], "local_rotation": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
         "values": "<float64 array, nodes>"},
        "<member 1, same keys>",
        "<member 2, same keys>"
      ],
      "cutters": [2],
      "regions": ["<R>"],
      "profile": {"scale": 0.5, "gx": "<P>", "gy": "<P>", "gz": "<P>"}
    }
  ]
}
```

`local_origin` and `local_rotation` are `null` when the interpolator has no
bounding box. A profile part `P` is one of:

```json
{"kind": "ones"}
{"kind": "zeros"}
{"kind": "smooth_peak"}
{"kind": "cubic", "w": [0, 0, 0, 1], "lim": [-1, 1], "min": -1e308, "max": 1e308, "value_points": [[0, 1], [1, 0]]}
{"kind": "composite", "positive": "<P>", "negative": "<P>"}
```

A region `R` is one of:

```json
{"kind": "positive", "fault": 2, "vector": [1, 0, 0], "point": [500, 500, 500]}
{"kind": "negative", "fault": 2, "vector": null, "point": null}
{"kind": "unconformity", "field": "<frame item>", "faults": [3, 4], "value": 12.5, "sign": true}
```

The spec of a `GridField` is one `frame` item.

#### D2c. Python side: `LoopStructural/modelling/features/_native_faults.py` (new)

```python
REQUIRED_FAULT_KERNEL_VERSION: int = 1

def native_enabled() -> bool
def native_threads() -> int
def chain_for(feature: BaseFeature) -> tuple[object, list] | None
def restore(feature: BaseFeature, points: np.ndarray, *, reverse: bool = False,
            block: dict | None = None, restored_through: Sequence = ()) -> np.ndarray | None
```

- `native_enabled()` (FR-006, FR-006a, FR-006b, C-002, C-008, C-009):
  1. If `LoopStructuralConfig.native_faults` is not None, use it.
  2. Else if the environment variable `LOOPSTRUCTURAL_NATIVE_FAULTS` is
     `0`, `false`, `off` or `no` (any case), the native path is off.
  3. Else on.
  If on, import `loop_interpolation_core`. If the import fails: off, with no
  warning. If `getattr(core, "FAULT_KERNEL_VERSION", None)` is not
  `REQUIRED_FAULT_KERNEL_VERSION`: off, and log one warning for each session
  (a module flag) that names both versions, for example
  `native fault kernel version None found, LoopStructural needs version 1; using the Python path`.
  The result of the import and the version check is kept in the module; the
  setting and the variable are read at each call.
- `native_threads()`: `LoopStructuralConfig.native_threads`, or 0 (all
  cores) if it is None (FR-012, C-006).
- `chain_for(feature)` returns `(FaultChain, faults)` or None (Python path).
  `faults` is `restore_order(feature.faults, owner=feature)` (D1a). It is
  None (not supported, FR-008, C-005, C-010) if any fault of the chain:
  - is not a `FaultSegment`, or its frame members are not three
    `GeologicalFeature` objects with a `DiscreteInterpolator` on a
    `StructuredGrid` with `_has_uniform_lookup()` (the same checks as
    `_shared_frame_members`, without the "no faults" and "same regions"
    checks);
  - has members whose fault lists are not the same objects in the same
    order (these become `cutters`);
  - has a region on a member that is not one of these (OP4): a
    `PositiveRegion` or `NegativeRegion` whose parent is a fault of the
    chain; an `UnconformityFeature` whose interpolator passes the checks of
    `field_for` (D3) and whose faults are all in the chain. Each region
    object is called with `restored_through` only if it has
    `takes_restored_through` (as `_calculate_mask` does); the native region
    follows the same rule;
  - has a `faultfunction` that is not a `FaultDisplacement`, or a profile
    part that is not `Ones`, `Zeros`, `CubicFunction`, `Composite` or the
    function `smooth_peak` (compare with `is`).
  The function calls `builder.up_to_date()` of each fault first, as
  `EvaluationContext.restore` does.
- **State token and cache (FR-009, FR-010).** The kernel is kept in
  `feature._native_chain = (token, refs, chain, faults)`. `token` is a tuple
  with, for each fault of `faults`: `id(fault)`, the displacement and steps
  (through `_hashable`), the `id` of each member interpolator `c` array, the
  grid geometry and bounding box values, the profile parameters (the spec
  dict of D2b without `values`, as a tuple), and the cutter ids. It also has
  `_fault_block.BLOCK_EXTENSION` and the other block settings. `refs` keeps a
  reference to each `c` array, so the ids stay valid (OP8). Each call makes
  the token again (small values only, no array copy). If it differs, a new
  `FaultChain` is made and the entry is replaced with one assignment. Each
  feature keeps its own entry, also when two features have the same faults.
  This is simple; the cost is one copy of the node values of the faults for
  each faulted feature.
- `restore(...)` returns None if `native_enabled()` is False or
  `chain_for(feature)` is None. Else it maps `block` to `sides` (faults not
  named: 2), `restored_through` to `restored` (by identity), uses the order of
  `restore_order(feature.faults, reverse=reverse, owner=feature)` with the
  faults of `restored_through` removed, and calls `FaultChain.restore`.
  Points with `dtype` object, or not shaped (N, 3), return None (Python
  path).

#### D2d. Where the Python path calls the native path

- `BaseFeature._apply_faults`: first
  `native = _native_faults.restore(self, evaluation_points, reverse=reverse, block=block, restored_through=restored_through)`;
  if it is not None, return it; else the present code with `owner=self`.
- `EvaluationContext.restore`: replace `restore_through(points, faults, reverse=reverse)`
  (both calls) with a private helper `_restore_direct(feature, points, reverse)`
  in the same module that tries `_native_faults.restore(feature, points, reverse=reverse)`
  and else calls `restore_through(points, faults, reverse=reverse, owner=feature)`.
  The scope keeps the result as today.
- No other call changes. `FaultSegment.apply_to_points` stays Python. It is
  used by `apply_to_vectors`, by the builders and by a direct user call.

#### D2e. Settings (`LoopStructural/__init__.py`)

Add to the dataclass `LoopStructuralConfig`:

```python
native_faults: bool | None = None   # None: use LOOPSTRUCTURAL_NATIVE_FAULTS, else on
native_threads: int | None = None   # None: all cores
```

Docstring entries in numpy style. A user sets the class attribute, as for
`nelements` (`LoopStructuralConfig.native_faults = False`).

#### D2f. Tests and benchmarks select a path (FR-013)

- The test suite reads the environment variable (D2c), so
  `LOOPSTRUCTURAL_NATIVE_FAULTS=0 pytest` runs the Python path.
- `benchmarks/fault_evaluation/run.py` gets the option
  `--native-faults {on,off}` (default `off`, so stored baselines keep their
  meaning). It sets `LoopStructuralConfig.native_faults`. The result file
  records the value in its `"suite"` object as `"native_faults"`.
- CI: `.github/workflows/native.yml` (it already builds the core) gets a step
  `pytest tests/fault_evaluation` with the native path on.
  `.github/workflows/tester.yml` keeps its run without the core; this is the
  "package not installed" run of SC-006.

### D3. Native field values and gradients (stage S3, FR-004)

New function in `_native_faults.py`:

```python
def field_for(feature: BaseFeature) -> object | None
```

It returns a `GridField` for the interpolator of `feature`, kept in
`feature._native_field = (token, refs, field)` with the token rules of D2c.
None if the native path is off or the interpolator is not a
`DiscreteInterpolator` on a `StructuredGrid` with the uniform lookup, or if
it overrides `evaluate_value`, `_evaluate_value_local` or
`_evaluate_gradient_local` (the checks of `_shared_frame_members`).

Changes in `GeologicalFeature` (`_geological_feature.py`):

- `evaluate_value`: if `field_for(self)` is not None, use
  `field.values(evaluation_points[mask, :], native_threads())` in place of
  `self.interpolator.evaluate_value(...)`. The interpolator first runs
  `_check_frame()` and `update()`, as `_evaluate_frame` does.
- `evaluate_gradient`, branch without faults: the same with
  `field.gradients(...)`.
- `evaluate_gradient`, branch with faults (OP5): if `field_for(self)` and
  `chain_for(self)` are not None, compute `element_scale_parameter` and the
  mask in Python as today, then call
  `chain.feature_gradient(field, pos, order, sides, restored, jumping, element_scale_parameter, native_threads())`.
  `order` is the restore order without the faults of `restored_through`;
  `jumping` is true for each fault of `self._unrestored_faults(restored_through)`
  that the block does not name. Set the rows outside the mask to NaN. If the
  refined count is above 0, log the one warning of the Python path
  (`"<n> points are too close to a fault. Refining the tetrahedron"`), once
  for the call. Else the present Python loop runs.

Masks and regions of the feature itself and `fillnan` of `evaluate_value`
stay in Python (C-003).

### D4. Packages (stage S4, FR-014)

- `.github/workflows/core-wheels.yml`: the smoke test also makes a
  `FaultChain` with one fault (constant profile: `ones` for gx, gy and gz),
  restores two points and checks the result to 1e-12 against the displacement
  along the slip direction. All other jobs (5 platforms, abi3-py39, sdist,
  publish through `pypi.yml`) stay.
- `packages/loop_interpolation_core/README.md` and `CHANGELOG.md` (by
  release-please) describe the fault kernel.

### Names that this plan uses (spec 003 branch, OP1)

`restore_through`, `with_cutters`, `cutters_first`, `_cutters`, `is_in`,
`force_side`, `BLOCK_SIDE_EPSILON`, `BLOCK_EXTENSION`, `BLOCK_EXTENSION_EXACT`,
`BLOCK_EXTENSION_LIMIT`, `BLOCK_EXTENSION_STEP` (`_fault_block.py`);
`EvaluationContext.restore`, `_hashable` (`_evaluation.py`);
`FaultSegment.apply_to_points`, `_evaluate_frame`, `_evaluate_frame_legacy`,
`_displacement_fraction`, `_shared_frame_members` (`_fault_segment.py`);
`FaultDisplacement`, `CubicFunction`, `Composite`, `Ones`, `Zeros`,
`smooth_peak`, `_saturate`, `_taylor`, `continue_profile`,
`block_displacement` (`_fault_function.py`); `BaseFeature._apply_faults`;
`GeologicalFeature.evaluate_value`, `evaluate_gradient`;
`StructuredGrid._has_uniform_lookup`, `_lookup_points`,
`_values_from_coefficients`, `_evaluate_gradient_from_lookup`.

### New names (each checked free with `git grep -n -w` on 2026-11-27)

`restore_order`, `_restore_order`, `_native_faults`, `_native_chain`,
`_native_field`, `native_enabled`, `native_threads`, `native_faults`,
`chain_for`, `field_for`, `_restore_direct`, `REQUIRED_FAULT_KERNEL_VERSION`,
`FAULT_KERNEL_VERSION`, `FaultChain`, `GridField`,
`LOOPSTRUCTURAL_NATIVE_FAULTS`, `KdTree`, `feature_gradient`,
`test_restore_order`, `test_native_faults`.

## Files

### Stage S1

- `LoopStructural/utils/_fault_block.py`: `restore_order`, keyword `owner`
  of `restore_through`.
- `LoopStructural/modelling/features/_base_geological_feature.py`:
  `_apply_faults` passes `owner=self`.
- `LoopStructural/modelling/features/_evaluation.py`: `owner=feature`.
- `LoopStructural/modelling/features/fault/_fault_segment.py`: D1b, only if
  OP3 needs it.
- `tests/fault_evaluation/test_restore_order.py` (new): the order is computed
  once (count calls of `cutters_first` with a mock), a new cutter changes
  the order (story 1 scenario 2), bit-identical values on `cutting_faults`
  and `f4_u3`.
- `specs/007-native-fault-kernel/profile-s1.md` (new): the measurement of
  S1-a, before and after.
- `ROADMAP.md`: the decision to add the native path and stage S1 (FR-015).

### Stage S2

- `packages/loop_interpolation_core/src/fault/mod.rs`, `grid_field.rs`,
  `profile.rs`, `region.rs`, `restore.rs`, `pool.rs` (new). `nearest.rs` and
  `gradient.rs` are stage S3 files.
- `packages/loop_interpolation_core/src/lib.rs`, `src/python.rs`.
- `packages/loop_interpolation_core/tests/fault_grid_field.rs`,
  `tests/fault_profile.rs`, `tests/fault_restore.rs` (new): grid lookup,
  profiles against values from Python, one fault with constant displacement,
  NaN input, empty input. (One file for each Rust module, so that tasks in
  one wave do not share a file.)
- `LoopStructural/modelling/features/_native_faults.py` (new).
- `LoopStructural/modelling/features/_base_geological_feature.py`,
  `LoopStructural/modelling/features/_evaluation.py` (D2d).
- `LoopStructural/__init__.py` (D2e).
- `tests/fault_evaluation/test_native_settings.py` (new): the setting, the
  variable, the version check and the thread count, with a mock kernel
  module (no build needed).
- `tests/fault_evaluation/test_native_safety.py` (new): the story 5
  scenarios (moved out of `test_native_faults.py`, so that two tasks can
  write tests at the same time).
- `tests/fault_evaluation/test_native_faults.py` (new): restore parity on
  each golden model (tolerance and NaN points, story 2), block, reverse and
  `restored_through` (story 2 scenario 2), switch and variable (story 4
  scenarios 2 and 3, with the real kernel). The settings tests are in
  `test_native_settings.py` and the story 5 tests in `test_native_safety.py`.
  Each parity test records the largest difference with
  `record_property("max_native_diff", value)` (FR-005c). Skip the module if
  the kernel is not installed.
- `benchmarks/fault_evaluation/run.py` (D2f).
- `.github/workflows/native.yml` (D2f).
- `specs/007-native-fault-kernel/verification.md` (new): the cause of the
  difference of the proof of concept (C-007), the largest native differences,
  SC-002, SC-003, SC-006 and SC-007 numbers.
- `ROADMAP.md`: stage S2.

### Stage S3

- `LoopStructural/modelling/features/_native_faults.py`: `field_for`.
- `packages/loop_interpolation_core/src/fault/nearest.rs`, `gradient.rs`
  (new), `src/python.rs` (`feature_gradient`).
- `packages/loop_interpolation_core/tests/fault_nearest.rs` (new): `KdTree`
  against a brute-force search.
- `packages/loop_interpolation_core/tests/fault_gradient.rs` (new): the
  tetrahedron and its gradient against values from Python.
- `tests/fault_evaluation/test_native_field.py` (new): the S3 parity tests
  (in place of more tests in `test_native_faults.py`).
- `LoopStructural/modelling/features/_geological_feature.py` (D3).
- (`test_native_field.py`, above): value and gradient parity,
  model ids with the unit-boundary rule (FR-005b), a straddling tetrahedron
  gives NaN (story 3 scenario 2).
- `specs/007-native-fault-kernel/verification.md`: SC-004 and SC-005.
- `ROADMAP.md`: stage S3.

### Stage S4

- `.github/workflows/core-wheels.yml` (D4).
- `packages/loop_interpolation_core/README.md`.
- `.github/workflows/native.yml` is changed in S2 only.
- `specs/007-native-fault-kernel/verification.md`: SC-008 (install on each
  platform, the golden tests with the native path on).
- `ROADMAP.md`: stage S4.

## Tasks of each stage (outline for spec-tasks)

- **S1-a** Measure the 94-fault model of spec 003 at 500 points (OP2): time of each cut
  fault and a profile. Write `profile-s1.md`.
- **S1-b** Tests of D1a (red). **S1-c** D1a. **S1-d** Measure again; D1b
  only if SC-001 fails (OP3). **S1-e** `ROADMAP.md`.
- **S2-a** Read the proof-of-concept crate (OP9) and the Python reference;
  record each difference in `verification.md`.
- **S2-b** Rust `grid_field.rs` and `profile.rs` with Rust tests.
- **S2-c** `LoopStructuralConfig` fields and `native_enabled`,
  `native_threads` with tests (no kernel needed: mock module).
- **S2-d** Rust `region.rs`, `restore.rs`, `pool.rs`, `mod.rs` with Rust
  tests.
- **S2-e** PyO3 binding and `FAULT_KERNEL_VERSION`.
- **S2-f** `chain_for`, `restore`, token and cache; D2d calls.
- **S2-g** `test_native_faults.py` for story 2, 4 and 5.
- **S2-h** Benchmarks: `--native-faults`, SC-002, SC-003, SC-007. CI step.
- **S2-i** `ROADMAP.md`, `verification.md`.
- **S3-a** Rust `nearest.rs` and `gradient.rs` with Rust tests; the
  `feature_gradient` binding.
- **S3-a2** `GridField` Python binding test, `field_for`, D3 changes, tests
  (value, gradient, refinement warning, unconformity region on a frame).
  **S3-b** SC-004, SC-005. **S3-c** `ROADMAP.md`.
- **S4-a** Wheel smoke test. **S4-b** SC-008 install check (the maintainer
  runs it on each platform; CI runs it on the 4 platforms with `smoke: true`).
  **S4-c** `ROADMAP.md`.

## Stages and rollback

| Stage | Merges alone | Rollback |
|---|---|---|
| S1 | Yes, Python only. | Revert the commit. |
| S2 | Yes, after S1 (uses `restore_order`). | `LoopStructuralConfig.native_faults = False` or `LOOPSTRUCTURAL_NATIVE_FAULTS=0`; or revert. |
| S3 | Yes, after S2. | As S2. |
| S4 | Yes, after S2. | Revert the workflow change. |

## Requirement coverage

| Item | Covered by |
|---|---|
| FR-001 | D1a, S1-b, S1-c |
| FR-002 | D1a key, `test_restore_order.py` (story 1 scenario 2) |
| FR-003 | D2a `restore.rs`, `profile.rs`, `region.rs`; S2-g |
| FR-004 | D3, `gradient.rs`, `nearest.rs` (OP5: fully native) |
| FR-005 | S2-g, S3-a parity tests with the spec 003 tolerance |
| FR-005a | S2-g, S3-a (NaN points equal) |
| FR-005b | S3-a (model ids, unit-boundary rule) |
| FR-005c | `record_property("max_native_diff")`, `verification.md` |
| FR-006 | D2c `native_enabled`, D2d |
| FR-006a | D2c version check and one warning |
| FR-006b | D2c, D2e |
| FR-006c | D2c (None returns the unchanged Python path); `tester.yml` run without the core; S1 bit-identical tests |
| FR-007 | No signature change; the existing test `tests/unit/test_public_api_contract.py` |
| FR-008 | D2c `chain_for` checks, D3 `field_for` checks, OP4 |
| FR-009 | D2b (data copied only in `__init__`), D2c cache; story 5 scenario 3 test |
| FR-010 | D2c token; story 5 scenario 4 test |
| FR-011 | D2b GIL release; story 5 scenario 5 test |
| FR-011a | Immutable kernel, single assignment of the cache; story 5 scenario 6 test |
| FR-012 | D2e `native_threads`, `pool.rs`; story 5 scenario 7 test |
| FR-013 | D2f |
| FR-014 | D4 (existing workflow: Linux x86_64 and aarch64, macOS arm64 and x86_64, Windows x64; abi3-py39) |
| FR-015 | S1-e, S2-i, S3-c, S4-c |
| SC-001 | S1-a, S1-d (OP2: case of spec 003; OP3) |
| SC-002 | S2-h |
| SC-003 | S2-h |
| SC-004 | S3-b (OP2: case of spec 003) |
| SC-005 | S3-b |
| SC-006 | D2f: native on (`native.yml`), off (variable), not installed (`tester.yml`) |
| SC-007 | S2-h (OP7) |
| SC-008 | S4-b |

Gaps:

- SC-001 and SC-004 need the 94-fault case of spec 003 (OP2). If spec 003
  merges without it, they cannot be measured.
- FR-004: the nearest-point fill uses the lowest index for equal distances.
  SciPy does not define the result for equal distances, so a tie can give a
  different side check. The parity tests show if this occurs.
- SC-008 on a clean computer: CI covers 4 platforms; macOS x86_64 is
  cross-compiled and not run in CI, so the maintainer checks it by hand.

## Handover

Other specs can use these names. All are private to LoopStructural except the
two settings.

- `LoopStructuralConfig.native_faults: bool | None`,
  `LoopStructuralConfig.native_threads: int | None`, environment variable
  `LOOPSTRUCTURAL_NATIVE_FAULTS`.
- `restore_order(faults, *, reverse=False, owner=None) -> list`
  (`LoopStructural/utils/_fault_block.py`) and the keyword `owner` of
  `restore_through`.
- `_native_faults.native_enabled()`, `native_threads()`, `chain_for(feature)`,
  `field_for(feature)`, `restore(feature, points, *, reverse, block, restored_through)`.
- `loop_interpolation_core.FAULT_KERNEL_VERSION`, `FaultChain(spec)` with
  `restore` and `frame_value`, `GridField(spec)` with `values` and
  `gradients`; the `FaultChain` spec dict of D2b. A change of this dict or of
  the semantics increases `FAULT_KERNEL_VERSION` and
  `REQUIRED_FAULT_KERNEL_VERSION` together.
- Spec 005 changes grid and transform code. If it changes
  `_lookup_points`, `_project_points` or the bounding box local frame, the
  same commit must change `grid_field.rs` and increase the kernel version.
