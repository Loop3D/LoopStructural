# Task: a native "fault snapshot" oracle for kiridashi (proof of concept)

You are working in the **kiridashi** repository (`/home/lgrose/dev/kiridashi`, Rust workspace).
Read its `AGENTS.md` / `CLAUDE.md` and follow the project rules (commands, formatting, lint,
specs). Do not change existing crates' behaviour. Add new code only.

## Goal

Today the field oracle that kiridashi meshes with is a Python object
(`crates/kiridashi-py/python/kiridashi/loop.py`, class `LoopOracle`) that calls
LoopStructural. Each call costs about 30–45 ms whatever the number of points, because of
Python overhead, and the projection phase makes about 2,000 one-point calls per horizon.

Build a **native implementation of the `FieldOracle` trait** (`crates/volmesh/src/field/oracle.rs`)
that evaluates a LoopStructural model from a **snapshot** exported to files. It must give the
same answers as the Python oracle (to the tolerances below) in microseconds per point.

Deliver:
1. A new crate (suggested name `crates/loop-snapshot`, workspace member) with:
   - a loader for the snapshot directory (`manifest.json` + raw `.bin` arrays),
   - a `SnapshotOracle` that implements `FieldOracle` (`values` with optional gradients,
     `fault_distance`, `tip`, `model`, `band`, `has_gradient`, `thread_safety`).
     It must be `Send + Sync` and reentrant (immutable after loading).
2. Validation tests against the reference results (below), one test group per phase.
3. A small benchmark (single point, 100 points, 2,000 points) of the snapshot oracle. Report
   the numbers; compare to the Python figures in "Context numbers".
4. A short `README.md` in the crate: what is supported, what is not, how to re-export.

Work in phases. **Do not start a phase until the previous phase's tests pass.** If a test
cannot pass because Python behaviour is odd, do not loosen the tolerance: write down the
discrepancy (input, Python value, your value, your explanation) in the README and stop that phase.

## Inputs

All inputs are in `/tmp/claude-1000/-home-lgrose-dev-LoopStructural/102a1e6e-6892-4ad1-9f14-903184d259cb/scratchpad/ls003merge/scripts/fault_snapshot/`. Copy the `out/`
directory to `fixtures/snapshot_poc/` in kiridashi (3.6 MB; if the repo has a size policy for
fixtures, ask rather than commit it).

* `out/<model>/manifest.json` — the snapshot. Models: `single_constant` (one fault, constant
  displacement), `single_cubic` (one fault, cubic profiles with slope at the fault, so the
  block continuation is exercised), `two_faults` (two independent faults `main`, `splay`).
* `out/<model>/data/*.bin` — node arrays: little-endian `float64`, C order, no header.
* `out/<model>/reference/index.json` + `*.bin` — results of the Python implementation at
  `reference/points.bin` (N×3, world coordinates; includes points outside the grid, points
  exactly on each fault, and points a hair either side). `index.json` gives each file's shape
  and meaning, and for each series the list of blocks (`series_<s>_blocks`).
* `export_snapshot.py` — the exporter (do not modify). `verify_export.py` — a 60-line numpy
  reader that reproduces the frame and series values from the files; use it as a worked example
  of the grid format. It agrees with the reference to 4e-16.

### Manifest (version 1)

* `constants`: `BLOCK_SIDE_EPSILON`, `BLOCK_EXTENSION` (`"smooth"`), `BLOCK_EXTENSION_EXACT`,
  `BLOCK_EXTENSION_LIMIT`, `BLOCK_EXTENSION_STEP`. Read these, never hard-code them.
* `gradient`: `{"method": "regular_tetrahedron", "scale": s}`.
* `faults[]`: `name`, `displacement` (model units, use as is), `steps`, `profile`
  (`scale`, `gx`, `gy`, `gz`; each is `ones`, `zeros`, `cubic{w,lim,limit_value,min_v,max_v}` or
  `composite{positive,negative}`), `frame` = three grids (gx, gy, gz), `frame_faults` (names of
  younger faults that displace this fault; may be empty), `has_tip`, `tip_ends`.
* `series[]`: `name`, `faults` (names, **youngest first** = restore order), `grid`.
* A grid: `origin`, `step`, `ncells` (3), `local_from_world` (`A` 3×3 row-major, `b`) so that
  `local = A·world + b`, `world_vector_from_local` (3×3 row-major) for gradient directions,
  `nodes` (relative path). Node index = `i + (nx+1)·(j + (ny+1)·k)`, `nx = ncells[0]`.

The exporter refuses anything else (regions/abutting faults, non-structured grids, other profile
types, a fault displaced by a fault that is not earlier in the series restore order). Your loader
should reject unknown `schema`/`version` and any structure it does not understand with a clear error.

## What to port (read these Python sources first; they are the specification)

Base path: `/tmp/claude-1000/-home-lgrose-dev-LoopStructural/102a1e6e-6892-4ad1-9f14-903184d259cb/scratchpad/ls003merge/`. The semantics below are the contract;
where this text and the Python differ, the Python wins — tell me.

1. **Grid evaluation** — `packages/loop_common/src/loop_common/geometry/_structured_grid_3d.py`
   (`position_to_cell_index`, `inside`, `_inside`, `position_to_local_coordinates`, `trilinear`,
   `evaluate_value`) and `.../supports/_3d_structured_grid.py` (`evaluate_gradient`,
   `_shape_gradients`). Points outside the grid give NaN (inside test has a tolerance
   `eps·32·max(|max−origin|,1)` per axis). Cell index is `clip(floor((l−o)/h), 0, n−1)`; local
   coordinates are clipped to [0,1]. The discrete interpolator wrapper
   (`packages/loop_interpolation/src/loop_interpolation/_discrete_interpolator.py`,
   `_evaluate_value_local`) returns **0.0 for a point with a NaN coordinate** and NaN outside the
   grid. Gradients are evaluated in the local frame and rotated with
   `world_vector_from_local`.
2. **Fault restore** — `LoopStructural/modelling/features/fault/_fault_segment.py`:
   `apply_to_points` (lines ~527–617), `_evaluate_frame`, `_displacement_fraction`,
   `evaluate_frame_value/_gradient`, `_evaluate_tip`/`evaluate_tip_value` (~467–513).
   Note: the mask of moved points is fixed from the first frame evaluation; gradients of gy are
   normalised and scaled by `(1/steps)·d`; NaN propagates exactly as in numpy; `block` side `0`
   returns the points unchanged.
3. **Profiles and the block extension** — `.../fault/_fault_function.py`: `Ones`, `Zeros`,
   `Composite` (note `v == 0` gives 0), `CubicFunction.__call__` (limits, floor, ceiling),
   `FaultDisplacement.__call__` (`scale·gx·gy·gz`), `_saturate`, `_taylor`, `continue_profile`,
   `block_displacement`; and `LoopStructural/utils/_fault_block.py` (`force_side`,
   `block_side`, constants). For a block that names a fault, points on the other side of that
   fault get the smooth continuation of the block's own profile (for a cubic: its own cubic at the
   saturated signed coordinate; for the constant profiles: the constant).
4. **Feature evaluation** — `LoopStructural/utils/_fault_block.py::restore_through` (each fault
   once, in series order, each told which faults are already done) and
   `LoopStructural/modelling/features/_geological_feature.py` (`evaluate_value`,
   `evaluate_gradient`). With a block that names every fault of the series (the oracle always does)
   the gradient is the P1 gradient over a regular tetrahedron of size `gradient.scale` centred on
   the point, its four vertices restored through the faults exactly like points
   (`LoopStructural/utils/maths.py`: `regular_tetraherdron_for_points`,
   `gradient_from_tetrahedron`). The "points too close to a fault" refinement loop only runs for
   faults *not* named in the block, which does not happen here.
5. **The oracle contract** — `crates/kiridashi-py/python/kiridashi/loop.py::LoopOracle`:
   `values(field, block, points, want_gradient)` uses only the faults of the field and passes
   `block[fault_index]` for them (+1 hanging wall, −1 footwall, 0 absent);
   `fault_distance` = `gx / |∇gx|` (NaN where the norm is 0); `tip` = `evaluate_tip_value`.
   Map `FaultSide::{HangingWall, Footwall, Absent}` to +1/−1/0 as `kiridashi-py/src/convert.rs`
   does.

## Phases and tests

Tolerances are absolute unless stated. NaN positions must match exactly (compare the NaN masks
first, then the values).

| Phase | What | Reference files | Tolerance |
|---|---|---|---|
| 0 | loader; grid value and gradient of every frame coordinate | `fault_<f>_frame_g{x,y,z}`, `fault_<f>_frame_grad_g{x,y,z}` | 1e-12 |
| 1 | `apply` and reverse, no block, single fault | `fault_<f>_apply`, `fault_<f>_apply_reverse` | 1e-12 (all models, each fault) |
| 2 | series restore chain and series value, no block | `series_<s>_restored_after_<j>`, `series_<s>_value` | 1e-12 |
| 3 | block values, **all** blocks listed in `series_<s>_blocks` | `series_<s>_block_<key>_value` | 1e-10 |
| 4 | block gradients (tetrahedron method), all blocks | `series_<s>_block_<key>_gradient` | 1e-8 |
| 5 | `fault_distance`, `tip` | `fault_<f>_distance`, `fault_<f>_tip` | 1e-12 (`inf` must match) |
| 6 | `FieldOracle` impl wired to `volmesh` types; `ModelInfo` built from the manifest | — | see below |

Phase 6: build `ModelInfo` from the manifest (domain, faults with `has_tip`, one field per
series with its fault ids). Add a test that calls the trait methods and compares with the same
references, including `band` behaviour and `has_gradient == true`. If a `volmesh` test helper
(`crates/testkit/src/oracles.rs`) shows how oracles are tested, follow it.

Extra tests (write them, small and fast):
* **Hand-made snapshot**: build a tiny snapshot in code with no Python: a planar fault
  (frame node arrays linear in x, y, z), constant profile, field linear in z. The restored value
  and gradient are known analytically; test the loader-free path.
* **Reentrancy**: evaluate from several threads with different blocks; results equal the serial ones.
* **Points on the fault**: every point in `points.bin` within 1e-3 of a fault is checked in all
  phases (they are in the reference); add a test that asserts the count is nonzero so nobody
  silently drops them.
* **Error paths**: wrong `schema`, wrong `version`, truncated `.bin`, node array of the wrong length.

Benchmark (`cargo bench` or an ignored test that prints numbers): `values` with gradient for 1,
100 and 2,000 points on `two_faults`, one block. Python reference figures for orientation:
a LoopStructural single-point value+gradient call is 31–47 ms in the real model and about
4–12 ms on these small fixtures; a prototype of a snapshot `apply` in Rust measured 2–3 µs per
moved point. Report whatever you measure; do not tune the tolerances to make it faster.

## Constraints

* Native Rust only for evaluation: no Python at test time (tests read the files). Use
  `serde`/`serde_json` only if the workspace already depends on them; otherwise a small parser
  and ask before adding a dependency.
* Evaluation must not allocate per point in the inner loops (allocate per call).
* Do not modify the exporter, the reference files or anything under LoopStructural; if you think
  the exporter is missing data, write the request in the README.
* Keep the schema version check strict. This is a proof of concept: prefer clear errors over
  partial support.
* When finished, report: phases passed, tolerances actually achieved (max abs error per
  reference file), benchmark numbers, and every place where you had to interpret the Python.

## Regenerating the inputs

```
cd /tmp/claude-1000/-home-lgrose-dev-LoopStructural/102a1e6e-6892-4ad1-9f14-903184d259cb/scratchpad/ls003merge
PYTHONPATH=. uv run --frozen --no-dev --with pytest python scripts/fault_snapshot/export_snapshot.py scripts/fault_snapshot/out
uv run --frozen --no-dev python scripts/fault_snapshot/verify_export.py scripts/fault_snapshot/out/two_faults
```
