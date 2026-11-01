# Research: Faster evaluation of faulted models

This file has the analysis of the present code and the baseline numbers for
[spec.md](spec.md). Line numbers are for commit `4a4bc499` (release 1.11.0).

## How a feature is evaluated today

`GeologicalFeature.evaluate_value`
(`LoopStructural/modelling/features/_geological_feature.py:116`) does these
steps for N points:

1. It builds the feature if it is not up to date (`self.builder.up_to_date()`).
2. It makes the mask: `_calculate_mask`
   (`_base_geological_feature.py:229`) calls each region with the points that
   are **not** restored.
3. It restores the points: `_apply_faults` (`_base_geological_feature.py:252`)
   calls `apply_to_points` of each fault of the feature, in order.
4. It evaluates the field at the restored points that are in the mask.

`FaultSegment.apply_to_points` (`features/fault/_fault_segment.py:315`):

- evaluates the three coordinates of the fault frame at all the points;
- then moves the points in `steps` steps (default 10, line 30). In each step
  it evaluates the three coordinates and the gradient of coordinate 1 at the
  points that the fault moves;
- each coordinate of the frame is itself a `GeologicalFeature`, so it applies
  its own faults (an older fault that cuts this fault) and its own regions.

Threads are off (`use_threads = False`, line 16), because of a known bug.

## Where the repeated work comes from

- `add_unconformity` (`core/_model_relationships.py:120`) makes an
  `UnconformityFeature` for the unconformity and adds it as a region of each
  older feature, until the next unconformity.
- `add_unconformity_above` (`core/_model_relationships.py:95`) adds each
  unconformity above a new feature as a region of the feature.
- `UnconformityFeature` (`features/_unconformity_feature.py`) has the faults
  of its parent and no regions. Its `evaluate` (line 62) calls
  `evaluate_value`, so it applies all the faults of the parent to the points
  again.
- When the parent of the unconformity has the same fault set as the feature
  (the usual case: the faults cut all the groups), the restored points are the
  same as in step 3. The work is done again for nothing.
- `GeologicalModel.evaluate_model` (`core/geological_model.py:1162`) and
  `evaluate_model_gradient` (line 1237) evaluate each group with a separate
  call, so the work of each group is done again for the next group.
- A fault that stops at an unconformity gets the unconformity as a region
  (`core/_model_feature_factory.py:772`). Then the groups above and below the
  unconformity have different fault sets.
- `evaluate_gradient` of a faulted feature (`_geological_feature.py:159`)
  makes a small tetrahedron around each point and applies the faults to its
  four corners, so it has four times more points to restore.

## Baseline 1: the plugin model

A model in the QGIS plugin: one group (`Group_0`), six faults, one
unconformity region on the group. Measured in QGIS 3.42 (Python 3.12,
NumPy 1.26.4, LoopStructural 1.11.0), 20 000 random points in the model box.

| What | Time |
|---|---|
| Field of `Group_0` only (`interpolator.evaluate_value`) | 0.01 s |
| `Group_0.evaluate_value` (with the faults) | 8.2 s |
| Coordinate 0 of one fault | 0.03 s |
| Displacement of one fault | 0.08 s |

Profile of `Group_0.evaluate_value` (20 000 points, cumulative time):

| Cumulative | Calls | Function |
|---|---|---|
| 8.80 s | 12 | `FaultSegment.apply_to_points` (6 faults, two times) |
| 3.81 s | 502 | `_calculate_mask` |
| 3.65 s | 1 | `UnconformityFeature.evaluate` (applies the 6 faults again) |
| 3.75 s | 120 | `GeologicalFeature.evaluate_gradient` (fault frame, in the steps) |
| 3.67 s | 382 | `StructuredGrid.evaluate_value` |
| 2.48 s | 1 244 | `position_to_cell_index` |
| 2.19 s | 120 | `get_element_gradient_for_location` |

Notes:

- More than 99 % of the time is in the faults.
- The 6 faults are applied two times: once for the group and once for the
  unconformity region. The two results are the same.
- Inside one fault application, most of the time is in the grid lookup of
  `loop_common` (`position_to_cell_index`, `position_to_cell_corners`,
  `get_element_gradient_for_location`). This is the cost of user story 5.

## Baseline 2: synthetic models

Script: [measure_baseline.py](measure_baseline.py). It builds models with F
vertical faults and U unconformities (U + 1 groups, all cut by all the faults,
3 000 elements for each feature), and evaluates 5 000 random points. It counts
the calls of `FaultSegment.apply_to_points`.

Computer: Intel Core i7-12700H, 32 GB, Windows 11, Python 3.14.4,
NumPy 2.4.6, SciPy 1.18.1, LoopStructural 1.11.0 (`solver='cg'`).

| F | U | Oldest group: time | Oldest group: fault applications | All groups: time | All groups: fault applications |
|---|---|---|---|---|---|
| 0 | 0 | 0.00 s | 0 | 0.00 s | 0 |
| 0 | 3 | 0.01 s | 0 | 0.02 s | 0 |
| 1 | 0 | 0.18 s | 1 | 0.21 s | 1 |
| 1 | 1 | 0.33 s | 2 | 0.64 s | 4 |
| 1 | 2 | 0.43 s | 3 | 1.38 s | 8 |
| 1 | 3 | 0.67 s | 4 | 2.05 s | 13 |
| 2 | 0 | 0.38 s | 2 | 0.41 s | 2 |
| 2 | 1 | 0.68 s | 4 | 1.44 s | 8 |
| 2 | 2 | 1.03 s | 6 | 2.80 s | 16 |
| 2 | 3 | 1.40 s | 8 | 4.54 s | 26 |
| 4 | 0 | 0.69 s | 4 | 0.67 s | 4 |
| 4 | 1 | 1.17 s | 8 | 2.62 s | 16 |
| 4 | 2 | 1.93 s | 12 | 5.12 s | 32 |
| 4 | 3 | 2.55 s | 16 | 8.65 s | 52 |

Notes:

- Without faults, unconformities cost almost nothing. The cost comes from the
  faults that each unconformity region applies again.
- Oldest group: the number of fault applications is F × (U + 1). It is F when
  each fault set is applied once.
- All groups: the number grows with the square of U (52 for F = 4, U = 3). It
  is F when each fault set is applied once, because all the groups have the
  same fault set in these models.
- One fault application costs about 0.17 s for 5 000 points, which is about
  34 µs for each point.
- The time is about proportional to the number of fault applications, so
  removing the repeated applications gives most of the gain of user story 2.

## The plugin model as a golden model

`data/plugin_model/` has the plugin model of baseline 1, as it is:

- `recipe.json`: the model box and the stratigraphic column, and the
  LoopStructural calls that the plugin made, in order. They were recorded from
  the plugin model manager on 2026-11-01: six `create_and_add_fault` calls
  (displacement 1000, dip 90), one `add_abutting_fault` (7439 abuts 3496), one
  `create_and_add_foliation` (Group_0) and one `add_unconformity`.
- `<feature>.csv`: the data of each feature.
- `fields.npz` and `fields.json`: the solved field and the grid of each
  interpolator (three for each fault). 1.2 MB.
- `qgis_check.npz`: values of the model in QGIS at 5000 fixed points.
- `plugin_model.py`: builds the model with the stored fields (no solve), or
  solves it again. Run it to compare with `qgis_check.npz`. The result is the
  same as in QGIS (difference 0, same NaN points).

Things that a golden model with stored fields must do, found when this model
was made:

- Use the box of the model with its local transform. The plugin gives a
  `BoundingBox` with `local_origin = 0`. `GeologicalModel(origin, maximum)`
  sets the local origin to the box origin, and then the points miss the grids
  of the faults.
- Set the grid of each fault. The build makes the grid of a fault larger (the
  fault buffer), so the grid before a build is smaller than the saved field.
- Give each feature the faults of its builder (`BaseBuilder.update` does this
  before a build), but only after all the calls: `add_unconformity` must see
  the feature before its first build.
- Mark each interpolator and each builder as up to date, so that an
  evaluation does not solve again.

## Notes for the plan

These notes are input for the plan. They are not requirements.

- The restored points depend only on the points and on the fault set. A
  feature and its unconformity regions can share them in one evaluation call
  when their fault sets are equal.
- The mask uses the points that are not restored, but each region restores
  the points with its own fault set. So the mask needs the restored points of
  the fault set of each region.
- `evaluate_model` can evaluate all the groups in one pass, and share the
  restored points of each fault set between the groups and their regions.
- A cache that lives longer than one evaluation call needs a reliable sign
  that the model changed (for example a version number that each builder
  changes when it builds). See C-003 in the spec.
- The present tests have no fixed reference for the model output, so the
  golden tests (user story 1) must come first.
