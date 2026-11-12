# Fault scenario gaps (input for a future spec)

Date: 2026-11-12. Source: the fault scenario examples in
`examples/6_fault_scenarios/` (13 examples), made from the fault meshing
benchmark suite (`fault_mesh_benchmarks/`).

This file is not a spec. It records the gaps that the examples found, with
the evidence and the code location, so that `spec-interview` can start from
it. Line numbers are for commit `f501197c` (origin/master on 2026-11-12).

## Summary

| ID | Gap | Kind | Examples |
|---|---|---|---|
| G1 | No named constant or one-sided displacement profile; profile classes are private | API | 6a, 6c, 6e, 6f, 6g, 6i, 6l, 6m |
| G2 | `FaultDisplacement(scale=0.5)` default halves a one-sided profile | API trap | 6a, 6c |
| G3 | Tip line is a rectangle (product of `gy` and `gz` profiles) | Model | 6b, 6c |
| G4 | Frame scaling differs: `gx` by `minor_axis / 2`, `gy`, `gz` by the radii | Model / docs | 6b, 6j |
| G5 | No separate decay distance for each wall | API | 6c, 6j |
| G6 | Double restoration through a younger fault when an older fault is displaced | Bug (design decision) | 6f, 6g, 6h, 6i |
| G7 | Slip direction of a displaced fault is slow (tetrahedron gradient) | Performance | 6f, 6g, 6h, 6i |
| G8 | `apply_to_points(reverse=True)` does not move points forward | Bug | 6f, 6g, 6h, 6i |
| G9 | No way to give two faults the same age (conjugate faults) | Model | 6g |
| G10 | An unconformity does not cut an older fault surface | Model | 6i |
| G11 | `evaluate_fault_displacements` returns the profile value, not the displacement | Bug / API | 6a |
| G12 | Isosurfaces are smeared across faults; no fault-aware surface extraction | Output | all |
| G13 | Narrow fault blocks need a fine output grid over the whole model (no local refinement) | Output / performance | 6l, 6m |
| G14 | No sampler for fault parameters (stochastic fault models) | Feature | 6k |
| G15 | `steps` (integration steps) is not a parameter of `create_and_add_fault` | API | all |
| G16 | `BaseFault` and `BaseFault3D` profile objects are shared between all faults | Bug risk | 6b, 6d, 6h, 6j, 6k |

## Details

### G1. Named displacement profiles

`FaultSegment.faultfunction` accepts only the strings `"BaseFault"` and
`"BaseFault3D"` (`LoopStructural/modelling/features/fault/_fault_segment.py:72-75`).
A rigid hanging wall (constant slip vector, a pure step) and a hanging-wall-only
profile with a tip line must be built from `FaultDisplacement`, `Composite`,
`Ones`, `Zeros`, `BaseFault.hw` and `smooth_peak`, which are in the private
module `_fault_function`.

The benchmark suite uses four profiles: `constant`, `tip` (one-sided with a
tip), `base` and `base3d`. Suggested: public named profiles (for example
`"Constant"`, `"HangingWallOnly"`, `"HangingWallOnly3D"`) and a public import
path for the profile classes.

### G2. Default `scale` of `FaultDisplacement`

`FaultDisplacement.__init__(..., scale=0.5)`
(`_fault_function.py:318`). The 0.5 is right for the antisymmetric
`BaseFault` profile (each wall moves half). A one-sided profile built with
the default moves the hanging wall only half of the displacement, with no
warning. Example 6a found this: the hanging wall moved 0.075 for a
displacement of 0.15. Suggested: derive the scale from the `gx` profile
(0.5 when both walls move, 1 when one wall moves), or warn.

### G3. Tip line shape

The displacement is `gx(gx) * gy(gy) * gz(gz)`. With `smooth_peak` in both
`gy` and `gz`, the tip line is the rectangle `|gy| = 1` or `|gz| = 1` in the
fault plane, with rounded displacement contours inside (figure in 6b).
Natural tip lines are close to ellipses. Suggested: an optional radial
profile, `f(sqrt(gy^2 + gz^2))`.

### G4. Frame scaling

Measured (benchmark README and 6b): `gz = d / major_axis`,
`gy = d / intermediate_axis`, but `gx = d / (minor_axis / 2)`. Major and
intermediate axes are radii; the minor axis is a full width. Suggested: make
the three consistent, or document the rule in the `create_and_add_fault`
docstring.

### G5. Decay distance for each wall

`BaseFault` decays over the same distance on both sides. A different decay
on the hanging wall and the footwall (for example drag on one side only, or
a rigid footwall with a broad hanging-wall rollover) needs a custom
`CubicFunction` (example 3c). Suggested: parameters such as
`hanging_wall_decay` and `footwall_decay`, or a profile factory.

### G6. Double restoration (design decision needed)

When a series lists faults `[A, B]` and A is built with `faults=[B]`, the
series restores a point through B, then evaluates the frame of A, which
restores the point through B again. Points in the hanging wall of B are
then tested against the wrong position of A.

Evidence (crossing faults, constant profiles, section y = 0.5, 10 201
points; `fault_mesh_benchmarks` probe, worktree code `f501197c`):

| Setup | Points with error > 0.02 | Median error |
|---|---|---|
| Plain LoopStructural | 5.99 % | 1.4e-5 |
| With `restore_once` workaround | 0.01 % | 9.3e-6 |

The examples use `fault_network_workaround(fault)`, which evaluates the frame
of the displaced fault with `faults_enabled = False` inside its
`apply_to_points`. A fix in core must change both the evaluation of a
faulted feature and the restoration of the data in the feature builder,
because both call `_apply_faults`. See also the memory note
"faulted-fault double restoration" and `fault_mesh_benchmarks/README.md`.

### G7. Slip direction of a displaced fault

`GeologicalFeature.evaluate_gradient` uses small tetrahedra that are moved
through the faults of the feature, and refines them while one crosses a
fault (`LoopStructural/modelling/features/_geological_feature.py:209-245`).
`apply_to_points` calls it at every integration step. Measured: restoring
2000 points through a displaced fault takes 15 to 22 s (about 8 to 11 ms a
point). The refinement loop has no limit: a point exactly on a fault never
resolves (the main checkout has an uncommitted fix with a limit of 20).

The workaround uses the gradient of `gy` at the restored point. This is
exact for a rigid younger fault and an approximation in its drag zone.
Related: spec 003 (fault evaluation performance).

### G8. Forward movement through a fault

`FaultSegment.apply_to_points(points, reverse=True)` reverses the
displacement only before the integration loop (`_fault_segment.py:357-359`);
inside the loop `d` is calculated again without the sign
(`_fault_segment.py:403`). So `reverse=True` restores points instead of
moving them forward. The main checkout has an uncommitted fix. Also, where
the frame gradient is NaN (outside the support), `g * d` makes the point NaN
(`_fault_segment.py:412`); the main checkout has an uncommitted fix for this
too.

The examples move restored points forward by fixed-point iteration
(`move_forward`). Suggested: a public forward-modelling API, for example
`GeologicalModel.forward(points, series)`, that is tested against
`restore(forward(p)) == p`.

### G9. Same-age faults

A constant profile moves the whole hanging-wall half-space. In a graben
(6g), each fault is in the hanging wall of the other, so one fault must be
declared younger, and the flanks also move. There is no way to state that
two faults are the same age and that each ends at their intersection (a
conjugate pair). Suggested: a fault relation "conjugate" or "coeval", with
each fault limited to its side of the other.

### G10. Unconformity and older faults

`FeatureRelationshipManager.add_unconformity` skips features of type
`FAULT` (`LoopStructural/modelling/core/_model_relationships.py:171`), and
`add_unconformity_above` returns early for faults (line 123). The surface of
a fault that is older than an unconformity therefore continues into the
younger rocks (6i). Kinematics are correct (a series is moved only by the
faults in its list). Suggested: add the unconformity as a region of the
older fault surface, for output and for abutting tests.

### G11. `evaluate_fault_displacements`

`GeologicalModel.evaluate_fault_displacements`
(`LoopStructural/modelling/core/geological_model.py:1287-1309`) sums
`displacementfeature.evaluate_value`, which is the profile value
(`scale * gx * gy * gz`), not multiplied by the fault displacement. For a
constant fault with displacement 0.15 it returned 0.5. Suggested: return the
displacement (or the movement vector), for each fault or summed.

### G12. Surfaces across faults

Marching cubes on a faulted scalar field smears each horizon across the
cells that a fault cuts (visible in the 3D views of 6b, 6d and 6j). The
benchmark suite solves this with one smooth field for each fault block
(`v_ext`) and clipping by the fault frame. Suggested: fault-block-aware
surface extraction in core, with cut-off lines on the fault surface.

### G13. Local refinement

The interpolation is not the problem. In 6m (a horst 0.022 to 0.058 wide)
the stratigraphy is interpolated in the restored position, where it is
smooth, so a support with cells of 0.05 gives an error below 0.001 in the
horst. A support with 150 000 elements (cell 0.026) took 109 s against
3.8 s and gave no better result.

The problem is the output: isosurfaces are extracted on the grid of
`model.bounding_box`, and a cell that contains the horst also contains both
faults. Only a fine grid over the whole model keeps the horst pieces of the
horizons. Suggested: refinement of the output grid near faults, or the
fault-block surface extraction of G12. For a narrow block in a fault frame
(a fault displaced by another fault), the support of the frame may still
need local refinement; this was not tested.

### G14. Fault parameter sampling

6k draws dip, displacement, minor axis and position by hand, and builds
each model again. Suggested: a sampler that takes distributions of fault
parameters and returns realisations, and that reuses interpolation work
where it can (for example, a fault frame that does not change).

### G15. Integration steps

`FaultSegment` has `steps=10` (`_fault_segment.py:30`), but
`create_and_add_fault` has no `steps` parameter. The examples set
`fault.steps` after the fault is built. One step is exact for a rigid planar
fault and is ten times faster. Suggested: a `steps` parameter, and a default
of 1 for a constant profile on a planar fault.

### G16. Shared profile objects

`BaseFault` and `BaseFault3D` define `hw`, `fw` and `fault_displacement` as
class attributes (`_fault_function.py:397-457`), so every fault that uses the
string `"BaseFault"` shares the same `FaultDisplacement` object. A change to
the profile of one fault changes all of them. Suggested: build a new object
for each fault.

## Examples that depend on the workarounds

| Example | Workaround |
|---|---|
| 6f crossing faults | `fault_network_workaround`, `move_forward` |
| 6g graben | `fault_network_workaround`, `move_forward` |
| 6h folded and faulted | `fault_network_workaround`, `move_forward` |
| 6i unconformity and faults | `fault_network_workaround`, `move_forward` |

When G6 to G8 are fixed, remove the workaround from these four examples and
check the result against the exact solution that 6f and 6g already print.
