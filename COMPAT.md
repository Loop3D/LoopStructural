# Compatibility changelog

Tracks deprecation shims added under the versioning policy in `ROADMAP.md`:
any module-path move/rename on the 1.x line gets a re-export shim with a
`DeprecationWarning`, kept for at least 2 minor releases, logged here.

| Old path | New path | Shim added | Shim first released in | Introduced by | Earliest removal version | Removed in | Owner |
|---|---|---|---|---|---|---|---|
| `LoopStructural.datatypes.BoundingBox` | `LoopStructural.geometry.BoundingBox` | 2026-07-24 | TBD (next release containing this shim) | `b5eb4742` (unshimmed at the time) | TBD (2 minor releases after first shim release) | Active | Core maintainers |
| `LoopStructural.datatypes.Surface` | `LoopStructural.geometry.Surface` | 2026-07-24 | TBD (next release containing this shim) | `b5eb4742` | TBD (2 minor releases after first shim release) | Active | Core maintainers |
| `LoopStructural.datatypes.ValuePoints` | `LoopStructural.geometry.ValuePoints` | 2026-07-24 | TBD (next release containing this shim) | `b5eb4742` | TBD (2 minor releases after first shim release) | Active | Core maintainers |
| `LoopStructural.datatypes.VectorPoints` | `LoopStructural.geometry.VectorPoints` | 2026-07-24 | TBD (next release containing this shim) | `b5eb4742` | TBD (2 minor releases after first shim release) | Active | Core maintainers |
| `LoopStructural.interpolators._finite_difference_interpolator` | `loop_interpolation._finite_difference_interpolator` | 2026-07-29 | 2026-07-29 | extracted-package integration | 2 minor releases after first shim release | Active | Core maintainers |
| `LoopStructural.interpolators._interpolator_builder` | `loop_interpolation._interpolator_builder` | 2026-07-29 | 2026-07-29 | extracted-package integration | 2 minor releases after first shim release | Active | Core maintainers |
| `LoopStructural.interpolators._interpolator_factory` | `loop_interpolation._interpolator_factory` | 2026-07-29 | 2026-07-29 | extracted-package integration | 2 minor releases after first shim release | Active | Core maintainers |

## Deprecation lifecycle

Every shim tracked in this file follows the same lifecycle:
1. **Announce:** document in release notes and this table.
2. **Warn:** emit runtime `DeprecationWarning` from the old path.
3. **Guard:** keep plugin-compat CI and API contract checks green.
4. **Schedule removal:** set `Earliest removal version` once the first
	shim-containing release is cut.
5. **Remove:** only after the earliest removal version is reached and known
	downstream consumers are migrated.

If a removal is postponed, update `Earliest removal version` with a short
reason in the release notes.

## Migration notices (not breaking, no shim needed)

Not a deprecation shim table entry — the old path keeps working
unchanged, but the new path is preferred going forward.

| Old (still works) | New (preferred) | Since |
|---|---|---|
| `LoopStructural.modelling.features._feature_converters.add_fold_to_feature` (private module, imported directly by the QGIS plugin) | `GeologicalModel.add_fold_to_feature` | 2026-07-24 |
| `LoopStructural.modelling.features._feature_converters.convert_feature_to_structural_frame` | `GeologicalModel.convert_feature_to_structural_frame` | 2026-07-24 |
| `GeologicalModel.create_and_add_fault(..., faults=[])` | `GeologicalModel.create_and_add_fault(..., faults=None)` (list built internally, same effective default) | 2026-07-27 |
| `GeologicalModel.create_and_add_intrusion(..., intrusion_frame_parameters={}, geometric_scaling_parameters={})` | `GeologicalModel.create_and_add_intrusion(..., intrusion_frame_parameters=None, geometric_scaling_parameters=None)` (dicts built internally, same effective default) | 2026-07-27 |
| `GeologicalModel.get_fault_surfaces(faults=[])` | `GeologicalModel.get_fault_surfaces(faults=None)` (list built internally, same effective default) | 2026-07-27 |
| `GeologicalModel.get_stratigraphic_surfaces(units=[])` | `GeologicalModel.get_stratigraphic_surfaces(units=None)` (list built internally, same effective default) | 2026-07-27 |
| `FaultBuilder.__init__`/`FoldedFeatureBuilder.__init__`/`StructuralFrameBuilder.__init__` `bounding_box` param annotated as `loop_common.geometry._bounding_box.BoundingBox` | annotated as `LoopStructural.geometry._bounding_box.BoundingBox` — `LoopStructural.geometry.BoundingBox` reverted to the local implementation (see `API.md`); accepted argument type is unchanged, only the class's canonical module path | 2026-07-30 |

## Compatibility debt summary

- Active shims: 7
- Oldest active shim added: 2026-07-24
- Next cleanup milestone: remove the interpolator-path shims after 2 minor releases from 2026-07-29
