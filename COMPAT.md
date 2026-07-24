# Compatibility changelog

Tracks deprecation shims added under the versioning policy in `ROADMAP.md`:
any module-path move/rename on the 1.x line gets a re-export shim with a
`DeprecationWarning`, kept for at least 2 minor releases, logged here.

| Old path | New path | Shim added | Introduced by | Remove after |
|---|---|---|---|---|
| `LoopStructural.datatypes.BoundingBox` | `LoopStructural.geometry.BoundingBox` | 2026-07-24 | `b5eb4742` (unshimmed at the time) | 2 minor releases after the shim ships |
| `LoopStructural.datatypes.Surface` | `LoopStructural.geometry.Surface` | 2026-07-24 | `b5eb4742` | 2 minor releases after the shim ships |
| `LoopStructural.datatypes.ValuePoints` | `LoopStructural.geometry.ValuePoints` | 2026-07-24 | `b5eb4742` | 2 minor releases after the shim ships |
| `LoopStructural.datatypes.VectorPoints` | `LoopStructural.geometry.VectorPoints` | 2026-07-24 | `b5eb4742` | 2 minor releases after the shim ships |

## Migration notices (not breaking, no shim needed)

Not a deprecation shim table entry — the old path keeps working
unchanged, but the new path is preferred going forward.

| Old (still works) | New (preferred) | Since |
|---|---|---|
| `LoopStructural.modelling.features._feature_converters.add_fold_to_feature` (private module, imported directly by the QGIS plugin) | `GeologicalModel.add_fold_to_feature` | 2026-07-24 |
| `LoopStructural.modelling.features._feature_converters.convert_feature_to_structural_frame` | `GeologicalModel.convert_feature_to_structural_frame` | 2026-07-24 |
