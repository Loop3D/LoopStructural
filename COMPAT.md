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

## Compatibility debt summary

- Active shims: 4
- Oldest active shim added: 2026-07-24
- Next cleanup milestone: set after first shim-containing release is tagged
