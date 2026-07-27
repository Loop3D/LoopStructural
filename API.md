# API contract

Companion to `ROADMAP.md` (staged plan) and `COMPAT.md` (deprecation-shim
changelog). This file defines *which* symbols are covered by the
versioning policy in `ROADMAP.md` ("1.x stays truly backward compatible"),
so that both contributors and the eventual Stage 5 graph-backend rewrite
know exactly what they must not break silently.

## Tiers

- **Stable** — covered by the `COMPAT.md` policy: any signature change,
  rename, or move requires a deprecation shim kept for >= 2 minor
  releases, logged in `COMPAT.md`. Checked in CI by
  `tests/unit/test_public_api_contract.py` against
  `tests/fixtures/api_surface_snapshot.json`.
- **Provisional** — works, used internally and/or by early adopters, but
  has not yet gone through a release cycle. No compatibility guarantee
  until promoted to stable.
- **Internal** — anything `_`-prefixed (methods, modules). No guarantee,
  can change or move at any time. The QGIS plugin currently imports one
  internal module directly (`modelling.features._feature_converters`) —
  see "Known internal-path consumers" below; this is being fixed by
  promoting that functionality to a provisional public method rather than
  changing the private module.

### Tier transition policy

- **Provisional -> Stable** requires all of the following:
  1. At least one public release containing the provisional symbol.
  2. No signature churn during that release cycle.
  3. Signature captured by the API snapshot test.
  4. At least one user-facing example or docs reference.
- **Stable -> Removed** requires deprecation handling per `COMPAT.md` and
  must never happen in 1.x without a compatibility shim period.

### What counts as a breaking API change

For **Stable** symbols, all of the following are treated as breaking unless
handled through documented deprecation + shims:
- Renaming, moving, or removing a symbol.
- Changing positional parameter order or required/optional status.
- Renaming parameters used by keyword arguments.
- Changing default argument values in a way that changes behavior.
- Changing return type/shape or units/coordinate conventions.
- Changing exception behavior that callers are expected to handle.

Marked with `@public_api(tier=...)` from `LoopStructural/utils/_api_registry.py`
— a no-op decorator at call time, it just records `(qualified_name,
signature, tier)` for the CI signature-snapshot check. It is deliberately
not a `Protocol`/ABC: forcing the future graph-backend engine (`ROADMAP.md`
Stage 5) to implement a fixed interface class is premature before that
stage's design work happens. A registry + signature-diff test catches
accidental breaks without pre-committing to a rigid shape now.

Contract-test scope note: the snapshot test protects symbol presence and
signatures, not full behavioral equivalence. Behavioral stability must be
covered by unit/integration/example tests for the relevant stable surface.

## Stable surface (as of this policy, 2026-07-24)

`LoopStructural.GeologicalModel`:
- Construction: `__init__`, `from_processor`, `from_file`
- Feature creation: `create_and_add_foliation`, `create_and_add_fold_frame`,
  `create_and_add_folded_foliation`, `create_and_add_folded_fold_frame`,
  `create_and_add_intrusion`, `create_and_add_domain_fault`,
  `create_and_add_fault`
- Unconformities: `add_unconformity`, `add_onlap_unconformity`
- Feature access: `__getitem__`, `__contains__`, `get_feature_by_name`,
  `feature_names`, `fault_names`, `faults` (property)
- Evaluation: `evaluate_model`, `evaluate_model_gradient`,
  `evaluate_feature_value`, `evaluate_feature_gradient`,
  `evaluate_fault_displacements`
- Solve: `update`
- Geometry: `scale`, `rescale`, `regular_grid`, `bounding_box` (property)
- Stratigraphy: `stratigraphic_column` (property), `stratigraphic_ids`
- Output: `get_fault_surfaces`, `get_stratigraphic_surfaces`,
  `get_block_model`, `save`, `to_file`, `to_dict`
- Data: `data` (property)

Other classes already depended on directly by the QGIS plugin, therefore
stable regardless of whether `GeologicalModel` alone would need them
exposed:
- `LoopStructural.StratigraphicColumn`, `LoopStructural.FaultTopology`
- `LoopStructural.modelling.core.fault_topology.FaultRelationshipType`
- `LoopStructural.modelling.core.stratigraphic_column.StratigraphicColumnElementType`
- `LoopStructural.modelling.features.FeatureType`
- `LoopStructural.modelling.features.StructuralFrame`
- `LoopStructural.modelling.features.fold.FoldFrame`
- `LoopStructural.modelling.features.builders.StructuralFrameBuilder`,
  `FaultBuilder`, `GeologicalFeatureBuilder`, `FoldedFeatureBuilder`
- `LoopStructural.geometry.BoundingBox` (+ `Surface`, `ValuePoints`,
  `VectorPoints`)
- `LoopStructural.getLogger`
- `LoopStructural.utils.observer.Observable`

## Provisional surface

- `GeologicalModel.create_and_add_feature(feature_type, name, **params)` —
  new generic dispatch entry point (see below). Existing
  `create_and_add_*` methods are stable wrappers around this; the generic
  method itself is provisional until it's shipped in a release and proven
  stable.
- `GeologicalModel.add_fold_to_feature`,
  `GeologicalModel.convert_feature_to_structural_frame` — promoted from
  the internal `_feature_converters` module (logic unchanged).
- Structured logging/timing (`ROADMAP.md` Stage 1b), all in
  `LoopStructural/utils/_log_sinks.py` and `_log_timing.py`, re-exported
  from `LoopStructural.utils` and top-level `LoopStructural`:
  - `LogSink` — ABC extension point for routing LoopStructural log records
    into a host application's own system (subclass and implement `emit`).
    A plain `Callable[[logging.LogRecord], None]` works too, without
    subclassing.
  - `StreamSink`, `FileSink`, `SqliteSink` — built-in `LogSink`
    implementations. `SqliteSink` stores structured fields (`stage`,
    `event`, `duration_s`, `run_id`) in dedicated columns and exposes
    `.query(...)` for querying run history.
  - `add_sink(sink)` / `remove_sink(handler)` — attach/detach a sink (or
    plain callable) to every current *and future* LoopStructural logger.
    This is the pattern the QGIS plugin's current "hook into the
    LoopStructural logger" approach should migrate to, though the old
    approach (calling `logging.getLogger(name).addHandler(...)` directly
    on a logger returned by `getLogger`) still works unchanged.
  - `timed_stage(logger, stage, **extra)` (context manager) /
    `timed(stage=None)` (decorator) — instrumentation helpers that log a
    `start`/`end` pair with a `duration_s` field around a block or
    function call. `GeologicalModel.update` is instrumented with
    `timed_stage(logger, "update", ...)` as the first real usage.
  - `getLogger` itself is promoted to enforced-**stable** (see below) as
    part of this work — its signature/behavior are unchanged, only
    tracked by the registry now.

## Known internal-path consumers

The QGIS plugin imports
`LoopStructural.modelling.features._feature_converters.add_fold_to_feature`
directly. This is not being broken — `_feature_converters` stays as-is —
but the plugin should migrate to `GeologicalModel.add_fold_to_feature`
when convenient, since that's now the supported, tested path.

Migration target: remove plugin dependence on `_feature_converters` private
imports before the Stage 5 (`2.0`) release candidate window opens.

## Extension mechanism: `FeatureBuilderRegistry`

`LoopStructural/modelling/core/_feature_registry.py` defines
`FeatureBuilderRegistry.register(feature_type: str, builder_factory)`.
Each of the 7 existing feature types is registered against a factory
function extracted unchanged from the corresponding `create_and_add_*`
method body. `GeologicalModel.create_and_add_feature(feature_type, name,
**params)` looks up the factory and calls it — this is the actual
"future-proof, allows for additions" mechanism: a new feature type (e.g.
the eventual intrusion-workflow rewrite, `ROADMAP.md` Stage 6) registers a
factory instead of requiring changes to `GeologicalModel`'s source.
