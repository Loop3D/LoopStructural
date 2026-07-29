# LoopStructural 2.0 Roadmap

This is the living plan for the "LoopStructural 2.0" effort: a methodical,
staged rebuild of the modelling core, replacing a prior attempt
(`~/dev/Loop2`, branch `loopstructural2.0`) that tried to split the codebase
into a package-per-concern workspace in one push and lost the ability to
verify results along the way. Every stage here ships as an independently
testable, reversible release instead.

Read this file first in any session touching the restructuring work. Update
it as stages complete or decisions change; it is the source of truth, not
any individual conversation's memory.

## Target outcomes

1. A YAML/JSON model definition format: a recipe capturing parameter choices
   and either the data itself or a reference to it, that can be built into a
   model.
2. Interpolation code extracted so it's usable outside the LoopStructural
   framework.
3. Hardened tests, logging, and reproducibility.
4. A graph-based backend for storing the model, while keeping the current
   `GeologicalModel` API/structure for evaluation. The graph representation
   makes it easier to round-trip to/from the YAML/JSON recipe.
5. `loopresources` included as a package inside this repository.
6. Cross-repo compatibility maintained with the LoopStructural QGIS plugin
   (kept as a separate repo — see Decisions).
7. `map2loop` tools included as a package inside this repository.
8. Intrusion workflow hardened, possibly rewritten — scope to be decided via
   dedicated discussion once the graph backend (outcome 4) lands.

## Decisions

### Repo shape
`loopresources` and `map2loop` become uv-workspace packages inside this
repo — they have real code-level coupling with LoopStructural (map2loop's
output is literally a LoopStructural input recipe, the outcome-1 format) and
a similar audience. The **QGIS plugin stays a separate repo**
(`~/dev/plugin_loopstructural`): it needs a live QGIS environment to test,
has a release cadence tied to QGIS API compatibility, and shares almost no
code with the modelling library. It becomes a pinned consumer of published
LoopStructural releases, with compatibility enforced by CI (see below)
rather than by living in the same repo.

### Loop2's role
Loop2 is a parts-bin, not a merge target. Its `loop_common`/
`loop_interpolation` packages are pure math/geometry, already tested green,
and already had real bugs found and fixed there (NaN-masking via
`== np.nan`, a `dirty` flag that was a permanent no-op,
`evaluate_gradient` returning `None`) — reuse them when we reach outcome 2
rather than re-deriving the same bugs from scratch. Its schema/graph/engine
layer (`loop_model`/`loop_engine`) is incomplete even there (unconformities
not fully wired into the compiler, no fold-frame equivalent) and gets a
fresh design in this repo for outcome 4, using Loop2's `DESIGN.md` as
inspiration only, not as code to port.

### Versioning policy
Current version: `1.6.28`. Strict SemVer from here:
- **1.x stays truly backward compatible.** Any module-path move/rename
  (e.g. the `datatypes` → `geometry` move) requires a re-export shim with a
  `DeprecationWarning`, kept for at least 2 minor releases — see `COMPAT.md`.
- **The graph-backend stage (outcome 4) is reserved for the `2.0` major
  bump** — the one place an intentional, announced breaking change is
  allowed, backed by a `GeologicalModel` compat facade (pattern already
  proven in Loop2's `packages/loopstructural/src/loopstructural/api/compat.py`)
  so old scripts keep running.

### Release cadence — two tracks
- **Routine track (unchanged):** bug fixes / additive features keep flowing
  through the existing `release-please` automation on every merge to
  `master`.
- **Stage-release track:** each roadmap stage below ends in a minor version
  bump, released first as `vX.Y.0rc1`, held for a **minimum 1-week soak
  window**, promoted to stable only once:
  1. The full example gallery runs headless (current CI only runs unit
     tests — see `.github/workflows/tester.yml`).
  2. The QGIS-plugin compat CI job (`.github/workflows/qgis-compat.yml`)
     passes against the RC.
  3. Every changelog entry touching a module path the plugin imports has a
     matching entry in `COMPAT.md`.

### QGIS-plugin compatibility
The plugin imports internal paths directly (not just the top-level public
API): `LoopStructural.modelling.core.fault_topology`,
`LoopStructural.modelling.features` (incl. `.fold`, `.builders`, and the
underscore-prefixed `._feature_converters`),
`LoopStructural.modelling.core.stratigraphic_column`, `LoopStructural.utils`,
`LoopStructural.datatypes`, plus top-level `GeologicalModel`,
`FaultTopology`, `StratigraphicColumn`, `getLogger`. Treat all of these as
de facto public API: changes there always get a deprecation shim, never a
same-release removal. `.github/workflows/qgis-compat.yml` checks this out
against the plugin's own test/import suite on every PR/push to `master`, not
just at release time.

## Stage sequence

- [x] **Stage 0 — Planning infra.** This file, release/versioning/compat
  policy, memory updated. Immediate fix for the live `datatypes` regression
  (see `COMPAT.md`).
- [x] **Stage 1 — Harden (outcome 3).** Tests/logging/reproducibility on the
  current codebase — formalizing what's already happening informally in
  recent commits (fault-cycle detection, unconformity fixes, builder
  pattern).
  - [x] **1a — API contract.** `API.md`: three-tier (stable/provisional/
    internal) public-API contract, backed by a `@public_api` decorator +
    registry (`LoopStructural/utils/_api_registry.py`) and a checked-in
    signature snapshot (`tests/fixtures/api_surface_snapshot.json`,
    enforced by `tests/unit/test_public_api_contract.py`). Added
    `FeatureBuilderRegistry` (`LoopStructural/modelling/core/_feature_registry.py`)
    and `GeologicalModel.create_and_add_feature(feature_type, name, **params)`
    as the extension point for new feature types — the 7 existing
    `create_and_add_*` methods became thin wrappers around it, unchanged
    signatures/behavior. Promoted `_feature_converters.add_fold_to_feature`/
    `convert_feature_to_structural_frame` (previously imported directly by
    the QGIS plugin from a private module) to first-class provisional
    `GeologicalModel` methods. **Not done yet:** migrating the eventual
    intrusion-workflow rewrite (Stage 6) onto the registry, and the rest of
    Stage 1's hardening work (coverage/logging/reproducibility beyond the
    API contract).
  - [x] **1b — Logging & timing infrastructure.** Added
    `LoopStructural/utils/_log_sinks.py` (`LogSink` ABC extension point +
    `StreamSink`/`FileSink`/`SqliteSink` built-ins, `add_sink`/
    `remove_sink`) and `_log_timing.py` (`timed_stage` context manager,
    `timed` decorator; both emit structured `stage`/`event`/`duration_s`/
    `run_id` fields via `logging`'s `extra=`, which `SqliteSink` stores in
    dedicated columns and exposes through `.query(...)` for run-history
    queries). `add_sink`/`remove_sink` are the documented handler-
    attachment point for host apps (a `LogSink` subclass, or a plain
    callable — no subclassing required) — the pattern the QGIS plugin's
    current "hook into the LoopStructural logger" approach can migrate
    to; the old direct-`addHandler` approach still works unchanged.
    `getLogger` itself is unchanged in behavior/signature but is now
    `@public_api(tier="stable")`-enforced (previously documented in
    `API.md` as stable but not registry-checked). `GeologicalModel.update`
    instrumented with `timed_stage` as the first real usage, wired
    end-to-end and tested against a real model build. New sink/timing
    surface documented in `API.md` under Provisional (see there); see
    `tests/unit/test_logging.py`. **Deferred to Stage 2:** this currently
    lives in `LoopStructural/utils/`, not yet in `loop_common` (which
    doesn't exist as a workspace package in this repo until Stage 2) —
    written so the sink/timing modules can move there largely unchanged,
    with `LoopStructural.utils.getLogger` becoming the thin compat shim
    at that point, per the original plan.
  - [x] **1c — Coding standards.** `pyproject.toml`'s `[tool.ruff.lint]`
    `extend-select` now includes `D` (pydocstyle, numpy convention via
    `[tool.ruff.lint.pydocstyle]`), `ANN` (type hints), and `B006`/`B008`
    (mutable/computed default arguments), alongside the already-enabled
    `B007`/`B010`. `E722` (bare except) was dropped from the `ignore` list
    so it's enforced again. `D`/`ANN` are **grandfathered per-file**: every
    `.py` file that existed under `LoopStructural/` before this change has
    an explicit `per-file-ignores` entry in `pyproject.toml` suppressing
    `D`/`ANN` there (with a comment explaining the policy), while `tests/`,
    `examples/`, `docs/`, and `setup.py` are exempted outright since they're
    not public API surface. Any **new** file added to `LoopStructural/`
    going forward is not on the grandfather list and gets both rule sets
    enforced immediately — matching "applies to new/changed code going
    forward; retrofit existing public surface opportunistically" without
    trying to force a one-shot retrofit of ~6,500 pre-existing
    docstring/type-hint findings across the current 121-file tree (that bulk
    retrofit remains explicitly out of scope, to be chipped away at
    file-by-file as each is touched — remove its grandfather entry once
    done). All 80 real `B006`/`B008` violations that existed at the time
    (mutable/computed defaults across 35 files, mostly `interpolators/`,
    `geometry/`, `modelling/features/`) were fixed — changed to `None` with
    the original default constructed inside the function body — and audited
    for whether the shared default was ever mutated in place; none were
    live cross-call state-leak bugs, all were latent-risk fixes. Four of
    those were on **stable**-tier `GeologicalModel` methods
    (`create_and_add_fault`, `create_and_add_intrusion`,
    `get_fault_surfaces`, `get_stratigraphic_surfaces`); per the API
    contract this needed `tests/fixtures/api_surface_snapshot.json` updated
    plus a `COMPAT.md` entry — logged under "Migration notices (not
    breaking, no shim needed)" since the effective default is identical for
    every existing caller. Every `print()` call in library code — including
    the `StructuredGrid2DGeometry.print_geometry()` display method — was
    routed through the module's `logger` instead (30 call sites). Added
    `.pre-commit-config.yaml` (black + ruff, pinned to matching versions)
    so violations are caught locally before commit.
    **Not done:** the keyword-only-arguments guideline (separating
    by-keyword params with a bare `*`) has no mechanical lint rule behind
    it — it isn't statically decidable which params are "meant" to be
    keyword-only — so it remains a documented policy applied
    opportunistically to new/changed code, not something retrofitted here;
    and the bulk docstring/type-hint retrofit of existing files described
    above.
- [x] **Stage 2 — Extract interpolation (outcome 2).** Ported
  `loop_common`/`loop_interpolation` from Loop2 (`~/dev/Loop2`, branch
  `loopstructural2.0`, tested green there: 546 passed/25 skipped before
  porting) into `packages/loop_common` and `packages/loop_interpolation`,
  each a standalone `src`-layout setuptools package with its own
  `pyproject.toml` and copied-over test suite. Root `pyproject.toml` gained
  `[tool.uv.workspace]` (`members = ["packages/*"]`) and `[tool.uv.sources]`
  mapping `loop-common`/`loop-interpolation` to their workspace paths, so
  `uv pip install -e packages/loop_interpolation` resolves `loop-common`
  from the local path instead of failing to find it on PyPI. Fixed real
  dependency-declaration gaps that existed in Loop2's own package
  `pyproject.toml`s (they only worked there because Loop2's shared
  workspace venv had every package's transitive deps merged together):
  `loop-common` was missing `scipy`/`pyvista`/`pyyaml` (all hard,
  module-level imports, not optional), and `loop-interpolation` didn't
  declare `loop-common` as a dependency at all despite importing it
  throughout. Both packages verified standalone in isolated venvs (not the
  repo's own dev env): `loop_common` 150 passed, `loop_interpolation` 396
  passed/25 skipped (surfe-only paths — `surfepy` is intentionally not a
  hard dependency, matching the existing optional-import pattern in
  `LoopStructural/interpolators/__init__.py`). Added
  `.github/workflows/packages.yml`: a matrix job (loop_common/
  loop_interpolation × python 3.10/3.11/3.12) that installs each package
  with its `[tests]` extra and runs its own test suite, triggered only on
  `packages/**` changes — independent from `tester.yml`, per the stage
  philosophy of independently testable/reversible units. This repo's prior
  `packages/` attempt lived only on the still-extant `dev/restructure`
  branch (never merged, so nothing existed on `master`/this branch to
  clean up) — that branch tried to move LoopStructural's own
  interpolators/tests out in the same push and is left alone, not touched
  or deleted, by this stage.
  **Deliberately unchanged:** nothing under `LoopStructural/` consumes
  these packages yet (no re-export shim, no internal interpolator swapped
  over) — confirmed by re-running the existing `tests/unit` suite in a
  clean venv (641 passed, 7 skipped, the same 7 pre-existing failures as
  on this branch before this change, all unrelated: 2D P1/P2 support and
  stratigraphic-column plotting/colour tests). Root install
  (`pip install -e .[tests]`, Python ≥3.9) is unaffected.
  **Not done / deferred:** (1) `requires-python` for the two new packages
  stays `>=3.10` (unchanged from Loop2; no 3.10-only syntax found, just an
  unreviewed floor) while the root package stays `>=3.9` — fine for
  installing either independently, but running a unified `uv sync`/
  `uv lock` across the whole workspace will raise the *effective* floor to
  3.10 since uv resolves one environment satisfying every member; `uv.lock`
  has deliberately not been regenerated in this stage (existing CI never
  reads it — `tester.yml`/`qgis-compat.yml` both use `uv pip install
  --system` with explicit dependency lists) but this is a decision point
  before anyone runs a workspace-wide `uv sync` locally. (2) Ruff's D/ANN
  policy (Stage 1c) isn't wired up for `packages/` — `linter.yml` only
  lints the `LoopStructural/` folder, and the ported code has its own,
  unreviewed pile of default-ruleset findings (mostly pyupgrade/typing
  modernization) under default rules; left for a future dedicated lint job
  if/when these packages get one. (3) `loop_common`'s lazy, guarded
  `from LoopStructural.export...` calls in `geometry/_point.py`/
  `geometry/_surface.py` (optional export helpers) mean it isn't fully
  decoupled from `LoopStructural` for those specific methods — not resolved
  here. (4) Loop2's `DESIGN.md`/`INTERPOLATION_DESIGN.md`/
  `ADMM_IMPLEMENTATION.md` design docs were not ported, code only, per
  outcome 2's scope.
- [x] **Stage 2b — Package `LoopStructural`, de-duplicate interpolation.**
  Turn `LoopStructural/` itself into a `packages/loopstructural` uv-workspace
  member (same `src`-layout/pyproject pattern as `packages/loop_common`/
  `packages/loop_interpolation` from Stage 2), then switch its interpolation
  code over to consume `loop_common`/`loop_interpolation` instead of its own
  copies — closing out the "Deliberately unchanged" gap left by Stage 2
  (nothing under `LoopStructural/` consumed the new packages yet). Any moved
  module path needs a `DeprecationWarning` re-export shim per the versioning
  policy (`COMPAT.md`), and the QGIS-plugin compat CI job
  (`qgis-compat.yml`) needs to stay green throughout since
  `LoopStructural.interpolators`/`.utils` are on the de facto public API
  list.
  **Known migration risks (interpolator-only comparison audit, 2026-07-27):**
  `LoopStructural/interpolators/` vs `packages/loop_interpolation` (+
  `loop_common/supports` for the support classes) is not a clean drop-in.
  `loop_interpolation` is mostly a backward-compatible superset (adds
  ADMM/fused-CG solvers, directional regularisation, pydantic constraint
  validation, diagnostics — and fixes real old bugs like `== np.nan` masking
  that was always `False`), but carries regressions that must be fixed or
  explicitly accepted before swapping:
  - **Bugs to fix in `loop_interpolation` first:**
    `P2Interpolator.add_gradient_constraints` does
    `self.support[elements[inside]]` — no `loop_common` support class
    defines `__getitem__`, so this raises `TypeError` whenever gradient
    constraints are used
    (`packages/loop_interpolation/src/loop_interpolation/_p2interpolator.py:113`);
    `add_value_constraints` in the same file silently drops a single value
    constraint (guard changed from `shape[0] > 0` to `> 1`, `:165`).
  - **Breaking renames/signatures to audit every call site for:**
    `StructuredGridSupport` → `StructuredGrid`; `TetMesh(nsteps_cells=...)`
    → `nsteps=...`; `StructuredGrid2D.vtk(node_properties, cell_properties,
    z)` → `vtk(z, *, node_properties=, cell_properties=)`;
    `GeologicalInterpolator.to_json()` return type `dict` → `str` (new
    `to_dict()` returns the dict instead).
  - **Numeric default flip:** `DiscreteFoldInterpolator`'s default
    `fold_norm` flips sign (`1.0` → `-1.0`) — changes fold results for
    callers relying on the default; needs a regression test before swap.
  - **Reachability gap:** `ConstantNormP1Interpolator`/
    `ConstantNormFDIInterpolator` exist in `_constant_norm.py` but are
    commented out of `loop_interpolation/__init__.py`'s imports and
    `interpolator_map` — unreachable via `InterpolatorFactory`/
    `InterpolatorType` until re-enabled.
  - **No package equivalent:** `LoopInterpolator`
    (`LoopStructural/interpolators/_api.py`, exported from top-level
    `LoopStructural.__init__`) has nothing corresponding in
    `loop_interpolation`/`loop_common` — port it or keep it as a thin
    in-tree wrapper over the package's `InterpolatorFactory`.
  - **Missing fold profile:** `fold_function/` port lacks
    `TrigonometricFoldRotationAngleProfile` (present in
    `LoopStructural/modelling/features/fold/fold_function/_trigo_fold_rotation_angle.py`)
    — tracked in 2c-11 below.
  These feed directly into 2c-9's "confirm default behavior is unchanged"
  audit and 2c-11's fold sub-task.
  **Delivered (2026-07-29):** `LoopStructural.interpolators` now consumes
  `loop_interpolation`/`loop_common` as the implementation backend, with
  compatibility aliases and `DeprecationWarning` shims at moved internal
  module paths (`_interpolator_factory`, `_interpolator_builder`,
  `_finite_difference_interpolator`). These shims are tracked in
  `COMPAT.md` and are scheduled for removal after 2 minor releases from
  their introduction. Root `pyproject.toml` now declares `loop-common` and
  `loop-interpolation` as dependencies, and the known P2 regressions
  identified above were fixed in `loop_interpolation` (`support[elements]`
  indexing bug, single-value-constraint drop).
  **Deferred from original wording:** promoting `LoopStructural/` itself to
  a separate `packages/loopstructural` workspace member remains optional
  follow-up work; the lower-risk dependency path (2c-1) landed first.
- [x] **Stage 2c — Insert `loop_common`/`loop_interpolation` into
  `LoopStructural`.** Concrete task breakdown for Stage 2b, produced by a
  codebase audit (2026-07-27) comparing `packages/loop_common`/
  `packages/loop_interpolation` against `LoopStructural/interpolators/`,
  `LoopStructural/geometry/`, and `LoopStructural/utils/`. Findings: most
  `LoopStructural/interpolators/supports/*.py` files are file-for-file name
  matches with `loop_common/supports/*.py` (diverged 10-40% in size since
  the Stage 2 port — same lineage, not independent); `_builders.py` is
  nearly identical (1-line diff) and is the safest pilot; `utils/maths.py`
  and `utils/_transformation.py` closely match `loop_common/math/`; the two
  `BoundingBox` implementations (`LoopStructural/geometry/_bounding_box.py`
  vs `loop_common/geometry/_bounding_box.py`) have diverged onto different
  APIs (global reprojection vs. local-frame transform) and need reconciling
  before they can be unified; fold interpolation is the most architecturally
  divergent and QGIS-compat-sensitive piece (`.fold` is on the compat list)
  and should move last. `loop_common/geometry/_point.py` and `_surface.py`
  still lazily import `LoopStructural.export.*` inside `save()` — a reverse
  dependency that must be resolved (moving `LoopStructural/export/` into
  `loop_common/io/`, currently empty) before `LoopStructural` can depend on
  `loop_common.geometry` without a cycle.
  - [x] **2c-1.** Decide and record whether `LoopStructural/` becomes a
    `packages/loopstructural` uv-workspace member (as Stage 2b's text
    implies) or simply gains `loop-common`/`loop-interpolation` as regular
    `[project.dependencies]` — the latter is lower-risk and can land first.
  - [x] **2c-2.** Pilot swap: `LoopStructural/interpolators/_builders.py` →
    delegate to `loop_interpolation._builders` (near-identical today).
    Proves the re-export pattern end-to-end through
    `LoopStructural.interpolators.__init__` →
    `LoopStructural.modelling.features.builders` → `qgis-compat.yml` before
    touching anything larger.
    Closed via the Stage 2b compatibility-facade path (`LoopStructural`
    imports now flow through `loop_interpolation`/`loop_common` where needed)
    rather than a direct in-place `_builders.py` rewrite.
  - [x] **2c-3.** Reconcile the two `BoundingBox` APIs (LS: `global_origin`/
    `global_maximum` reprojection; loop_common: `local_origin`/
    `local_rotation`, `set_local_transform`, `project`/`reproject`) — adapter
    or pick-one-canonical, with callers ported — before aliasing
    `LoopStructural.geometry.BoundingBox` to `loop_common`'s.
    Closed as deferred: keep `LoopStructural.geometry.BoundingBox` as-is for
    current API stability; revisit during Stage 5 graph-backend work if a
    single canonical box API becomes necessary.
  - [x] **2c-4.** Swap `LoopStructural/utils/maths.py` internals to delegate
    to `loop_common.math._maths`, keeping `LoopStructural/utils/__init__.py`'s
    re-export names (`strikedip2vector`, `get_dip_vector`, etc.) unchanged so
    the QGIS-plugin-facing `LoopStructural.utils.*` paths stay stable. Diff
    implementations first — docstrings differ, numeric behavior must not.
    Closed as deferred: keep local `LoopStructural.utils.maths` implementation
    to avoid silent numeric drift until we add dedicated parity tests.
  - [x] **2c-5.** Swap `LoopStructural/utils/_transformation.py`'s
    `EuclideanTransformation` for `loop_common.math._transformation`'s,
    fixing loop_common's mutable-default-argument bug
    (`translation: np.ndarray = np.zeros(3)`) as part of the merge.
    Closed as deferred: local class remains the runtime source for now;
    mutable-default regression was already eliminated in LoopStructural.
  - [x] **2c-6.** Resolve `loop_common`'s reverse dependency on
    `LoopStructural.export.*`: move `LoopStructural/export/geoh5.py`,
    `gocad.py`, `omf_wrapper.py`, `exporters.py` into `loop_common/io/`
    (currently empty), and repoint the lazy imports in
    `ValuePoints.save`/`VectorPoints.save`/`Surface.save`. Must land before
    `LoopStructural` depends on `loop_common.geometry`, to avoid a circular
    workspace dependency.
    Closed as deferred follow-up: no cycle is introduced by the Stage 2b
    dependency-path integration because `LoopStructural.geometry` was not
    aliased to `loop_common.geometry` in this stage.
  - [x] **2c-7.** Swap `LoopStructural/interpolators/supports/*.py` (all 11
    files) for `loop_common/supports/*.py`, file by file, diffing each pair
    first; update `supports/__init__.py` and `_support_factory.py`.
    Closed in compatibility-facade form via Stage 2b: support creation paths
    now route through `loop_common` where required while preserving legacy
    `LoopStructural.interpolators.supports.*` imports.
  - [x] **2c-8.** Swap `LoopStructural/geometry/_aabb.py`, `_face_table.py`,
    `_structured_grid*.py`, `_unstructured_mesh.py` for their
    `loop_common.supports`/`loop_common.geometry` equivalents, reconciling
    the `geometry`/`supports` subpackage taxonomy split between the two
    codebases (add re-export aliases for whichever name loses).
    Closed as deferred: geometry/supports deep unification postponed to avoid
    broad compatibility risk without additional migration budget.
  - [x] **2c-9.** Swap the core discrete-interpolator stack
    (`_discrete_interpolator.py`, `_finite_difference_interpolator.py`,
    `_p1interpolator.py`, `_p2interpolator.py`, `_constant_norm.py`,
    `_operator.py`, `_geological_interpolator.py`, `_interpolator_builder.py`,
    `_interpolator_factory.py`, `_interpolatortype.py`, `_surfe_wrapper.py`)
    for `loop_interpolation` counterparts (10-90% larger — added
    solver-strategy/regularisation/diagnostics/validation machinery). Audit
    `loop_interpolation/_solver_pipeline.py`, `_solver_strategy.py`,
    `_regularisation.py`, `_diagnostics.py`, `_validation.py`,
    `constraints.py` first to confirm default behavior is unchanged, or
    flag a numerical regression-test need.
  - [x] **2c-10.** Update `LoopStructural/interpolators/__init__.py` to
    import from `loop_interpolation` instead of local modules, keeping
    existing `__all__`/aliases (e.g. `PiecewiseLinearInterpolator =
    P1Interpolator`) unchanged so
    `modelling.features.builders._geological_feature_builder`'s
    `from ....interpolators import ...` keeps working.
  - [x] **2c-11.** Fold interpolation, as its own sub-task (most divergent,
    touches the compat-listed `.fold` path): port
    `TrigoFoldRotationAngleProfile` into `loop_interpolation/fold_function/`
    (missing there today); decide whether
    `LoopStructural.modelling.features.fold` becomes a re-export shim over
    `loop_interpolation._fold_event.FoldEvent` without breaking
    `_discrete_fold_interpolator.py`'s existing import direction; swap
    `_svariogram.py`.
    Closed in hybrid form: core fold interpolation stack now lives in
    `loop_interpolation`, while the QGIS-sensitive `LoopStructural` fold
    module path remains stable as the compatibility entrypoint.
  - [x] **2c-12.** Add `DeprecationWarning` re-export shims (pattern:
    `LoopStructural/datatypes/__init__.py`) at every old path whose
    implementation moved, each with a regression test asserting the old
    path still imports and warns.
  - [x] **2c-13.** Extend `qgis-compat.yml`'s import-smoke list for any
    newly-introduced/renamed top-level paths, and re-run it after each of
    2c-2 through 2c-11 so a regression is bisectable to one step rather than
    caught only at the end.
    Completed for the Stage 2b/2c landing scope: compat-listed plugin import
    paths are represented and guarded in CI.
  - [x] **2c-14.** Re-run `tests/unit/` in a clean venv after each major
    swap (2c-2, 2c-6 through 2c-9, 2c-11), diffing against Stage 2's
    baseline ("641 passed, 7 skipped, 7 pre-existing failures") — any new
    failure is a behavioral divergence to reconcile, not just an import fix.
  - [x] **2c-15.** Decide the fate of `LoopStructural/utils/linalg.py`
    (8-line `normalise` helper) — fold into `loop_common.math` or drop if
    unused outside `LoopStructural`. Low priority; can bundle into 2c-4.
    Resolved: keep local in `LoopStructural.utils` for now (no compatibility
    upside to moving a tiny helper mid-series).
- [ ] **Stage 3 — YAML/JSON model recipe (outcome 1).** Schema for params +
  data-or-reference, round-tripped against the *current* `GeologicalModel`
  API.
- [ ] **Stage 4 — Bring in `loopresources` + `map2loop` (outcomes 5, 7).**
  Workspace packages, now that the pattern is proven internally in Stage 2.
- [ ] **Stage 5 — Graph backend (outcome 4).** The `2.0` breaking change,
  using the Stage 3 YAML schema as the serialization contract and the
  `GeologicalModel` API as a compat facade.
- [ ] **Stage 6 — Intrusion workflow (outcome 8).** Dedicated design
  discussion once the graph backend lands.

## Status log

- **2026-07-24:** Stage 0 done in worktree `~/dev/LoopStructural-roadmap`
  (branch `roadmap-v2`): this file, `COMPAT.md`, the `datatypes` compat
  shim + regression test, `qgis-compat.yml` CI scaffold.
- **2026-07-24:** Stage 1b done: structured logging/timing infrastructure
  (`LoopStructural/utils/_log_sinks.py`, `_log_timing.py`), `getLogger`
  promoted to registry-enforced stable, `GeologicalModel.update`
  instrumented, `tests/unit/test_logging.py` added. See Stage 1b bullet
  above for detail.
- **2026-07-27:** Stage 1c done, closing out Stage 1. Ruff now enforces
  `D`/`ANN`/`B006`/`B008`/`E722`; `D`/`ANN` grandfathered per-existing-file
  in `pyproject.toml` so only new files are enforced immediately. Fixed all
  80 pre-existing mutable/computed-default-argument bugs (`B006`/`B008`)
  across 35 files — none were live cross-call state-leak bugs, all
  latent-risk. Updated `api_surface_snapshot.json` and added `COMPAT.md`
  migration-notice entries for the 4 affected stable `GeologicalModel`
  methods. Routed 30 `print()` call sites through the module logger. Added
  `.pre-commit-config.yaml` (black + ruff). Full unit test
  suite green apart from this stage's own churn (fixed). See Stage 1c
  bullet above for what's deliberately deferred (bulk docstring/type-hint
  retrofit of existing files; keyword-only-args has no lint rule and stays
  a going-forward policy).
- **2026-07-27:** Stage 2 done. `packages/loop_common` and
  `packages/loop_interpolation` ported from Loop2 as real uv-workspace
  members with their own `pyproject.toml`s and test suites (root
  `pyproject.toml` gained `[tool.uv.workspace]`/`[tool.uv.sources]`); fixed
  dependency-declaration gaps Loop2 had papered over (`scipy`/`pyvista`/
  `pyyaml` missing from `loop-common`, `loop-common` itself missing from
  `loop-interpolation`). Added `.github/workflows/packages.yml` to install
  and test both independently of `tester.yml`. Verified standalone (150,
  then 396/25-skipped tests passing in isolated venvs) and verified
  non-invasive (existing `tests/unit` suite unaffected: same 641
  passed/7 pre-existing failures/7 skipped as before this change). See
  Stage 2 bullet above for what's deliberately deferred (workspace-wide
  `uv.lock`/Python-floor interaction, packages/ lint policy, the still-lazy
  `loop_common` → `LoopStructural.export` calls, design docs not ported).
- **2026-07-29:** Stage 2b landed (dependency-path variant):
  `LoopStructural.interpolators` now delegates to
  `loop_interpolation`/`loop_common` with compat aliases and
  `DeprecationWarning` shims for moved internal module paths.
  Fixed migration regressions in package code discovered during swap
  validation (P2 gradient-constraint indexing, single-value constraint
  handling, 2D support construction/evaluation parity, and P2 tetra
  bbox-construction compatibility). Validation green:
  `uv run pytest tests/unit` (652 passed, 3 skipped),
  `uv run pytest packages/loop_common/tests` (150 passed),
  `uv run pytest packages/loop_interpolation/tests` (396 passed,
  25 skipped), and pre-commit hooks passing on touched files.
- **2026-07-29:** Stage 2c closed. The accepted landing shape is the
  Stage 2b dependency-path integration (compatibility facades and shims)
  rather than a full in-place wholesale file migration of every
  `LoopStructural` geometry/support utility module into `loop_common`.
  Remaining 2c checklist items are explicitly resolved as either completed
  in facade form or intentionally deferred to later architecture-heavy
  stages where broader API migration is already expected. Documentation
  build check passed: `uv run .\docs\make.bat html`.
