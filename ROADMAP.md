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
- [ ] **Stage 1 — Harden (outcome 3).** Tests/logging/reproducibility on the
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
  - [ ] **1c — Coding standards.** All functions must have docstrings.
    All function arguments meant to be passed by keyword must be
    keyword-only, separated from positional arguments with a bare `*` in
    the signature, to prevent positional contamination (callers passing
    by position and silently breaking when parameter order changes).
    Applies to new/changed code going forward; retrofit existing public
    surface opportunistically, but changing a **stable** (`API.md`)
    signature from positional-or-keyword to keyword-only is itself a
    breaking change per the API contract and needs a deprecation shim,
    not a silent edit. Also:
    - **Enforce docstrings via ruff's `D` (pydocstyle) rules.**
      `pyproject.toml` already has a `[tool.pydocstyle]` numpy-convention
      block, but it isn't wired into `ruff.lint.extend-select` so nothing
      currently checks it — add `D` to `extend-select` so missing/malformed
      docstrings fail CI instead of the config sitting unused.
    - **Type hints on public signatures.** Parameters and return values on
      public (`API.md` stable/provisional) functions must be typed; add
      mypy or ruff's `ANN` rules to check it.
    - **No mutable default arguments.** Enable ruff `B006`/`B008` to catch
      mutable defaults and function-call defaults.
    - **Re-enable bare-except lint.** Drop the `E722` entry from the
      `ignore` list in `[tool.ruff.lint]` (currently marked "temporary")
      so broad bare excepts get flagged again.
    - **No `print()` for diagnostics.** Route through the logger instead —
      depends on the 1b logging infrastructure landing first.
    - **Pre-commit hook.** Add `.pre-commit-config.yaml` running
      black/ruff locally, so violations are caught before commit instead
      of only after push via the auto-fix-PR bot in `linter.yml`.
- [ ] **Stage 2 — Extract interpolation (outcome 2).** Port
  `loop_common`/`loop_interpolation` from Loop2 into a real uv workspace
  under `packages/`, with CI that actually installs and tests it (this
  repo's `packages/` directory currently exists but is empty/orphaned from
  an earlier abandoned attempt on branch `dev/restructure` — redo it
  properly).
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
