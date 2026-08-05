# Intrusions Module — Review, User Guide & Hardening Plan

This document is the "dedicated discussion" input for **ROADMAP.md's Stage 6
(outcome 8: intrusion workflow hardened, possibly rewritten)**. It has four
parts:

1. **How the module works today** — architecture summary.
2. **User guide** — what data and parameters are actually required to build
   an intrusion, worked from the one working example in the codebase.
3. **Review findings** — concrete bugs, dead code, and design smells, each
   with a file:line.
4. **Simplification & hardening plan** — a phased proposal for Stage 6.

Everything here was verified against the code as of commit `160c6e2a`
(2026-08-05), not against docs or memory of older versions.

**Status (2026-08-05): Phases A-C below are done.** The fixes, the data-
contract validation, and new regression tests (including one that exercises
`marginal_faults` end-to-end for the first time, and one that pins the newly-
discovered `intrusion_steps` breakage — see finding 1) are all in
`tests/unit/modelling/intrusions/test_intrusions.py`. Findings are left
in place below as the historical record of what was found and fixed;
each fixed item is marked **[FIXED]**. Phase D (the larger
simplification/rewrite) is still the open item for the dedicated Stage 6
discussion.

---

## 1. Architecture summary

Files in `LoopStructural/modelling/intrusions/`:

| File | Role |
|---|---|
| `intrusion_frame.py` | `IntrusionFrame(StructuralFrame)` — marker subclass, no added behaviour. |
| `intrusion_frame_builder.py` | `IntrusionFrameBuilder(StructuralFrameBuilder)` — builds the curvilinear coordinate system (`c0`=growth away from the reference contact, `c1`=position along strike, `c2`=lateral/side distance), optionally conditioned by fault "steps" and "marginal faults". |
| `intrusion_builder.py` | `IntrusionBuilder(BaseBuilder)` — given the frame, prepares per-side/per-contact data and fits RBF interpolators over `(c1, c2)` for lateral thresholds and `(c1, c2)` for vertical thresholds, constrained by a "conceptual model" function. |
| `intrusion_feature.py` | `IntrusionFeature(BaseFeature)` — evaluates a signed-distance-like scalar field to the intrusion contact from the frame + the fitted thresholds. |
| `geom_conceptual_models.py` | Three conceptual-geometry functions: `ellipse_function`, `constant_function`, `obliquecone_function`. |
| `geometric_scaling_functions.py` | Empirical length→thickness scaling from the literature (pluton/laccolith/sill scaling laws) — see bugs, this is **not wired up**. |
| `intrusion_support_functions.py` | Fast-marching "shortest path" helpers. **Dead code** — see §3. |

Entry point: `GeologicalModel.create_and_add_intrusion(...)` →
`_build_intrusion(...)` (`LoopStructural/modelling/core/geological_model.py:1373-1525`),
registered in the `FeatureBuilderRegistry` under `"intrusion"`.

Build sequence inside `_build_intrusion`:

```
IntrusionFrameBuilder.set_intrusion_frame_parameters(intrusion_data, params)
IntrusionFrameBuilder.create_constraints_for_c0()      # synthesise c0=0 points/gradients
IntrusionFrameBuilder.set_intrusion_frame_data(frame_data)
IntrusionFrameBuilder.build(...)                        # -> IntrusionFrame
IntrusionBuilder(frame, lateral_extent_model=..., vertical_extent_model=...)
IntrusionBuilder.set_data_for_extent_calculation(intrusion_data)
IntrusionBuilder.update_build_arguments({"geometric_scaling_parameters": ...})
# lazy: IntrusionBuilder.build() runs on first evaluate_value() via up_to_date()
```

---

## 2. User guide — what's required to build an intrusion

As of 2026-08-05, `tests/unit/modelling/intrusions/test_intrusions.py`
covers: a tabular intrusion built from a single roof-or-floor contact plus
one stratigraphic conformable feature as the inflation-gradient proxy
(`load_tabular_intrusion()`, the original path), and — new — one
`marginal_faults` example built from scratch
(`test_intrusion_marginal_faults`/`_build_marginal_fault_model`) showing a
sill offset against a single bounding fault. `intrusion_steps` is *not*
demonstrated as working, because it currently can't be — see §3, finding 1;
`test_intrusion_steps_broken_with_current_stratigraphic_column` pins the
failure instead. The shortest-path network method remains undemonstrated
and is now dead code (§3, finding 5) rather than just untested.

### 2.1 Two blocks of input data

Both blocks live in `model.data` (one `pandas.DataFrame`, distinguished by
`feature_name`), matching `intrusion_name` and `intrusion_frame_name` passed
to `create_and_add_intrusion`.

**A. Intrusion frame data** (`feature_name == intrusion_frame_name`) — this
is an ordinary structural-frame dataset, the same schema used for faults and
fold frames: rows tagged `coord` 0/1/2, each with either a `val` (point lies
on that isovalue) or a gradient (`nx, ny, nz`, vector the coordinate's
gradient should be parallel to). You do **not** need to supply `coord=0`
data yourselves — `IntrusionFrameBuilder.set_intrusion_frame_data` (`intrusion_frame_builder.py:935-978`)
synthesises coord-0 point and gradient constraints from the intrusion
contact data below and appends them for you. You only need to supply
`coord=1` and `coord=2` constraints (an origin point + two orthogonal
direction vectors is enough — see the example).

**B. Intrusion contact data** (`feature_name == intrusion_name`) — points
on the intrusion margin:

| Column | Required? | Meaning |
|---|---|---|
| `X, Y, Z` | yes | point location |
| `intrusion_contact_type` | yes | `"roof"`/`"top"` or `"floor"`/`"base"` — which margin this point is on |
| `intrusion_side` | yes for lateral extent | boolean; `True` marks points used to constrain the lateral (lengthwise) margins. Split into "min side"/"max side" by the sign of `c2` at that point. |
| `intrusion_anisotropy` | only for steps/marginal-faults/shortest-path | name of the host-rock series or fault feature the point is associated with |

None of these columns are validated: a missing `intrusion_contact_type`
column raises a bare pandas `KeyError` deep inside `set_intrusion_frame_c0_data`
with no indication of what's wrong (see §3, finding 3).

### 2.2 `intrusion_frame_parameters` dict

Passed to `create_and_add_intrusion(..., intrusion_frame_parameters={...})`.

| Key | Default | Notes |
|---|---|---|
| `contact` | `"floor"` | which `intrusion_contact_type` value is the *reference* contact used to build `c0=0` |
| `contact_anisotropies` | `None` | **de facto required** — a list of series-type features; `create_constraints_for_c0` unconditionally indexes `[0]` (`intrusion_frame_builder.py:890`) to get an inflation-gradient proxy. Omitting it, or passing `[]`, crashes with `TypeError`/`IndexError`, not a helpful error. |
| `g_w` | `None` → 100 pts | weight/count controlling how many synthetic gradient constraints get added for `c0` |
| `intrusion_steps` | `None` | **currently broken, don't use** — see §3, finding 1. Was meant to be a dict of step definitions, each needing `structure` (fault), `unit_from`/`unit_to` (stratigraphic unit names), `series_from`/`series_to` (series features) |
| `marginal_faults` | `None` | dict of fault definitions bounding the intrusion; each needs `structure` (a `FaultSegment`/`FaultSegment`-like object, indexed as `structure[0]` for its coordinate-0 feature — same convention as `intrusion_steps`' `structure`), `block` (`"hanging wall"`/`"foot wall"`), `series` (a series feature). Working example: `test_intrusion_marginal_faults` in `test_intrusions.py`. |
| `delta_c`, `delta_f` | `[1]*n` | multiplies the std-dev band used to detect points near a contact/fault when synthesising `c0` constraints |

### 2.3 Conceptual model functions

`intrusion_lateral_extent_model` and `intrusion_vertical_extent_model` are
plain functions from `geom_conceptual_models.py` (or custom ones matching
the same **dual-arity calling convention** — see §3, finding 8):
`ellipse_function` (lateral) and `constant_function` or `obliquecone_function`
(vertical) are the built-ins.

### 2.4 Worked example (from the passing test)

```python
from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_tabular_intrusion
from LoopStructural.modelling.intrusions import ellipse_function, constant_function

data, bounding_box = load_tabular_intrusion()

model = GeologicalModel(bounding_box[0, :], bounding_box[1, :])
# NOT `model.data = data` -- this dataset has no nx/ny/nz/tx/ty/tz columns,
# and `prepare_data` is what normalises those in; skipping it builds a model
# that crashes with a bare KeyError the moment anything calls
# `evaluate_value()` on it. See §3, finding 1b.
model.data = model.prepare_data(data)

stratigraphy = model.create_and_add_foliation("stratigraphy")

intrusion = model.create_and_add_intrusion(
    "tabular_intrusion",
    "tabular_intrusion_frame",
    intrusion_frame_parameters={
        "contact": "roof",
        "contact_anisotropies": [stratigraphy],
    },
    intrusion_lateral_extent_model=ellipse_function,
    intrusion_vertical_extent_model=constant_function,
)
```

`intrusion.evaluate_value(xyz)` then returns a signed pseudo-distance field
to the intrusion contact (negated internally so it can double as an
unconformity — `intrusion_feature.py:329-331`). Note `evaluate_gradient` is
`NotImplementedError` (`intrusion_feature.py:233-235`) — anything in the
model stack that needs gradients of a feature (fold axes, some fault
throw calculations, vector-field visualisation) cannot be pointed at an
`IntrusionFeature`.

---

## 3. Review findings

Ordered roughly most → least impactful. **[FIXED]** markers were added
2026-08-05 after implementing Phases A-C (see §4); the rest of the text is
left as originally written for the record.

1. **`intrusion_steps` is not just untested — it's provably broken against
   the current `StratigraphicColumn`.** Discovered while writing the Phase C
   regression test (below), this is worse than originally written here.
   `IntrusionFrameBuilder.set_intrusion_steps_parameters` reads
   `self.model.stratigraphic_column[series_from_name.name][unit_from_name].get("id"/"min"/"max")`
   — i.e. it expects `model.stratigraphic_column` to still be the old
   nested `{group_name: {unit_name: {...}}}` dict. It no longer is:
   `GeologicalModel.stratigraphic_column` is now a `StratigraphicColumn`
   object (`modelling/core/stratigraphic_column.py`) whose `__getitem__`
   looks up a single element by `uuid`, not by group/series name, and
   returns a `StratigraphicColumnElement` (no `.get()`). Worse,
   `GeologicalModel.set_stratigraphic_column` — the old dict-based setter —
   now unconditionally `raise DeprecationWarning(...)` after logging, so
   there is no longer any way to put the model into the shape
   `intrusion_steps` expects at all. This is a real regression from an
   unrelated stratigraphic-column refactor that was never propagated to the
   intrusions module. **Confirmed by** `test_intrusion_steps_broken_with_current_stratigraphic_column`
   in `test_intrusions.py`, which pins the exact current failure
   (`KeyError: 'No element found with uuid: ...'`). Fixing this for real
   needs a design decision (what should `unit_from`/`series_from` resolve
   against now?) — that's Stage 6 work, not something patched in Phase A-C.

1b. **[FIXED, as a documentation/test fix]** The packaged
   `datasets/data/tabular_intrusion.csv` — the *only* dataset this module
   has — has no `nx`/`ny`/`nz`/`tx`/`ty`/`tz` columns, only `gx`/`gy`/`gz`.
   `GeologicalModel.data = df` does not normalise columns (unlike
   `GeologicalModel.prepare_data(df)`, which adds any of `all_heading()`
   missing as `NaN`). `GeologicalFeatureBuilder.add_data_to_interpolator`
   unconditionally indexes `normal_vec_names()`/`tangent_vec_names()` on a
   builder's own data, so **any full end-to-end build+`evaluate_value()`
   call through the real `create_and_add_intrusion` API** (as opposed to
   the two pre-existing tests, which manually construct
   `IntrusionFrameBuilder`/`IntrusionBuilder` and pass already-`prepare_data`-d
   slices, or never call `evaluate_value` at all) crashes with a bare
   pandas `KeyError` on the missing columns. This is a general
   `GeologicalModel`/`GeologicalFeatureBuilder` sharp edge, not specific to
   intrusions, but intrusions is where it surfaced, precisely because (per
   finding 1) nothing had ever driven this dataset through the real,
   public, end-to-end path before. The fix used throughout the new tests is
   simply `model.data = model.prepare_data(df)` instead of `model.data = df`
   — worth calling out prominently in any future user guide/docs for this
   dataset, since the obvious `model.data = data` (as the original
   `test_intrusion_frame_builder`/`test_intrusion_builder` imply, and as
   this document's own §2.4 example originally showed) silently sets up a
   model that will crash the moment anyone actually evaluates it.

2. **[FIXED] `geometric_scaling_parameters` is fully non-functional.**
   `IntrusionBuilder.create_geometry_using_geometric_scaling`
   (`intrusion_builder.py:113-154`) always raised `NotImplementedError`,
   even on the branch where `thickness` *is* provided — the
   `raise NotImplementedError("Not implemented")` at line 146 was
   unconditional, after the early-return-shaped `if estimated_thickness is
   None` check. The only way to reach this function is if the *other*
   contact has zero data points **and** `geometric_scaling_parameters` is
   non-empty, so passing this parameter to `create_and_add_intrusion` for a
   single-contact intrusion always crashed. `geometric_scaling_functions.py`'s
   `contact_pts_using_geometric_scaling` (which this should call) is fully
   implemented and unit-testable, but dead — its only caller was commented
   out. **Fix applied:** the function now raises immediately, with a clear
   message naming what's missing, instead of after several lines that look
   like partial progress; the dead commented-out code and the now-pointless
   wildcard import of `geometric_scaling_functions` were removed. Pinned by
   `test_intrusion_geometric_scaling_not_implemented`.

3. **[FIXED] A typo silently drops a user-supplied build weight.**
   `GeologicalModel._build_intrusion` called
   `intrusion_frame_builder.build(nelements=..., w2=weights[0],
   w1=weights[1], gxygz=weights[2])`. `StructuralFrameBuilder.build`
   only recognises the deprecated alias `gyxgz` for `w3` — `gxygz` (letters
   transposed) fell into `**kwargs` and was silently ignored. Any caller
   passing `gyxgz=...` to `create_and_add_intrusion` to weight the frame's
   third-coordinate orthogonality constraint had no effect; `w3` was always
   `1.0`. This was the only call site in the codebase using this kwarg name
   for a frame build — faults/foliations don't have this bug, which is
   consistent with the intrusion path being far less exercised. **Fix
   applied:** `gxygz` → `gyxgz`. Pinned by
   `test_intrusion_gyxgz_weight_reaches_frame_build`.

4. **[FIXED] No data-contract validation.** Missing `intrusion_contact_type`,
   `intrusion_side`, or `intrusion_anisotropy` columns failed with bare
   `KeyError`/`AttributeError` deep in helper methods, not at the
   `create_and_add_intrusion` boundary where the user could get a useful
   message. Likewise `contact_anisotropies` being `None`/empty crashed with
   no hint of what's missing. **Fix applied:**
   `GeologicalModel._validate_intrusion_inputs` now runs at the top of
   `_build_intrusion` and raises a clear `ValueError` naming the missing
   feature/column/parameter. Pinned by
   `test_intrusion_missing_data_raises_clear_error`,
   `test_intrusion_missing_contact_type_column_raises_clear_error`,
   `test_intrusion_missing_contact_anisotropies_raises_clear_error`.

5. **Dead code — mostly [FIXED]:**
   - **[FIXED, deleted]** `intrusion_support_functions.py` (shortest-path
     grid/graph helpers, ~390 lines) — not imported anywhere outside itself.
   - **[FIXED, deleted]** `IntrusionFeature.evaluate_value_test` — an entire
     parallel implementation of `evaluate_value`, unused.
   - **[FIXED, deleted]** `IntrusionFrameBuilder.update()` was an exact
     duplicate of the parent `StructuralFrameBuilder.update()` — redundant
     override.
   - **[FIXED, deleted]** the two discarded
     `self.marginal_faults[fault_i].get("emplacement_mechanism")` calls.
   - **[FIXED]** `IntrusionFeature.evaluate_value`'s dead `if/else` with
     identical bodies (the `marginal_faults is not None` branch) collapsed
     to one unconditional assignment.
   - **Left as-is (not removed):** `IntrusionFeature.add_assisting_faults`/
     `self.assisting_faults` is still set but never populated by any caller
     in the codebase, so the "asymmetry weight" branch in `evaluate_value`
     stays unreachable in practice — removing a public method felt like a
     bigger call than Phase A's "no behaviour change" scope; flagging again
     for Phase D. `intrusion_builder.py`'s comment
     `# intrusion_frame_builder.post_intrusion_faults = faults  # LG unused?`
     (your own prior TODO) is also still there, unaddressed.

6. **Reproducibility, corrected.** Originally written up here as "KMeans
   hardcodes `random_state=0` instead of going through the shared,
   seedable `rng`" — that framing turned out to be backwards.
   `loop_common.utils.rng` (`from ...utils import rng`, used e.g. for
   `rng.shuffle`) is `np.random.default_rng()` created **unseeded**, once,
   at import time — it's not actually reproducible across process runs
   either, and there's no `set_seed`-style entry point anywhere in the
   codebase. Routing `KMeans` through it would have made cluster labels
   (and therefore the synthesised `c0` constraints for steps/marginal
   faults) vary run-to-run instead of being the one deterministic piece.
   **Fix applied instead:** kept the fixed seed (it's the right call), but
   replaced the two duplicated `random_state=0` literals and the stale
   `# TODO create global loopstructural random state variable` comment
   (which described a solution — the shared `rng` — that doesn't actually
   solve this) with one named module constant, `_KMEANS_RANDOM_STATE`.

7. **Scikit-learn dependency, corrected.** Originally flagged here as "a
   hard dependency on scikit-learn just for the steps/marginal-faults code
   path." On checking `pyproject.toml`, `scikit-learn` is already a
   top-level `dependencies` entry for the whole package (also used in
   `LoopStructural/utils/helper.py` and `_transformation.py`), so this
   isn't an extra install cost specific to intrusions — retracted as a
   finding. The module-scope `try/except ImportError: ... raise` in
   `intrusion_frame_builder.py` is still slightly unusual style (failing at
   import time rather than at first use), but not a real problem worth
   spending Phase A/B effort on.

8. **Conceptual-model functions overload call arity to mean two different
   things.** `ellipse_function`/`constant_function`/`obliquecone_function`
   are called once with **no arguments** (during
   `set_conceptual_models_parameters`, to fetch bounds:
   `intrusion_builder.py:277-278`) and later called **with data** to
   evaluate the model at points. This works only because every branch
   checks `if <data>.empty` and returns a different tuple shape. It's an
   implicit, undocumented protocol — a custom conceptual model has to
   reverse-engineer this from the three built-ins rather than from any
   docstring or interface. A named two-method interface (e.g. `.bounds()`
   / `.evaluate(data, ...)`, whether a `Protocol` or small ABC) would make
   custom conceptual models actually writable by a user.

9. **`IntrusionFrameBuilder` is a ~20-attribute god object** mixing frame
   geometry, steps, marginal faults, deprecated shortest-path indicator
   functions (`IFf`, `IFc`, single/double-letter names), and sill-splitting
   in one class. `create_constraints_for_c0` alone is ~200 lines covering
   four largely-independent concerns (steps / marginal faults / gradient
   constraints / sill-splitting) in one method with deep nesting.

10. **Copy-pasted min/max-side blocks.** `IntrusionBuilder.set_data_for_lateral_thresholds`
    (`intrusion_builder.py:344-563`, ~220 lines) repeats the same
    conceptual/residual computation twice — once for the `L<0` side, once
    for `L>0` — with only sign/index flips between them. Same shape of
    duplication in `interpolate_lateral_thresholds`
    (`intrusion_feature.py:60-137`). Both are natural candidates for a
    single side-parameterised helper.

11. **Fragile cross-feature reach-through.** Sill-splitting
    (`create_constraints_for_c0`, `intrusion_frame_builder.py:807-815`) does
    `self.model.__getitem__(splits_from_sill_name).intrusion_frame.builder.intrusion_steps`
    — reaching through a private-ish attribute chain on another feature
    entirely, with no cycle detection if two sills reference each other.

12. **[FIXED] Wildcard import.** `intrusion_builder.py:6` had
    `from .geometric_scaling_functions import *`. Removed as part of the
    finding-2 fix, once `create_geometry_using_geometric_scaling` no longer
    called anything from that module.

13. **Naming/typo debt** (low severity, but adds up for anyone trying to
    read this code for the first time): `strigraphic` (`intrusion_frame_builder.py:306,474`),
    `framce` (`:307,475`), `aftecting` (`:904`), `yo may increase` (`:822`),
    `indentify` (multiple), single-letter/cryptic names (`If`, `Ic`, `IFf`,
    `IFc`, `s`, `d`/`e`/`f` as boolean masks in `evaluate_value`).

---

## 4. Simplification & hardening plan

Proposed as the concrete content for **ROADMAP.md Stage 6**. Phased so each
step is independently shippable and testable, following the same
small-stages philosophy as the rest of the v2 roadmap.

**Phase A — stop the bleeding (no behaviour change to the working path) —
DONE 2026-08-05:**
- Fixed the `gxygz`/`gyxgz` typo (finding 3); regression test
  `test_intrusion_gyxgz_weight_reaches_frame_build`.
- Deleted confirmed-dead code: `intrusion_support_functions.py`,
  `evaluate_value_test`, the unused `emplacement_mechanism` `.get()` calls,
  the redundant `IntrusionFrameBuilder.update()` override, the now-pointless
  wildcard import.
- `geometric_scaling_parameters` now fails immediately with a clear message
  instead of after several lines that looked like partial progress
  (finding 2); regression test `test_intrusion_geometric_scaling_not_implemented`.
  `add_assisting_faults` was deliberately **not** touched — removing a
  public method felt out of scope for "no behaviour change"; still flagged
  for Phase D.
- Replaced the two hardcoded `KMeans(random_state=0)` with one named
  constant and corrected the misleading comment about routing through the
  shared `rng` (finding 6 — that would have made results *less*
  reproducible, not more; see the corrected finding for why).

**Phase B — make the failure modes legible — DONE 2026-08-05:**
- `GeologicalModel._validate_intrusion_inputs` now runs at the top of
  `_build_intrusion`: checks both `feature_name`s have data, required
  columns (`intrusion_contact_type`, `intrusion_side`) are present, and
  `contact_anisotropies` is non-empty — clear `ValueError`s instead of a
  `KeyError` several calls deep. Regression tests
  `test_intrusion_missing_data_raises_clear_error`,
  `test_intrusion_missing_contact_type_column_raises_clear_error`,
  `test_intrusion_missing_contact_anisotropies_raises_clear_error`.
- The "make scikit-learn optional" item was **retracted**: scikit-learn is
  already a top-level dependency of the whole package (finding 7,
  corrected), so there was nothing to fix here.

**Phase C — cover what's untested before touching it further — DONE
2026-08-05, with a significant finding:**
- Added `test_intrusion_marginal_faults` (`_build_marginal_fault_model`): a
  from-scratch synthetic sill offset by one bounding fault, built and
  evaluated end-to-end through the real `create_and_add_intrusion` public
  API — the first time this code path has ever been exercised.
- Attempting the equivalent for `intrusion_steps` surfaced finding 1:
  it doesn't work at all against the current `StratigraphicColumn`, not
  just "untested." `test_intrusion_steps_broken_with_current_stratigraphic_column`
  pins the current failure instead of demonstrating success. **This means
  Phase D's scope for `intrusion_steps` is bigger than originally framed
  here**: it's not "simplify working code," it's "decide what this feature
  should even look like against the new stratigraphic column API, or drop
  it." That decision belongs in the dedicated Stage 6 discussion.
- Also surfaced finding 1b (the packaged `tabular_intrusion.csv` dataset
  needs `model.prepare_data()`, not a bare `model.data = df`, to survive an
  end-to-end `evaluate_value()` call) — fixed in the tests and in this
  document's own worked example (§2.4), which had the same bug.

**Phase D — the actual simplification (larger, still open — needs the
Stage 6 discussion):**
- Split `IntrusionFrameBuilder` along its four concerns (frame geometry /
  steps / marginal faults / deprecated shortest-path) into
  composable pieces, so the common tabular case doesn't carry ~20 unused
  attributes.
- Replace the conceptual-model dual-arity calling convention (finding 8)
  with an explicit two-method interface; migrate the three built-ins.
- Collapse the min/max-side and lateral/vertical duplicated blocks
  (findings 10) into single side-parameterised helpers.
- Decide whether the shortest-path network method is kept (and given the
  same test/doc treatment as the rest) or removed outright — it's now fully
  dead code (finding 5/1), not half-present: `intrusion_network_type` is
  initialised to `None` and never set to `"shortest path"` anywhere reachable.
- Decide what `intrusion_steps` should even look like against the current
  `StratigraphicColumn` object (finding 1) — the old nested-dict lookup it
  was written against no longer exists, and there's no bridge back to it.
  This is bigger than a refactor: it needs a real design decision before
  any of the "simplify the god object" work below can safely touch it.

Phases A-C landed 2026-08-05 (see the "DONE" notes above) — they didn't
need to wait for the Stage-5 graph backend. Phase D is what still needs the
"dedicated discussion" ROADMAP.md defers to post-Stage-5, since a graph
backend may change how frame/feature dependencies (sill-splitting, step
faults) are expressed anyway, and the `intrusion_steps`/`StratigraphicColumn`
question needs a decision either way.
