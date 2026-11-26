# Implementation plan: loop_link, a Python link to geoviewer

**Feature**: `006-loop-link` | **Date**: 2026-11-06 (S5 added 2026-11-07, S4 added 2026-11-08) | **Spec**: [spec.md](spec.md)

<!-- A plan says HOW. Written in ASD-STE100 Simplified Technical English.
Scope of this plan: stories 1 and 5 (stages S1 to S3, built), and stories 2 to 4
(stage S5, planned 2026-11-07, see "Design of stage S5"), and stories 6 and 7 (stage S4,
planned 2026-11-08, see "Design of stage S4"). -->

## Open points

Each point has a default. A point that has no owner is for the maintainer.

- **OP1**: The `[NEEDS CLARIFICATION]` mark of story 7, scenario 7 (FR-021). Answered
  2026-11-08 (see "Answers"): all colour settings go as follow-up commands on the tool
  connection, in the same call. The `scene.upsert` data array can carry `colour_map`
  (`name`, `range`, `clip`, `scale`), but it has no field for the property that colours
  the object or for the plain colour, so the follow-up commands are needed anyway. One
  path is simpler than two. The user sees one call.
- **OP2**: The default `crs` of `hello`. `Hello.crs` is a required string in the
  viewer. The viewer may refuse an empty text. Default: `Viewer(crs="")`, and the
  library keeps the CRS that the viewer answers in `Viewer.scene_crs`. Task T-S1
  for `viewer.py` reads `crates/geoviewer/src` hello handling in the geoviewer
  repository first. If the viewer refuses `""`, the library raises
  `InvalidInputError` that says "give crs=..." for an empty scene.
- **OP3**: FR-010 says that pyvista input works for all seven types. A pyvista
  mesh cannot hold a collar, a survey and intervals. Default: `from_pyvista`
  raises `InvalidInputError` for a Drillhole, and the message names `Drillhole`
  and the pandas adapter. The coverage table shows this gap until the spec owner
  answers. Blocks S2 only for this one combination.
- **OP4**: FR-017 sets Python 3.10 for this package. The QGIS plugin runs on the
  Python of QGIS, which can be 3.9 on some systems. Owner: plugin repository.
  Default: keep 3.10. This blocks no task here.
- **OP5**: The root `pyproject.toml` has `requires-python = '>=3.9'` and the uv
  workspace takes all `packages/*`. A member with `>=3.10` can make `uv lock`
  refuse the workspace, or can raise the floor of the lock. Default: add the
  member, and run `uv lock`. If uv refuses, or if the lock drops 3.9 for the other
  members, add `packages/loop_link` to `[tool.uv.workspace] exclude` (as
  `packages/loop_interpolation_core` is) and make the CI job install it with
  `pip install -e packages/loop_link[tests]`.
- **OP6**: Stories 6 and 7 (filter and colour, FR-018 to FR-021, SC-007) are planned in
  stage S4 (see "Design of stage S4"). They use the tool connection of S5.
- **OP7**: Stories 2, 3 and 4 (select, read selection, fetch) are planned in stage S5
  (see "Design of stage S5"). They use geoviewer spec 067-selection-and-fetch
  (protocol 1.18), merged in the geoviewer repository. The text of this spec and of
  the earlier parts of this plan calls it "066"; read "067". The spec owner should
  change the number in `spec.md`. The names of the messages are final (they are in
  `docs/protocol.md` of the geoviewer repository).
- **OP8**: SC-004 needs the folder `fixtures/protocol/python` of the geoviewer
  repository. CI of this repository cannot see it. Default: the tests read the
  folder from the environment variable `LOOP_LINK_GEOVIEWER_FIXTURES`. When it
  is not set, they skip, with the reason in the skip text. A copy of the frames
  is not stored here (one source). The maintainer runs the tests with the
  variable before a release.
- **OP9**: The structured grid of the spec is the viewer `BlockModel` (a regular
  axis-aligned grid). Default: a non-regular pyvista `StructuredGrid` or a
  `RectilinearGrid` with uneven spacing gives `InvalidInputError`. The viewer has
  no type for them.

- **OP10**: A fetch of a Drillhole that the viewer made by an import can have log
  sets, sample sets and per-hole values that the `Drillhole` dataclass cannot hold
  (it has one set of intervals). Default: `fetch` keeps the first interval set and
  the survey values, and drops the rest with a `UserWarning` that names each dropped
  property. The other choice is an error, which would make most imported holes
  impossible to fetch. Owner: spec owner. Blocks no task.
- **OP11**: The tool fixtures of geoviewer are in `fixtures/protocol/tool`, next to
  `fixtures/protocol/python`. Default: the tests read them from
  `<LOOP_LINK_GEOVIEWER_FIXTURES>/../tool` and skip when the folder is missing (OP8
  rule). No new variable.
- **OP12**: Which object types have a selection table is a viewer fact that the
  protocol text does not list. Default for the fake viewer: `Points` has nodes only;
  `Curve`, `Surface`, `BlockModel` and `TetraMesh` have nodes and cells; `Drillhole`
  and `Group` have no table. Manual check M3 repeats this with the real viewer. If the
  real viewer differs, change the fake viewer (one task).
- **OP13**: The spec entity "Selection" is a list of object ids. Spec 067 gives one
  object for each selection (one table). Default: `Selection` has one `object_id`
  (or `None`). The spec owner should change the entity text.
- **OP14**: The plan of geoviewer spec 067 says that a part has an `offset` and that
  the parts follow the order of `arrays`. It does not say that the parts of one array
  arrive in order with no gap. Default: the assembler requires this and gives an error
  otherwise. M3 checks it with the real viewer. If the real viewer sends the parts in
  another order, relax the rule to "each byte is filled one time" (one task).

- **OP15**: ANSWERED 2026-11-09 (see "Answers"). FR-019 asks the call to wait for the end of a
  filter build. Geoviewer spec 072 (protocol 1.20) makes `filter.edit` a long command with
  progress and an end result. The repeat of `set_on` of the first default is not used.
- **OP16**: The viewer does not say if `add` for an array that has a filter is an error, or
  if a new filter is on after `add` and `set_range`. Default: the library sends `add`,
  goes on if it is refused, and sends `set_on` true at the end. M4 step 9 checks it.
- **OP17**: A value filter takes the `u32` ids of the value map. Story 6 scenario 2 says "the
  values to show" and does not say if the user gives names. Default: S4 takes ids
  (integers). The library keeps no value map of the objects it sent. A later stage can add
  `names=`. Owner: spec owner. Blocks no task.
- **OP18**: The type rules (a filter on a text property, a map on a category) are assumed to
  give `refused` with a text. If the real viewer accepts them without a change, scenario 4
  of story 6 and scenario 5 of story 7 fail against the real viewer. M4 steps 4 and 6 check
  it. The fix would then be a type check in the library, which needs the types of the
  properties (a fetch), or a viewer change.
- **OP19**: `send(..., colour=...)` waits 0.5 s (`flush`) before the colour commands. A
  large object can take longer to enter the scene, and the first colour command then gets
  an error for an unknown object. M4 step 8 checks it. If it fails, replace the fixed wait
  by a retry for a time that the user sets. Owner: maintainer. Blocks no task.
- **OP20**: ANSWERED 2026-11-09 (see "Answers"). The result of `filter.edit` of protocol 1.20
  has `passing` and `total`. The filter calls return them in a `FilterResult`.

## Answers

- 2026-11-06 (orchestrator, wave 10 review): pyvista adapter. The adapter skips the arrays that
  pyvista itself adds (`Normals`, `Texture Coordinates`) without an error, so that `pv.Sphere()`
  and `pv.Plane()` can be sent as they are. Any other 2-D array is still refused with the array
  name. The core object rule (a 2-D property is refused) does not change.
- 2026-11-06 (orchestrator, end review, verified against the geoviewer Rust codec): a Drillhole
  `interval_properties` array goes out with the association `cell` (one value for each interval,
  protocol 1.4 and later), not `interval`. The association `interval` needs a `log_sets` entry in the
  header, and `loop_link` sends none. Where this plan says `interval`, read `cell`.
- 2026-11-06 (orchestrator, end review): `push_objects` takes items `(id, obj)` or
  `(id, obj, options)`; `options` is a dict with the keys `name` and `parent`. The removal step of
  `replace_prefix` never removes an object of type `Group`.
- 2026-11-06 (orchestrator, wave 10 review): `Viewer.inventory()` shows the known scene: the list
  from the hello frame, updated by `send` (add or update the entry), `remove` (drop it) and
  `clear` (drop all). `push_objects` with `replace_prefix` keeps the ids of the objects that were
  sent or skipped in this call, and does no removal when the push was cancelled.
- 2026-11-07 (orchestrator, manual check M2 with the real window viewer; this replaces the
  2026-11-06 entry about `cell`): a Drillhole `interval_properties` array goes out with the
  association `interval`, the field `log_set` = `"intervals"` in its description, and the header
  field `log_sets` = `[{"name": "intervals", "kind": "interval"}]`. The viewer offers
  drillhole interval values for colour only from such arrays; with `cell` it accepted the object
  but showed no property. A Drillhole with no interval properties has no `log_sets`. A tetra
  mesh has no Data button in the viewer; its properties are used by the isosurface and filter.
- 2026-11-07 (maintainer): the geoviewer spec for stories 2 to 4 is 067-selection-and-fetch
  (the text calls it 066).
- 2026-11-07 (orchestrator, wave 11): `Viewer.select` always returns both keys, `{"cell": n, "vertex": m}`.
  A kind that the reply does not list has the count 0. This replaces "only the kinds that the reply has" in
  the text of "Viewer (`viewer.py`)" of stage S5. The tests of T048 check this shape.
- 2026-11-08 (planner, S4): OP1 is answered. `send(..., colour=ColourSettings(...))` sends the
  object, waits with `flush()`, and then sends the colour commands (`colour.set_active`,
  `colour.set_map`, `colour.set_range`, `colour.set_scale`, `object.set_colour`) on the tool
  connection. The `colour_map` field of the data array is not used. The commands use the
  protocol version 1.11 of the fixtures.
- 2026-11-09 (spec owner: geoviewer spec 072-filter-result is merged, protocol 1.20): OP15 and OP20
  are answered. `filter.edit` is a long command: the viewer sends `command.progress` while it builds and
  one `command.result` when the build ends. The result is `{"status": "done", "passing": n, "total": m}`,
  or `{"status": "done"}` (no filter, or the `cancel` edit), or `{"status": "done", "state": "unbuilt"}`
  (no build for the new filters: the object is hidden, the draw data is not ready, the edit repeats a
  cancelled or failed edit, the object has a new revision, or another build runs), or
  `{"status": "done", "state": "building"}` (the viewer stops to wait after 60 s; the build goes on).
  `command.cancel` stops a build: the request gets `request.error` `cancelled` and the filter stays as it was.
  The `commands` entry of `filter.edit` has `"long": true` in 1.20. Decision for `loop_link` (stage S4d):
  (1) The tool protocol version of the client stays 1.18; the viewer sends progress to older clients.
  (2) `loop_link.commands.FilterResult(passing, total, state)` is new (frozen dataclass, each field `None` when
  the viewer gives none). `filter_range`, `filter_values`, `filter_enable` and `filter_remove` return it.
  (3) `wait=True` (default) waits for the result of each edit, so the call returns when the build ends. The
  arguments `cancel` (an event that sends `command.cancel`, error `CancelledError`), `timeout` (the time to wait
  for the next frame, as in `ToolConnection.request`) and the new `progress` (a callable that gets each
  fraction) are used. A `state` of `building` is returned, not raised.
  (4) `wait=False` sends the last `set_on` with the new `ToolConnection.post` (no read) and returns an empty
  `FilterResult()`. A refusal of that edit is not seen.
  (5) A viewer of 1.19 or less answers at once with no counts. The call then returns an empty `FilterResult()`.
  Manual check M4 step 7 needs a viewer of 1.20. The fake viewer gets a long `filter.edit` for
  `protocol_minor >= 20` (tasks T088 and T089).

## Summary

`packages/loop_link` is a new Python package. It has its own codec for protocol
1.14 (numpy and the standard library only), a stdlib WebSocket client (ported from
the plugin), discovery of the viewer and the token, the seven object types as
dataclasses, one `Viewer` class that holds the one data connection (with a
keepalive thread), optional adapters (pyvista, LoopStructural, pandas), a plugin
push API with no Qt, and a fake viewer in `loop_link.testing`. The package
registers in `ROADMAP.md`, the release-please files, the CI matrix and the Ruff
config.

The work has three built stages: S1 (send, story 1), S2 (adapters, FR-010), S3
(plugin API, story 5). Stage S5 (planned 2026-11-07) adds the tool connection, the
selection calls and the fetch (stories 2 to 4). Stage S4 (filter and colour) is for
a later run.

## Technical context

- **Packages changed**: new `packages/loop_link`. Root files: `pyproject.toml`,
  `uv.lock`, `ROADMAP.md`, `.release-please-config.json`,
  `.release-please-manifest.json`, `.github/workflows/packages.yml`. No change
  in `LoopStructural/`, `loop_common`, or `packages/loop_interpolation_core`.
- **Languages and versions**: Python `>=3.10` for the package (FR-017). Root
  Ruff and Black target py39. So the package has its own `ruff.toml` with
  `target-version = "py310"` (see "Project structure"). Every module starts with
  `from __future__ import annotations` (Ruff rule `FA`). The code uses no
  `match` statement and no parenthesised `with`, because `black .` from the root
  parses with target py39.
- **New dependencies**: required: `numpy`. Extras: `pyvista` (`pyvista`),
  `loopstructural` (`loop-common`), `pandas` (`pandas`), `tests` (`pytest`).
  There is no WebSocket package. The plugin client uses the standard library
  because QGIS does not ship a WebSocket package for each platform. The spec says
  "a WebSocket client", and a stdlib client meets this with no dependency
  (principle VII of the geoviewer constitution has the same intent).
- **Public API and import changes**: new import path `loop_link`. No existing
  import changes. The public names are in "Functions and public API".
- **Saved-data format changes**: none. The package writes no file. It reads
  `link.json` and `instances/<pid>.json` of the viewer.
- **Coordinate contract**: the library sends the coordinates that the caller
  gives, as world coordinates (spec 005: public data is in world coordinates). It
  subtracts the origin of the scene (from the `hello` answer of the viewer) in
  `f64`, then writes `f32`, or `f64` when a value is more than 4000 from the
  origin. No local frame or geological frame crosses the boundary. The LoopStructural
  adapter reads the public attributes `origin`, `step_vector`, `nsteps`,
  `vertices` (world coordinates) and applies no other transform.
- **Performance target**: SC-001 (1,000,000 cells, one float property, under 5 s
  to the inventory) and SC-007 (S4). Measure SC-001 with the fake viewer and
  with the real viewer by hand (manual check M1).
- **Platforms and backends**: Windows, Linux, macOS; Python 3.10 to 3.12. The CI
  matrix variable `PYTHON_VERSIONS` can have 3.9. The package job skips 3.9. The
  optional extras may be missing. A missing extra gives `MissingExtraError` when
  the adapter runs, and never at `import loop_link`.

## Constitution check

| Project rule | Result | How the plan agrees |
|---|---|---|
| Work stays within the existing library and package responsibilities | Pass | The new package is a client of the geoviewer protocol only. It does not change `LoopStructural/` or another package. The extra `loopstructural` reads the public classes of `loop_common.geometry`. |
| Development commands come from .dev-tools/project.json and the constitution | Pass | The commands are `pytest`, `ruff check LoopStructural packages` and `black .`. The package tests run with `uv run --package loop-link pytest packages/loop_link/tests` (the form of `packages.yml`). |
| Both pytest and ruff check LoopStructural packages must pass before merge | Pass | Each stage lists both commands. The new package is under `packages`, so Ruff checks it. |
| Python changes use numpy-style docstrings and type hints; Ruff targets py39 | Pass, with one deliberate difference | New files get no `D`/`ANN` exemption in `src`. Tests get the same exemption as the other package tests. The package target is py310 by `packages/loop_link/ruff.toml` (FR-017 asks the plan to set it). The repository default stays py39. |
| Black uses line length 100 and preserves string quoting | Pass | The package uses the root Black config. The code avoids syntax that py39 cannot parse. |
| ROADMAP.md owns architecture stages and decisions; update it when either changes | Pass | S1 adds the package to `ROADMAP.md` (FR-016, SC-006). |
| Stages can be tested and reversed independently | Pass | See "Stages and rollback". The new package only adds files. Rollback of each stage is a revert of that stage. |
| The public GeologicalModel API stays compatible; a break needs a ROADMAP.md decision and a major-version note | Not applicable | No change in `GeologicalModel` or `LoopStructural/`. |
| The QGIS plugin stays in its separate repository | Pass | The plugin change and the removal of the old client stay in the plugin repository (spec, "Out of scope"). S3 only gives the API. |
| Feature documents stay under specs/NNN-name/ | Pass | This file is `specs/006-loop-link/plan.md`. |
| Releases use release-please with per-package paths | Pass | S1 adds `packages/loop_link` to the config and manifest, and adds it to `exclude-paths` of the root entry. |
| A constitution amendment requires a pull request and maintainer acceptance | Not applicable | The plan does not change the constitution. |

## Design

### Types and data

All modules are in `packages/loop_link/src/loop_link/`. Names are checked free
(`git grep -n -w` on 2026-11-06 found no use of `loop_link`, `LinkError`,
`ViewerConnection` or `FakeViewer` outside `specs/`).

**Errors** (`errors.py`). Each type has its own class and a text that says what
to do (FR-008).

```python
class LoopLinkError(Exception): ...
class NoViewerError(LoopLinkError): ...         # no instance file, no listener
class TokenError(LoopLinkError): ...            # no token, or the viewer refused it
class BusyError(LoopLinkError): ...             # another program holds the data connection
class MissingFeatureError(LoopLinkError):       # viewer too old
    feature: str
class InvalidInputError(LoopLinkError): ...
class ObjectTooLargeError(LoopLinkError):
    object_id: str; size: int; limit: int
class ViewerClosedError(LoopLinkError): ...     # the viewer closed in the middle of a call
class ViewerRefusedError(LoopLinkError):
    refusals: list[Refusal]
class MissingExtraError(LoopLinkError):         # an optional extra is not installed
    extra: str
@dataclass(frozen=True)
class Refusal: object_id: str; code: str; message: str   # code as in protocol "error"
```

**Codec** (`protocol.py`, no Qt, no pyvista):

```python
PROTOCOL_VERSION = {"major": 1, "minor": 14}
MAX_FRAME_BYTES = 256 * 1024 * 1024
MAX_HEADER_BYTES = 1024 * 1024
MAX_OBJECT_BYTES = 2 * 1024 ** 3
F32_LIMIT = 4000.0
@dataclass(frozen=True)
class Frame:
    header: dict[str, Any]
    buffers: dict[str, np.ndarray]     # by name, shape as in the header, little-endian
def encode_frame(header: dict, buffers: Sequence[tuple[str, np.ndarray]] = ()) -> bytes
def decode_frame(data: bytes) -> Frame
def decode_header(data: bytes) -> dict
class ProtocolError(LoopLinkError): ...
```

`encode_frame` adds `version` and `buffers`, writes `json.dumps(...,
separators=(",", ":"), allow_nan=False)`, and refuses a header over
`MAX_HEADER_BYTES` or a frame over `MAX_FRAME_BYTES` (`ObjectTooLargeError` is
raised by the caller, which knows the id). The dtype names are `f32`, `f64`,
`u32`, `i32`, `u8`. `decode_frame` checks that `4 + header_len + sum(buffer
sizes) == len(data)`, and gives `ProtocolError` otherwise.

**Objects** (`objects.py`). All are frozen dataclasses. The arrays are numpy
arrays after `__post_init__` (`np.asarray`, no copy when possible). A shape error
gives `InvalidInputError` with the field name. `properties` has one value per
vertex (or per node for a `BlockModel`). `cell_properties` has one value per cell.
A value is a 1-D numpy array, or `Categorical`.

```python
ArrayLike = Union[np.ndarray, Sequence[float]]
@dataclass(frozen=True)
class Categorical:                       # protocol data type "referenced"
    codes: np.ndarray                    # (n,) integers >= 0, below 2**32
    names: Sequence[str]                 # name of value i is names[i]
    colours: Optional[Sequence[Sequence[int]]] = None   # RGBA per name; default DEFAULT_CATEGORY_COLOURS cycle
DEFAULT_CATEGORY_COLOURS: tuple[tuple[int, int, int, int], ...]   # 10 colours of the tab10 set

@dataclass(frozen=True)
class Points:
    vertices: ArrayLike                  # (n, 3) world coordinates
    properties: Mapping[str, Property] = field(default_factory=dict)
@dataclass(frozen=True)
class Vectors:                           # a Points object with a property group
    vertices: ArrayLike                  # (n, 3)
    orientation: Optional[tuple[ArrayLike, ArrayLike]] = None   # (dip_direction, dip) degrees, each (n,)
    lineation: Optional[tuple[ArrayLike, ArrayLike]] = None     # (trend, plunge) degrees, each (n,)
    properties: Mapping[str, Property] = field(default_factory=dict)
    @classmethod
    def from_directions(cls, vertices, directions, *, kind: str, properties=None) -> "Vectors"
        # kind: "normal" (plane normal -> orientation) or "lineation" (line direction)
@dataclass(frozen=True)
class Lines:
    vertices: ArrayLike                  # (n, 3)
    segments: ArrayLike                  # (m, 2) vertex indices
    properties: ...; cell_properties: ...   # cell = segment
@dataclass(frozen=True)
class Surface:
    vertices: ArrayLike; triangles: ArrayLike     # (n, 3), (m, 3)
    properties: ...; cell_properties: ...
@dataclass(frozen=True)
class BlockModel:                        # the "structured grid"
    origin: ArrayLike                    # (3,) world coordinates of the lowest corner
    cell_sizes: ArrayLike                # (3,) > 0
    cells: ArrayLike                     # (3,) integers >= 1: nx, ny, nz
    properties: ...                      # nodes, (nx+1)(ny+1)(nz+1), x fastest
    cell_properties: ...                 # cells, nx*ny*nz, x fastest
@dataclass(frozen=True)
class TetraMesh:
    vertices: ArrayLike; tetrahedra: ArrayLike    # (n, 3), (m, 4)
    properties: ...; cell_properties: ...
@dataclass(frozen=True)
class Drillhole:
    collar: ArrayLike                    # (3,) world coordinates
    survey: ArrayLike                    # (k, 3): depth, azimuth degrees, dip degrees (positive down)
    intervals: ArrayLike                 # (j, 2): from, to
    end_of_hole: Optional[float] = None
    properties: ...                      # one value per survey station (association "vertex")
    interval_properties: ...             # one value per interval (association "interval")
SceneObject = Union[Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole]
```

Validity rules (all checked in `__post_init__`, so a bad object never reaches the
codec): the vertex array has finite values; index arrays are integers with
`0 <= i < n`; an array property has the length of its association; a 2-D property
is refused; `cells` and `cell_sizes` are positive; a `Points` or `Vectors` object
has no `cell_properties`.

**Property to protocol data type** (`encode.py`, function `_data_array`):

| Python value | `data_type` | Buffer |
|---|---|---|
| `float32` | `float` | `f32` |
| other float dtype | `float` | `f64` |
| `bool` | `boolean` | `u8` |
| integer dtype, all values in the `i32` range | `integer` | `i32` |
| integer dtype with a value outside `i32` | refused | `InvalidInputError` |
| numpy string dtype (`U`, `S`) | `text` | `u8` bytes plus `u32` offsets |
| `Categorical` | `referenced` | `u32`, with `value_map` (`id`, `name`, `colour`) |
| `object` dtype, other | refused | `InvalidInputError` that names the accepted dtypes |

The data array `id` and `name` are the property name. The buffer name is
`data/<name>`, and the offsets buffer name is `data/<name>/offsets`. A name that
is not unique after this rule gives `InvalidInputError`.

**Vector conversion** (`orientation.py`). The directions use x east, y north, z up.
The result is in degrees.

```python
def normals_to_orientation(normals: ArrayLike) -> tuple[np.ndarray, np.ndarray]
    # n = normalise(normal); if n_z < 0: n = -n
    # dip = degrees(arccos(clip(n_z, -1, 1)));  dip_direction = degrees(arctan2(n_x, n_y)) % 360
    # n_x = n_y = 0: dip_direction = 0.0
def directions_to_lineation(directions: ArrayLike) -> tuple[np.ndarray, np.ndarray]
    # d = normalise(direction); if d_z > 0: d = -d     (a lineation is an axis)
    # plunge = degrees(arcsin(clip(-d_z, -1, 1)));  trend = degrees(arctan2(d_x, d_y)) % 360
    # d_x = d_y = 0: trend = 0.0
def orientation_to_normals(dip_direction, dip) -> np.ndarray      # helper for S5 (fetch); inverse of the first
def lineation_to_directions(trend, plunge) -> np.ndarray          # inverse of the second
```

A vector of length 0 or with a NaN gives `InvalidInputError`. The upward normal
of a plane that dips to azimuth `a` with dip `d` is
`(sin d sin a, sin d cos a, cos d)`. The unit tests use this as the independent
expected result. `Vectors.from_directions(kind="normal")` uses the first function.
`kind="lineation"` uses the second. Another `kind` gives `InvalidInputError`.

**Property groups of `Vectors`**: group id and name `orientation` (role
`orientation`, members `dip_direction` and `dip`) or `lineation` (role `lineation`,
members `trend` and `plunge`). The data arrays are `float` `f64`, association
`vertex`, ids `orientation/dip_direction`, `orientation/dip`,
`lineation/trend`, `lineation/plunge`. A `Vectors` object with both groups has both.

**Connection state** (`viewer.py`):

```python
@dataclass(frozen=True)
class InventoryEntry: id: str; object_type: str; revision: int
```

### Functions and public API

`loop_link/__init__.py` exports: `Viewer`, `PushResult`, `push_objects`,
`Points`, `Vectors`, `Lines`, `Surface`, `BlockModel`, `TetraMesh`, `Drillhole`,
`Categorical`, `InventoryEntry`, `Refusal`, all error classes, and `__version__`.
`import loop_link` imports numpy and the standard library only (SC-005). A test
checks that `sys.modules` has no key that starts with `PyQt`, `PySide`, `qgis`,
`pyvista`, `pandas`, `loop_common` or `LoopStructural` after the import.

```python
# viewer.py
class Viewer:
    def __init__(self, *, port: Optional[int] = None, token: Optional[str] = None,
                 host: str = "127.0.0.1", crs: str = "", origin: Sequence[float] = (0.0, 0.0, 0.0),
                 settings_dir: Optional[Path] = None, timeout: float = 5.0,
                 heartbeat_interval: float = 5.0) -> None
    def connect(self) -> "Viewer"
    def send(self, obj: Any, *, id: str, name: Optional[str] = None,
             parent: Optional[str] = None) -> int
    def send_group(self, id: str, name: str, parent: Optional[str] = None) -> int
    def remove(self, id: str) -> None
    def clear(self) -> None
    def inventory(self) -> list[InventoryEntry]
    def flush(self, settle: float = 0.5, *, raise_on_refusal: bool = True) -> list[Refusal]
    def release_data(self) -> None
    def close(self) -> None
    @property
    def revisions(self) -> Mapping[str, int]       # read-only copy
    @property
    def scene_origin(self) -> tuple[float, float, float]
    @property
    def scene_crs(self) -> str
    @property
    def connected(self) -> bool
    def __enter__(self) -> "Viewer"; def __exit__(self, *exc) -> None
```

Behaviour of each method:

- `send` connects when it is not connected (`connect`). `obj` is a `SceneObject`,
  or an object that an adapter accepts (see S2). `id` is required. `name` defaults
  to `id`. It returns the new revision. The revision of an id is `previous + 1`,
  where `previous` is the larger of the revision that this `Viewer` sent and the
  revision in the inventory of the last `hello` answer; for a new id it is 1
  (FR-003, scenario 2). It raises `InvalidInputError` for an unaccepted input
  (the text names the accepted inputs of FR-010), `ObjectTooLargeError` before any
  byte goes out, and `ViewerClosedError` when the socket closes.
  After each frame the method drains the socket without a wait and adds each
  `error` frame that has an `object` to an internal list of `Refusal`.
- `flush` waits up to `settle` seconds for late `error` frames, returns the list
  (and clears it), and raises `ViewerRefusedError` when the list is not empty and
  `raise_on_refusal` is true.
- `release_data` closes the data connection and stops the keepalive thread. It
  keeps `revisions`. A later `send` opens a new data connection (FR-014). If the
  viewer answers `busy`, `send` and `connect` raise `BusyError` ("Another program
  holds the data connection of the viewer. It can call `release_data()`."). The
  objects of the other program do not change, because the viewer closes only the
  second connection.
- `close` calls `release_data`. `Viewer` is not thread safe. One thread uses it at
  a time. A lock inside protects the socket against the keepalive thread only.
- `clear` sends `scene.clear`, `remove` sends `scene.remove`; both drop the
  revisions of the removed ids. `inventory` returns the entries of the last
  inventory. It does not poll the viewer (the viewer sends the inventory after
  `hello` only).

```python
# discovery.py
@dataclass(frozen=True)
class ViewerAddress: host: str; port: int; token: str
DEFAULT_PORT = 47615
def settings_dir() -> Path          # %APPDATA%\geoviewer on Windows, else ~/.geoviewer
def find_viewer(*, port: Optional[int] = None, token: Optional[str] = None,
                host: str = "127.0.0.1", settings: Optional[Path] = None,
                timeout: float = 0.5) -> ViewerAddress
```

`find_viewer` order for the token: the argument, the environment variable
`GEOVIEWER_TOKEN` (empty is not a token), the field `token` of `link.json` in the
settings folder. No token gives `TokenError` ("Set GEOVIEWER_TOKEN, or save a
token in the link popup of the viewer."). Order for the port: the argument; else
the instance files `instances/<pid>.json` (fields `version`, `pid`, `port`,
`kind`, `started`, `focused`): drop the files whose port does not accept a TCP
connection within `timeout`, prefer `kind == "window"` over `"headless"`, then
the largest `focused`; damaged files are ignored; else `DEFAULT_PORT` if it
accepts a connection; else `NoViewerError`. The library never starts a viewer
(scenario 3). A wrong token gives HTTP 401 at the handshake, and the transport
maps it to `TokenError`.

```python
# _websocket.py  (port of loopstructural/main/geoviewer/websocket_client.py in the plugin)
class WebSocketClient:
    @classmethod
    def connect(cls, host: str, port: int, *, token: str = "", timeout: float = 5.0,
                path: str = "/") -> "WebSocketClient"
    def send_binary(self, data: bytes) -> None
    def receive(self, timeout: float) -> Optional[bytes]
    def close(self) -> None
class WebSocketError(LoopLinkError); class WebSocketClosed(WebSocketError)
MAX_INCOMING_BYTES = MAX_FRAME_BYTES + MAX_HEADER_BYTES   # the plugin limit of 16 MiB is too small for S5
```

The port changes three things: the base class of the errors, the incoming limit
(above), and `connect` maps an HTTP 401 to `TokenError` and a refused connection to
`NoViewerError`. The handshake header is `Authorization: Bearer <token>`. Data
clients use the path `/`. The path `/tool` is for S4.

```python
# encode.py
def encode_object(obj: SceneObject, *, id: str, name: str, parent: Optional[str],
                  revision: int, origin: Sequence[float]) -> bytes
def hello_frame(crs: str, origin: Sequence[float], version: str) -> bytes
def heartbeat_frame(seq: int) -> bytes
def remove_frame(id: str) -> bytes
def clear_frame() -> bytes
def group_frame(id: str, name: str, parent: Optional[str], revision: int) -> bytes
```

`encode_object` gives the `scene.upsert` header with these fields: `type`, `id`,
`object_type`, `name`, `parent`, `revision`, `edit` (`[]`), `replaces` (`null`),
`geometry`, `params`, `data`, `property_groups`, `drape` (`null`), `style` (`{}`).
Mapping of the types:

| Class | `object_type` | `geometry` roles and buffers | `params` |
|---|---|---|---|
| `Points`, `Vectors` | `Points` | `vertices` | none |
| `Lines` | `Curve` | `vertices`, `segments` (`u32`) | none |
| `Surface` | `Surface` | `vertices`, `triangles` (`u32`) | none |
| `TetraMesh` | `TetraMesh` | `vertices`, `tetrahedra` (`u32`) | none |
| `BlockModel` | `BlockModel` | none | `origin` (3, minus scene origin), `cell_sizes` (3), `cells` (3) as JSON lists |
| `Drillhole` | `Drillhole` | `collar` (3, minus origin), `survey` `f64` (k, 3), `intervals` `f64` (j, 2) | `end_of_hole` as `[value]` when set |

Coordinate buffers (`vertices`, `collar`) are `f32` when the largest absolute
value after the origin shift is at most `F32_LIMIT`, and `f64` otherwise. Before
the frame goes out, `encode_object` sums the declared bytes of all buffers (and
`entries * 8 + name bytes` for each `value_map`). A sum over `MAX_OBJECT_BYTES`
or a frame over `MAX_FRAME_BYTES` raises `ObjectTooLargeError` (edge case, FR-008).

### Messages and saved data

No new message and no saved data. The library writes: `hello` (`application`
`{"name": "loop_link", "version": __version__}`, `crs`, `origin`), `heartbeat`
(`seq`), `scene.upsert`, `scene.remove`, `scene.clear`. It reads: `hello`,
`scene.inventory`, `error`, `heartbeat`. Example `hello` header:

```json
{"type":"hello","application":{"name":"loop_link","version":"0.1.0"},
 "crs":"EPSG:28350","origin":[500000.0,7000000.0,0.0],
 "version":{"major":1,"minor":14},"buffers":[]}
```

A reply with an unknown type is ignored. An `error` with the code
`unsupported_version` raises `MissingFeatureError("protocol 1.x")` with the
`supported` field in the text.

### Behaviour

**Connect** (`Viewer.connect`): `find_viewer`; `WebSocketClient.connect`; send
`hello`; read until both `hello` and `scene.inventory` arrived (timeout `timeout`,
then `ViewerClosedError`). An `error` frame with the code `busy` raises
`BusyError` and closes the socket. Another `error` frame raises `ViewerRefusedError`
with a `Refusal` that has an empty `object_id`. The scene origin is the `origin`
of the viewer `hello`; when it is `null` the origin that the caller proposed
applies. The library keeps both `scene_origin` and `scene_crs`. When the scene
origin differs from the proposed origin, all frames use the scene origin
(the viewer rule; the client converts). The revisions of the inventory fill
the internal table `max(local, inventory)`.

**Keepalive**: the viewer closes a data connection that sends no frame for 15 s
(spec 003 of geoviewer). After `connect`, a daemon thread `loop-link-heartbeat`
tries every `heartbeat_interval` seconds (default 5.0) to take the socket lock
without a wait. When it gets the lock and the last send was more than
`heartbeat_interval` ago, it sends `heartbeat_frame(seq)`. When the lock is not
free, a send is running, which counts as traffic. A socket error stops the thread
and marks the connection closed. `release_data` and `close` stop the thread
(`threading.Event`) and join it with a timeout of `2 * heartbeat_interval`.

**Plugin API** (S3, `push.py`):

```python
@dataclass
class PushResult:
    sent: int = 0
    removed: int = 0
    refused: list[Refusal] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)   # (id, reason)
    cancelled: bool = False
    def summary(self) -> str
def push_objects(viewer: Viewer, objects: Iterable[tuple[str, Any]], *,
                 replace_prefix: Optional[str] = None,
                 progress: Optional[Callable[[str], None]] = None,
                 cancel: Optional[threading.Event] = None) -> PushResult
```

`objects` yields `(id, object)` pairs. The function connects the viewer when it
is not connected. For each pair it checks `cancel` (stop with `cancelled=True`),
calls `progress("Sending <id>...")`, and calls `viewer.send(obj, id=id)`. An
`InvalidInputError` or `ObjectTooLargeError` for one object goes to `skipped` and
the push goes on (the plugin rule). A `ViewerClosedError` or `BusyError` stops the
push and goes up. After the loop, when `replace_prefix` is set, it sends
`viewer.remove(id)` for each id of the inventory that starts with the prefix and
was not sent in this call. Then it calls `viewer.flush(raise_on_refusal=False)` and
puts the refusals in `refused`. The function uses no Qt and no module state. It
calls `progress` on the thread that calls it. The plugin runs it on a worker thread.

### Stages and rollback

| Stage | Files and behaviour | Acceptance check | Rollback target and compatibility check |
|---|---|---|---|
| S1 | New package: codec, errors, WebSocket, discovery, objects, orientation, encode, `Viewer`, fake viewer; repository registration (ROADMAP, release-please, CI, Ruff, uv). Inputs: the dataclasses and numpy arrays only. | `uv run --package loop-link pytest packages/loop_link/tests`; `ruff check LoopStructural packages`; root `pytest`; SC-001, SC-004 (with the variable), SC-005, SC-006; manual check M1 | Revert the stage: the files are new, and the root files change by added lines only. The root `pytest` result is the same as before. |
| S2 | `adapters/` for pyvista, LoopStructural, pandas; `Viewer.send` calls `adapters.convert`. | The same commands. Adapter tests (pyvista, LoopStructural, pandas extras installed in the test environment). | Revert S2. S1 still sends the dataclasses. |
| S3 | `push.py` and its tests; frame comparison with the plugin frames; the slow-viewer and cancel tests. | The same commands; SC-004 frames of the plugin. | Revert S3. S1 and S2 stay. |
| S4 | Planned in "Design of stage S4": S4a (frames, commands module, fake commands), S4b (filter and colour members), S4c (colour in send), S4d (filter wait), S4e (real codec checks, manual check M4, docs). | See "Stages and rollback of S4". | See "Stages and rollback of S4". |
| S5 | Planned in "Design of stage S5": S5a (tool connection and selection), S5b (fetch), S5c (real codec checks, manual check M3, docs). | See "Stages and rollback of S5". | See "Stages and rollback of S5". |

## Design of stage S5 (stories 2 to 4)

Planned 2026-11-07, after S1 to S3 were built (PR 406) and geoviewer spec 067 merged
(protocol 1.18). S5 adds the tool connection, the selection calls and the fetch. S4
(filter and colour) is not in S5, but S4 reuses the tool connection of S5.

### Facts that S5 uses (geoviewer spec 067, `docs/protocol.md` version 1.18)

- A tool connects to the same server on the path `/tool`, with the same bearer token.
  It sends no `hello`. The viewer accepts 8 tool connections at the same time. A tool
  connection does not take the one data slot, so it works while another program holds
  the data connection (FR-014).
- Each tool message has the field `request` (a number that the tool chooses). Replies
  have the same number. The new messages: `selection.set`, `selection.get`,
  `selection.result`, `object.fetch`, `object.header`, `object.part`, `object.end`.
  Fixtures of these messages are in `fixtures/protocol/tool/` of the geoviewer
  repository. They use the version `{"major":1,"minor":18}`. `command.cancel` and
  `command.progress` come from spec 035 (the fixture of `command.cancel` has minor 11).
- `selection.set`: fields `object`, optional `revision`, `entries` (list of
  `{"kind": "cells" | "nodes", "indices": <buffer name>}`), one `u32` buffer for each
  entry. An empty `entries` list clears the selection of the object. A kind that is
  not in the list is cleared. The viewer counts a repeated index one time. The reply
  is `selection.result` with `object`, `revision` and `entries` that have `kind` and
  `count` only.
- `selection.get`: no fields. The reply is `selection.result` with one entry for each
  kind that has a selection, and a `u32` buffer for each entry (`indices` gives its
  name). The selection is the selection of the table that the user (or a tool)
  changed last, so it holds one object. With no table the reply has no `object`, no
  `revision` and `"entries": []`.
- `object.fetch`: field `id`. Replies, in this order: `object.header` (fields
  `origin`, `crs`, `object` = the `scene.upsert` fields of the object with no buffers,
  and `arrays` = the list of `{name, dtype, shape}` that the parts fill), then
  `object.part` messages (fields `array`, `offset`, `bytes`; one `u8` buffer named
  `part`; at most 64 MiB; `offset` is a multiple of 8), then `object.end` (field
  `parts`). Between parts the viewer sends `command.progress` with the same
  `request`. A tool cancels with `command.cancel` (field `request`); the reply is
  `request.error` with the code `cancelled`. A fetch needs no window.
- All coordinates of a fetch are relative to `origin` (f32 or f64 as stored).
- `request.error` has `error` = `{code, message, field}`. The codes that S5 meets:
  `invalid_input` (`field` is `object`, `id`, `kind`, `indices` or `revision`),
  `needs_window`, `cancelled`, `busy`, `failed`, `too_large`.
- A viewer older than 1.18 answers a tool message that it does not know with an
  `error` frame that has the code `unknown_message`. The tool connection has no
  `hello`, so this reply is the only way to detect an old viewer (FR-013).

### Names and errors

All new names below are free (checked with `git grep -n -w` on 2026-11-07).

New errors in `errors.py` (all derive from `LoopLinkError`):

```python
class StaleRevisionError(LoopLinkError):
    # attributes: object_id: str, revision: Optional[int]
    def __init__(self, object_id: str, revision: Optional[int] = None) -> None
    # text: "The object '<id>' is not at that revision in the viewer. Read the
    #        selection or fetch the object again, then try again."
class CancelledError(LoopLinkError):
    def __init__(self, message: str = "The request was cancelled.") -> None
class UnsupportedObjectError(LoopLinkError):
    # attributes: object_id: str, object_type: str
    def __init__(self, object_id: str, object_type: str) -> None
    # text: "The object '<id>' has the type '<type>'. loop_link can fetch Points, Curve,
    #        Surface, BlockModel, TetraMesh and Drillhole objects only."
```

`loop_link/__init__.py` exports them, and `Selection`, `SelectionWatcher`,
`FetchResult`.

Mapping of a tool failure to a library error (function `raise_for_request_error(error:
Mapping[str, Any], *, object_id: Optional[str]) -> NoReturn` in `tool.py`):

| Tool reply | Library error |
|---|---|
| `request.error`, code `invalid_input`, field `revision` | `StaleRevisionError(object_id)` |
| `request.error`, code `invalid_input`, other field | `InvalidInputError("The viewer refused the request: <message>")` |
| `request.error`, code `needs_window` | `MissingFeatureError("a viewer window (this viewer has no window)")` |
| `request.error`, code `cancelled` | `CancelledError()` |
| `request.error`, code `busy` | `BusyError()` |
| `request.error`, any other code | `ViewerRefusedError([Refusal(object_id or "", code, message)])` |
| `error` frame, code `unknown_message` | `MissingFeatureError(feature)`, where `feature` is the argument of `request` (for S5: `"selection and fetch (geoviewer protocol 1.18)"`) |
| `error` frame, code `unsupported_version` | `MissingFeatureError("protocol " + str(supported))` |
| other `error` frame | `ViewerRefusedError([Refusal("", code, message)])` |
| HTTP 401 at connect | `TokenError()` (the client already does this) |
| socket closed, or no frame in time | `ViewerClosedError` |

### Protocol changes (`protocol.py`, `encode.py`)

```python
TOOL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 18}
CANCEL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}
def encode_frame(header, buffers=(), *, version: Optional[dict[str, int]] = None) -> bytes
    # version=None keeps PROTOCOL_VERSION (1.14). Else the header gets this version.
```

New frame builders in `encode.py` (each returns the bytes of one frame):

```python
def selection_set_frame(request: int, object_id: str, cells: Optional[np.ndarray],
                        vertices: Optional[np.ndarray], revision: Optional[int]) -> bytes
def selection_get_frame(request: int) -> bytes
def fetch_frame(request: int, id: str) -> bytes
def cancel_frame(request: int) -> bytes     # command.cancel, version CANCEL_PROTOCOL_VERSION
```

`selection_set_frame`: the header is `{"type": "selection.set", "request": n, "object":
id, ["revision": r,] "entries": [...]}`. An entry is added for `cells` (kind `cells`,
buffer name `cells`) and for `vertices` (kind `nodes`, buffer name `nodes`) only when
the array is not empty. Each buffer is `u32` with shape `[k]`. With no entry, there is
no buffer. The `revision` key is left out when `revision` is `None`. The version is
`TOOL_PROTOCOL_VERSION`. `selection_get_frame` and `fetch_frame` use the same version.

### Tool connection (`tool.py`)

```python
class ToolConnection:
    def __init__(self, *, port: Optional[int] = None, token: Optional[str] = None,
                 host: str = "127.0.0.1", settings_dir: Optional[Path] = None,
                 timeout: float = 5.0) -> None
    @property
    def connected(self) -> bool
    def connect(self) -> ToolConnection          # no hello; path "/tool"
    def request(self, build: Callable[[int], bytes], *, until: Sequence[str],
                feature: str, object_id: Optional[str] = None,
                on_frame: Optional[Callable[[Frame], None]] = None,
                progress: Optional[Callable[[float], None]] = None,
                cancel: Optional[threading.Event] = None,
                timeout: Optional[float] = None) -> Frame
    def close(self) -> None
```

Behaviour of `request`:

1. It takes the lock `self._lock` (one request at a time on one connection). It
   connects when there is no socket (`discovery.find_viewer`, then
   `WebSocketClient.connect(host, port, token=..., timeout=..., path="/tool")`). It
   sends no `hello`.
2. It makes the number `n` (starts at 1, increases by 1 for each request of this
   connection), sends `build(n)`, and notes the time.
3. It reads frames with `client.receive(0.1)` in a loop. It decodes each with
   `protocol.decode_frame`. A `ProtocolError` closes the connection and raises
   `ViewerClosedError("The viewer sent a frame that is not valid. ...")`.
4. A frame with `type == "error"` goes to `raise_for_request_error` (the `error` rows
   of the table). A frame with another `request` number is dropped (a late frame of an
   older request). For a frame with the number `n`: `request.error` goes to
   `raise_for_request_error`; `command.progress` calls `progress(fraction)` when
   `fraction` is a number and `progress` is not `None`; a frame whose type is in
   `until` ends the loop and is returned (before that, `on_frame` has seen it too);
   any other frame goes to `on_frame` when it is not `None`, else it is dropped.
   `on_frame` sees every frame of the request except `command.progress` and
   `request.error`.
5. The wait time restarts at each frame of the request. When `timeout` seconds (default
   `self._timeout`) pass with no frame, the function closes the connection and raises
   `ViewerClosedError("The viewer did not answer in time. ...")`.
6. When `cancel` is set, the function sends `cancel_frame(n)` one time and goes on
   reading until the `cancelled` error arrives (then `CancelledError` goes up), or the
   `until` frame arrives (the request ended first, so it is returned), or the time
   runs out.
7. A socket error sets the socket to `None` and raises `ViewerClosedError`. The next
   call of `request` connects again.

`connect` raises `NoViewerError` and `TokenError`. A tool connection has no keepalive
thread (the viewer sets no idle time for it).

### Selection (`selection.py`)

```python
@dataclass(frozen=True, eq=False)
class Selection:
    object_id: Optional[str]        # None when nothing is selected
    revision: Optional[int]
    cells: np.ndarray               # uint32, 1-D, empty when none
    vertices: np.ndarray            # uint32, 1-D, empty when none
    @property
    def empty(self) -> bool         # no cells and no vertices
    def same_as(self, other: Selection) -> bool
        # same object_id, and the sorted cells and the sorted vertices are equal.
        # The revision is not compared.

def parse_selection_result(frame: Frame) -> Selection
    # kind "cells" gives cells, kind "nodes" gives vertices. The array is
    # frame.buffers[entry["indices"]]. An entry with a count but no buffer, a buffer
    # that is not uint32, or a length that differs from count gives ProtocolError.

class SelectionWatcher:
    @property
    def running(self) -> bool
    @property
    def error(self) -> Optional[BaseException]    # the error that stopped the thread
    def stop(self) -> None                         # stop, then join (timeout 2 * interval + 1 s)
    def __enter__(self) -> SelectionWatcher
    def __exit__(self, *exc) -> None
```

Only `Viewer.watch_selection` makes a `SelectionWatcher`. It has a daemon thread
`loop-link-selection-poll`. The thread holds a `weakref` to the `Viewer`, so a dropped
`Viewer` ends it. Loop of the thread: wait `interval` seconds on a `threading.Event`
(so `stop` is fast), call `viewer.selection()`, and when the result is not `same_as`
the last one, call `callback(selection)` one time and keep it as the last one. An
exception of the poll or of the callback goes to `error`, and the thread ends. No
thread runs before `watch_selection` is called (FR-005).

### Fetch (`fetch.py`)

```python
@dataclass(frozen=True, eq=False)
class FetchResult:
    id: str
    object: SceneObject          # Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole
    object_type: str             # the viewer type name, for example "Curve"
    revision: int
    name: Optional[str]
    parent: Optional[str]
    origin: Optional[tuple[float, float, float]]
    crs: Optional[str]

def fetch_object(tool: ToolConnection, id: str, *,
                 progress: Optional[Callable[[float], None]] = None,
                 cancel: Optional[threading.Event] = None,
                 timeout: float = 30.0) -> FetchResult

def decode_object(header: Mapping[str, Any], arrays: Mapping[str, np.ndarray],
                  origin: Optional[Sequence[float]]) -> SceneObject
```

`fetch_object` calls `tool.request(lambda n: encode.fetch_frame(n, id), until=("object.end",),
feature="selection and fetch (geoviewer protocol 1.18)", object_id=id,
on_frame=assembler.feed, progress=progress, cancel=cancel, timeout=timeout)`. The
private class `_Assembler` keeps the state and the rules:

1. The first frame must be `object.header`. Any other first frame gives `ProtocolError`.
2. Header checks. Each entry of `arrays` has a unique `name`, a `dtype` in `f32 f64 u32
   i32 u8`, and a `shape` of non-negative integers. The sum of the bytes of all arrays
   must be at most `protocol.MAX_OBJECT_BYTES`, else `ObjectTooLargeError(id, size,
   limit)`. Each buffer name in `object.geometry`, in `object.data[*].buffer`, and in the
   `offsets` fields must be in `arrays`, else `ProtocolError`. Then the assembler makes
   one `numpy.empty(nbytes, dtype=numpy.uint8)` for each array.
3. `object.part`: the array must be in `arrays`; `offset % 8 == 0`; the buffer `part` must
   be `u8` with `bytes` entries; `offset` must equal the bytes already filled in that
   array (parts of one array arrive in order, with no gap and no overlap); `offset +
   bytes` must be at most the size of the array. A break gives `ProtocolError`. The
   assembler copies the part into the array.
4. `object.end`: `parts` must equal the number of parts received, and each array must be
   full. A break gives `ProtocolError`.
5. The assembler views each array with its dtype and shape (`numpy.ndarray.view`,
   `reshape`) and calls `decode_object(header["object"], arrays, header["origin"])`. The
   `revision` of the result is `header["object"]["revision"]`.

`ProtocolError` is the existing `loop_link.protocol.ProtocolError`. `fetch_object`
catches it, closes the tool connection (the stream is out of step), and raises
`ViewerClosedError("The viewer sent an object that is not valid: <text>")`.
`CancelledError` and the other library errors go up unchanged.

Rules of `decode_object` (world coordinates, no other change of values):

| `object_type` | Result |
|---|---|
| `Points` | `Vectors` when `property_groups` has an entry with the role `orientation` or `lineation`; else `Points`. |
| `Curve` | `Lines(vertices, segments, ...)` |
| `Surface` | `Surface(vertices, triangles, ...)` |
| `BlockModel` | `BlockModel(origin = params["origin"] + origin, cell_sizes = params["cell_sizes"], cells = params["cells"], ...)` |
| `TetraMesh` | `TetraMesh(vertices, tetrahedra, ...)` |
| `Drillhole` | `Drillhole(collar + origin, survey, intervals, end_of_hole = params["end_of_hole"][0] when present, ...)` |
| any other (`Grid2D`, `GeoImage`, `Group`) | `UnsupportedObjectError(id, object_type)` |

- `origin` is `None` or three numbers; `None` means `(0, 0, 0)`. The `vertices` and
  `collar` are computed in float64: `array.astype(float64) + origin`. The result is
  float64. This is the only dtype change: the viewer holds f32 relative to the origin,
  and the sum in f64 is the world coordinate.
- Property rules, for each entry of `object["data"]` (the name is `name`; the values are
  the array of `buffer`):
  - Association `vertex` goes to `properties`. Association `cell` goes to
    `cell_properties`, except for a Drillhole, where `cell` and `interval` (with no
    `log_set`, or with the first interval log set) go to `interval_properties`. For
    `Points` or `Vectors`, an association `cell` gives `UnsupportedObjectError`.
  - A property of another log set of a Drillhole, a property with the association `sample`
    or `object`, and any other association are dropped. The function calls
    `warnings.warn("loop_link dropped the property '<name>' of object '<id>': <reason>",
    UserWarning, stacklevel=2)` one time for each (OP10).
  - Types: `data_type` `float` keeps the stored dtype (f32 stays f32, f64 stays f64);
    `integer` gives int32; `boolean` gives `numpy.bool_` (`values != 0`); `text` gives a
    numpy `U` array (decode each item as UTF-8 from the data and the `offsets` arrays);
    `referenced` gives `Categorical(codes, names, colours)`: `names` are the `name` fields
    of `value_map` ordered by `id`; the codes are renumbered with a lookup when the ids
    are not `0..n-1`; a code that is not in `value_map` gives `InvalidInputError`;
    `colours` are the `colour` lists as 4-tuples (alpha 255 when a list has 3 values).
  - Vectors: for each entry of `property_groups` with the role `orientation` (members
    `dip_direction`, `dip`) or `lineation` (members `trend`, `plunge`), the member arrays
    (as float64) form the tuple `orientation` or `lineation`. The member arrays are not
    in `properties`.

### Viewer (`viewer.py`)

New members (the code of the data connection does not change):

```python
class Viewer:
    @property
    def tool_connected(self) -> bool
    def select(self, object_id: str, *, cells: Optional[ArrayLike] = None,
               vertices: Optional[ArrayLike] = None,
               revision: Optional[int] = None) -> dict[str, int]
    def selection(self) -> Selection
    def watch_selection(self, callback: Callable[[Selection], None], *,
                        interval: float = 0.5, initial: bool = False) -> SelectionWatcher
    def fetch(self, id: str, *, progress: Optional[Callable[[float], None]] = None,
              cancel: Optional[threading.Event] = None) -> FetchResult
```

- The `Viewer` makes its `ToolConnection` at the first call of `select`, `selection`,
  `watch_selection` or `fetch`, with the same `port`, `token`, `host`, `settings_dir` and
  `timeout`. These calls do not open the data connection. So a program that gets
  `BusyError` from `connect` or `send` can still select, read and fetch (FR-014).
  `release_data` does not close the tool connection. `close` stops all watchers, closes
  the tool connection and marks the viewer closed. A later tool call raises
  `ViewerClosedError`, as `send` does.
- `select(object_id, cells, vertices, revision)`: `object_id` must be a non-empty `str`.
  At least one of `cells`, `vertices` must not be `None`, else `InvalidInputError("Give
  cells or vertices. Use an empty list to clear the selection.")`. Each given array must
  be a 1-D sequence of integers from 0 to 2**32 - 1 (no float, no bool, no negative
  value), else `InvalidInputError` that names the argument. The call sends one
  `selection.set` frame (all indices in one message, SC-003). It returns
  `{"cell": n, "vertex": m}` with the counts of the reply (only the kinds that the reply
  has). Empty lists clear the selection. When the viewer refuses the call, the selection
  does not change.
- `selection()` sends `selection.get` and returns a `Selection`.
- `watch_selection(callback, interval, initial)`: `interval` must be more than 0, else
  `InvalidInputError`. The call first reads the selection in the calling thread, so
  `MissingFeatureError` and `ViewerClosedError` go to the caller. When `initial` is true
  and this first selection is not empty, the callback is called once with it. Then the
  call starts the `SelectionWatcher` thread. The `Viewer` keeps a list of its watchers.
- `fetch(id, progress, cancel)` calls `fetch.fetch_object(...)` on the tool connection.
  `progress` gets a float from 0 to 1. `cancel` is a `threading.Event`; when it is set,
  the fetch stops and raises `CancelledError`. An unknown id gives `InvalidInputError`
  (the viewer reply `invalid_input`, field `id`).

### Fake viewer changes (`testing/fake_viewer.py`)

The tool endpoint answers the 1.18 messages. New constructor arguments (all have
defaults, so no S1 to S3 test changes): `protocol_minor: int = 18` (a value below 18
makes each selection or fetch message get an `error` frame with the code
`unknown_message`), `part_bytes: int = 64 * 1024 * 1024` (the largest part; it must be a
multiple of 8), `fetch_delay: float = 0.0` (the pause after each part; it makes a cancel
test possible). `kind="headless"` (it exists) makes `selection.set` and `selection.get`
answer `request.error` with `needs_window`. New members:

- `.selection` (property): the current selection as `{"object": str, "revision": int,
  "cells": list[int], "nodes": list[int]}`, or `None`.
- `.set_user_selection(object_id: str, *, cells=None, nodes=None) -> None`: it simulates a
  selection that the user made in the window. It uses the checks of `selection.set`.
- `.tool_requests` (list of `Frame`): each tool frame that arrived, for tests.

Rules of the fake tool endpoint (they follow spec 067; the real viewer is the
authority, and manual check M3 repeats them):

- `selection.set`: an unknown object gives `invalid_input` (field `object`). A `revision`
  that is not the current revision gives `invalid_input` (field `revision`). An object
  with no table (the types `Drillhole` and `Group`, see OP12), or a kind that the object
  lacks (`Points` has `nodes` only), gives `invalid_input` (field `kind`). An index at or
  above the element count gives `invalid_input` (field `indices`). A buffer that is not
  `u32`, or is not in the frame, gives `invalid_input` (field `indices`). Counts: nodes
  are the vertex count (a BlockModel has `(nx+1)(ny+1)(nz+1)`); cells are the segments,
  triangles, tetrahedra, or `nx*ny*nz`. A repeated index counts one time. The new
  selection replaces the old one, also when the old one is on another object. An empty
  `entries` list clears it. The reply has `object`, `revision`, and `entries` with `kind`
  and `count`.
- `selection.get`: the reply has `object`, `revision`, and one entry for each kind that
  has indices, with a `u32` buffer named `cells` or `nodes`. With no selection the reply
  has `"entries": []`.
- `object.fetch`: an unknown id gives `invalid_input` (field `id`). Else the fake sends
  `object.header` (`origin` and `crs` of the fake, the stored header of the object with
  `buffers: []`, and `arrays` = the stored buffers), then `object.part` frames of at
  most `part_bytes` for each array in order, then `object.end`. When there is more than
  one part, a `command.progress` frame (`fraction`, `text`) follows each part. After
  each part the fake reads for at most 0.001 s. When it finds a `command.cancel` with the
  same `request`, it sends `request.error` with the code `cancelled` and stops. A fetch
  does not change the selection or the objects.

### Files (S5)

- `packages/loop_link/src/loop_link/errors.py` - Change. The three new errors.
- `packages/loop_link/src/loop_link/protocol.py` - Change. `TOOL_PROTOCOL_VERSION`,
  `CANCEL_PROTOCOL_VERSION`, `encode_frame(..., version=None)`.
- `packages/loop_link/src/loop_link/encode.py` - Change. The four frame builders.
- `packages/loop_link/src/loop_link/tool.py` - Create. `ToolConnection`,
  `raise_for_request_error`.
- `packages/loop_link/src/loop_link/selection.py` - Create. `Selection`,
  `parse_selection_result`, `SelectionWatcher`.
- `packages/loop_link/src/loop_link/fetch.py` - Create. `FetchResult`, `fetch_object`,
  `decode_object`, `_Assembler`.
- `packages/loop_link/src/loop_link/viewer.py` - Change. The new members.
- `packages/loop_link/src/loop_link/__init__.py` - Change. The new exports.
- `packages/loop_link/src/loop_link/testing/fake_viewer.py` - Change. The tool endpoint.
- `packages/loop_link/tests/test_tool.py` - Create. `ToolConnection` and each row of the
  error table; request numbers; a late frame of an older request is dropped; progress;
  cancel; the connection closes in the middle of a request; reconnect at the next
  request; wrong token; no viewer; `unknown_message` with `protocol_minor=14`.
- `packages/loop_link/tests/test_selection.py` - Create. Story 2 scenarios 1 to 4; story 3
  scenarios 1 to 3; the clear rule and the duplicate rule; `cells` and `vertices`
  together; `StaleRevisionError`; `needs_window` for `kind="headless"`; a selection of
  100,000 cells is one frame (SC-003: `fake.tool_requests` has one `selection.set` frame
  with a `u32` buffer of shape `[100000]`); a selection that the user made
  (`set_user_selection`); `watch_selection` (the callback runs once for a change and not
  for no change, `initial`, `stop`, no thread before the call, the error of a callback
  is in `.error`, a dropped `Viewer` ends the thread); FR-014 (program A holds the data
  connection; program B gets `BusyError` from `send` but selects, reads and fetches; after
  A calls `release_data`, B sends and the objects of A stay).
- `packages/loop_link/tests/test_fetch.py` - Create. `decode_object` for each row of the
  table, each property rule, and each drop with its warning. The cases are built by hand
  as a viewer writes them: an orientation group, a value map with ids that are not
  `0..n-1`, text, boolean, a Drillhole with two log sets and a sample set, a non-zero
  origin, `UnsupportedObjectError` for a `Group`. `_Assembler` with hostile frames: the
  first frame is not a header; an array too large for the 2 GiB limit (a header only);
  a part for an array that the header does not list; an offset that is not a multiple of
  8; a part past the end; a gap; an overlap; a missing part; a wrong `parts` count; a
  `part` buffer that is not `u8`. A fetch of several parts (`part_bytes=64`) with
  `progress` calls that rise. A cancel with `fetch_delay=0.05` gives `CancelledError`
  and the scene stays the same. An unknown id gives `InvalidInputError`. A headless fake
  still fetches. The revision is in the result.
- `packages/loop_link/tests/test_roundtrip.py` - Create. SC-002 and story 4 scenario 2:
  for each of the seven types (and a `Vectors` made with `from_directions`), `send` then
  `fetch`; compare with `numpy.testing.assert_allclose(rtol=1e-6)` for floats and
  `array_equal` for integers, booleans, text, `Categorical` (codes, names, colours) and
  index arrays; the dtype of each property stays; the coordinates are float64 world
  coordinates for a scene with a non-zero origin and for a point 10 km from it.
- `packages/loop_link/tests/test_codec_real.py` - Change. Add the cases `selection.set`
  (cells and nodes, with and without `revision`, with no entry), `selection.get`,
  `object.fetch` and `command.cancel`, written by the new builders, and check that the
  real codec accepts them (the helper and the negative control stay).
- `packages/loop_link/tests/test_fixtures.py` - Change. When `LOOP_LINK_GEOVIEWER_FIXTURES`
  is set, read `<variable>/../tool/{selection.set,selection.get,object.fetch,command.cancel}.json`
  (OP11). For each, build the frame with the `request`, `object` or `id`, `revision` and
  indices of the fixture, and assert that the header equals the fixture (including
  `version` and `buffers`). Skip with a reason when the variable or the folder is missing.
- `packages/loop_link/tests/test_fake_viewer.py` - Change. Add the rules of the fake tool
  endpoint (one test for each row of the list above).
- `packages/loop_link/tests/manual_checks.md` - Change. Add M3 (see "Test strategy of S5").
- `packages/loop_link/README.md` - Change. A short section "Select and fetch" with an
  example of about 10 lines.
- `ROADMAP.md` - Change (last task of S5). Update the Stage 7 entry: S5 done.

### Stages and rollback of S5

| Stage | Files and behaviour | Acceptance check | Rollback target and compatibility check |
|---|---|---|---|
| S5a | `errors.py`, `protocol.py`, `encode.py`, `tool.py`, `selection.py`, `viewer.py` (select, selection, watch_selection), `__init__.py`, fake viewer (selection part), `test_tool.py`, `test_selection.py`, `test_fake_viewer.py` | `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests`; `ruff check LoopStructural packages`; `black --check packages/loop_link`; FR-004, FR-005, FR-013 (selection part), FR-014 (tool part), SC-003 | Revert S5a. S1 to S3 do not use the new files. The new argument of `encode_frame` has a default that keeps the old frames. |
| S5b | `fetch.py`, `viewer.py` (fetch), fake viewer (fetch part), `test_fetch.py`, `test_roundtrip.py` | The same commands; FR-006, FR-013 (fetch part), SC-002 | Revert S5b. S5a stays. |
| S5c | `test_codec_real.py`, `test_fixtures.py`, `manual_checks.md` (M3), `README.md`, `ROADMAP.md` | The same commands with `LOOP_LINK_GEOVIEWER_REPO` and `LOOP_LINK_GEOVIEWER_FIXTURES` set; M3 by hand with the 067 viewer build | Revert S5c. |

### Test strategy of S5

- Unit and regression tests: the files above. The fake viewer applies the rules of spec
  067, as the S1 fake applies the rules of the codec.
- Integration and compatibility: the real geoviewer codec checks the frames that the
  library writes (`test_codec_real.py`, needs `LOOP_LINK_GEOVIEWER_REPO`). The fixtures
  of `fixtures/protocol/tool` check the headers (`test_fixtures.py`).
- Manual check M3 (real window viewer built from geoviewer spec 067):
  1. Send a `BlockModel` of 10 x 10 x 10 cells. Select cells 12, 13 and 14 with
     `Viewer.select`. Check that the viewer outlines three cells and that
     `Viewer.selection()` returns them.
  2. Select cells in the window with the mouse. Read them with `Viewer.selection()`.
  3. `select` with an index of 5000 gives `InvalidInputError`, and the selection stays.
  4. Send, replace, then `select` with the old revision. It gives `StaleRevisionError`.
  5. Fetch each of the seven types back and compare with the original (the script of M2
     plus `fetch`).
  6. Fetch a drillhole that the viewer imported. Check the warning for the log sets that
     the library drops.
  7. Start `watch_selection`. Change the selection in the window. Check that the callback
     runs once.
  8. Start a second program that sends. Check that it gets `BusyError` while the first
     program still selects and fetches.
  M3 checks the points that the fake viewer assumes (OP12, OP13), and the order of the
  parts of a fetch (OP14). Record the time of a fetch of 1,000,000 cells (no limit).
- Backend coverage and cross-language fixtures: as in the main "Test strategy".
- Performance: SC-002 has no time limit.
- Merge gates: as in the main "Test strategy".

### Constitution check for S5

All rows of the table above hold for S5. The differences: the new files are under
`packages/loop_link` only; docstrings and type hints follow the same rules (`D` and
`ANN` apply to `src`); the ROADMAP entry is updated in the last task; each stage is a
separate revert.

## Design of stage S4 (stories 6 and 7)

Planned 2026-11-08, after S5 was built. S4 adds the filter calls and the colour calls. It
reuses the tool connection of S5 (`tool.py`). It needs no viewer change, except that the
list of colour map names (FR-020) needs a viewer of protocol 1.18 (geoviewer spec 067 gives
the names in the schema of `colour.set_map`; the text of this spec calls it "066").

### Facts that S4 uses (geoviewer `ui-fixes` build, read 2026-11-08)

The facts come from `crates/geoviewer/src/tool/input_colour.rs`, `tool/result.rs`,
`tool/host.rs`, `tool/input_object.rs`, `fixtures/protocol/tool/` and `docs/protocol.md`
of the geoviewer repository. The real viewer is the authority. Manual check M4 repeats the
points marked "assumed".

- A command goes on the tool connection (path `/tool`) in a `command.run` message:
  `{"type": "command.run", "version": {"major": 1, "minor": 11}, "request": n,
  "command": name, "input": {...}, "buffers": []}`. The fixtures use minor 11 for it. The
  reply is `command.result` with `"result": {"status": ...}`, or `request.error`.
- `commands.list` (field `request` only) has the reply `commands`: `{"type": "commands",
  "request": n, "viewer": "window" | "headless", "commands": [{"name", "description",
  "needs_window", "tool_available", "unavailable_reason", "input_schema",
  "result_schema"}]}`.
- The seven commands that S4 uses, with their input. Each input rejects an unknown field.

| Command | Input | Note |
|---|---|---|
| `filter.edit` | `{"object": id, "edit": {"kind": ..., ...}}` | The kinds are in the next list. |
| `colour.set_active` | `{"id": id, "array": name}`; leave `array` out for no array | |
| `colour.set_map` | `{"id": id, "name": text}` | In a viewer of protocol 1.18 the schema of `name` is `{"type": "string", "enum": [...]}`. The names are `viridis`, `grey` and `coolwarm`. |
| `colour.set_range` | `{"id": id, "min": text, "max": text}` | The numbers are texts. |
| `colour.clear_range` | `{"id": id}` | |
| `colour.set_scale` | `{"object": id, "scale": "linear" or "log"}` | The key is `object`, not `id`. It needs a window. It does not wait when the scene is busy (the reply is `busy` at once). |
| `object.set_colour` | `{"id": id, "colour": [r, g, b] or null}` | `null` removes the plain colour. |

- The edits of `filter.edit` (field `kind`): `add` (`array`), `set_range` (`array`, `min`
  text, `max` text, `keep_nan` bool), `set_values` (`array`, `values`: list of `u32`, the
  `id` values of the `value_map` of the array), `set_on` (`array`, `on` bool),
  `set_exclude` (`array`, `exclude` bool), `remove` (`array`), `cancel` (no other field;
  it cancels the running job).
- The result of each of these commands is `{"status": "done"}` or `{"status":
  "refused", "text": ...}`. A filter edit that fails (for example an array that does not
  exist) gives `refused` with a text. The result has no count of elements that pass.
- A command that "waits when busy" waits at most 5 s for a scene that is not busy, then
  gets `request.error` with the code `busy`. All commands of the list above wait, except
  `colour.set_scale`.
- `filter.edit` is not a long command: it sends no `command.progress`, and it answers when the
  edit is accepted. The viewer builds the filter after that. No message tells the end of
  the build (OP15).
- Assumed (M4 checks): a filter on a text property, a filter on a property that does not
  exist, and a map, a range or a scale on a category or boolean property are `refused`
  with a text. A new revision of an object keeps its filters and its colour settings. A
  command for an unknown object gives `request.error` `invalid_input` or `refused`.
- The data array in `scene.upsert` can carry `colour_map` (`name`, `range`, `clip`,
  `scale`). It has no field for the array that colours the object or for the plain
  colour. See the answer of OP1 for the S4 decision.

### Names and errors

All new names below are free (checked with `git grep -n -w` on 2026-11-08).

No new error class. S4 uses the errors of S5 and S1:

| Reply | Library error |
|---|---|
| `command.result`, status `done` | none |
| `command.result`, status `refused`, with `text` | `InvalidInputError("The viewer refused the command '<command>': <text>")` |
| `command.result`, any other status | `ProtocolError("The viewer sent an unknown result status '<status>'.")` |
| `request.error` and `error` frames | `raise_for_request_error` (the table of S5). `needs_window` gives `MissingFeatureError`. `busy` gives `BusyError`. `invalid_input` gives `InvalidInputError`. |
| `error` frame, code `unknown_message` for `commands.list` or `command.run` | `MissingFeatureError("filter and colour commands (geoviewer protocol 1.11)")` |
| A map name that is not in the list of the viewer | `InvalidInputError("The viewer has no colour map '<name>'. The names are: <list>.")`, before any frame goes out |

### Protocol changes (`protocol.py`, `encode.py`)

```python
COMMAND_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}   # protocol.py
def commands_list_frame(request: int) -> bytes                          # encode.py
def command_frame(request: int, command: str, input: Mapping[str, Any]) -> bytes
```

`commands_list_frame` writes `{"type": "commands.list", "request": n, "buffers": []}`.
`command_frame` writes `{"type": "command.run", "request": n, "command": command,
"input": dict(input), "buffers": []}`. Both use the version `COMMAND_PROTOCOL_VERSION`, so
a viewer of protocol 1.11 to 1.17 accepts them (FR-018: no viewer change). The header must
equal `fixtures/protocol/tool/commands.list.json` and `command.run.json` when the
arguments are the same.

### Commands module (`commands.py`)

```python
@dataclass(frozen=True)
class ColourSettings:
    array: Optional[str] = None                      # property that colours the object
    colour_map: Optional[str] = None                 # name of a map of the viewer
    range: Optional[tuple[float, float]] = None      # (minimum, maximum)
    scale: Optional[str] = None                      # "linear" or "log"
    plain: Optional[tuple[int, int, int]] = None     # plain colour of the object
    # A field that is None is not sent. __post_init__ calls the check functions below.

def run_command(tool: ToolConnection, command: str, input: Mapping[str, Any], *,
                object_id: Optional[str] = None, timeout: Optional[float] = None,
                cancel: Optional[threading.Event] = None) -> Mapping[str, Any]
def list_colour_map_names(tool: ToolConnection) -> Optional[tuple[str, ...]]
def parse_colour_map_names(frame: Frame) -> Optional[tuple[str, ...]]
def raise_for_command_result(result: Mapping[str, Any], *, command: str) -> None
def format_number(value: float) -> str
def check_number(value: object, name: str) -> float
def check_text(value: object, name: str) -> str
def check_rgb(value: object, name: str) -> tuple[int, int, int]
def check_codes(values: object, name: str) -> list[int]
```

- `run_command` calls `tool.request(lambda n: encode.command_frame(n, command, input),
  until=("command.result",), feature="filter and colour commands (geoviewer protocol 1.11)",
  object_id=object_id, ...)`. It reads `frame.header["result"]`, calls
  `raise_for_command_result`, and returns the result mapping. A result that is not a
  mapping with a `status` text gives `ProtocolError`.
- `list_colour_map_names` calls `commands_list_frame` with `until=("commands",)` and returns
  `parse_colour_map_names(frame)`.
- `parse_colour_map_names` finds the entry with `name == "colour.set_map"`. It reads
  `input_schema["properties"]["name"]["enum"]`. It returns the tuple of names. It returns
  `None` when the entry, the property or the `enum` is missing (a viewer older than 1.18). A
  malformed field (the `enum` is not a list of text) gives `ProtocolError`.
- `format_number(value)` returns `repr(float(value))` (for example `"0.5"` or `"1e-05"`).
- `check_number` accepts `int` and `float` (and numpy numbers) that are not `bool` and not
  NaN or infinite, and returns `float`. Other input gives `InvalidInputError` that names
  `name`. `check_text` accepts a non-empty `str`. `check_rgb` accepts three integers from 0
  to 255 (not `bool`). `check_codes` accepts a sequence or array of integers from 0 to
  2**32 - 1 (not `bool`, not float) and returns a `list[int]`.

### Viewer (`viewer.py`)

New members. They use the tool connection of S5 (`self._tool()`), so they work while another
program holds the data connection (FR-014). They do not open the data connection. Each one
raises `InvalidInputError` for a wrong argument before any frame goes out.

```python
class Viewer:
    def filter_range(self, object_id: str, array: str, minimum: float, maximum: float, *,
                     keep_nan: bool = False, wait: bool = True,
                     cancel: Optional[threading.Event] = None, timeout: float = 60.0) -> None
    def filter_values(self, object_id: str, array: str, values: ArrayLike, *,
                      exclude: bool = False, wait: bool = True,
                      cancel: Optional[threading.Event] = None, timeout: float = 60.0) -> None
    def filter_enable(self, object_id: str, array: str, on: bool = True) -> None
    def filter_remove(self, object_id: str, array: str) -> None
    def cancel_filter(self, object_id: str) -> None
    def colour_map_names(self) -> tuple[str, ...]
    def set_colour_property(self, object_id: str, array: Optional[str]) -> None
    def set_colour_map(self, object_id: str, name: str) -> None
    def set_colour_range(self, object_id: str, minimum: float, maximum: float) -> None
    def clear_colour_range(self, object_id: str) -> None
    def set_colour_scale(self, object_id: str, scale: str) -> None
    def set_object_colour(self, object_id: str, colour: Optional[tuple[int, int, int]]) -> None
    def apply_colour(self, object_id: str, settings: ColourSettings) -> None
    # send() gets one new keyword-only argument, before **hints:
    #   colour: Optional[ColourSettings] = None
```

Filter rules (each call sends the edits in the order shown, with `run_command`; the first
refusal stops the call and raises):

- `filter_range`: `add`, `set_range` (`min` and `max` are `format_number` texts, `keep_nan`),
  `set_on` with `on` true. `minimum` must not be more than `maximum`, else
  `InvalidInputError`. If `add` is refused, the call goes on with `set_range` (OP16: the
  filter can exist already). If `set_range` is refused, its text goes up.
- `filter_values`: `add` (the same refusal rule), `set_values` (the codes from
  `check_codes`), `set_exclude` (`exclude`), `set_on` true. The codes are the `id` values of
  the `value_map` of the property. A `Categorical` that the library sent has the ids
  `0..n-1` in the order of its names (OP17: names as input are a later stage).
- `filter_enable`: one `set_on` edit. `filter_remove`: one `remove` edit. `cancel_filter`:
  one `cancel` edit.
- `wait` is true by default and is stage S4d. In S4a to S4c the argument exists, and a
  true value is accepted and does nothing more than return when the viewer accepted the
  edits. S4d adds the wait (OP15). `cancel` and `timeout` are used only by the wait.

Colour rules:

- `colour_map_names()` returns the names from `list_colour_map_names`. It keeps the tuple
  for the life of the tool connection. If the result is `None`, it raises
  `MissingFeatureError("colour map names (geoviewer protocol 1.18)")`.
- `set_colour_map(object_id, name)`: `name` must be a non-empty `str`. The method reads the
  names. When the viewer gives names and `name` is not one of them, it raises the error of
  the table above. When the viewer gives no names (older viewer), it sends the name and
  the viewer decides. It sends `colour.set_map` `{"id", "name"}`.
- `set_colour_property`: `colour.set_active` `{"id", "array"}`. With `array=None` the key
  `array` is left out (no property).
- `set_colour_range`: `colour.set_range` `{"id", "min", "max"}` with `format_number` texts.
  `minimum` must be less than `maximum`. `clear_colour_range`: `colour.clear_range`.
- `set_colour_scale`: `scale` must be `"linear"` or `"log"`, else `InvalidInputError`.
  It sends `colour.set_scale` `{"object", "scale"}`.
- `set_object_colour`: it sends `object.set_colour` `{"id", "colour"}` with the list
  `[r, g, b]` (from `check_rgb`) or `null`.
- `apply_colour(object_id, settings)`: it sends, in this order, `colour.set_active`,
  `colour.set_map`, `colour.set_range`, `colour.set_scale`, `object.set_colour`, each one
  only when its field is not `None`. It checks all fields before it sends any frame.
- `send(..., colour=settings)`: the call sends the object as before. When `colour` is not
  `None`, it then calls `self.flush()` (this waits 0.5 s for a late refusal of the object
  and raises it), and then `self.apply_colour(id, settings)`. If a colour command fails,
  the object is already in the viewer, and the error goes up. The user can call the
  `set_colour_*` methods later. The send returns the revision as before.

### Fake viewer changes (`testing/fake_viewer.py`)

New constructor arguments, all with defaults, so no earlier test changes:
`filter_delay: float = 0.0` (the seconds that the scene stays busy after a `filter.edit`
that changes a filter) and `busy_wait: float = 5.0` (the seconds that a waiting command
waits before `busy`; tests set it to 0.05). New members: `.filters` (property): a
`dict[str, dict[str, dict]]`, by object id and array name, with the keys `kind` (`"range"`
or `"values"`), `on`, `minimum`, `maximum`, `keep_nan`, `values`, `exclude`; `.colouring`
(property): a `dict[str, dict]` by object id with the keys `array`, `colour_map`,
`range`, `scale`, `plain`. Both return copies.

Rules of the fake tool endpoint (they follow the facts above; the real viewer is the
authority and M4 repeats them):

- `commands.list` answers `commands` with `viewer` = `self.kind` and one entry for each
  command of the table (and no other). `needs_window` is true only for `colour.set_scale`.
  `input_schema` of `colour.set_map` has `properties.name.enum` = `["viridis", "grey",
  "coolwarm"]` when `protocol_minor >= 18`, and `properties.name.type` = `"string"` only
  below 18. With `protocol_minor < 11` the endpoint answers `unknown_message`.
- `command.run` of an unknown command: `request.error` `invalid_input`, field `command`.
  An input with an unknown or a missing field, or a field of a wrong type (the key and type
  rules of the table above), gives `request.error` `invalid_input` with that field. An
  object that does not exist gives `request.error` `invalid_input`, field `id` (or
  `object` for the two commands that use that key). `kind="headless"` gives `needs_window`
  for `colour.set_scale` only.
- A command that waits: while the scene is busy (the time `filter_delay` after a filter
  edit), each command except `colour.set_scale` waits up to `busy_wait` seconds. Then it
  runs, or it gets `request.error` `busy` when the scene is still busy. `colour.set_scale`
  gets `busy` at once. `filter.edit` with `cancel` ends the busy time.
- Property rules, from the header of the stored object (`data` entries: `name`,
  `data_type`, `association`): a filter on an array that does not exist, or on a `text`
  array or an `object` association array, is `refused` with a text that names the array. A
  range edit on a `referenced` or `boolean` array, and a values edit on a number array,
  are `refused`. `colour.set_active` with an unknown array is `refused`. `colour.set_map`,
  `colour.set_range` and `colour.set_scale` on an object whose active array is `referenced`
  or `boolean` are `refused`. These rules are assumed from the spec (OP18).
- Effects: `add` creates an off filter. `set_range`, `set_values`, `set_exclude` and
  `set_on` change it (a missing filter is `refused`). `remove` deletes it. A
  `scene.upsert` of an object that already exists keeps its `filters` and `colouring`
  (story 6, scenario 6). When the new revision has no array of that name, the fake drops
  the filter of that array.
- The reply is `command.result` `{"status": "done"}` or `{"status": "refused", "text":
  ...}`, with the same `request` number.

### Files (S4)

- `packages/loop_link/src/loop_link/protocol.py` - Change. `COMMAND_PROTOCOL_VERSION`.
- `packages/loop_link/src/loop_link/encode.py` - Change. `commands_list_frame`, `command_frame`.
- `packages/loop_link/src/loop_link/commands.py` - Create. The names of "Commands module".
- `packages/loop_link/src/loop_link/viewer.py` - Change. The new members and the `colour`
  argument of `send`.
- `packages/loop_link/src/loop_link/__init__.py` - Change. Export `ColourSettings`.
- `packages/loop_link/src/loop_link/testing/fake_viewer.py` - Change. The commands endpoint.
- `packages/loop_link/tests/test_commands.py` - Create. The frame builders (header
  equal to the two fixtures when the arguments are the same), `check_*` functions,
  `raise_for_command_result` for each row of the table, `parse_colour_map_names` (with
  names, without names, malformed), `ColourSettings`.
- `packages/loop_link/tests/test_filter.py` - Create. Story 6 scenarios 1 to 4 and 6:
  range, values with exclude, off, on and remove (data does not change: the stored object
  stays the same), a missing array, a text array, `minimum > maximum`, a new revision keeps
  the filters, `add` refused and the call goes on, `cancel_filter`; the viewer without
  window; an old viewer (`protocol_minor=10`) gives `MissingFeatureError`; the tool part of
  FR-014 (program A holds the data connection; program B filters and colours). S4d adds the
  tests of the wait.
- `packages/loop_link/tests/test_colour.py` - Create. Story 7 scenarios 1 to 7: property,
  no property, map, range, clear, scale, plain colour, a category and a boolean property
  (`InvalidInputError` with the viewer text), an unknown map name (the error lists the
  names, and the fake received no `colour.set_map`), an older viewer without names,
  `colour_map_names`, `send(..., colour=...)` (the object is in the fake before the
  commands, the commands arrive in the order of `apply_colour`, a failed command leaves the
  object in the viewer), the headless viewer for `set_colour_scale`. SC-007: with an object
  of 1,000,000 cells, each of the commands of the table gets its reply in less than 250 ms.
- `packages/loop_link/tests/test_fake_viewer.py` - Change. One test for each rule of the
  fake commands endpoint.
- `packages/loop_link/tests/test_codec_real.py` - Change. Add `commands.list` and
  `command.run` for `filter.edit` (each kind) and each colour command; the real codec
  accepts them.
- `packages/loop_link/tests/test_fixtures.py` - Change. When `LOOP_LINK_GEOVIEWER_FIXTURES`
  is set, read `<variable>/../tool/{commands.list,command.run}.json` (OP11). Build the
  frames with the `request`, `command` and `input` of the fixtures, and assert that the
  headers are equal. Skip with a reason when the folder is missing.
- `packages/loop_link/tests/manual_checks.md` - Change. Add M4 (see "Test strategy of S4").
- `packages/loop_link/README.md` - Change. A section "Filter and colour", about 10 lines.
- `ROADMAP.md` - Change (last task of S4). Update the Stage 7 entry: S4 done.

### Stages and rollback of S4

| Stage | Files and behaviour | Acceptance check | Rollback target and compatibility check |
|---|---|---|---|
| S4a | `protocol.py`, `encode.py`, `commands.py`, fake viewer (commands endpoint), `test_commands.py`, `test_fake_viewer.py` | `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests`; `ruff check LoopStructural packages`; `black --check packages/loop_link` | Revert S4a. Nothing else uses the new names. |
| S4b | `viewer.py` (filter and colour members, no wait, no `send` argument), `__init__.py`, `test_filter.py`, `test_colour.py` | The same commands; FR-018, FR-020, FR-014 (tool part), SC-007 | Revert S4b. S4a stays. |
| S4c | `viewer.py` (`colour` argument of `send`, `apply_colour` use), tests in `test_colour.py` | The same commands; FR-021 | Revert S4c. The default `colour=None` keeps the old `send`. |
| S4d | `commands.py` (`FilterResult`), `viewer.py` (the wait and the result of the filter calls), `tool.py` (`post`), `testing/fake_viewer.py` (long `filter.edit`), tests, M4 step 7, `README.md` | The same commands; FR-019 | Revert S4d. The filter calls return `None` and do not wait. |
| S4e | `test_codec_real.py`, `test_fixtures.py`, `manual_checks.md` (M4), `README.md`, `ROADMAP.md` | The same commands with `LOOP_LINK_GEOVIEWER_REPO` and `LOOP_LINK_GEOVIEWER_FIXTURES` set; M4 by hand | Revert S4e. |

### Test strategy of S4

- Unit and regression tests: the files above. The fake viewer applies the rules above, as the
  S5 fake applies the rules of spec 067.
- Integration and compatibility: the real geoviewer codec checks the frames
  (`test_codec_real.py`, needs `LOOP_LINK_GEOVIEWER_REPO`). The fixtures check the two
  headers that exist (`test_fixtures.py`).
- Manual check M4 (real window viewer of protocol 1.18). Use a fresh scene (`Viewer.clear`)
  and place test objects near the origin of the scene (M3 found that far objects give wrong
  pictures):
  1. Send a `BlockModel` with a number property `grade`. Call `filter_range` for `grade`
     from 0.5 to 1.0. Check that the viewer shows only the cells in the range.
  2. Send an object with a category property. Call `filter_values` with two codes, then with
     `exclude=True`. Check the picture each time.
  3. Call `filter_enable(False)` and then `filter_remove`. Check that all cells show again.
  4. Call `filter_range` for a text property and for a property that does not exist. Check
     the `InvalidInputError` and its text.
  5. Call `set_colour_property`, `set_colour_map` (`coolwarm`), `set_colour_range`,
     `set_colour_scale("log")` and `set_object_colour`. Check the picture after each call.
  6. Call the colour map, range and scale methods for a category property. Check the error.
     Call `colour_map_names()` and check that it gives the names of the viewer.
  7. Send the same id again (a new revision). Check that the filter and the colours stay.
     Start a filter that takes time (an object with 1,000,000 cells), call `filter_range`
     with the default `wait=True`, and check that it returns after the build with
     `FilterResult.passing` and `total` (viewer 1.20). Call it with `wait=False` and call
     `cancel_filter` during the build.
  8. Send an object of 1,000,000 cells with `colour=ColourSettings(...)`. Check that the
     colour settings apply (OP19: the 0.5 s of `flush` can be too short).
  9. Check the questions "assumed" in the facts: does `add` for an existing filter get
     `refused` (OP16), and is a new filter on after `add` and `set_range`.
  M4 records the time of one `set_range` command on the object of 1,000,000 cells.
- Backend coverage and cross-language fixtures: as in the main "Test strategy".
- Performance: SC-007 is a test with the fake viewer (a limit of 250 ms per command).
- Merge gates: as in the main "Test strategy".

### Constitution check for S4

All rows of the table above hold for S4. The differences: the new files are under
`packages/loop_link` only; docstrings and type hints follow the same rules; the ROADMAP
entry is updated in the last task; each stage is a separate revert.

## Test strategy

- **Unit and regression tests** (`packages/loop_link/tests/`):
  - `test_protocol.py`: `encode_frame` then `decode_frame` for every dtype; a frame
    with a wrong size, a header over 1 MiB, a frame over 256 MiB each give
    `ProtocolError`/`ObjectTooLargeError`. Expected values come from the frames
    in `fixtures/protocol/expected/*.json` of the geoviewer repository when
    `LOOP_LINK_GEOVIEWER_FIXTURES` is set (OP8), else from literal bytes in the test.
  - `test_orientation.py`: a table of planes (`dip_direction`, `dip`) turned into
    normals by the formula `(sin d sin a, sin d cos a, cos d)` written in the test,
    then converted back with the library; the same for lineations; the edge cases
    (vertical plane, horizontal plane, vertical line, zero vector, NaN).
    Tolerance `1e-9` degrees.
  - `test_objects.py`: each validity rule, each of the seven types, the property
    type table (all rows), the refusal of an `object` array and of a 2-D property.
  - `test_encode.py`: header fields and buffers of the seven types; `f32` below
    4000 from the origin and `f64` above; the origin shift in `f64`;
    `ObjectTooLargeError` for 2 GiB + 1 declared bytes (a test that builds the
    header sizes with `numpy.broadcast_to`, with no real allocation).
  - `test_discovery.py`: token order (argument, environment, `link.json`), a
    damaged `link.json`, instance files (window before headless, focus order, dead
    port dropped, damaged file ignored), no viewer, default port.
  - `test_viewer.py` (with `FakeViewer`): story 1 scenarios 1 to 7 (the two-program
    scenario uses two `Viewer` objects on one fake viewer), revision rules, `busy`
    with the first program unchanged, `release_data` then send from the second
    program, the keepalive (the fake closes after `heartbeat_timeout=0.6`, the
    viewer sends with `heartbeat_interval=0.1` and stays up for 2 s), the viewer
    that closes in the middle of a send gives `ViewerClosedError`, the wrong token
    gives `TokenError`, the `unsupported_version` error.
  - `test_import.py`: SC-005 (`sys.modules` check in a subprocess run with
    `python -I -c`).
- **Integration and compatibility tests**: `test_fixtures.py` reads every
  `*.bin` in `LOOP_LINK_GEOVIEWER_FIXTURES/` (`fixtures/protocol/python` of the
  geoviewer repository) with `FakeViewer.accept_frame`, and fails for a frame
  that the fake refuses (SC-004). It also encodes the plugin test objects and
  compares header and buffers to the `.bin` frames. It skips with a reason when
  the variable is missing (OP8). No public API of `GeologicalModel` is touched.
- **Round-trip tests**: S1 to S3 compare the arrays that the fake viewer stored with
  the input (`test_viewer.py`). S5 adds `test_roundtrip.py` (send, then fetch).
- **Backend coverage**: the library has one backend (stdlib sockets). S2 tests need
  pyvista, `loop_common` and pandas; each test file uses `pytest.importorskip`,
  and the CI `tests` extra installs them (`tests = ["pytest", "pyvista",
  "loop-common", "pandas"]`), so they do not skip in CI.
- **Cross-language fixtures**: the frames of `fixtures/protocol/python` and
  `expected` (OP8). The Rust side is in the geoviewer repository.
- **Grid and display checks**: `BlockModel` node and cell counts and order. The
  test builds an analytic property `f = x + 10 y + 100 z` on the nodes of a
  3 x 4 x 5 grid in x-fastest order, sends it, and reads the stored buffer from
  the fake viewer. The LoopStructural adapter test repeats this with a
  `loop_common.geometry.StructuredGrid` (x varies fastest in this class, so no
  reorder is made; the test proves it).
- **Performance check**: SC-001. `tests/test_performance.py` is marked
  `@pytest.mark.slow` (the marker is registered in `packages/loop_link/pyproject.toml`;
  the default run does not exclude it). It sends a
  `BlockModel` of 100 x 100 x 100 cells with one `float32` cell property to the
  fake viewer, and asserts that the inventory has the object in under 5 s. Manual
  check M1 repeats it with the real headless viewer
  (`geoviewer-headless --port 0`) and records the time.
- **Merge gates**: `pytest` and `ruff check LoopStructural packages` must pass. The
  package gate is `uv run --package loop-link pytest packages/loop_link/tests`
  (the line of `.github/workflows/packages.yml`, which S1 extends). The root
  `pytest` does not run the package tests. Black: `black --check packages/loop_link`.

## Project structure

### New package `packages/loop_link`

**Responsibility**: the Python client of the geoviewer protocol: codec, transport,
discovery, object model, `Viewer`, adapters, plugin push API, fake viewer.

**Files** (S1 unless the text says S2 or S3):

- `packages/loop_link/pyproject.toml` - Create. Name `loop-link`, version `0.1.0`,
  `requires-python = ">=3.10"`, `dependencies = ["numpy"]`, extras `pyvista`,
  `loopstructural = ["loop-common"]`, `pandas`, `tests`, setuptools
  `where = ["src"]`, `include = ["loop_link", "loop_link.*"]`, pytest marker
  `slow`.
- `packages/loop_link/README.md` - Create. Install, the token, a 10-line example.
- `packages/loop_link/LICENSE` - Create. Copy of `packages/loop_common/LICENSE`.
- `packages/loop_link/ruff.toml` - Create. `extend = "../../pyproject.toml"` and
  `target-version = "py310"`.
- `packages/loop_link/src/loop_link/__init__.py` - Create. Exports (see "Functions
  and public API").
- `packages/loop_link/src/loop_link/errors.py` - Create. The error classes and `Refusal`.
- `packages/loop_link/src/loop_link/protocol.py` - Create. The codec.
- `packages/loop_link/src/loop_link/_websocket.py` - Create. Port of the plugin client.
- `packages/loop_link/src/loop_link/discovery.py` - Create. `find_viewer`, `settings_dir`.
- `packages/loop_link/src/loop_link/objects.py` - Create. The dataclasses and `Categorical`.
- `packages/loop_link/src/loop_link/orientation.py` - Create. The four conversion functions.
- `packages/loop_link/src/loop_link/encode.py` - Create. `encode_object` and the frame builders.
- `packages/loop_link/src/loop_link/viewer.py` - Create. `Viewer`, `InventoryEntry`, keepalive thread.
- `packages/loop_link/src/loop_link/testing/__init__.py` - Create. Exports `FakeViewer`.
- `packages/loop_link/src/loop_link/testing/fake_viewer.py` - Create (FR-015).
  `FakeViewer(token="secret", crs="EPSG:28350", origin=(0.0, 0.0, 0.0), kind="window",
  heartbeat_timeout=15.0, reply_errors=None, read_delay=0.0)` with `.port`,
  `.objects` (id to `StoredObject(header, buffers)`), `.received` (list of
  `Frame`), `.inventory()`, `.write_instance_file(settings: Path) -> Path`,
  `.accept_frame(data: bytes) -> list[Frame]` (apply one frame, return the
  replies), `.close()`, and the context manager. It listens on `127.0.0.1` port 0
  in a thread; it accepts any number of connections; it checks the bearer token
  (HTTP 401 if wrong); it answers `hello` with `hello` and `scene.inventory`; it
  gives `busy` to a second data connection and closes it; it closes a data
  connection that sends no frame for `heartbeat_timeout` seconds; it applies the
  protocol rules of `scene.upsert` (roles, index range, array length, finite
  coordinates, object size limit) and answers `invalid_object` with `object` when
  a rule fails; `reply_errors` maps an id to a forced error text. It accepts the
  path `/tool` and answers every tool message with `unknown_message` in S1.
- `packages/loop_link/src/loop_link/adapters/__init__.py` - Create (S2).
  `convert(obj, **hints) -> SceneObject`: return `obj` when it is a `SceneObject`;
  else choose the adapter by `type(obj).__module__` prefix (`pyvista`,
  `loop_common`), without an import of the extra; else `InvalidInputError` with
  the text "Accepted inputs: loop_link objects; pyvista meshes and grids;
  numpy arrays in a loop_link object; LoopStructural Surface and StructuredGrid;
  pandas tables for drillholes (use loop_link.adapters.drillholes_from_pandas)".
- `packages/loop_link/src/loop_link/adapters/pyvista.py` - Create (S2).
  `from_pyvista(mesh, *, kind: Optional[str] = None, vectors: Optional[str] = None, vectors_kind: str = "normal") -> SceneObject`.
  Rules: `PolyData` with only vertex cells gives `Points` (`Vectors` when `vectors`
  names a point array of shape (n, 3), converted by `Vectors.from_directions`);
  `PolyData` with `lines` gives `Lines`; `PolyData` with faces gives `Surface`
  (`triangulate()` first when a face is not a triangle);
  `ImageData` or a `RectilinearGrid` with equal spacing on each axis gives
  `BlockModel`; `UnstructuredGrid` whose cells are all type `VTK_TETRA` (10) gives
  `TetraMesh`; any other mesh, and `kind="Drillhole"`, give `InvalidInputError`
  (OP3, OP9). `point_data` goes to `properties`, `cell_data` to `cell_properties`;
  a 2-D array is refused. `kind` overrides the rule for `Points` versus `Vectors`
  only.
- `packages/loop_link/src/loop_link/adapters/loopstructural.py` - Create (S2).
  `from_loopstructural(obj) -> SceneObject`: `loop_common.geometry.Surface` gives
  `Surface(vertices, triangles, properties=obj.properties or {}, cell_properties=obj.cell_properties or {})`;
  `loop_common.geometry.StructuredGrid` gives `BlockModel(origin=obj.origin,
  cell_sizes=obj.step_vector, cells=obj.nsteps_cells, properties=obj.properties,
  cell_properties=obj.cell_properties)`. The adapter finds the class with
  `sys.modules.get("loop_common.geometry")`, so it never imports the module when the
  caller has no LoopStructural object.
- `packages/loop_link/src/loop_link/adapters/pandas.py` - Create (S2).
  `drillholes_from_pandas(collar, survey, logs=None, *, hole_id="hole_id", x="x", y="y", z="z", depth="depth", azimuth="azimuth", dip="dip", from_="from", to="to", dip_positive_down=True) -> dict[str, Drillhole]`.
  It groups each table by `hole_id`; the result maps a hole id to a `Drillhole`;
  extra columns of `logs` become `interval_properties` (a text column of strings
  becomes a numpy `U` array; any other `object` column gives `InvalidInputError`);
  `dip_positive_down=False` changes the sign of `dip` (the `geoh5` convention). A
  hole in `survey` or `logs` that is not in `collar` gives `InvalidInputError`.
  `Viewer.send` accepts one `Drillhole`. The user sends the dict with a loop of
  `viewer.send(hole, id=hole_id)`.
- `packages/loop_link/src/loop_link/push.py` - Create (S3). `PushResult`, `push_objects`.
- `packages/loop_link/tests/conftest.py` - Create. Fixture `fake_viewer` (starts
  `FakeViewer`, writes `GEOVIEWER_TOKEN` with `monkeypatch`, sets `settings_dir` to
  `tmp_path`), and the helper `geoviewer_fixtures()` (OP8).
- `packages/loop_link/tests/test_protocol.py`, `test_orientation.py`, `test_objects.py`,
  `test_encode.py`, `test_discovery.py`, `test_viewer.py`, `test_import.py`,
  `test_performance.py`, `test_fixtures.py` - Create (S1).
- `packages/loop_link/tests/test_adapter_pyvista.py`,
  `test_adapter_loopstructural.py`, `test_adapter_pandas.py` - Create (S2).
- `packages/loop_link/tests/test_push.py` - Create (S3). Scenarios of story 5:
  `push_objects` frames equal `Viewer.send` frames; `replace_prefix` removes old
  ids; one bad object goes to `skipped`; `cancel` stops after the current object;
  a `FakeViewer(read_delay=0.05)` push in a worker thread leaves the main thread
  free (the main thread counts to 20 ticks of 0.01 s while the push runs).
- `packages/loop_link/tests/manual_checks.md` - Create (S1). M1: performance with the
  real headless viewer; M2: a send of the seven types, then a visual check in the
  window viewer.

### Repository registration (S1)

**Responsibility**: FR-016, SC-006, and the repository settings for the new package.

**Files**:

- `ROADMAP.md` - Change. Add a stage entry for `packages/loop_link` after the
  last entry of "Stage sequence" (the task reads the end of the list first and
  uses the next stage number), and add one status log entry. Both entries name
  `loop_link` (SC-006, first entry).
- `.release-please-config.json` - Change. Add `"packages/loop_link": {"release-type":
  "python", "component": "loop-link"}` and add `"packages/loop_link"` to `exclude-paths` of
  the root entry (SC-006, second entry; the fix of 2026-09-27 keeps sub-package
  commits out of the root bump).
- `.release-please-manifest.json` - Change. Add `"packages/loop_link": "0.1.0"`.
- `pyproject.toml` - Change. Add `loop-link = { workspace = true }` to
  `[tool.uv.sources]` (or the `exclude` entry of OP5); add
  `"packages/loop_link/tests/*" = ["D", "ANN"]` to the Ruff `per-file-ignores`. No
  entry for `packages/loop_link/src`.
- `uv.lock` - Change. Run `uv lock` (OP5).
- `.github/workflows/packages.yml` - Change. Add `{ name: loop-link, dir: loop_link }`
  to the matrix, and add the condition `if: ${{ !(matrix.package.name == 'loop-link' &&
  startsWith(matrix.python-version, '3.9')) }}` to the two steps that install and test.
  The paths filter `packages/**` already covers the new folder.

## Requirement coverage

| Requirement | Design and files | Verification | Gap |
|---|---|---|---|
| FR-001 | "Functions and public API", `discovery.py` | `test_discovery.py`; story 1 scenario 1 in `test_viewer.py` | None |
| FR-002 | `objects.py`, `encode.py` (type table) | `test_encode.py`, `test_viewer.py` (seven types) | None |
| FR-003 | `Viewer.send` revision rule | `test_viewer.py` scenario 2 | None |
| FR-004 | `Viewer.select`, `encode.selection_set_frame` (S5a) | `test_selection.py`, `test_codec_real.py`, `test_fixtures.py` | None (OP12, OP13 are viewer facts that M3 checks) |
| FR-005 | `Viewer.selection`, `Viewer.watch_selection`, `SelectionWatcher` (S5a) | `test_selection.py` (poll once for a change, interval, no thread before the call) | None |
| FR-006 | `Viewer.fetch`, `fetch.py` (S5b) | `test_fetch.py`, `test_roundtrip.py` | OP10 (a drillhole with several log sets loses the extra sets, with a warning), OP14 |
| FR-007 | `import loop_link` has no Qt; `push.py` has no Qt | `test_import.py` | None |
| FR-008 | `errors.py`; each edge case maps to one class (token: `TokenError`; no viewer: `NoViewerError`; too large: `ObjectTooLargeError`; closed: `ViewerClosedError`; bad property: `InvalidInputError`; busy: `BusyError`; older viewer: `MissingFeatureError`; origin differs: handled by conversion, not an error) | `test_viewer.py`, `test_encode.py`, `test_discovery.py` | None for S1 to S3. The older-viewer case for select and read is S5 (OP7). |
| FR-009 | `push.py` (`push_objects`, `PushResult`) | `test_push.py` | None. The plugin change is in the plugin repository. |
| FR-010 | `adapters/` | `test_adapter_*.py` | OP3 (pyvista to Drillhole); OP9 (non-regular grids) |
| FR-011 | `orientation.py`, `Vectors` | `test_orientation.py`; story 1 scenario 6 | None |
| FR-012 | One codec in `protocol.py` and `encode.py`; core imports numpy and the standard library only; adapters are extras | `test_import.py`; `test_protocol.py` | None |
| FR-013 | `unknown_message` of the tool connection gives `MissingFeatureError` (`raise_for_request_error`) | `test_tool.py`, `test_selection.py` (`protocol_minor=14`) | None |
| FR-014 | `Viewer.release_data`, `BusyError`, keepalive thread; the tool connection (S5) for select, read and fetch while another program holds the data connection. Filter and colour use the same tool connection (S4b). | `test_viewer.py` scenario 7 (send part); `test_selection.py` (tool part); `test_filter.py` (filter and colour while another program holds the data connection) | None |
| FR-015 | `testing/fake_viewer.py` | `test_viewer.py`, `test_fixtures.py` | None |
| FR-016 | "Repository registration" | `git grep -n loop_link ROADMAP.md .release-please-config.json` gives 2 entries or more | None |
| FR-017 | `pyproject.toml` (`>=3.10`), `ruff.toml` (py310) | `ruff check` and CI on 3.10 to 3.12 | OP4 (QGIS Python), OP5 (uv) |
| FR-018 | `Viewer.filter_range`, `filter_values`, `filter_enable`, `filter_remove`, `commands.run_command` (S4b) | `test_filter.py`, `test_codec_real.py`, `test_fixtures.py`; M4 steps 1 to 4 | OP16, OP17, OP18, OP20 |
| FR-019 | `Viewer.cancel_filter` (S4b); the `wait`, `cancel` and `progress` arguments (S4d) | `test_filter.py`; M4 step 7 | OP15 answered (spec 072) |
| FR-020 | `Viewer.set_colour_property`, `set_colour_map`, `set_colour_range`, `clear_colour_range`, `set_colour_scale`, `set_object_colour`, `colour_map_names` (S4b) | `test_colour.py`, `test_commands.py`; M4 steps 5 and 6 | OP18 |
| FR-021 | `ColourSettings`, `Viewer.apply_colour`, the `colour` argument of `Viewer.send` (S4c) | `test_colour.py`; M4 step 8 | OP19 |
| SC-001 | `test_performance.py`; M1 | assert under 5 s | None |
| SC-002 | `fetch.py`, `Viewer.fetch` | `test_roundtrip.py` (rtol 1e-6, seven types) | None |
| SC-003 | `encode.selection_set_frame` (one frame, one `u32` buffer) | `test_selection.py` (100,000 cells in one frame) | None |
| SC-004 | `test_fixtures.py` | frames of `fixtures/protocol/python`; the geoviewer tests read the frames of the library in the geoviewer repository | OP8 |
| SC-005 | `test_import.py`; CI matrix | `sys.modules` check | None |
| SC-006 | "Repository registration" | `git grep` of two files | None |
| SC-007 | `run_command` on the tool connection (S4b) | `test_colour.py` (each command under 250 ms with an object of 1,000,000 cells in the fake) | None |
| Story 1 scenarios 1 to 8 | `Viewer`, `discovery.py`, `test_viewer.py`, `test_import.py` | as above | None |
| Story 2 scenarios 1 to 4 | `Viewer.select`, errors table | `test_selection.py` | None |
| Story 3 scenarios 1 to 3 | `Viewer.selection`, `watch_selection`, `parse_selection_result` | `test_selection.py` | None |
| Story 4 scenarios 1 to 3 | `Viewer.fetch`, `decode_object` | `test_fetch.py`, `test_roundtrip.py` | None |
| Story 6 scenarios 1 to 6 | `Viewer.filter_*`, `cancel_filter`, the fake `filters` | `test_filter.py` (scenario 5 in S4d); M4 | OP15, OP16, OP17, OP18, OP20 |
| Story 7 scenarios 1 to 7 | `Viewer.set_colour_*`, `apply_colour`, `send(colour=)` | `test_colour.py`; M4 | OP18, OP19 |
| Story 5 scenarios 1 and 2 | `push.py`, `test_push.py` | frame comparison with plugin frames (OP8); worker thread test | OP8 |

## Handover

Other specs and repositories use these contracts.

- **Python API**: the names in "Functions and public API", `loop_link.objects`,
  `loop_link.errors`, and `loop_link.testing.FakeViewer`. The `0.x` version means
  a breaking change can come in a minor version, and release-please records it.
- **Plugin repository**: `push_objects(viewer, objects, replace_prefix=, progress=,
  cancel=)` replaces `GeoViewerLink.push`. `PushResult.summary()` keeps the text of
  the old `PushResult.summary()`. The plugin builds `objects` as `(id, SceneObject)`
  pairs and runs the call on its worker thread. The plugin owns the change and the
  removal of the old client (FR-009).
- **S4** (filter and colour) gives: `Viewer.filter_range`, `filter_values`, `filter_enable`,
  `filter_remove`, `cancel_filter`, `colour_map_names`, `set_colour_property`,
  `set_colour_map`, `set_colour_range`, `clear_colour_range`, `set_colour_scale`,
  `set_object_colour`, `apply_colour`, the `colour` argument of `Viewer.send`,
  `loop_link.commands.ColourSettings`, and the fake viewer members `filters` and
  `colouring`. It uses the tool messages `commands.list`, `command.run`, `command.result`
  and the commands of "Facts that S4 uses". The result of a viewer of 1.20 has the count of
  elements that pass (OP20, answered).
- **S5** (select, read selection, fetch) gives: `Viewer.select`, `Viewer.selection`,
  `Viewer.watch_selection`, `Viewer.fetch`, `Selection`, `SelectionWatcher`,
  `FetchResult`, the errors `StaleRevisionError`, `CancelledError` and
  `UnsupportedObjectError`, and the fake viewer members `selection`,
  `set_user_selection`, `tool_requests`. S4 (filter and colour) reuses
  `loop_link.tool.ToolConnection`: `request(build, until=("command.result",),
  feature=..., object_id=...)` with a frame builder for `command.run`, and
  `raise_for_request_error` for the failures. The plugin can use `select`,
  `selection` and `watch_selection` on a worker thread (the tool connection takes
  one request at a time).
- **geoviewer repository**: the frames that `encode_object` writes must pass the
  geoviewer protocol tests. The coordinator copies the frames to
  `fixtures/protocol/python` by hand when the fixtures change (OP8).
- **Shared files**: `ROADMAP.md`, `pyproject.toml`, `uv.lock`,
  `.release-please-config.json`, `.release-please-manifest.json`,
  `.github/workflows/packages.yml`. Specs 003 and 005 also work in these
  repository areas. S1 changes only the lines named above. A merge conflict in
  `uv.lock` is solved by `uv lock` again.

## Complexity tracking

None. Every applicable constitution check passes. The py310 target of the package
is a decision that FR-017 asks of this plan, and it does not break a project rule
(the repository target py39 stays for all other code).
- 2026-11-09 (orchestrator, T089): the fake gives `passing` and `total` for `set_on` always, and for `set_range`, `set_values` and `set_exclude` only when the filter is on after the edit (else plain `{"status": "done"}`). This is a fake assumption; the real viewer rule is checked in M4.
- 2026-11-09 (orchestrator, T081): the fake viewer keeps a non-cancel frame that arrives during a long build and handles it after the build (as a viewer that does one command at a time). This changed `testing/fake_viewer.py` in the T081 round (a file of T089).
