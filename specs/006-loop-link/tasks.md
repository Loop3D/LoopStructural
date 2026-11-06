# Tasks: loop_link, a Python link to geoviewer

**Feature**: `006-loop-link` | **Plan**: [plan.md](plan.md)

<!--
Scope: stages S1 to S3 of the plan (stories 1 and 5, built, T001 to T038) and stage S5
(stories 2 to 4, T039 to T063). Stage S4 (stories 6, 7; OP6) gets tasks after its own planner run.
Package test command: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/<file>`.
The root `pytest` does not run these tests (testpaths = ["tests"]).
OP5 RESULT (2026-11-06): the uv workspace fallback was used (`packages/loop_link` is in the
workspace `exclude` list; `uv.lock` is unchanged). Test command now: `packages/loop_link/.venv/Scripts/python.exe -m pytest <file>`.
Make the venv once: `uv venv --python 3.12 packages/loop_link/.venv` then `uv pip install --python packages/loop_link/.venv/Scripts/python.exe -e "packages/loop_link[tests]" -e packages/loop_common`.
`uv lock --check` fails on the unchanged master with the local uv 0.11.8 (lock revision 5), so it is not a gate here.
If OP5 forces the pip fallback, use `pytest packages/loop_link/tests/<file>` after
`pip install -e packages/loop_link[tests]`.
Common rules for every Python task: `from __future__ import annotations` first;
numpy-style docstrings; type hints on all public names; no `match` statement; line
length 100; keep string quotes as written (Black `skip-string-normalization`).
-->

## Orchestrator notes

- Run this file with the spec-orchestrator skill. It has the rules for the
  worktree, the sub agents, the checks and the commits.
- Waves: 19. Tiers: [S] = simple model, [M] = standard model.
- Do not start a task with "Blocked by".

## Phase 1: Setup

- [x] T001 [Setup] [M] Create the package skeleton `packages/loop_link`.
  - Files: `packages/loop_link/pyproject.toml`, `packages/loop_link/README.md`, `packages/loop_link/LICENSE`, `packages/loop_link/ruff.toml`, `packages/loop_link/src/loop_link/__init__.py`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import tomllib,pathlib;tomllib.loads(pathlib.Path('packages/loop_link/pyproject.toml').read_text())"` passes and `ruff check packages/loop_link` passes.
  - Context: Model on `packages/loop_common/pyproject.toml` (setuptools, `where = ["src"]`). `[project]` name `loop-link`, version `0.1.0`, `requires-python = ">=3.10"`, `dependencies = ["numpy"]`. Extras: `pyvista = ["pyvista"]`, `loopstructural = ["loop-common"]`, `pandas = ["pandas"]`, `tests = ["pytest", "pyvista", "loop-common", "pandas"]`. Setuptools `include = ["loop_link", "loop_link.*"]`. Register the pytest marker `slow` under `[tool.pytest.ini_options]`. `ruff.toml` has exactly `extend = "../../pyproject.toml"` and `target-version = "py310"`. `LICENSE` is a copy of `packages/loop_common/LICENSE`. `README.md`: install, the token (`GEOVIEWER_TOKEN` or the saved token of the viewer), a 10-line example with `Viewer`, `Points`. `__init__.py` holds a docstring and `__version__ = "0.1.0"` only (T023 adds the exports).

- [x] T002 [P] [Setup] [S] Register `loop_link` for release-please.
  - Files: `.release-please-config.json`, `.release-please-manifest.json`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import json;c=json.load(open('.release-please-config.json'));m=json.load(open('.release-please-manifest.json'));assert 'packages/loop_link' in c['packages'] and 'packages/loop_link' in c['packages']['.']['exclude-paths'] and m['packages/loop_link']=='0.1.0'"` passes.
  - Context: In the config, add the key `"packages/loop_link"` under `packages` with the value `{"release-type": "python", "component": "loop-link"}`. Also append `"packages/loop_link"` to the `exclude-paths` list of the `"."` entry (keep the list sorted by path as it is). In the manifest, add `"packages/loop_link": "0.1.0"`. Keep the 2-space JSON indent.

- [x] T003 [Setup] [M] Add `loop-link` to the uv workspace and to the Ruff config (OP5).
  - Files: `pyproject.toml`, `uv.lock`
  - Needs: T001
  - Wave: 2
  - Done when: `uv lock --check` passes and `ruff check LoopStructural packages` passes.
  - Context: In root `pyproject.toml` add `loop-link = { workspace = true }` to `[tool.uv.sources]`. Add `"packages/loop_link/tests/*" = ["D", "ANN"]` to `[tool.ruff.lint.per-file-ignores]` next to the other `packages/*/tests/*` lines. Do NOT add an entry for `packages/loop_link/src`. Run `uv lock`. OP5 fallback: if uv refuses the workspace, or the lock drops Python 3.9 for the other members, remove the source line, add `"packages/loop_link"` to `[tool.uv.workspace] exclude`, run `uv lock` again, and write one line in the commit message that says the fallback was used. Do not change another dependency in `uv.lock`.

- [x] T004 [P] [Setup] [S] Add `loop-link` to the package CI matrix.
  - Files: `.github/workflows/packages.yml`
  - Needs: none
  - Wave: 1
  - Done when: `python -c "import yaml;yaml.safe_load(open('.github/workflows/packages.yml'))"` passes and `git grep -n "loop-link" .github/workflows/packages.yml` shows 3 lines.
  - Context: In `jobs.packages-test.strategy.matrix.package` add the line `- { name: loop-link, dir: loop_link }`. Add this line to the two steps that install and test (`Install package with test extras` and `pytest`): `if: ${{ !(matrix.package.name == 'loop-link' && startsWith(matrix.python-version, '3.9')) }}`. The package needs Python 3.10 or later. Do not change another step.

- [x] T005 [P] [Setup] [M] Add `loop_link` to `ROADMAP.md`.
  - Files: `ROADMAP.md`
  - Needs: none
  - Wave: 1
  - Done when: `git grep -n -c loop_link ROADMAP.md` gives 2 or more.
  - Context: Read the end of the section "Stage sequence" and the end of "Status log" first. Add one unchecked stage entry for `packages/loop_link` (spec `specs/006-loop-link`) with the next free stage number: the Python link to geoviewer (send, adapters, plugin push API; select, fetch, filter and colour wait for geoviewer spec 066). Add one status log entry dated 2026-11-06 that says the spec and plan exist and the package is new. Use the style of the neighbouring entries. Principle: `ROADMAP.md` owns the stages (constitution).

## Phase 2: Foundation (blocks all user stories)

- [x] T006 [Found] [S] Errors module with its test.
  - Files: `packages/loop_link/src/loop_link/errors.py`, `packages/loop_link/tests/test_errors.py`
  - Needs: T001
  - Wave: 2
  - Done when: write the test first and see it fail; then `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_errors.py` passes.
  - Context: Classes, all subclasses of `LoopLinkError(Exception)`: `NoViewerError`, `TokenError`, `BusyError`, `InvalidInputError`, `ViewerClosedError`, plus `MissingFeatureError(feature: str)` (attribute `feature`; text names the feature), `ObjectTooLargeError(object_id: str, size: int, limit: int)` (attributes with the same names; text has all three), `ViewerRefusedError(refusals: list[Refusal])` (attribute `refusals`), `MissingExtraError(extra: str)` (attribute `extra`; text says `pip install loop-link[<extra>]`). `@dataclass(frozen=True) class Refusal: object_id: str; code: str; message: str`. Each default text says what to do (FR-008). Test: each class is a `LoopLinkError`; the attributes are kept; `str(error)` is not empty; `Refusal` is hashable.

- [x] T007 [P] [Found] [M] Test: the codec in `test_protocol.py`.
  - Files: `packages/loop_link/tests/test_protocol.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_protocol.py` fails with ImportError of `loop_link.protocol`, then passes after T008.
  - Context: API under test: `encode_frame(header: dict, buffers: Sequence[tuple[str, np.ndarray]] = ()) -> bytes`, `decode_frame(data: bytes) -> Frame` (`Frame(header: dict, buffers: dict[str, np.ndarray])`), `decode_header(data: bytes) -> dict`, `ProtocolError`, constants `MAX_FRAME_BYTES=256*1024*1024`, `MAX_HEADER_BYTES=1024*1024`, `PROTOCOL_VERSION={"major":1,"minor":14}`. Frame layout: `[u32 LE header_len][JSON header][buffers]`, little-endian. `encode_frame` adds `version` and `buffers` (`[{"name","dtype","shape"}]`). dtype names `f32 f64 u32 i32 u8`. Tests: a round trip for each of the 5 dtypes with shape (2,3); a zero-element array is valid; a frame whose length differs from `4+header_len+sum(sizes)` gives `ProtocolError`; a header over 1 MiB gives `ProtocolError` on encode; a frame under 4 bytes gives `ProtocolError`; NaN in the header gives an error (`allow_nan=False`); an unsupported dtype (`complex128`) gives `ProtocolError`. A byte-literal test: `encode_frame({"type":"heartbeat","seq":1})` starts with `struct.pack('<I', len(header_json))` and the JSON has no spaces. If the environment variable `LOOP_LINK_GEOVIEWER_FIXTURES` is set, also decode each `*.bin` file in it and assert `decode_frame` does not raise (skip with a reason when not set).

- [x] T008 [Found] [M] Implement the codec `protocol.py`.
  - Files: `packages/loop_link/src/loop_link/protocol.py`
  - Needs: T007
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_protocol.py` passes and `ruff check packages/loop_link` passes.
  - Context: Port the idea of `encode_frame` and `decode_header` from `C:\Users\lachl\dev\plugin_loopstructural-geoviewer\loopstructural\main\geoviewer\protocol.py` (read it). Public names: `PROTOCOL_VERSION`, `MAX_FRAME_BYTES`, `MAX_HEADER_BYTES`, `MAX_OBJECT_BYTES = 2*1024**3`, `F32_LIMIT = 4000.0`, `Frame` (frozen dataclass), `encode_frame`, `decode_frame`, `decode_header`, `ProtocolError(LoopLinkError)`. The codec raises `ProtocolError` for size problems. A frame over `MAX_FRAME_BYTES` raises `ProtocolError` in `encode_frame` (the object-level `ObjectTooLargeError` comes from `encode.py`, T014). JSON: `json.dumps(full, separators=(",", ":"), allow_nan=False)`. Buffers: `np.ascontiguousarray`, little-endian `tobytes`. `decode_frame` returns arrays with the shape of the header (`np.frombuffer(...).reshape(shape)`, a copy is not needed). Numpy and the standard library only. Read `MAX_FRAME_BYTES` as a module global at call time (a test patches it); do not bind it as a default argument.

- [x] T009 [P] [Found] [M] Test: vector conversion in `test_orientation.py`.
  - Files: `packages/loop_link/tests/test_orientation.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_orientation.py` fails with ImportError, then passes after T010.
  - Context: Under test: `normals_to_orientation(normals) -> (dip_direction, dip)`, `directions_to_lineation(directions) -> (trend, plunge)`, `orientation_to_normals(dip_direction, dip) -> (n,3)`, `lineation_to_directions(trend, plunge) -> (n,3)`. Degrees. x east, y north, z up. Independent expected value in the test: a plane that dips to azimuth `a` with dip `d` has the upward normal `(sin d sin a, sin d cos a, cos d)` (write this formula in the test, not by calling the library). Tests: 12 planes (`a` in 0, 45, 90, 135, 180, 270, 359; `d` in 10, 45, 80) convert to `(a, d)` within 1e-9 degrees; a downward normal gives the same result as its upward twin; a vertical plane (n_z=0) gives dip 90; a horizontal plane gives dip 0 and dip direction 0.0; lineations: `(trend, plunge)` to direction `(cos p sin t, cos p cos t, -sin p)` and back; a direction with z>0 flips (axis); vertical line gives trend 0.0 and plunge 90; the inverse functions round trip; a zero vector and a NaN vector each raise `InvalidInputError`. Dip direction and trend are in [0, 360).

- [x] T010 [Found] [M] Implement `orientation.py`.
  - Files: `packages/loop_link/src/loop_link/orientation.py`
  - Needs: T009
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_orientation.py` passes.
  - Context: Signatures take array-like (n,3) or (3,) input and return numpy arrays: `normals_to_orientation(normals) -> tuple[np.ndarray, np.ndarray]`; `directions_to_lineation(directions) -> tuple[np.ndarray, np.ndarray]`; `orientation_to_normals(dip_direction, dip) -> np.ndarray`; `lineation_to_directions(trend, plunge) -> np.ndarray`. Rules: `n = v/|v|`; for normals, if `n_z < 0` use `-n`; `dip = degrees(arccos(clip(n_z,-1,1)))`; `dip_direction = degrees(arctan2(n_x, n_y)) % 360`; if `n_x == n_y == 0` the dip direction is 0.0. For lineations, if `d_z > 0` use `-d`; `plunge = degrees(arcsin(clip(-d_z,-1,1)))`; `trend = degrees(arctan2(d_x, d_y)) % 360`; if `d_x == d_y == 0` the trend is 0.0. Inverse of normals: `(sin d sin a, sin d cos a, cos d)`. Inverse of lineations: `(cos p sin t, cos p cos t, -sin p)`. A zero-length or non-finite vector raises `InvalidInputError` (from `loop_link.errors`).

- [x] T011 [P] [Found] [M] Test: the object model in `test_objects.py`.
  - Files: `packages/loop_link/tests/test_objects.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_objects.py` fails with ImportError of `loop_link.objects`, then passes after T012.
  - Context: Frozen dataclasses in `loop_link.objects`: `Points(vertices, properties={})`; `Vectors(vertices, orientation=None, lineation=None, properties={})` with `Vectors.from_directions(vertices, directions, *, kind, properties=None)` (`kind` is `"normal"` or `"lineation"`, else `InvalidInputError`); `Lines(vertices, segments, properties, cell_properties)`; `Surface(vertices, triangles, properties, cell_properties)`; `BlockModel(origin, cell_sizes, cells, properties, cell_properties)`; `TetraMesh(vertices, tetrahedra, properties, cell_properties)`; `Drillhole(collar, survey, intervals, end_of_hole=None, properties, interval_properties)`; `Categorical(codes, names, colours=None)`; `DEFAULT_CATEGORY_COLOURS` (10 RGBA tuples). Rules to test (each raises `InvalidInputError` that names the field): vertices not (n,3) or not finite; an index array with a negative value or a value >= n; a property of the wrong length (vertex: n; cell: m; `BlockModel` nodes: (nx+1)(ny+1)(nz+1), cells: nx*ny*nz, x fastest; Drillhole `properties`: k stations, `interval_properties`: j rows); a 2-D property; `cell_sizes` <= 0 or `cells` < 1 or non-integer; `Points` or `Vectors` with `cell_properties` (they have no such field, so the test checks `TypeError`); `Categorical` with a code >= `len(names)`. Valid input of all 7 types is accepted and arrays become numpy arrays.

- [x] T012 [Found] [M] Implement `objects.py`.
  - Files: `packages/loop_link/src/loop_link/objects.py`
  - Needs: T011, T010
  - Wave: 5
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_objects.py` passes.
  - Context: Implement the dataclasses that T011 names, with the same fields and defaults (`field(default_factory=dict)` for the mappings). Each class does its checks in `__post_init__` and converts arrays with `object.__setattr__(self, name, np.asarray(...))` (frozen). Vertex and index arrays: `float64` and integer. `Vectors.from_directions(..., kind="normal")` calls `normals_to_orientation` and `kind="lineation"` calls `directions_to_lineation` (both in `loop_link.orientation`) and stores the pair in `orientation` or `lineation`. `Vectors` accepts both `orientation` and `lineation`. Add `SceneObject = Union[Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole]`, `Property = Union[np.ndarray, Sequence, Categorical]` and `ArrayLike`. `DEFAULT_CATEGORY_COLOURS` = the 10 tab10 colours as RGBA with alpha 255: (31,119,180), (255,127,14), (44,160,44), (214,39,40), (148,103,189), (140,86,75), (227,119,194), (127,127,127), (188,189,34), (23,190,207). A category without a given colour uses the colour at index `i % 10`. Do not import pyvista, pandas or Qt. Keep a `Categorical` property as the same object (the test checks `is`). `orientation` is a tuple of two arrays `(dip_direction, dip)` and `lineation` is `(trend, plunge)`, in degrees. Index arrays must be integer dtype (float indexes raise `InvalidInputError`). A `Categorical` code must be in `[0, min(len(names), 2**32))`.

- [x] T013 [Found] [M] Test: object encoding in `test_encode.py`.
  - Files: `packages/loop_link/tests/test_encode.py`
  - Needs: T012, T008
  - Wave: 6
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_encode.py` fails with ImportError of `loop_link.encode`, then passes after T014.
  - Context: Under test: `encode_object(obj, *, id, name, parent, revision, origin) -> bytes` (decode with `decode_frame`), `hello_frame(crs, origin, version)`, `heartbeat_frame(seq)`, `remove_frame(id)`, `clear_frame()`, `group_frame(id, name, parent, revision)`. Header fields of `scene.upsert`: `type id object_type name parent revision edit([]) replaces(null) geometry params data property_groups drape(null) style({})`. Type map: Points and Vectors -> `Points` (geometry `vertices`); Lines -> `Curve` (`vertices`, `segments` u32); Surface -> `Surface` (`vertices`, `triangles` u32); TetraMesh -> `TetraMesh` (`vertices`, `tetrahedra` u32); BlockModel -> `BlockModel` (no geometry; `params` `origin`, `cell_sizes`, `cells` as JSON lists, origin minus scene origin); Drillhole -> `Drillhole` (`collar` minus origin, `survey` f64 (k,3), `intervals` f64 (j,2); `params.end_of_hole=[v]` only when set). Vertices are f32 when max abs after the shift <= 4000, else f64; the shift is done in f64 (test with origin 500000, 7000000, 0 and vertices near it: f32 buffer, values equal to input minus origin within 1e-3). Data types: float32 -> `float`/f32; float64 -> `float`/f64; bool -> `boolean`/u8; int64 inside i32 -> `integer`/i32; int64 outside i32 -> `InvalidInputError`; `U` strings -> `text` (u8 bytes + u32 offsets named `data/<name>/offsets`, header field `offsets`); `Categorical` -> `referenced`/u32 with `value_map` entries `{"id","name","colour"}`; `object` dtype -> `InvalidInputError`; two properties that give the same buffer name -> `InvalidInputError`. Data array: `id`=`name`=property name, `buffer`=`data/<name>`, `association` `vertex` or `cell` (Drillhole `interval_properties`: `interval`). Vectors: `property_groups` `[{"id":"orientation","name":"orientation","role":"orientation","members":[{"role":"dip_direction","data":"orientation/dip_direction"},{"role":"dip","data":"orientation/dip"}]}]` and the f64 vertex arrays; lineation uses `lineation`, `trend`, `plunge`. Too large: build an object whose declared size exceeds 2 GiB with `np.broadcast_to` views (no real allocation) and expect `ObjectTooLargeError` with the id, size and limit.

- [x] T014 [Found] [M] Implement `encode.py`.
  - Files: `packages/loop_link/src/loop_link/encode.py`
  - Needs: T013
  - Wave: 7
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_encode.py` passes.
  - Context: Public: `encode_object`, `hello_frame(crs: str, origin: Sequence[float], version: str) -> bytes` (header `type hello`, `application {"name":"loop_link","version":version}`, `crs`, `origin` as floats), `heartbeat_frame(seq: int)`, `remove_frame(id)` (`scene.remove`), `clear_frame()` (`scene.clear`), `group_frame(id, name, parent, revision)` (a `scene.upsert` of `object_type` `Group`). Private `_data_array(name, value, association) -> (desc: dict, buffers: list[tuple[str, np.ndarray]])` follows the type table in the test of T013. Use `loop_link.protocol.encode_frame`. Before building the frame, sum the declared bytes of all buffers plus `entries*8+name bytes` of each `value_map`; if the sum is over `MAX_OBJECT_BYTES`, or the frame is over `MAX_FRAME_BYTES` (catch `ProtocolError` from `encode_frame`), raise `ObjectTooLargeError(id, size, limit)`. Use broadcast views without copying when computing sizes (read `array.size * itemsize`, not `tobytes`). Coordinates: reject NaN and infinity with `InvalidInputError` (objects already check this). Principle IV: subtract the origin in float64 before any cast to f32. Check for NaN/inf in `encode_object` itself (the tests build objects that skip the object checks). Compute the declared size of the object from array shapes and itemsize BEFORE any copy or shift of the vertices (a test uses a 2**28 x 3 broadcast view that must not be allocated). The test expects `ObjectTooLargeError.limit == 2*1024**3` and `size > limit`; `geometry` is `{}` for BlockModel; a lineation group has role `lineation` with members `trend`, `plunge` in that order; data entries have `id`, `name`, `buffer`, `association`, `data_type` (plus `offsets` for text and `value_map` for referenced).

- [x] T015 [Found] [M] Port the stdlib WebSocket client to `_websocket.py`.
  - Files: `packages/loop_link/src/loop_link/_websocket.py`
  - Needs: T008
  - Wave: 5
  - Done when: `ruff check packages/loop_link` passes and `packages/loop_link/.venv/Scripts/python.exe -c "import loop_link._websocket"` passes. Behaviour is tested through `test_viewer.py` (T020).
  - Context: Read and port `C:\Users\lachl\dev\plugin_loopstructural-geoviewer\loopstructural\main\geoviewer\websocket_client.py` (245 lines, standard library only). Keep the class `WebSocketClient` with `connect(cls, host, port, *, token="", timeout=5.0, path="/")`, `send_binary(data: bytes)`, `receive(timeout: float) -> Optional[bytes]`, `close()`. Changes: (1) `WebSocketError(LoopLinkError)` and `WebSocketClosed(WebSocketError)`; (2) `MAX_INCOMING_BYTES = MAX_FRAME_BYTES + MAX_HEADER_BYTES` (import from `loop_link.protocol`; the old limit was 16 MiB); (3) in `connect`, an HTTP 401 answer raises `TokenError`, and a refused or timed-out TCP connection raises `NoViewerError` (both from `loop_link.errors`); (4) add numpy-style docstrings and type hints to every public method (Ruff `D` and `ANN` apply to new files); (5) the handshake header stays `Authorization: Bearer <token>`. Keep the thread rule: not thread safe.

- [x] T016 [P] [Found] [M] Test: discovery in `test_discovery.py`.
  - Files: `packages/loop_link/tests/test_discovery.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_discovery.py` fails with ImportError, then passes after T017.
  - Context: Under test: `find_viewer(*, port=None, token=None, host="127.0.0.1", settings=None, timeout=0.5) -> ViewerAddress` (frozen dataclass `host, port, token`), `settings_dir() -> Path` (`%APPDATA%\geoviewer` on Windows, else `~/.geoviewer`), `DEFAULT_PORT = 47615`. Use `tmp_path` as `settings` and `monkeypatch` for `GEOVIEWER_TOKEN`. Token order: the argument; then `GEOVIEWER_TOKEN` (an empty value is not a token); then the field `token` of `settings/link.json`; none gives `TokenError`; a damaged `link.json` (not JSON) gives `TokenError`. Instance files: `settings/instances/<pid>.json` with fields `version (1), pid, port, kind ("window" or "headless"), started, focused` (ms). To make a live port, open a `socket` listener on port 0 in the test. Rules: a file whose port has no listener is dropped; `window` before `headless`; then the largest `focused`; a damaged file is ignored; an explicit `port` skips the files; no file and no listener on the default port gives `NoViewerError`. The function never starts a program. `settings_dir()` test: monkeypatch `sys.platform`, `APPDATA` and `HOME`.

- [x] T017 [Found] [M] Implement `discovery.py`.
  - Files: `packages/loop_link/src/loop_link/discovery.py`
  - Needs: T016
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_discovery.py` passes.
  - Context: Public: `DEFAULT_PORT = 47615`, `ViewerAddress`, `settings_dir()`, `find_viewer(...)`. Token text on failure: "Set GEOVIEWER_TOKEN, or save a token in the link popup of the viewer." (`TokenError`). Alive check: `socket.create_connection((host, port), timeout)` then close. Sort key of live instances: `(kind != "window", -focused)`. Read JSON with `json.loads`; catch `OSError` and `ValueError` per file. The `NoViewerError` text says: "No viewer found. Start geoviewer and turn the link on. This library never starts a viewer." Use only the standard library. An explicit `port` argument is returned as given: do no alive check and do not read the instance files (the tests of T016 assume this). The test file imports `from loop_link import discovery`.

- [x] T019 [Found] [M] Implement the fake viewer.
  - Files: `packages/loop_link/src/loop_link/testing/__init__.py`, `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T008, T006
  - Wave: 5
  - Done when: `ruff check packages/loop_link` passes and `packages/loop_link/.venv/Scripts/python.exe -c "from loop_link.testing import FakeViewer"` passes. Behaviour is tested by T020 and T025.
  - Context: Base: `C:\Users\lachl\dev\plugin_loopstructural-geoviewer\tests\unit\geoviewer_fake_server.py` (179 lines, read it). New class `FakeViewer(token="secret", crs="EPSG:28350", origin=(0.0, 0.0, 0.0), kind="window", heartbeat_timeout=15.0, reply_errors=None, read_delay=0.0)`. Members: `.port`, `.objects: dict[str, StoredObject]` (`StoredObject(header: dict, buffers: dict[str, np.ndarray])`), `.received: list[Frame]`, `.inventory() -> list[dict]` (`id`, `object_type`, `revision`), `.write_instance_file(settings: Path) -> Path` (writes `instances/<pid>.json` and `link.json` is not written), `.accept_frame(data: bytes) -> list[Frame]` (apply one data-client frame and return the reply frames), `.close()`, context manager. Rules: listen on `127.0.0.1` port 0 in a thread; many connections; bearer token check (HTTP 401 when wrong); answer `hello` with `hello` (`application {"name":"fake-viewer"}`, `crs`, `origin`) then `scene.inventory` (`entries`); the second data connection (path not `/tool`) gets an `error` frame code `busy` then is closed; close a data connection with no frame for `heartbeat_timeout` seconds; the path `/tool` is accepted and every message gets `error` code `unknown_message`; `scene.upsert` rules (geometry roles present for the type, index < vertex count, data array length from association, finite coordinates, buffer sizes match, object size <= 2 GiB) else `error` code `invalid_object` with field `object`; `reply_errors` maps id to a forced error text; `scene.remove` and `scene.clear` change `objects`; `read_delay` waits that many seconds after each frame. The object keeps the data until a `scene.remove`. Frames use `loop_link.protocol`. Add docstrings and type hints.

- [x] T018 [Found] [M] Test support: `conftest.py`.
  - Files: `packages/loop_link/tests/conftest.py`
  - Needs: T019, T017
  - Wave: 6
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests --collect-only -q` passes.
  - Context: Fixture `fake_viewer` (function scope): start `FakeViewer()`, write its instance file into `tmp_path / "geoviewer"`, `monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")`, yield an object with `.viewer` (the `FakeViewer`) and `.settings` (the settings `Path`), and close at the end. Fixture `viewer_factory`: returns a function that makes `Viewer(settings_dir=<settings>, heartbeat_interval=0.1)` objects and closes all of them at teardown. Helper `geoviewer_fixtures() -> Optional[Path]`: the path in `LOOP_LINK_GEOVIEWER_FIXTURES` when it is a folder, else `None`. No test code in this file.

- [x] T020 [Found] [M] Test: `Viewer` in `test_viewer.py` (story 1).
  - Files: `packages/loop_link/tests/test_viewer.py`
  - Needs: T018, T012
  - Wave: 7
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_viewer.py` fails with ImportError of `loop_link.viewer`, then passes after T021.
  - Context: Use fixtures `fake_viewer` and `viewer_factory` of `conftest.py`. API: `Viewer(*, port=None, token=None, host="127.0.0.1", crs="", origin=(0,0,0), settings_dir=None, timeout=5.0, heartbeat_interval=5.0)`; `connect()`, `send(obj, *, id, name=None, parent=None) -> int`, `send_group`, `remove`, `clear`, `inventory() -> list[InventoryEntry]`, `flush(settle=0.5, *, raise_on_refusal=True) -> list[Refusal]`, `release_data()`, `close()`, properties `revisions`, `scene_origin`, `scene_crs`, `connected`. Tests (one per scenario of story 1): 1 discovery with no port or token argument and the object in `fake_viewer.viewer.objects`; 2 same id twice gives revisions 1 then 2 and the stored header has revision 2; 3 no viewer gives `NoViewerError` and no process starts; 4 properties are the same names and values in the stored buffers; 5 one object of each of the 7 types gives 7 objects with the right `object_type`; 6 a `Vectors.from_directions` object gives stored `dip_direction` and `dip` equal to the expected values; 7 a second `Viewer` gets `BusyError` and the objects of the first do not change, then after `release_data()` of the first the second sends and the objects of the first stay. Also: a wrong token gives `TokenError`; `reply_errors` makes `flush()` raise `ViewerRefusedError` with a `Refusal` that has the id; the scene origin of the viewer differs from the proposed origin: coordinates use the scene origin; a viewer that closes in the middle of a send (`FakeViewer.close()` in a thread) gives `ViewerClosedError`; keepalive: `FakeViewer(heartbeat_timeout=0.6)` with `heartbeat_interval=0.1` stays connected for 2 s and a later `send` works; an `error` code `unsupported_version` gives `MissingFeatureError`.

- [x] T021 [Found] [M] Implement `Viewer` in `viewer.py`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T020, T014, T015, T017
  - Wave: 8
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_viewer.py` passes.
  - Context: Class and method signatures as in T020 (also `@dataclass(frozen=True) class InventoryEntry: id: str; object_type: str; revision: int`). `connect()`: `find_viewer(...)`, `WebSocketClient.connect(host, port, token=token, timeout=timeout, path="/")`, send `hello_frame(crs, origin, __version__)`, read until both `hello` and `scene.inventory` arrive (else `ViewerClosedError`); an `error` frame with code `busy` raises `BusyError` and closes the socket; `unsupported_version` raises `MissingFeatureError`; other `error` raises `ViewerRefusedError`. Scene origin: the viewer `hello` origin when not null, else the proposed origin; keep `scene_crs` from the viewer `hello`. All frames use the scene origin. Revision of an id: `max(local, inventory) + 1`, 1 for a new id. `send(obj, ...)` calls `adapters.convert` only when `obj` is not a `SceneObject` (T027 adds the hook; before that, a non-`SceneObject` raises `InvalidInputError` that names the accepted inputs). After each frame, drain the socket with a 1 ms wait and keep each `error` frame that has an `object` as a `Refusal`. `flush` waits up to `settle` for late frames. Keepalive: daemon thread `loop-link-heartbeat`; every `heartbeat_interval` it tries a non-blocking `lock.acquire`; if the last send was older than the interval, it sends `heartbeat_frame(seq)`; an `OSError` stops the thread and marks the connection closed. `release_data` and `close` set a `threading.Event`, join with timeout `2*heartbeat_interval`, close the socket, and keep `revisions`. OP2: default `crs=""`; first read the hello handling of the geoviewer repository (`C:\Users\lachl\dev\geoviewer\crates`, search `Hello`) and, if the viewer refuses an empty CRS, raise `InvalidInputError("give crs=...")` for an empty scene. Convert `OSError` and `WebSocketClosed` on send to `ViewerClosedError`. Test assumptions of T020 (satisfy them): `viewer.revisions` is a mapping indexable by id (id -> int); `InventoryEntry` has attributes `id`, `object_type`, `revision`; `scene_origin` is a tuple of 3 floats; `flush()` returns a list of `Refusal`; `connect()` raises `BusyError`, `TokenError`, `NoViewerError` and `MissingFeatureError` (for an `unsupported_version` error); the `NoViewerError` test monkeypatches `subprocess.Popen` and passes `token="secret"` with an empty settings dir.

- [x] T023 [US1] [S] Export the public names in `__init__.py`.
  - Files: `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "import loop_link as l; [getattr(l, n) for n in l.__all__]"` passes and `ruff check packages/loop_link` passes.
  - Context: Keep `__version__ = "0.1.0"`. Import and list in `__all__`: from `.viewer`: `Viewer`, `InventoryEntry`; from `.objects`: `Points`, `Vectors`, `Lines`, `Surface`, `BlockModel`, `TetraMesh`, `Drillhole`, `Categorical`; from `.errors`: `LoopLinkError`, `NoViewerError`, `TokenError`, `BusyError`, `MissingFeatureError`, `InvalidInputError`, `ObjectTooLargeError`, `ViewerClosedError`, `ViewerRefusedError`, `MissingExtraError`, `Refusal`. Do not export `PushResult` or `push_objects` yet (T036). Do not import pyvista, pandas, Qt or LoopStructural here.

- [x] T024 [P] [US1] [S] Test: SC-001 speed in `test_performance.py`.
  - Files: `packages/loop_link/tests/test_performance.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_performance.py` passes.
  - Context: Mark the test `@pytest.mark.slow`. Use fixtures `fake_viewer` and `viewer_factory`. Build `BlockModel(origin=(0,0,0), cell_sizes=(1,1,1), cells=(100,100,100), cell_properties={"grade": np.random.default_rng(0).random(1_000_000, dtype=np.float32)})`. Time with `time.perf_counter()` from the start of `viewer.send(obj, id="grid")` to the moment that `fake_viewer.viewer.inventory()` has an entry with `id == "grid"` (poll every 10 ms, give up after 10 s). Assert the time is below 5.0 s (SC-001). The stored buffer must equal the input (`np.array_equal`).

- [x] T025 [P] [US1] [M] Test: frames of the geoviewer repository (SC-004).
  - Files: `packages/loop_link/tests/test_fixtures.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fixtures.py -rs` passes, or skips with the reason "set LOOP_LINK_GEOVIEWER_FIXTURES"; with the variable set to `C:\Users\lachl\dev\geoviewer\fixtures\protocol\python` it passes with no skip.
  - Context: Use the helper `geoviewer_fixtures()` of `conftest.py` (OP8). `pytest.skip("set LOOP_LINK_GEOVIEWER_FIXTURES to the folder fixtures/protocol/python of the geoviewer repository")` when it is `None`. Test 1: for each `*.bin` file in the folder, `FakeViewer.accept_frame(bytes)` returns no `error` reply (use a fresh `FakeViewer` and send a `hello_frame` first). Test 2: for the file `scene_upsert_drillhole_end_of_hole.bin`, decode it with `decode_frame`, rebuild the same object with `Drillhole` and `encode_object` (read the header for `id`, `name`, collar, survey, intervals, `end_of_hole`), decode the result, and assert that the headers are equal except `version` and except the buffer NAMES (the real frame names its buffers `c`, `s`, `i`; `encode_object` names them `geometry/collar`, `geometry/survey`, `geometry/intervals`; the viewer takes any name). Normalise: replace each buffer name by its role through `header['geometry']` in both `geometry` and the `buffers` list, then compare; the buffers are equal by role with `np.array_equal`. Do not change `encode.py` for this. Other plugin frames: compare the header keys `type`, `object_type`, `geometry` only.

- [x] T027 [US1] [M] Add the adapter dispatch `adapters/__init__.py` and hook it into `Viewer.send`.
  - Files: `packages/loop_link/src/loop_link/adapters/__init__.py`, `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_viewer.py packages/loop_link/tests/test_import.py` passes (no regression) and `ruff check packages/loop_link` passes.
  - Context: `adapters/__init__.py`: `convert(obj: Any, **hints: Any) -> SceneObject`. Return `obj` when it is a `SceneObject`. Else choose by the first part of `type(obj).__module__`: `pyvista` -> `adapters.pyvista.from_pyvista(obj, **hints)`; `loop_common` -> `adapters.loopstructural.from_loopstructural(obj)`; import the adapter module lazily inside the branch. Any other type raises `InvalidInputError("Accepted inputs: loop_link objects (Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole); pyvista meshes and grids; LoopStructural Surface and StructuredGrid; pandas tables for drillholes (use loop_link.adapters.pandas.drillholes_from_pandas).")`. The adapter modules are created by T029 and T031; `convert` imports them lazily, so a missing module raises `ImportError` only when used. Do not create stubs for them. In `viewer.py`, `Viewer.send` calls `convert(obj)` first. Add `**hints` to `Viewer.send` (`send(obj, *, id, name=None, parent=None, **hints)`) and pass them on.

- [x] T034 [P] [US5] [M] Test: plugin push API in `test_push.py` (story 5).
  - Files: `packages/loop_link/tests/test_push.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_push.py` fails with ImportError of `loop_link.push`, then passes after T035.
  - Context: Under test: `push_objects(viewer, objects, *, replace_prefix=None, progress=None, cancel=None) -> PushResult` and `PushResult(sent=0, removed=0, refused=[], skipped=[], cancelled=False)` with `summary() -> str`. `objects` is an iterable of `(id, object)`. Use `fake_viewer` and `viewer_factory`. Tests: (1) frames equal: push 3 objects and send the same 3 with `Viewer.send` to a second `FakeViewer`; the stored headers (without `version`) and buffers are equal (story 5, scenario 1); (2) `replace_prefix="model/"`: the fake starts with ids `model/old` and `other/keep` (send them first); after the push of `model/new` the ids `model/old` is removed, `other/keep` stays, `removed == 1`; (3) an object that raises `InvalidInputError` (for example an unsupported type `42`) goes to `skipped` as `(id, reason)` and the push goes on; (4) `cancel` set after the first object: `cancelled is True`, `sent == 1`; (5) `progress` gets the text `Sending <id>...` before each object, on the calling thread; (6) a slow viewer: `FakeViewer(read_delay=0.05)`; run `push_objects` in a `threading.Thread` for 5 objects; in the main thread count 20 ticks of 0.01 s with `time.sleep`; all 20 ticks finish while the push still runs or within 0.5 s of its start (the main thread is not blocked; story 5, scenario 2); (7) `summary()` text: "3 objects sent to GeoViewer", plus ", 1 old objects removed" and ", 1 objects skipped" when set, and ". The viewer refused N messages" when `refused` is not empty, ending with ".".

- [x] T037 [P] [Finish] [S] Manual checks document.
  - Files: `packages/loop_link/tests/manual_checks.md`
  - Needs: T021
  - Wave: 9
  - Done when: `git grep -n -c "^## M" packages/loop_link/tests/manual_checks.md` gives 2.
  - Context: Two sections. `## M1: speed with the real headless viewer (SC-001)`: start `geoviewer-headless --port 0` with `GEOVIEWER_TOKEN=secret` (it prints the port); run a script that sends a `BlockModel` of 100 x 100 x 100 cells with one `float32` cell property to it and prints the time until the object is in `Viewer.inventory()` after a reconnect; pass if under 5 s; record the time. `## M2: seven types in the window viewer`: start the viewer with a window, turn the link on, save a token, send one object of each of the 7 types with one property each, check by eye that 7 objects with the right types show, and that the property names show in the property table. Keep each section to 8 steps or fewer, in ASD-STE100 English.

## Phase 3: Stage S2 adapters (FR-010)

- [x] T028 [P] [US1] [M] Test: the pyvista adapter.
  - Files: `packages/loop_link/tests/test_adapter_pyvista.py`
  - Needs: T027
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pyvista.py` fails with ImportError of `loop_link.adapters.pyvista`, then passes after T029.
  - Context: `pv = pytest.importorskip("pyvista")`. Under test: `from_pyvista(mesh, *, kind=None, vectors=None, vectors_kind="normal") -> SceneObject` and `Viewer.send(mesh, id=...)`. Cases: `pv.PolyData(points)` with only vertices -> `Points`, `point_data` -> `properties`; with a point array `"n"` of shape (n,3) and `vectors="n"` -> `Vectors` (check `dip` against the formula arccos(n_z) for upward normals); `pv.lines_from_points` PolyData -> `Lines` with segments (m,2); `pv.Sphere().triangulate()` -> `Surface` (vertices, triangles (m,3), cell data -> `cell_properties`); a quad `pv.Plane()` face set is triangulated -> `Surface`; `pv.ImageData(dimensions=(4,5,6), spacing=(1,2,3), origin=(10,20,30))` -> `BlockModel(origin=(10,20,30), cell_sizes=(1,2,3), cells=(3,4,5))` with point data on nodes and cell data on cells; `pv.RectilinearGrid` with uneven spacing -> `InvalidInputError` (OP9); `pv.StructuredGrid` -> `InvalidInputError`; an `UnstructuredGrid` of all tetra cells (`pv.UnstructuredGrid(cells, [pv.CellType.TETRA], points)`) -> `TetraMesh` with tetrahedra (m,4); a mixed-cell `UnstructuredGrid` -> `InvalidInputError`; `kind="Drillhole"` -> `InvalidInputError` whose text names `Drillhole` and `drillholes_from_pandas` (OP3); a 2-D array in `point_data` of shape (n,2) -> `InvalidInputError`. Analytic grid test: `f = x + 10 y + 100 z` on the nodes of the ImageData grid; send it and read the stored buffer from the `FakeViewer`; it equals `f` at the node order x fastest.

- [x] T029 [US1] [M] Implement the pyvista adapter.
  - Files: `packages/loop_link/src/loop_link/adapters/pyvista.py`
  - Needs: T028
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pyvista.py` passes.
  - Context: `from_pyvista(mesh, *, kind: Optional[str] = None, vectors: Optional[str] = None, vectors_kind: str = "normal") -> SceneObject`. Import `pyvista` inside the module; if it is missing raise `MissingExtraError("pyvista")` at the first call of `from_pyvista`, not at `import loop_link`. Rules: `PolyData` with `n_lines > 0` -> `Lines` (segments from the line cell array: each line of k points gives k-1 segments); `PolyData` with faces -> `Surface` (call `mesh.triangulate()` when `not mesh.is_all_triangles`; triangles `faces.reshape(-1, 4)[:, 1:]`); `PolyData` of vertices only -> `Points`, or `Vectors.from_directions(points, mesh.point_data[vectors], kind=vectors_kind)` when `vectors` is set; `ImageData` -> `BlockModel(origin=mesh.origin, cell_sizes=mesh.spacing, cells=np.array(mesh.dimensions) - 1)`; `RectilinearGrid` with equal spacing on each axis (`np.allclose(np.diff(axis), step)`) -> `BlockModel`, else `InvalidInputError` (OP9); `UnstructuredGrid` with all cells `pyvista.CellType.TETRA` -> `TetraMesh` (`cells_dict[pyvista.CellType.TETRA]`); any other type or `kind="Drillhole"` -> `InvalidInputError` that names the accepted inputs (OP3 text: "A pyvista mesh cannot hold a drillhole. Use loop_link.Drillhole or loop_link.adapters.pandas.drillholes_from_pandas."). `point_data` -> `properties`, `cell_data` -> `cell_properties`, as numpy arrays; a 2-D array raises `InvalidInputError` that names the array. Order of pyvista image data is x fastest: no reorder. The arrays that pyvista adds itself, named `Normals` and `Texture Coordinates`, are skipped without an error (a plain `pv.Sphere()` or `pv.Plane()` must be accepted); any other 2-D array still raises `InvalidInputError` that names the array. Read the `vectors` array before the skip. Add docstrings and type hints.

- [x] T030 [P] [US1] [M] Test: the LoopStructural adapter.
  - Files: `packages/loop_link/tests/test_adapter_loopstructural.py`
  - Needs: T027
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_loopstructural.py` fails with ImportError of `loop_link.adapters.loopstructural`, then passes after T031.
  - Context: `geometry = pytest.importorskip("loop_common.geometry")`. Classes: `loop_common.geometry.Surface` (fields `vertices`, `triangles`, `properties: dict|None`, `cell_properties: dict|None`, `name`) and `loop_common.geometry.StructuredGrid` (fields `origin`, `step_vector`, `nsteps` (node counts), `properties` (per node, x fastest), `cell_properties` (per cell, x fastest), property `nsteps_cells`). Under test: `from_loopstructural(obj)` and `Viewer.send(obj, id=...)`. Cases: `Surface` with 4 vertices, 2 triangles, `properties={"a": ...}`, `cell_properties={"b": ...}` -> `Surface` with the same arrays; `properties=None` -> empty mapping; `StructuredGrid(origin=[1,2,3], step_vector=[1,2,3], nsteps=[4,5,6])` -> `BlockModel(origin=(1,2,3), cell_sizes=(1,2,3), cells=(3,4,5))`; analytic check: set `properties["f"]` to `grid.nodes @ [1, 10, 100]` and send it; the stored buffer from the `FakeViewer` equals it (so no reorder is needed); any other class from `loop_common.geometry` (for example `ValuePoints`) -> `InvalidInputError` that names the accepted inputs.

- [x] T031 [US1] [M] Implement the LoopStructural adapter.
  - Files: `packages/loop_link/src/loop_link/adapters/loopstructural.py`
  - Needs: T030
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_loopstructural.py` passes.
  - Context: `from_loopstructural(obj: Any) -> SceneObject`. Do not import `loop_common` at module level. Find the classes with `module = sys.modules.get("loop_common.geometry")`; if `module is None` or `obj` is not an instance of `module.Surface` or `module.StructuredGrid`, raise `InvalidInputError` that names the accepted inputs. If `loop_common` cannot be imported when the extra is needed, raise `MissingExtraError("loopstructural")`. `Surface` -> `Surface(vertices=obj.vertices, triangles=obj.triangles, properties=obj.properties or {}, cell_properties=obj.cell_properties or {})`. `StructuredGrid` -> `BlockModel(origin=obj.origin, cell_sizes=obj.step_vector, cells=obj.nsteps_cells, properties=obj.properties, cell_properties=obj.cell_properties)`. Coordinates are the public attributes (world coordinates, spec 005). No other transform. Add docstrings and type hints.

- [x] T032 [P] [US1] [M] Test: the pandas adapter.
  - Files: `packages/loop_link/tests/test_adapter_pandas.py`
  - Needs: T027
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pandas.py` fails with ImportError of `loop_link.adapters.pandas`, then passes after T033.
  - Context: `pd = pytest.importorskip("pandas")`. Under test: `drillholes_from_pandas(collar, survey, logs=None, *, hole_id="hole_id", x="x", y="y", z="z", depth="depth", azimuth="azimuth", dip="dip", from_="from", to="to", dip_positive_down=True) -> dict[str, Drillhole]`. Data: collar `hole_id,x,y,z` for `DH1`, `DH2`; survey `hole_id,depth,azimuth,dip` (2 rows each); logs `hole_id,from,to,grade,lith` (grade float, lith strings). Cases: result keys `DH1`, `DH2`; `Drillhole.collar` equals the row; `survey` shape (2,3) in the order depth, azimuth, dip; `intervals` shape (j,2); `interval_properties["grade"]` float array and `["lith"]` a numpy `U` array; `dip_positive_down=False` changes the sign of `dip`; a custom column name through `x="east"` works; a hole in `survey` that is not in `collar` -> `InvalidInputError`; an `object` column with a non-string value in `logs` -> `InvalidInputError`; `logs=None` gives holes with 0 intervals (`intervals` shape (0,2)); the result sends with `Viewer.send(hole, id=hole_id)` to the `FakeViewer` and the stored `object_type` is `Drillhole`.

- [x] T033 [US1] [M] Implement the pandas adapter.
  - Files: `packages/loop_link/src/loop_link/adapters/pandas.py`
  - Needs: T032
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pandas.py` passes.
  - Context: Signature as in T032. Import `pandas` inside the function body of the first use; if it is missing raise `MissingExtraError("pandas")`. Group the three tables by `hole_id` (`groupby`, `sort=False`); sort the survey rows by `depth`; `Drillhole(collar=[x,y,z], survey=array of [depth, azimuth, dip], intervals=array of [from, to], properties={}, interval_properties={extra columns of logs})`. `dip_positive_down=False` multiplies `dip` by -1 (the `geoh5` sign is the opposite of the viewer, which has positive down). A column of strings (all `str`) becomes `np.asarray(values, dtype="U")`; a numeric or bool column keeps its dtype; any other `object` column raises `InvalidInputError` that names the column. A hole id in `survey` or `logs` that is not in `collar` raises `InvalidInputError` that names the id. Return a `dict` that keeps the order of `collar`. Add docstrings and type hints.

## Phase 4: Stage S3 plugin API (story 5)

- [x] T035 [US5] [M] Implement `push.py`.
  - Files: `packages/loop_link/src/loop_link/push.py`
  - Needs: T034
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_push.py` passes.
  - Context: `@dataclass class PushResult: sent: int = 0; removed: int = 0; refused: list[Refusal] = field(default_factory=list); skipped: list[tuple[str, str]] = field(default_factory=list); cancelled: bool = False` and `summary(self) -> str` (text as in the test of T034: "<sent> objects sent to GeoViewer", ", <n> old objects removed", ", <n> objects skipped", ". The viewer refused <n> messages", final "."). `push_objects(viewer: Viewer, objects: Iterable[tuple[str, Any]], *, replace_prefix: Optional[str] = None, progress: Optional[Callable[[str], None]] = None, cancel: Optional[threading.Event] = None) -> PushResult`. Steps: `viewer.connect()` if `not viewer.connected`; for each `(id, obj)`: if `cancel` is set, set `cancelled=True` and stop; call `progress(f"Sending {id}...")`; `viewer.send(obj, id=id)` and `sent += 1`; `InvalidInputError` or `ObjectTooLargeError` -> `skipped.append((id, str(error)))` and go on; `ViewerClosedError` and `BusyError` are not caught. After the loop, if `replace_prefix`, for each id in `viewer.inventory()` that starts with the prefix and was not sent in this call: `viewer.remove(id)` and `removed += 1`. Then `refused = viewer.flush(raise_on_refusal=False)`. No Qt, no module state; `progress` runs on the calling thread. Docstrings and type hints.

- [x] T036 [US5] [S] Export `PushResult` and `push_objects`.
  - Files: `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T035, T023
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "from loop_link import PushResult, push_objects"` passes and `ruff check packages/loop_link` passes.
  - Context: Add `from .push import PushResult, push_objects` to `__init__.py` and add both names to `__all__`. Keep all other lines.

- [x] T022 [Finish] [S] Test: the core import needs no Qt or optional extras (SC-005).
  - Files: `packages/loop_link/tests/test_import.py`
  - Needs: T036
  - Wave: 12
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_import.py` passes.
  - Context: Run in a subprocess: `subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)` where `code` imports `loop_link`, then prints `sorted(m for m in sys.modules if m.split(".")[0] in {"PyQt5","PyQt6","PySide2","PySide6","qgis","pyvista","pandas","loop_common","LoopStructural"})`. The test asserts the printed list is `[]` and the return code is 0. Add a second test: `loop_link.__all__` names all exist (`getattr`).

## Phase 5: Finish

- [x] T038 [Finish] [M] Run the merge gates and fix the findings.
  - Files: `packages/loop_link/`
  - Needs: T002, T003, T004, T005, T022, T024, T025, T029, T031, T033, T037
  - Wave: 13
  - Done when: `ruff check LoopStructural packages && pytest` passes, `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests` passes, and `black --check packages/loop_link` passes.
  - Context: Constitution merge gate: `pytest` and `ruff check LoopStructural packages` must both pass. The root `pytest` does not run the package tests, so run the package command too. Run `black packages/loop_link` to fix format (line length 100, quotes kept). Fix only problems in `packages/loop_link/`. Do not change a gate or add a Ruff ignore for `packages/loop_link/src`. Run `git grep -n loop_link ROADMAP.md .release-please-config.json .release-please-manifest.json` and check 2 files or more show the name (SC-006). Report the result of each command in the commit message.

## Phase 6: Stage S5a, tool connection and selection (stories 2 and 3)

Test command in this phase: `packages/loop_link/.venv/Scripts/python.exe -m pytest <file>`. Geoviewer spec 067 (protocol 1.18) is the viewer side. The text of the spec calls it "066" (plan OP7).

- [x] T039 [P] [Found] [S] Three new errors, with a test.
  - Files: `packages/loop_link/src/loop_link/errors.py`, `packages/loop_link/tests/test_errors.py`
  - Needs: none
  - Wave: 1
  - Done when: write the test first and see it fail; then `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_errors.py` passes.
  - Context: Add three classes in the style of `ObjectTooLargeError`, all subclasses of `LoopLinkError`. `StaleRevisionError(object_id: str, revision: Optional[int] = None)`: attributes `object_id`, `revision`; text "The object '<id>' is not at that revision in the viewer. Read the selection or fetch the object again, then try again." `CancelledError(message: str = "The request was cancelled.")`. `UnsupportedObjectError(object_id: str, object_type: str)`: attributes with the same names; text "The object '<id>' has the type '<type>'. loop_link can fetch Points, Curve, Surface, BlockModel, TetraMesh and Drillhole objects only." Add the names to `__all__` if the module has one. Test: each is a `LoopLinkError`; attributes are kept; `str(error)` holds the id (and the type); `CancelledError().args` has the default text.

- [x] T040 [P] [Found] [S] Tool protocol versions and the `version` argument of `encode_frame`.
  - Files: `packages/loop_link/src/loop_link/protocol.py`, `packages/loop_link/tests/test_protocol.py`
  - Needs: none
  - Wave: 1
  - Done when: write the test first and see it fail; then `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_protocol.py` passes (the old tests stay green).
  - Context: Add `TOOL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 18}` and `CANCEL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}` next to `PROTOCOL_VERSION`. New signature: `encode_frame(header, buffers=(), *, version: Optional[dict[str, int]] = None) -> bytes`. `version=None` keeps `PROTOCOL_VERSION` (1.14), so every old frame is the same. Else the header gets `dict(version)`. Update the numpy-style docstring (new parameter). `decode_frame` does not change. Tests: default header version is 1.14; `version=TOOL_PROTOCOL_VERSION` gives `{"major":1,"minor":18}`; `version=CANCEL_PROTOCOL_VERSION` gives minor 11; the argument is keyword-only (a third positional argument raises `TypeError`); the dict that the caller gives is not changed.

- [x] T041 [Found] [M] Test: the four tool frame builders in `test_encode.py`.
  - Files: `packages/loop_link/tests/test_encode.py`
  - Needs: T040
  - Wave: 2
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_encode.py` fails with AttributeError for the new builders, then passes after T042.
  - Context: Add tests for `loop_link.encode`; decode each result with `decode_frame`. Signatures: `selection_set_frame(request: int, object_id: str, cells: Optional[np.ndarray], vertices: Optional[np.ndarray], revision: Optional[int]) -> bytes`; `selection_get_frame(request: int)`; `fetch_frame(request: int, id: str)`; `cancel_frame(request: int)`. Expected headers (add `version`, `buffers`): `{"type":"selection.set","request":3,"object":"grid","revision":2,"entries":[{"kind":"cells","indices":"cells"},{"kind":"nodes","indices":"nodes"}]}` with two `u32` buffers `cells` and `nodes`, shape `[k]`; `{"type":"selection.get","request":7}`; `{"type":"object.fetch","request":2,"id":"grid"}`; `{"type":"command.cancel","request":2}`. Version of the first three: `{"major":1,"minor":18}`. Version of `cancel_frame`: `{"major":1,"minor":11}`. Rules to test: an empty or `None` array adds no entry and no buffer; both empty gives `"entries": []` and `buffers: []`; `revision=None` leaves out the key `revision`; only `vertices` gives one entry with kind `nodes` and buffer `nodes`; list input becomes `uint32`; the entry order is cells, then nodes. Old tests stay green.

- [x] T042 [Found] [M] Implement the four tool frame builders in `encode.py`.
  - Files: `packages/loop_link/src/loop_link/encode.py`
  - Needs: T041
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_encode.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add `selection_set_frame`, `selection_get_frame`, `fetch_frame`, `cancel_frame` with the signatures and headers of T041, and add the names to `__all__`. Use `protocol.encode_frame(header, buffers, version=protocol.TOOL_PROTOCOL_VERSION)`; `cancel_frame` uses `protocol.CANCEL_PROTOCOL_VERSION`. Buffers: `np.ascontiguousarray(np.asarray(x, dtype=np.uint32))`, 1-D, names `cells` and `nodes`. Do not check the index range here (`Viewer.select` does). Numpy-style docstrings and type hints.

- [x] T043 [Found] [M] Test: the tool endpoint of the fake viewer (selection and fetch rules).
  - Files: `packages/loop_link/tests/test_fake_viewer.py`
  - Needs: T042
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py` fails for the new tests (the tool endpoint still answers `unknown_message`), then passes after T045.
  - Context: Talk to the fake with `WebSocketClient.connect("127.0.0.1", fake.port, token="secret", path="/tool")`, `send_binary(frame)`, `receive(timeout)` and `decode_frame`. Send objects first with `fake.accept_frame(encode_object(...))`. New `FakeViewer` arguments: `protocol_minor=18` (below 18: each selection or fetch message gets an `error` frame code `unknown_message`), `part_bytes=64*1024*1024` (multiple of 8), `fetch_delay=0.0`; `kind="headless"` makes `selection.set` and `selection.get` answer `request.error` code `needs_window`. Members: `.selection` (`{"object","revision","cells","nodes"}` or `None`), `.set_user_selection(object_id, *, cells=None, nodes=None)`, `.tool_requests` (list of `Frame`). One test for each rule. `selection.set` replies `selection.result` (`request`, `object`, `revision`, `entries` with `kind`, `count`) or `request.error` (`request`, `error`={code, message, field}) with code `invalid_input` and field: unknown object `object`; wrong `revision` `revision`; kind the object lacks `kind` (Points: nodes only; Curve, Surface, BlockModel, TetraMesh: nodes and cells; Drillhole, Group: no table); index at or above the count `indices`; buffer not `u32` or missing `indices`. Nodes: vertex count (BlockModel `(nx+1)(ny+1)(nz+1)`); cells: segments, triangles, tetrahedra or `nx*ny*nz`. A repeated index counts one time. A new selection replaces the old one (also on another object). Empty `entries` clears it. `selection.get`: reply has `object`, `revision` and one entry per kind with a `u32` buffer named `cells` or `nodes`; no selection gives `"entries": []` and no `object`. `object.fetch`: unknown id gives `invalid_input` field `id`; else `object.header` (`origin`, `crs`, `object` = stored header with `buffers: []`, `arrays` = `[{name,dtype,shape}]`), then `object.part` (`array`, `offset`, `bytes`, one `u8` buffer `part`, at most `part_bytes`, offset multiple of 8) per array in order, then `object.end` (`parts` = count). With more than one part, `command.progress` (`fraction`, `text`) follows each part. A `command.cancel` with the same `request` while `fetch_delay=0.05` gives `request.error` code `cancelled`. A fetch does not change objects or selection. All replies carry the `request` of the message. Replies use version 1.18.

- [x] T044 [Found] [M] Fake viewer: tool endpoint, selection part.
  - Files: `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T043, T019
  - Wave: 6
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py -k "selection or protocol_minor or headless or user_selection or tool_requests"` passes, and the whole old test file still passes.
  - Context: Change `FakeViewer.__init__` (new keyword arguments with defaults `protocol_minor=18`, `part_bytes=64*1024*1024`, `fetch_delay=0.0`, as in T043) and `_run` (the `tool=True` branch now calls a new method `tool_frame(data: bytes) -> list[Frame]` and appends each frame to `.tool_requests`; keep `read_delay`). Implement `selection.set` and `selection.get` with the rules in T043; replies are `Frame(header, buffers)` with `version` 1.18 and the `request` of the message (look at `_error` and `_encode` for the frame helpers). Messages of another type and all messages when `protocol_minor < 18` get `error` code `unknown_message`. Add `.selection`, `.set_user_selection`, `.tool_requests`. The data path (`accept_frame`) does not change. Add docstrings and type hints.

- [x] T045 [Found] [M] Fake viewer: tool endpoint, fetch part.
  - Files: `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T044
  - Wave: 7
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add `object.fetch` to `tool_frame` (T044) with the rules in T043. The tool loop must send several reply frames and read between parts: after each `object.part` the fake reads for at most 0.001 s on the connection (use a short socket timeout, restore it after) and, on a `command.cancel` with the same `request`, sends `request.error` code `cancelled` and stops. So `_run` needs a streaming path for fetch, not a list of replies; keep `tool_frame` for the one-reply messages. `fetch_delay` is a `time.sleep` after each part. Parts: for each stored buffer in header order, slices of at most `part_bytes` bytes (`array.tobytes()`), `offset` = bytes already sent for that array; `arrays` shape from the stored array. A stored buffer of zero bytes still is in `arrays` and gets no part. `parts` in `object.end` is the total part count. Origin and crs in the header are those of the fake.

- [x] T046 [Found] [M] Test: `ToolConnection` in `test_tool.py`.
  - Files: `packages/loop_link/tests/test_tool.py`
  - Needs: T045, T039
  - Wave: 8
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_tool.py` fails with ImportError of `loop_link.tool`, then passes after T047.
  - Context: Use fixtures `fake_viewer` (`.viewer` is the `FakeViewer`, `.settings` the settings folder) and `viewer_factory` of `conftest.py`. API: `ToolConnection(*, port=None, token=None, host="127.0.0.1", settings_dir=None, timeout=5.0)` with `connected`, `connect()`, `request(build, *, until, feature, object_id=None, on_frame=None, progress=None, cancel=None, timeout=None) -> Frame`, `close()`; `raise_for_request_error(error, *, object_id)`. Frames come from `loop_link.encode` (`selection_get_frame`, `fetch_frame`, ...). Tests: request numbers start at 1 and rise by 1 for each request; the reply of `selection.get` is returned for `until=("selection.result",)`; a frame with another request number is dropped (send a stale frame with a raw server or by calling twice); `progress` gets rising floats for a fetch with `part_bytes=64`; `cancel` set gives `CancelledError` (with `fetch_delay=0.05`) and the fake sees one `command.cancel`; the socket closes in the middle of a request gives `ViewerClosedError`, and the next request connects again; wrong token gives `TokenError`; no viewer gives `NoViewerError`; no frame within `timeout=0.2` gives `ViewerClosedError`. Error table (call `raise_for_request_error` directly with `{"code","message","field"}`): `invalid_input`+field `revision` -> `StaleRevisionError(object_id)`; `invalid_input` other field -> `InvalidInputError` (text starts "The viewer refused the request:"); `needs_window` -> `MissingFeatureError`; `cancelled` -> `CancelledError`; `busy` -> `BusyError`; other code -> `ViewerRefusedError` with one `Refusal(object_id or "", code, message)`. An `error` frame code `unknown_message` (fake with `protocol_minor=14`) -> `MissingFeatureError` whose `feature` is the `feature` argument; code `unsupported_version` -> `MissingFeatureError`.

- [ ] T047 [Found] [M] Implement `tool.py`.
  - Files: `packages/loop_link/src/loop_link/tool.py`
  - Needs: T046
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_tool.py` passes and `ruff check packages/loop_link` passes.
  - Context: Create `ToolConnection` and `raise_for_request_error(error: Mapping[str, Any], *, object_id: Optional[str]) -> NoReturn` (signatures and error table in T046; the `error` frame rows take `code`, `message` from the frame header). `request`: take `self._lock`; connect when `self._client is None` (`discovery.find_viewer(...)`, then `WebSocketClient.connect(host, port, token=..., timeout=..., path="/tool")`; no `hello`; no keepalive thread). Number `n` starts at 1. Send `build(n)`. Loop `client.receive(0.1)`; `protocol.decode_frame` (a `ProtocolError` closes the socket and raises `ViewerClosedError("The viewer sent a frame that is not valid. ...")`). `type == "error"` -> `raise_for_request_error`; a frame with another `request` is dropped; for `request == n`: `request.error` -> `raise_for_request_error(header["error"], object_id=object_id)` (for `unknown_message` of an `error` frame the feature text is the `feature` argument); `command.progress` -> `progress(fraction)` when `fraction` is a number; `on_frame(frame)` sees every other frame of the request; a type in `until` ends the loop and is returned. The wait restarts at each frame; after `timeout` (default `self._timeout`) with no frame close and raise `ViewerClosedError("The viewer did not answer in time. ...")`. When `cancel` is set: send `cancel_frame(n)` one time and keep reading until `cancelled` (raise `CancelledError`) or an `until` frame (return it). A socket error (`OSError`, `WebSocketClosed`, `WebSocketError`) sets the socket to `None` and raises `ViewerClosedError`. Docstrings and type hints.

- [ ] T048 [US2] [M] Test: select and read the selection in `test_selection.py` (stories 2 and 3).
  - Files: `packages/loop_link/tests/test_selection.py`
  - Needs: T047, T045
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_selection.py` fails with ImportError of `loop_link.selection`, then passes after T050.
  - Context: Use `fake_viewer` and `viewer_factory`. Send a `BlockModel` (cells 4x4x4, id `grid`) first. API on `Viewer`: `select(object_id, *, cells=None, vertices=None, revision=None) -> dict[str, int]` (returns `{"cell": n, "vertex": m}`), `selection() -> Selection`, `watch_selection(callback, *, interval=0.5, initial=False) -> SelectionWatcher`, `tool_connected`. `loop_link.selection`: `Selection(object_id, revision, cells, vertices)` (frozen, `uint32` arrays, `.empty`, `.same_as(other)` compares object id and sorted cells and vertices, not the revision), `parse_selection_result(frame)` (an entry with a count but no buffer, a buffer not `uint32`, or a length not equal to `count` gives `ProtocolError`), `SelectionWatcher` (`running`, `error`, `stop()`, context manager). Tests: story 2 scenarios 1 to 4 (select cells 12 and 13: `fake.selection` holds them; index 64 gives `InvalidInputError` and the selection stays; `StaleRevisionError` for an old `revision`; `FakeViewer(protocol_minor=14)` gives `MissingFeatureError`); clear with `cells=[]`; `select()` with neither argument, a float index, a bool, a negative value or a value above 2**32-1 each give `InvalidInputError` that names the argument; `cells` and `vertices` together; duplicates count one time; 100,000 cells are one `selection.set` frame with a `u32` buffer of shape `[100000]` (`fake.tool_requests`, SC-003); `kind="headless"` gives `MissingFeatureError`; story 3: `selection()` after `set_user_selection` returns object id, cells, vertices; a 100,000-cell selection returns a `uint32` array; empty selection has `object_id is None`; `watch_selection` calls the callback once per change at the next poll (`interval=0.05`), not for no change, `initial=True` calls once with a non-empty first selection, `interval<=0` gives `InvalidInputError`, `stop()` ends the thread, no thread named `loop-link-selection-poll` exists before the call, an exception in the callback sets `.error` and ends the thread, a dropped `Viewer` ends the thread. A `Viewer` that only selects does not open the data connection (`viewer.connected` is false). `close()` stops watchers and a later `selection()` raises `ViewerClosedError`.

- [ ] T049 [US3] [M] Implement `selection.py`.
  - Files: `packages/loop_link/src/loop_link/selection.py`
  - Needs: T048
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_selection.py -k "parse or same_as or empty or dtype"` passes and `ruff check packages/loop_link` passes. The watcher tests pass after T050.
  - Context: Create `Selection` (`@dataclass(frozen=True, eq=False)`: `object_id: Optional[str]`, `revision: Optional[int]`, `cells: np.ndarray`, `vertices: np.ndarray`; property `empty`; `same_as`), `parse_selection_result(frame: Frame) -> Selection` (kind `cells` -> `cells`, kind `nodes` -> `vertices`; array is `frame.buffers[entry["indices"]]`; checks as in T048; `ProtocolError` from `loop_link.protocol`) and `SelectionWatcher`. The watcher is made only by `Viewer.watch_selection`: `SelectionWatcher(viewer_ref: weakref.ref, callback, interval, last: Selection)`, a daemon thread named `loop-link-selection-poll`. Loop: wait `interval` on a `threading.Event` (so `stop` is fast); get the viewer from the weak reference (end the thread if it is gone); call `viewer.selection()`; when the result is not `same_as` the last one, call `callback(selection)` one time and keep it. An exception of the poll or of the callback goes to `error` and the thread ends. `stop()` sets the event and joins with timeout `2 * interval + 1`. `__enter__` returns self, `__exit__` calls `stop`. No thread before `start`. Docstrings and type hints.

- [ ] T050 [US2] [M] Add `select`, `selection`, `watch_selection` and the tool connection to `Viewer`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T049, T047, T042
  - Wave: 12
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_selection.py packages/loop_link/tests/test_viewer.py packages/loop_link/tests/test_push.py` passes and `ruff check packages/loop_link` passes.
  - Context: New members as in T048. The `Viewer` makes its `ToolConnection` (from `loop_link.tool`) at the first tool call, with its own `port`, `token`, `host`, `settings_dir`, `timeout`; tool calls never open the data connection. `release_data` keeps the tool connection. `close` stops all watchers, closes the tool connection and sets `_dead`; a later tool call raises `ViewerClosedError`. `select`: `object_id` must be a non-empty `str`; at least one of `cells`, `vertices` not `None`, else `InvalidInputError("Give cells or vertices. Use an empty list to clear the selection.")`; each given array is a 1-D sequence of integers from 0 to 2**32-1 (no float, no bool, no negative, check `np.asarray`, dtype kind `i` or `u`, and bool dtype refused), else `InvalidInputError` that names the argument. One `selection.set` frame (`encode.selection_set_frame`), `until=("selection.result",)`, `feature="selection and fetch (geoviewer protocol 1.18)"`; return `{"cell": n, "vertex": m}` for the kinds in the reply. `selection()`: `selection.get` -> `parse_selection_result` -> `Selection` with `object_id`/`revision` from the header (`None` when absent). `watch_selection(callback, *, interval=0.5, initial=False)`: `interval <= 0` gives `InvalidInputError`; read the selection in the calling thread first (errors go to the caller); call the callback once when `initial` and the selection is not empty; then start `SelectionWatcher` with a `weakref.ref(self)`; keep it in a list. Do not change the data-connection code.

## Phase 7: Stage S5b, fetch (story 4)

- [x] T051 [P] [US4] [M] Test: `decode_object` in `test_fetch.py`.
  - Files: `packages/loop_link/tests/test_fetch.py`
  - Needs: T039
  - Wave: 2
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py` fails with ImportError of `loop_link.fetch`, then passes after T052.
  - Context: Under test: `loop_link.fetch.decode_object(header: Mapping, arrays: Mapping[str, np.ndarray], origin: Optional[Sequence[float]]) -> SceneObject`. `header` is the `object` field of a viewer `object.header`: the `scene.upsert` fields (`id`, `object_type`, `name`, `revision`, `geometry` {role: buffer name}, `params`, `data` [{id,name,buffer,data_type,association,offsets?,value_map?}], `property_groups`). Build cases by hand as the viewer writes them (use `encode_object` and `decode_frame` to get a header and buffers, then change them). Rules: `Points` -> `Vectors` when a group has role `orientation` or `lineation` (member arrays float64 become the tuples; they are not in `properties`), else `Points`; `Curve` -> `Lines`; `Surface`; `TetraMesh`; `BlockModel` (origin = `params["origin"] + origin`); `Drillhole` (collar + origin; `end_of_hole = params["end_of_hole"][0]` when present); `Group`, `Grid2D`, `GeoImage` -> `UnsupportedObjectError(id, type)`. `vertices` and `collar` = `array.astype(float64) + origin`; `origin` None means (0,0,0); result float64. Properties: association `vertex` -> `properties`; `cell` -> `cell_properties` (Drillhole: `interval_properties`); `interval` with no `log_set` or the first interval log set -> Drillhole `interval_properties`; `cell` on Points/Vectors -> `UnsupportedObjectError`; other log set, `sample`, `object` -> dropped with `warnings.warn("loop_link dropped the property '<name>' of object '<id>': <reason>", UserWarning)` one time each (use `pytest.warns`). Types: `float` keeps f32 or f64; `integer` -> int32; `boolean` -> numpy bool; `text` -> numpy `U` array from data bytes (u8) and `offsets` (u32); `referenced` -> `Categorical(codes, names, colours)` with names ordered by `id`, codes renumbered when ids are not 0..n-1, 3-value colours get alpha 255, a code not in `value_map` -> `InvalidInputError`.

- [x] T052 [US4] [M] Implement `FetchResult` and `decode_object` in `fetch.py`.
  - Files: `packages/loop_link/src/loop_link/fetch.py`
  - Needs: T051
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py` passes and `ruff check packages/loop_link` passes.
  - Context: Create the module with `FetchResult` (`@dataclass(frozen=True, eq=False)`: `id: str`, `object: SceneObject`, `object_type: str`, `revision: int`, `name: Optional[str]`, `parent: Optional[str]`, `origin: Optional[tuple[float, float, float]]`, `crs: Optional[str]`) and `decode_object` with the rules in T051. Use the classes of `loop_link.objects`; for the Vectors groups use members `dip_direction`, `dip` (role `orientation`) and `trend`, `plunge` (role `lineation`), each from the data array named in `members[].data`. Data array lookup: `arrays[entry["buffer"]]`; text uses `arrays[entry["buffer"]]` bytes and `arrays[entry["offsets"]]`. Do not import pyvista, pandas or Qt. `_Assembler` and `fetch_object` come in T054. Docstrings and type hints.

- [ ] T053 [US4] [M] Test: assembler, `fetch_object` and `Viewer.fetch` in `test_fetch.py`.
  - Files: `packages/loop_link/tests/test_fetch.py`
  - Needs: T052, T050
  - Wave: 13
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py` fails for the new tests (AttributeError on `fetch_object`, `_Assembler` and `Viewer.fetch`), then passes after T055.
  - Context: Names fixed by this task list: `loop_link.fetch._Assembler(object_id: str)` with `feed(frame: Frame) -> None` and `finish() -> FetchResult`; `fetch_object(tool, id, *, progress=None, cancel=None, timeout=30.0) -> FetchResult`; `Viewer.fetch(id, *, progress=None, cancel=None) -> FetchResult`. Name the Viewer tests `test_viewer_fetch_*`. Assembler tests feed hand-built frames (`Frame(header, buffers)`): first frame not `object.header` -> `ProtocolError`; arrays with the sum of bytes over 2 GiB (header only, no allocation) -> `ObjectTooLargeError`; part for an array not in `arrays`; `offset` not a multiple of 8; part past the end; a gap; an overlap; a missing part at `object.end`; wrong `parts` count; `part` buffer not `u8` -> each `ProtocolError`; a buffer name in `geometry`, `data[*].buffer` or `offsets` that is not in `arrays` -> `ProtocolError`; a good stream gives a `FetchResult` with the revision of `header["object"]["revision"]`. `fetch_object` and `Viewer.fetch` with the fake (`fake_viewer`, `viewer_factory`): fetch of a sent `Surface` equals the input; `FakeViewer(part_bytes=64)` gives several parts and `progress` calls that rise to 1; `FakeViewer(fetch_delay=0.05, part_bytes=64)` with `cancel` set gives `CancelledError` and `fake.objects` is the same; unknown id -> `InvalidInputError`; `kind="headless"` fake still fetches; a stream error closes the tool connection and raises `ViewerClosedError("The viewer sent an object that is not valid: ...")`; the revision is in the result; FR-014: program A holds the data connection (`send`), program B (`viewer_factory()`) gets `BusyError` from `send`, but `select`, `selection` and `fetch` work; after A calls `release_data`, B sends and the objects of A stay (`fake.objects`).

- [ ] T054 [US4] [M] Implement `_Assembler` and `fetch_object` in `fetch.py`.
  - Files: `packages/loop_link/src/loop_link/fetch.py`
  - Needs: T053
  - Wave: 14
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py -k "not test_viewer_fetch"` passes and `ruff check packages/loop_link` passes.
  - Context: `fetch_object` calls `tool.request(lambda n: encode.fetch_frame(n, id), until=("object.end",), feature="selection and fetch (geoviewer protocol 1.18)", object_id=id, on_frame=assembler.feed, progress=progress, cancel=cancel, timeout=timeout)`, then returns `assembler.finish()`. It catches `ProtocolError`, calls `tool.close()` and raises `ViewerClosedError("The viewer sent an object that is not valid: <text>")`; other library errors go up. `_Assembler` rules: (1) first frame must be `object.header`; (2) header checks: each `arrays` entry has a unique `name`, a `dtype` in `f32 f64 u32 i32 u8`, a `shape` of non-negative integers; the sum of bytes over `protocol.MAX_OBJECT_BYTES` -> `ObjectTooLargeError(id, size, limit)`; each buffer name in `object.geometry`, `object.data[*].buffer` and `offsets` is in `arrays`; then one `numpy.empty(nbytes, dtype=numpy.uint8)` per array; (3) `object.part`: array known; `offset % 8 == 0`; buffer `part` is `u8`; `offset` equals the bytes already filled in that array (in order, no gap, no overlap); `offset + bytes` at most the array size; copy in; (4) `object.end`: `parts` equals the count received and each array is full; (5) `finish`: view each array with its dtype and shape (`.view`, `.reshape`), call `decode_object(header["object"], arrays, header["origin"])`; fill `FetchResult` (`origin` as tuple or None, `crs`, `name`, `parent`, `revision`, `object_type`). Every break raises `ProtocolError`.

- [ ] T055 [US4] [M] Add `Viewer.fetch`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T054, T050
  - Wave: 15
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py packages/loop_link/tests/test_selection.py packages/loop_link/tests/test_viewer.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add `Viewer.fetch(self, id: str, *, progress: Optional[Callable[[float], None]] = None, cancel: Optional[threading.Event] = None) -> FetchResult`. It uses the tool connection of T050 (same lazy creation, same closed rule) and calls `fetch.fetch_object(tool, id, progress=progress, cancel=cancel)`. It never opens the data connection. An unknown id gives `InvalidInputError` (the tool layer maps `invalid_input`). `progress` gets floats from 0 to 1. Import `FetchResult` and `fetch_object` from `loop_link.fetch`. Numpy-style docstring with Parameters, Returns, Raises (`InvalidInputError`, `CancelledError`, `MissingFeatureError`, `ViewerClosedError`, `UnsupportedObjectError`).

- [ ] T056 [US4] [M] Test: send then fetch for the seven types in `test_roundtrip.py` (SC-002).
  - Files: `packages/loop_link/tests/test_roundtrip.py`, `packages/loop_link/src/loop_link/fetch.py`
  - Needs: T055
  - Wave: 16
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_roundtrip.py` passes.
  - Context: Use `fake_viewer` and `viewer_factory` (`FakeViewer(origin=(500000.0, 7000000.0, 0.0))` for the origin cases: build it in the test and write its instance file as `conftest.py` does). For each of `Points`, `Vectors` (made with `Vectors.from_directions`), `Lines`, `Surface`, `BlockModel`, `TetraMesh`, `Drillhole`: `viewer.send(obj, id=...)`, then `viewer.fetch(id).object`. Compare: floats with `np.testing.assert_allclose(rtol=1e-6)`; integers, booleans, text, `Categorical` (codes, names, colours) and index arrays with `np.array_equal`; the dtype of each property stays (float32 stays float32); coordinates are float64 world coordinates, also for a point 10 km from the origin. Vectors: compare `orientation` angles to the input within 1e-6 degrees (the viewer keeps angles; an x,y,z input is equal up to the conversion). The test file is the task. Change `fetch.py` only if a case shows a bug in `decode_object` or the assembler; then say so in the task result.

## Phase 8: Stage S5c, real codec, fixtures, manual check, docs

- [ ] T057 [Found] [S] Export the S5 public names.
  - Files: `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T055, T049, T039
  - Wave: 16
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "import loop_link as l; [getattr(l, n) for n in l.__all__]"` passes and `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_import.py` passes.
  - Context: Add imports and `__all__` entries (keep `__all__` sorted as it is): from `.errors`: `CancelledError`, `StaleRevisionError`, `UnsupportedObjectError`; from `.selection`: `Selection`, `SelectionWatcher`; from `.fetch`: `FetchResult`. Keep all old lines. Do not import pyvista, pandas, Qt or LoopStructural (SC-005, `test_import.py`).

- [x] T058 [P] [US2] [M] Test: tool frames against the real geoviewer codec.
  - Files: `packages/loop_link/tests/test_codec_real.py`
  - Needs: T042
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_codec_real.py -rs` passes or skips with the reason that the variable is not set; with `LOOP_LINK_GEOVIEWER_REPO=C:\Users\lachl\dev\geoviewer` and `cargo` on the PATH it passes with no skip.
  - Context: In `_cases()` add the frames: `selection_set` (cells and nodes, revision 2), `selection_set_no_revision`, `selection_set_cells_only`, `selection_set_clear` (both `None`), `selection_get`, `object_fetch`, `command_cancel`. Build them with `encode.selection_set_frame(3, "grid", np.array([1, 2, 3]), np.array([4]), 2)`, `encode.selection_get_frame(4)`, `encode.fetch_frame(5, "grid")`, `encode.cancel_frame(6)`. The existing helper, the validator program and the negative control stay. Do not change `codec_validator/`. If the validator reports FAIL for a type it does not know, make the task fail with the validator output as the reason; do not drop the case.

- [x] T059 [P] [US2] [M] Test: tool fixtures of the geoviewer repository (OP11).
  - Files: `packages/loop_link/tests/test_fixtures.py`
  - Needs: T042
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fixtures.py -rs` passes or skips with a reason; with `LOOP_LINK_GEOVIEWER_FIXTURES=C:\Users\lachl\dev\geoviewer\fixtures\protocol\python` it passes with no skip for the new tests.
  - Context: Tool fixtures are in `<LOOP_LINK_GEOVIEWER_FIXTURES>/../tool/` as JSON files `selection.set.json`, `selection.get.json`, `object.fetch.json`, `command.cancel.json`. Use `geoviewer_fixtures()` of `conftest.py`; skip with a reason when the variable or the `tool` folder is missing (same style as `_folder()`). Read each fixture (look at the file to see its form: header fields `type`, `request`, `object`/`id`, `revision`, `entries`, `version`, `buffers`). Build the frame with `encode.selection_set_frame`, `selection_get_frame`, `fetch_frame`, `cancel_frame` from the `request`, `object`, `id`, `revision` and index counts of the fixture (use `np.arange(count)` for indices), decode it, and assert the header equals the fixture header including `version` and `buffers`. A fixture with a different form than the plan expects: make the test fail with a clear message.

- [x] T060 [P] [Finish] [S] Manual check M3 in `manual_checks.md`.
  - Files: `packages/loop_link/tests/manual_checks.md`
  - Needs: none
  - Wave: 1
  - Done when: `git grep -n -c "^## M3" packages/loop_link/tests/manual_checks.md` gives 1.
  - Context: Add a section `## M3: select, read and fetch with the real window viewer (stories 2 to 4)` before `## Record`, in ASD-STE100 English, 8 steps or fewer in the style of M2. Steps: (1) send a `BlockModel` of 10x10x10 cells, select cells 12, 13, 14 with `Viewer.select`, check the viewer outlines three cells and `Viewer.selection()` returns them; (2) select cells in the window with the mouse and read them with `Viewer.selection()`; (3) `select` with index 5000 gives `InvalidInputError` and the selection stays; (4) send, replace, then `select` with the old revision gives `StaleRevisionError`; (5) fetch each of the seven types (script of M2 plus `fetch`) and compare with the original; (6) fetch a drillhole that the viewer imported and check the warning for the dropped log sets; (7) start `watch_selection`, change the selection in the window, check the callback runs once; (8) start a second program that sends and check it gets `BusyError` while the first program still selects and fetches. Add a closing paragraph: M3 checks the fake viewer assumptions OP12 (which types have a selection table), OP13 (one object per selection) and OP14 (parts arrive in order); record the time of a fetch of 1,000,000 cells. Pass line: all 8 steps give the expected result.

- [ ] T061 [P] [Finish] [S] README section "Select and fetch".
  - Files: `packages/loop_link/README.md`
  - Needs: T057
  - Wave: 17
  - Done when: `git grep -n -c "^## Select and fetch" packages/loop_link/README.md` gives 1.
  - Context: Add `## Select and fetch` before `## Test against the real codec`, with an example of about 10 lines: `from loop_link import Viewer, BlockModel`; `viewer = Viewer()`; `viewer.send(grid, id="grid")`; `viewer.select("grid", cells=[12, 13])`; `sel = viewer.selection()` then `sel.object_id`, `sel.cells`; `watcher = viewer.watch_selection(print, interval=0.5)` and `watcher.stop()`; `result = viewer.fetch("grid")` and `result.object`. Add 3 short sentences: select, read and fetch use the tool connection, so they work while another program holds the data connection; a viewer without geoviewer protocol 1.18 gives `MissingFeatureError`; a fetch of a drillhole keeps the first interval set and warns about dropped properties. ASD-STE100 English.

## Phase 9: Finish

- [ ] T062 [Finish] [S] Update `ROADMAP.md`: stage S5 done.
  - Files: `ROADMAP.md`
  - Needs: T056, T057, T058, T059, T060, T061
  - Wave: 18
  - Done when: `git grep -n -c "loop_link" ROADMAP.md` gives 3 or more and `git grep -n "stage S5" ROADMAP.md` shows 1 or more lines.
  - Context: In the entry `- [ ] **Stage 7 — `loop_link`: Python link to geoviewer.**` replace the sentence "Select, fetch, filter and colour wait for geoviewer spec 066." with "Stage S5 (select, read the selection and fetch; geoviewer spec 067) is done. Filter and colour (stage S4) are not planned yet." Add one status log entry at the end of "Status log", in the style of the neighbours, dated with the date of the run (YYYY-MM-DD): "Stage 7 stage S5 implemented: tool connection, select, read selection, fetch and round trip in `packages/loop_link`." Do not tick the Stage 7 box (S4 is open).

- [ ] T063 [Finish] [M] Run the merge gates for S5 and fix the findings.
  - Files: `packages/loop_link/`
  - Needs: T062, T038
  - Wave: 19
  - Done when: `ruff check LoopStructural packages && pytest` passes, `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests` passes, and `black --check packages/loop_link` passes.
  - Context: Constitution merge gate: `pytest` and `ruff check LoopStructural packages` must both pass. The root `pytest` does not run the package tests, so run the package command too. Run `black packages/loop_link` to fix the format (line length 100, quotes kept). Fix only problems in `packages/loop_link/`. Do not change a gate or add a Ruff ignore for `packages/loop_link/src`. Also run `packages/loop_link/tests/test_import.py` alone (SC-005: no Qt or extras at import). If `LOOP_LINK_GEOVIEWER_REPO` and `LOOP_LINK_GEOVIEWER_FIXTURES` are set on this machine, run `test_codec_real.py` and `test_fixtures.py` with them. Report the result of each command in the commit message.
## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001, T002, T004, T005, T039, T040, T060 | 5 | 2 |
| 2 | T003, T006, T041, T051 | 1 | 3 |
| 3 | T007, T009, T011, T016, T042, T052 | 0 | 6 |
| 4 | T008, T010, T017, T043, T058, T059 | 0 | 6 |
| 5 | T012, T015, T019 | 0 | 3 |
| 6 | T013, T018, T044 | 0 | 3 |
| 7 | T014, T020, T045 | 0 | 3 |
| 8 | T021, T046 | 0 | 2 |
| 9 | T023, T024, T025, T027, T034, T037, T047 | 3 | 4 |
| 10 | T028, T030, T032, T035, T048 | 0 | 5 |
| 11 | T029, T031, T033, T036, T049 | 1 | 4 |
| 12 | T022, T050 | 1 | 1 |
| 13 | T038, T053 | 0 | 2 |
| 14 | T054 | 0 | 1 |
| 15 | T055 | 0 | 1 |
| 16 | T056, T057 | 1 | 1 |
| 17 | T061 | 1 | 0 |
| 18 | T062 | 1 | 0 |
| 19 | T063 | 0 | 1 |

## Critical path

T001 -> T006 -> T007 -> T008 -> T019 -> T044 -> T045 -> T046 -> T047 -> T048 -> T049 -> T050 -> T053 -> T054 -> T055 -> T057 -> T061 -> T062 -> T063

## Dependencies

T001 to T038 (stages S1 to S3) are done. Stage S5 adds T039 to T063.

- S5a: errors (T039) and tool versions (T040) start at once. Frame builders follow (T041, T042). The fake
  viewer tool endpoint is tested (T043) and built in two steps, selection then fetch (T044, T045; both
  change `fake_viewer.py`, so they run in a chain after T019). `ToolConnection` is tested and built
  (T046, T047). Then selection tests, `selection.py` and the `Viewer` members (T048, T049, T050).
- S5b: the pure `decode_object` runs early and apart (T051, T052). The assembler tests, `fetch_object`
  and `Viewer.fetch` need the tool connection (T053, T054, T055). The round trip test closes it (T056).
- S5c: real codec cases (T058) and fixture checks (T059) need only the frame builders (T042), so they
  run in wave 4. The manual check M3 (T060) has no dependency. Exports (T057), README (T061),
  ROADMAP (T062) and the gate (T063) close the work.
- Open points: none blocks a task. OP1 and OP6 block S4, which has no tasks yet.
