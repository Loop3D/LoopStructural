# Tasks: loop_link, a Python link to geoviewer

**Feature**: `006-loop-link` | **Plan**: [plan.md](plan.md)

<!--
Scope: stages S1 to S3 of the plan (stories 1 and 5, built, T001 to T038), stage S5
(stories 2 to 4, T039 to T063) and stage S4 (stories 6 and 7, T064 to T092).
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
- Waves: 23. Tiers: [S] = simple model, [M] = standard model.
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

- [x] T047 [Found] [M] Implement `tool.py`.
  - Files: `packages/loop_link/src/loop_link/tool.py`
  - Needs: T046
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_tool.py` passes and `ruff check packages/loop_link` passes.
  - Context: Create `ToolConnection` and `raise_for_request_error(error: Mapping[str, Any], *, object_id: Optional[str]) -> NoReturn` (signatures and error table in T046; the `error` frame rows take `code`, `message` from the frame header). `request`: take `self._lock`; connect when `self._client is None` (`discovery.find_viewer(...)`, then `WebSocketClient.connect(host, port, token=..., timeout=..., path="/tool")`; no `hello`; no keepalive thread). Number `n` starts at 1. Send `build(n)`. Loop `client.receive(0.1)`; `protocol.decode_frame` (a `ProtocolError` closes the socket and raises `ViewerClosedError("The viewer sent a frame that is not valid. ...")`). `type == "error"` -> `raise_for_request_error`; a frame with another `request` is dropped; for `request == n`: `request.error` -> `raise_for_request_error(header["error"], object_id=object_id)` (for an `error` frame code `unknown_message` the feature text is the `feature` argument; for code `unsupported_version` the feature text is `"protocol " + str(header.get("supported"))`); `command.progress` -> `progress(fraction)` when `fraction` is a number; `on_frame(frame)` sees every other frame of the request; a type in `until` ends the loop and is returned. The wait restarts at each frame; after `timeout` (default `self._timeout`) with no frame close and raise `ViewerClosedError("The viewer did not answer in time. ...")`. When `cancel` is set: send `cancel_frame(n)` one time and keep reading until `cancelled` (raise `CancelledError`) or an `until` frame (return it). A socket error (`OSError`, `WebSocketClosed`, `WebSocketError`) sets the socket to `None` and raises `ViewerClosedError`. Docstrings and type hints.

- [x] T048 [US2] [M] Test: select and read the selection in `test_selection.py` (stories 2 and 3).
  - Files: `packages/loop_link/tests/test_selection.py`
  - Needs: T047, T045
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_selection.py` fails with ImportError of `loop_link.selection`, then passes after T050.
  - Context: Use `fake_viewer` and `viewer_factory`. Send a `BlockModel` (cells 4x4x4, id `grid`) first. API on `Viewer`: `select(object_id, *, cells=None, vertices=None, revision=None) -> dict[str, int]` (returns `{"cell": n, "vertex": m}`), `selection() -> Selection`, `watch_selection(callback, *, interval=0.5, initial=False) -> SelectionWatcher`, `tool_connected`. `loop_link.selection`: `Selection(object_id, revision, cells, vertices)` (frozen, `uint32` arrays, `.empty`, `.same_as(other)` compares object id and sorted cells and vertices, not the revision), `parse_selection_result(frame)` (an entry with a count but no buffer, a buffer not `uint32`, or a length not equal to `count` gives `ProtocolError`), `SelectionWatcher` (`running`, `error`, `stop()`, context manager). Tests: story 2 scenarios 1 to 4 (select cells 12 and 13: `fake.selection` holds them; index 64 gives `InvalidInputError` and the selection stays; `StaleRevisionError` for an old `revision`; `FakeViewer(protocol_minor=14)` gives `MissingFeatureError`); clear with `cells=[]`; `select()` with neither argument, a float index, a bool, a negative value or a value above 2**32-1 each give `InvalidInputError` that names the argument; `cells` and `vertices` together; duplicates count one time; 100,000 cells are one `selection.set` frame with a `u32` buffer of shape `[100000]` (`fake.tool_requests`, SC-003); `kind="headless"` gives `MissingFeatureError`; story 3: `selection()` after `set_user_selection` returns object id, cells, vertices; a 100,000-cell selection returns a `uint32` array; empty selection has `object_id is None`; `watch_selection` calls the callback once per change at the next poll (`interval=0.05`), not for no change, `initial=True` calls once with a non-empty first selection, `interval<=0` gives `InvalidInputError`, `stop()` ends the thread, no thread named `loop-link-selection-poll` exists before the call, an exception in the callback sets `.error` and ends the thread, a dropped `Viewer` ends the thread. A `Viewer` that only selects does not open the data connection (`viewer.connected` is false). `close()` stops watchers and a later `selection()` raises `ViewerClosedError`.

- [x] T049 [US3] [M] Implement `selection.py`.
  - Files: `packages/loop_link/src/loop_link/selection.py`
  - Needs: T048
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_selection.py -k "parse or same_as or empty or dtype"` passes and `ruff check packages/loop_link` passes. The watcher tests pass after T050.
  - Context: Create `Selection` (`@dataclass(frozen=True, eq=False)`: `object_id: Optional[str]`, `revision: Optional[int]`, `cells: np.ndarray`, `vertices: np.ndarray`; property `empty`; `same_as`), `parse_selection_result(frame: Frame) -> Selection` (kind `cells` -> `cells`, kind `nodes` -> `vertices`; array is `frame.buffers[entry["indices"]]`; checks as in T048; `ProtocolError` from `loop_link.protocol`) and `SelectionWatcher`. The watcher is made only by `Viewer.watch_selection`: `SelectionWatcher(viewer_ref: weakref.ref, callback, interval, last: Selection)`, a daemon thread named `loop-link-selection-poll`. Loop: wait `interval` on a `threading.Event` (so `stop` is fast); get the viewer from the weak reference (end the thread if it is gone); call `viewer.selection()`; when the result is not `same_as` the last one, call `callback(selection)` one time and keep it. An exception of the poll or of the callback goes to `error` and the thread ends. `stop()` sets the event and joins with timeout `2 * interval + 1`. `__enter__` returns self, `__exit__` calls `stop`. No thread before `start`. Docstrings and type hints.

- [x] T050 [US2] [M] Add `select`, `selection`, `watch_selection` and the tool connection to `Viewer`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T049, T047, T042
  - Wave: 12
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_selection.py packages/loop_link/tests/test_viewer.py packages/loop_link/tests/test_push.py` passes and `ruff check packages/loop_link` passes.
  - Context: New members as in T048. The `Viewer` makes its `ToolConnection` (from `loop_link.tool`) at the first tool call, with its own `port`, `token`, `host`, `settings_dir`, `timeout`; tool calls never open the data connection. `release_data` keeps the tool connection. `close` stops all watchers, closes the tool connection and sets `_dead`; a later tool call raises `ViewerClosedError`. `select`: `object_id` must be a non-empty `str`; at least one of `cells`, `vertices` not `None`, else `InvalidInputError("Give cells or vertices. Use an empty list to clear the selection.")`; each given array is a 1-D sequence of integers from 0 to 2**32-1 (no float, no bool, no negative, check `np.asarray`, dtype kind `i` or `u`, and bool dtype refused), else `InvalidInputError` that names the argument. One `selection.set` frame (`encode.selection_set_frame`), `until=("selection.result",)`, `feature="selection and fetch (geoviewer protocol 1.18)"`; return `{"cell": n, "vertex": m}` with the counts of the reply (both keys are always present; 0 for a kind that the reply does not list). `selection()`: `selection.get` -> `parse_selection_result` -> `Selection` with `object_id`/`revision` from the header (`None` when absent). `watch_selection(callback, *, interval=0.5, initial=False)`: `interval <= 0` gives `InvalidInputError`; read the selection in the calling thread first (errors go to the caller); call the callback once when `initial` and the selection is not empty; then start `SelectionWatcher` with a `weakref.ref(self)`; keep it in a list. Do not change the data-connection code.

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

- [x] T053 [US4] [M] Test: assembler, `fetch_object` and `Viewer.fetch` in `test_fetch.py`.
  - Files: `packages/loop_link/tests/test_fetch.py`
  - Needs: T052, T050
  - Wave: 13
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py` fails for the new tests (AttributeError on `fetch_object`, `_Assembler` and `Viewer.fetch`), then passes after T055.
  - Context: Names fixed by this task list: `loop_link.fetch._Assembler(object_id: str)` with `feed(frame: Frame) -> None` and `finish() -> FetchResult`; `fetch_object(tool, id, *, progress=None, cancel=None, timeout=30.0) -> FetchResult`; `Viewer.fetch(id, *, progress=None, cancel=None) -> FetchResult`. Name the Viewer tests `test_viewer_fetch_*`. Assembler tests feed hand-built frames (`Frame(header, buffers)`): first frame not `object.header` -> `ProtocolError`; arrays with the sum of bytes over 2 GiB (header only, no allocation) -> `ObjectTooLargeError`; part for an array not in `arrays`; `offset` not a multiple of 8; part past the end; a gap; an overlap; a missing part at `object.end`; wrong `parts` count; `part` buffer not `u8` -> each `ProtocolError`; a buffer name in `geometry`, `data[*].buffer` or `offsets` that is not in `arrays` -> `ProtocolError`; a good stream gives a `FetchResult` with the revision of `header["object"]["revision"]`. `fetch_object` and `Viewer.fetch` with the fake (`fake_viewer`, `viewer_factory`): fetch of a sent `Surface` equals the input; `FakeViewer(part_bytes=64)` gives several parts and `progress` calls that rise to 1; `FakeViewer(fetch_delay=0.05, part_bytes=64)` with `cancel` set gives `CancelledError` and `fake.objects` is the same; unknown id -> `InvalidInputError`; `kind="headless"` fake still fetches; a stream error closes the tool connection and raises `ViewerClosedError("The viewer sent an object that is not valid: ...")`; the revision is in the result; FR-014: program A holds the data connection (`send`), program B (`viewer_factory()`) gets `BusyError` from `send`, but `select`, `selection` and `fetch` work; after A calls `release_data`, B sends and the objects of A stay (`fake.objects`).

- [x] T054 [US4] [M] Implement `_Assembler` and `fetch_object` in `fetch.py`.
  - Files: `packages/loop_link/src/loop_link/fetch.py`
  - Needs: T053
  - Wave: 14
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py -k "not test_viewer_fetch"` passes and `ruff check packages/loop_link` passes.
  - Context: `fetch_object` calls `tool.request(lambda n: encode.fetch_frame(n, id), until=("object.end",), feature="selection and fetch (geoviewer protocol 1.18)", object_id=id, on_frame=assembler.feed, progress=progress, cancel=cancel, timeout=timeout)`, then returns `assembler.finish()`. It catches `ProtocolError`, calls `tool.close()` and raises `ViewerClosedError("The viewer sent an object that is not valid: <text>")`; other library errors go up. `_Assembler` rules: (1) first frame must be `object.header`; (2) header checks: each `arrays` entry has a unique `name`, a `dtype` in `f32 f64 u32 i32 u8`, a `shape` of non-negative integers; the sum of bytes over `protocol.MAX_OBJECT_BYTES` -> `ObjectTooLargeError(id, size, limit)`; each buffer name in `object.geometry`, `object.data[*].buffer` and `offsets` is in `arrays`; then one `numpy.empty(nbytes, dtype=numpy.uint8)` per array; (3) `object.part`: array known; `offset % 8 == 0`; buffer `part` is `u8`; `offset` equals the bytes already filled in that array (in order, no gap, no overlap); `offset + bytes` at most the array size; copy in; (4) `object.end`: `parts` equals the count received and each array is full; (5) `finish`: view each array with its dtype and shape (`.view`, `.reshape`), call `decode_object(header["object"], arrays, header["origin"])`; fill `FetchResult` (`origin` as tuple or None, `crs`, `name`, `parent`, `revision`, `object_type`). Every break raises `ProtocolError`.

- [x] T055 [US4] [M] Add `Viewer.fetch`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T054, T050
  - Wave: 15
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fetch.py packages/loop_link/tests/test_selection.py packages/loop_link/tests/test_viewer.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add `Viewer.fetch(self, id: str, *, progress: Optional[Callable[[float], None]] = None, cancel: Optional[threading.Event] = None) -> FetchResult`. It uses the tool connection of T050 (same lazy creation, same closed rule) and calls `fetch.fetch_object(tool, id, progress=progress, cancel=cancel)`. It never opens the data connection. An unknown id gives `InvalidInputError` (the tool layer maps `invalid_input`). `progress` gets floats from 0 to 1. Import `FetchResult` and `fetch_object` from `loop_link.fetch`. Numpy-style docstring with Parameters, Returns, Raises (`InvalidInputError`, `CancelledError`, `MissingFeatureError`, `ViewerClosedError`, `UnsupportedObjectError`).

- [x] T056 [US4] [M] Test: send then fetch for the seven types in `test_roundtrip.py` (SC-002).
  - Files: `packages/loop_link/tests/test_roundtrip.py`, `packages/loop_link/src/loop_link/fetch.py`
  - Needs: T055
  - Wave: 16
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_roundtrip.py` passes.
  - Context: Use `fake_viewer` and `viewer_factory` (`FakeViewer(origin=(500000.0, 7000000.0, 0.0))` for the origin cases: build it in the test and write its instance file as `conftest.py` does). For each of `Points`, `Vectors` (made with `Vectors.from_directions`), `Lines`, `Surface`, `BlockModel`, `TetraMesh`, `Drillhole`: `viewer.send(obj, id=...)`, then `viewer.fetch(id).object`. Compare: floats with `np.testing.assert_allclose(rtol=1e-6)`; integers, booleans, text, `Categorical` (codes, names, colours) and index arrays with `np.array_equal`; the dtype of each property stays (float32 stays float32); coordinates are float64 world coordinates, also for a point 10 km from the origin. Vectors: compare `orientation` angles to the input within 1e-6 degrees (the viewer keeps angles; an x,y,z input is equal up to the conversion). The test file is the task. Change `fetch.py` only if a case shows a bug in `decode_object` or the assembler; then say so in the task result.

## Phase 8: Stage S5c, real codec, fixtures, manual check, docs

- [x] T057 [Found] [S] Export the S5 public names.
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

- [x] T061 [P] [Finish] [S] README section "Select and fetch".
  - Files: `packages/loop_link/README.md`
  - Needs: T057
  - Wave: 17
  - Done when: `git grep -n -c "^## Select and fetch" packages/loop_link/README.md` gives 1.
  - Context: Add `## Select and fetch` before `## Test against the real codec`, with an example of about 10 lines: `from loop_link import Viewer, BlockModel`; `viewer = Viewer()`; `viewer.send(grid, id="grid")`; `viewer.select("grid", cells=[12, 13])`; `sel = viewer.selection()` then `sel.object_id`, `sel.cells`; `watcher = viewer.watch_selection(print, interval=0.5)` and `watcher.stop()`; `result = viewer.fetch("grid")` and `result.object`. Add 3 short sentences: select, read and fetch use the tool connection, so they work while another program holds the data connection; a viewer without geoviewer protocol 1.18 gives `MissingFeatureError`; a fetch of a drillhole keeps the first interval set and warns about dropped properties. ASD-STE100 English.

## Phase 9: Finish

- [x] T062 [Finish] [S] Update `ROADMAP.md`: stage S5 done.
  - Files: `ROADMAP.md`
  - Needs: T056, T057, T058, T059, T060, T061
  - Wave: 18
  - Done when: `git grep -n -c "loop_link" ROADMAP.md` gives 3 or more and `git grep -n "stage S5" ROADMAP.md` shows 1 or more lines.
  - Context: In the entry `- [ ] **Stage 7 — `loop_link`: Python link to geoviewer.**` replace the sentence "Select, fetch, filter and colour wait for geoviewer spec 066." with "Stage S5 (select, read the selection and fetch; geoviewer spec 067) is done. Filter and colour (stage S4) are not planned yet." Add one status log entry at the end of "Status log", in the style of the neighbours, dated with the date of the run (YYYY-MM-DD): "Stage 7 stage S5 implemented: tool connection, select, read selection, fetch and round trip in `packages/loop_link`." Do not tick the Stage 7 box (S4 is open).

- [!] T063 [Finish] [M] Run the merge gates for S5 and fix the findings. Failed: the package gates pass (ruff, black, 660 package tests twice, import test, real codec and fixtures 45 passed) but the root `pytest` has no result: tests/integration/test_fold_models.py is extremely slow here and a tests/unit run was stopped by the system for low memory.
  - Files: `packages/loop_link/`
  - Needs: T062, T038
  - Wave: 19
  - Done when: `ruff check LoopStructural packages && pytest` passes, `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests` passes, and `black --check packages/loop_link` passes.
  - Context: Constitution merge gate: `pytest` and `ruff check LoopStructural packages` must both pass. The root `pytest` does not run the package tests, so run the package command too. Run `black packages/loop_link` to fix the format (line length 100, quotes kept). Fix only problems in `packages/loop_link/`. Do not change a gate or add a Ruff ignore for `packages/loop_link/src`. Also run `packages/loop_link/tests/test_import.py` alone (SC-005: no Qt or extras at import). If `LOOP_LINK_GEOVIEWER_REPO` and `LOOP_LINK_GEOVIEWER_FIXTURES` are set on this machine, run `test_codec_real.py` and `test_fixtures.py` with them. Report the result of each command in the commit message.

## Phase 10: Stage S4, filter and colour (stories 6 and 7)

Test command for this phase: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/<file>`. All S4 tasks use the tool connection of S5 (`loop_link.tool.ToolConnection`, method `request`). Waves follow the shared graph: a task that changes a file that an earlier task changed `Needs` that task.

### Stage S4a: frames, commands module, fake commands endpoint

- [x] T064 [Found] [S] Add `COMMAND_PROTOCOL_VERSION` to `protocol.py`.
  - Files: `packages/loop_link/src/loop_link/protocol.py`
  - Needs: T040
  - Wave: 2
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "from loop_link import protocol; assert protocol.COMMAND_PROTOCOL_VERSION == {'major': 1, 'minor': 11}"` passes and `ruff check packages/loop_link` passes.
  - Context: In `protocol.py`, after the line `CANCEL_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}`, add the line `COMMAND_PROTOCOL_VERSION: dict[str, int] = {"major": 1, "minor": 11}` and a one-line comment above it: the version of the `commands.list` and `command.run` frames (the fixtures of the geoviewer repository use 1.11, so a viewer of 1.11 to 1.17 accepts them). Change nothing else.

- [x] T065 [Found] [S] Test: the frames `commands_list_frame` and `command_frame` in `test_commands.py`.
  - Files: `packages/loop_link/tests/test_commands.py`
  - Needs: T064
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_commands.py` fails with ImportError of `commands_list_frame`, then passes after T066.
  - Context: Create the file (module docstring, `from __future__ import annotations`). Under test, in `loop_link.encode`: `commands_list_frame(request: int) -> bytes` and `command_frame(request: int, command: str, input: Mapping[str, Any]) -> bytes`. Decode with `loop_link.protocol.decode_frame`. Expected headers (the geoviewer fixtures, decoded): `commands_list_frame(1)` gives `{"type": "commands.list", "version": {"major": 1, "minor": 11}, "request": 1, "buffers": []}`; `command_frame(2, "object.set_visible", {"id": "p1", "visible": False})` gives `{"type": "command.run", "version": {"major": 1, "minor": 11}, "request": 2, "command": "object.set_visible", "input": {"id": "p1", "visible": False}, "buffers": []}`. Tests: the two header equalities; the frame has no buffers (`frame.buffers == {}`); `command_frame` does not change the input mapping (pass a dict, compare after); a nested input (`{"object": "g", "edit": {"kind": "set_values", "array": "a", "values": [0, 2]}}`) round trips; a NaN in the input raises `ProtocolError` (`allow_nan=False`).

- [x] T066 [Found] [S] Implement `commands_list_frame` and `command_frame` in `encode.py`.
  - Files: `packages/loop_link/src/loop_link/encode.py`
  - Needs: T065, T042
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_commands.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add after `cancel_frame` (same style, numpy-style docstring with Parameters and Returns): `def commands_list_frame(request: int) -> bytes` returning `protocol.encode_frame({"type": "commands.list", "request": request}, version=protocol.COMMAND_PROTOCOL_VERSION)`, and `def command_frame(request: int, command: str, input: Mapping[str, Any]) -> bytes` returning `protocol.encode_frame({"type": "command.run", "request": request, "command": command, "input": dict(input)}, version=protocol.COMMAND_PROTOCOL_VERSION)`. `encode_frame` adds `buffers`. Add both names to `__all__` of `encode.py` (keep it sorted). Import `Mapping` from `collections.abc` if it is not imported (the file has `from __future__ import annotations`).

- [x] T067 [Found] [M] Test: the commands module in `test_commands.py`.
  - Files: `packages/loop_link/tests/test_commands.py`
  - Needs: T066
  - Wave: 5
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_commands.py` fails with ImportError of `loop_link.commands`, then passes after T068.
  - Context: Add tests to the file of T065. Under test, in `loop_link.commands`: `ColourSettings(array=None, colour_map=None, range=None, scale=None, plain=None)` (frozen dataclass; `__post_init__` checks each field that is not None: `array` and `colour_map` by `check_text`; `range` a pair of numbers with minimum < maximum; `scale` `"linear"` or `"log"`; `plain` by `check_rgb`; each error is `InvalidInputError` that names the field); `format_number(value) -> str` (`repr(float(value))`: `0.5` gives `"0.5"`, `1e-05` gives `"1e-05"`, `3` gives `"3.0"`); `check_number(value, name) -> float` (accepts int, float, numpy numbers; refuses `bool`, NaN, infinity, a str, None); `check_text(value, name) -> str` (non-empty str); `check_rgb(value, name) -> tuple[int, int, int]` (three ints 0 to 255, not bool; refuses 256, -1, 2 values, floats); `check_codes(values, name) -> list[int]` (sequence or array of ints 0 to 2**32-1; refuses bool, float, negative, 2**32; returns a `list[int]`); `raise_for_command_result(result, *, command)` (`{"status": "done"}` returns None; `{"status": "refused", "text": "x"}` raises `InvalidInputError` with the text "The viewer refused the command 'filter.edit': x"; another status raises `ProtocolError` with "The viewer sent an unknown result status 'weird'."); `parse_colour_map_names(frame) -> tuple[str, ...] | None` with a `Frame(header, {})` of type `commands`: an entry `{"name": "colour.set_map", "input_schema": {"properties": {"name": {"type": "string", "enum": ["viridis", "grey", "coolwarm"]}}}}` gives the 3 names; no `enum` gives None; no entry gives None; an `enum` that is not a list of str gives `ProtocolError`; `run_command(tool, command, input, *, object_id=None, timeout=None, cancel=None) -> Mapping`: use a small stub class with a method `request(build, *, until, feature, object_id=None, on_frame=None, progress=None, cancel=None, timeout=None)` that records its arguments and returns `Frame({"type": "command.result", "request": 1, "result": {"status": "done"}}, {})`; check `until == ("command.result",)`, `feature == "filter and colour commands (geoviewer protocol 1.11)"`, that `decode_frame(build(7)).header` has `request` 7, `command` and `input`, that the result mapping is returned, that a `refused` result raises, and that a result that is not a mapping with a `status` str raises `ProtocolError`; `list_colour_map_names(tool)`: the stub returns a `commands` frame; check `until == ("commands",)` and the names.

- [x] T068 [Found] [M] Implement `commands.py`.
  - Files: `packages/loop_link/src/loop_link/commands.py`
  - Needs: T067
  - Wave: 6
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_commands.py` passes and `ruff check packages/loop_link` passes.
  - Context: Create the module with the names and rules of T067 (docstrings and type hints on all public names; `__all__`). `ColourSettings` is `@dataclass(frozen=True)` with fields `array: Optional[str] = None`, `colour_map: Optional[str] = None`, `range: Optional[tuple[float, float]] = None`, `scale: Optional[str] = None`, `plain: Optional[tuple[int, int, int]] = None`; `__post_init__` stores the checked values with `object.__setattr__` (`range` as a tuple of two floats). `run_command` calls `tool.request(lambda n: encode.command_frame(n, command, input), until=("command.result",), feature="filter and colour commands (geoviewer protocol 1.11)", object_id=object_id, cancel=cancel, timeout=timeout)`, reads `frame.header.get("result")`, raises `ProtocolError` when it is not a mapping with a str `status`, calls `raise_for_command_result`, returns the mapping. `list_colour_map_names(tool)` calls `tool.request(lambda n: encode.commands_list_frame(n), until=("commands",), feature="filter and colour commands (geoviewer protocol 1.11)")` and returns `parse_colour_map_names(frame)`. `parse_colour_map_names` looks in `frame.header["commands"]` (a list) for `name == "colour.set_map"` and reads `input_schema["properties"]["name"]["enum"]`. `raise_for_command_result`: `refused` gives `InvalidInputError(f"The viewer refused the command '{command}': {text}")`; another status that is not `done` gives `ProtocolError(f"The viewer sent an unknown result status '{status}'.")`. Imports: `loop_link.encode`, `loop_link.errors.InvalidInputError`, `loop_link.protocol.Frame` and `ProtocolError`, `loop_link.tool.ToolConnection` (type only). No pyvista, pandas or Qt.

- [x] T069 [P] [Found] [M] Test: the commands endpoint of the fake viewer in `test_fake_viewer.py`.
  - Files: `packages/loop_link/tests/test_fake_viewer.py`
  - Needs: T066, T043
  - Wave: 5
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py -k "fake_command or fake_filter or fake_colour"` fails (TypeError on the new constructor arguments or missing replies), then passes after T071. The old tests of the file still pass.
  - Context: Test the rules below through `FakeViewer.tool_frame(bytes) -> list[Frame]` with frames from `encode.command_frame` and `encode.commands_list_frame`, as the selection tests of this file do (read them first, and `_upsert` in this file for how to put an object in the fake). Name the tests `test_fake_command_*`, `test_fake_filter_*`, `test_fake_colour_*`. New constructor arguments: `filter_delay: float = 0.0` (seconds that the scene stays busy after a `filter.edit` that changes a filter) and `busy_wait: float = 5.0` (seconds that a waiting command waits before `request.error` `busy`; tests use 0.05). New read-only properties that return copies: `.filters` (`dict[object id, dict[array name, dict]]`, keys `kind` `"range"` or `"values"`, `on`, `minimum`, `maximum`, `keep_nan`, `values`, `exclude`) and `.colouring` (`dict[object id, dict]`, keys `array`, `colour_map`, `range`, `scale`, `plain`). Rules: `commands.list` answers `commands` with `viewer` = the kind and one entry for each of the 7 commands `filter.edit`, `colour.set_active`, `colour.set_map`, `colour.set_range`, `colour.clear_range`, `colour.set_scale`, `object.set_colour` (and no other; entry keys `name`, `description`, `needs_window` true only for `colour.set_scale`, `tool_available`, `unavailable_reason`, `input_schema`, `result_schema`); `input_schema["properties"]["name"]["enum"]` of `colour.set_map` is `["viridis", "grey", "coolwarm"]` when `protocol_minor >= 18`, and the property has only `type: "string"` below 18; `protocol_minor < 11` answers `error` `unknown_message` to both messages. `command.run`: an unknown command, an unknown or missing input field, or a wrong field type gives `request.error` `invalid_input` with that `field`; an unknown object gives `invalid_input` with field `id` (`object` for `colour.set_scale` and `filter.edit`); `kind="headless"` gives `needs_window` for `colour.set_scale` only. Reply `{"type": "command.result", "request": n, "result": {"status": "done"}}`. Inputs: `filter.edit` `{"object", "edit": {"kind": ...}}` with `add` (`array`), `set_range` (`array`, `min` str, `max` str, `keep_nan` bool), `set_values` (`array`, `values` list of u32), `set_on` (`array`, `on`), `set_exclude` (`array`, `exclude`), `remove` (`array`), `cancel` (nothing else); `colour.set_active` `{"id", "array"?}`; `colour.set_map` `{"id", "name"}`; `colour.set_range` `{"id", "min" str, "max" str}`; `colour.clear_range` `{"id"}`; `colour.set_scale` `{"object", "scale"}` (`linear` or `log`); `object.set_colour` `{"id", "colour"}` (list of 3 ints or null). Effects: `add` makes an off filter; `set_*` change it (a missing filter gives `refused`); `remove` deletes it; the colour commands fill `.colouring`; a later `scene.upsert` of the same id keeps `.filters` and `.colouring`, and drops the filter of an array that the new revision does not have. `refused` (`{"status": "refused", "text": ...}` with the array name in the text): a filter on an array that is missing, `text`, or of association `object`; a range edit on a `referenced` or `boolean` array; a values edit on a number array; `colour.set_active` with an unknown array; `colour.set_map`, `colour.set_range` and `colour.set_scale` when the active array is `referenced` or `boolean`. Busy: while `filter_delay` has not passed since a filter change, each command except `colour.set_scale` waits up to `busy_wait` s, then runs or gets `request.error` `busy`; `colour.set_scale` gets `busy` at once; `filter.edit` with `kind: "cancel"` ends the busy time.

- [x] T070 [Found] [M] Fake viewer: commands endpoint, list and colour part.
  - Files: `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T069, T045
  - Wave: 8
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py -k "fake_command or fake_colour"` passes and `ruff check packages/loop_link` passes.
  - Context: Add the constructor arguments `filter_delay` and `busy_wait` (defaults and meaning in T069; add them to the class docstring), the properties `.filters` and `.colouring` (copies, read under `self._state_lock`), and the commands endpoint inside `_tool_reply` (the types `command.run` and `commands.list`; only when `protocol_minor >= 11`, else the existing `unknown_message` reply). Use `_request_error(request, code, message, field_name)` for errors. Build the replies as `Frame(header={"type": "command.result", "request": request, "result": {...}}, buffers={})` and `Frame({"type": "commands", "request": request, "viewer": self.kind, "commands": [...]}, {})`. Implement here: `commands.list`; the validation of the input of all 7 commands (unknown command, unknown or missing field, wrong type, unknown object); the `needs_window` rule; the 5 colour commands with their `refused` rules and their effects in `.colouring`; the keep rule for `scene.upsert` (store the colouring and the filters apart from `self.objects`; in `_upsert` drop only the filters whose array is not in the new header `data`). For `filter.edit`, call a method `_filter_edit(object_id, edit) -> dict` that returns `{"status": "done"}` and does nothing yet (T071 fills it). Put the busy rule of T069 in one helper `_wait_if_busy(command) -> bool` that T071 also uses (it returns False when the reply must be `request.error` `busy`). Data entries have `id`, `name`, `data_type`, `association`. Docstrings and type hints.

- [x] T071 [Found] [M] Fake viewer: commands endpoint, filter part.
  - Files: `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T070
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py` passes and `ruff check packages/loop_link` passes.
  - Context: Fill `_filter_edit` of T070 with the rules of T069 for the edits `add`, `set_range`, `set_values`, `set_on`, `set_exclude`, `remove` and `cancel`, and their effects in `.filters` (a filter is a dict with keys `kind`, `on`, `minimum`, `maximum`, `keep_nan`, `values`, `exclude`; `add` makes `{"kind": None, "on": False, ...}`; the first `set_range` or `set_values` sets the kind; `minimum` and `maximum` are stored as floats made from the texts `min` and `max`). The property rules (refusals with a text that names the array) are in T069. Start the busy time (`time.monotonic() + filter_delay`) after each edit that changes a filter; `cancel` clears it. Use `_wait_if_busy` of T070 at the start of each command. All old tests of the file must still pass.

### Stage S4b: filter and colour members of `Viewer`

- [x] T072 [P] [US1] [M] Test: filter calls in `test_filter.py` (story 6).
  - Files: `packages/loop_link/tests/test_filter.py`
  - Needs: T071
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_filter.py` fails with AttributeError on `filter_range`, then passes after T074.
  - Context: Use `fake_viewer` and `viewer_factory` (`fake_viewer.viewer` is the `FakeViewer`; use its `.filters`, `.objects`, `.tool_requests`). First send a `BlockModel` (cells 4x4x4, id `grid`) with `cell_properties={"grade": float array, "unit": Categorical(...)}` and `properties={"t": np.array(["a"] * 125)}`. API on `Viewer`: `filter_range(object_id, array, minimum, maximum, *, keep_nan=False, wait=True, cancel=None, timeout=60.0)`, `filter_values(object_id, array, values, *, exclude=False, wait=True, cancel=None, timeout=60.0)`, `filter_enable(object_id, array, on=True)`, `filter_remove(object_id, array)`, `cancel_filter(object_id)`. In this task `wait` has no effect (S4d). Tests: story 6 scenario 1 (`filter_range("grid", "grade", 0.2, 0.8, keep_nan=True)` gives in `fake.filters["grid"]["grade"]` kind `range`, `on` True, minimum 0.2, maximum 0.8, `keep_nan` True; the edits arrive in the order `add`, `set_range`, `set_on` in `fake.tool_requests`; `min` and `max` are the texts `"0.2"` and `"0.8"`); scenario 2 (`filter_values("grid", "unit", [0, 2], exclude=True)` gives kind `values`, `values` [0, 2], `exclude` True, `on` True; edits `add`, `set_values`, `set_exclude`, `set_on`); scenario 3 (`filter_enable(..., False)` gives `on` False; `filter_remove` deletes the filter; the buffers of `fake.objects["grid"]` equal the sent buffers before and after: the data does not change); scenario 4 (an array that does not exist, and `"t"` (text), each give `InvalidInputError` whose text has the array name; `filter_range` on `unit` and `filter_values` on `grade` give `InvalidInputError`; `minimum > maximum` gives `InvalidInputError` and `fake.tool_requests` has no `command.run`; values with a float, a bool, a negative number or 2**32 give `InvalidInputError`); scenario 6 (send `grid` again as revision 2: `fake.filters` is unchanged); an `add` that is refused for an existing filter: monkeypatch `FakeViewer._filter_edit` so that kind `add` returns `{"status": "refused", "text": "exists"}`, then `filter_range` still succeeds; `cancel_filter("grid")` sends one `filter.edit` with edit `{"kind": "cancel"}`; `FakeViewer(kind="headless")` still filters; `FakeViewer(protocol_minor=10)` gives `MissingFeatureError`; FR-014 tool part: program A (`viewer_factory()`) sends, program B (`viewer_factory()`) gets `BusyError` from `send`, but B can call `filter_range` and `filter_remove`; none of these opens the data connection (`B.connected` is false).

- [x] T073 [P] [US1] [M] Test: colour calls in `test_colour.py` (story 7).
  - Files: `packages/loop_link/tests/test_colour.py`
  - Needs: T071
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_colour.py` fails with AttributeError on `set_colour_property`, then passes after T075.
  - Context: Use `fake_viewer` and `viewer_factory`. First send a `BlockModel` (cells 4x4x4, id `grid`) with `cell_properties={"grade": float array, "unit": Categorical(...), "flag": bool array}`. API on `Viewer` (and `loop_link.commands.ColourSettings` for `apply_colour`): `colour_map_names() -> tuple[str, ...]`, `set_colour_property(object_id, array: Optional[str])`, `set_colour_map(object_id, name)`, `set_colour_range(object_id, minimum, maximum)`, `clear_colour_range(object_id)`, `set_colour_scale(object_id, scale)`, `set_object_colour(object_id, colour: Optional[tuple[int, int, int]])`, `apply_colour(object_id, settings)`. Check the state in `fake.colouring["grid"]` (keys `array`, `colour_map`, `range`, `scale`, `plain`). Tests: scenario 1 (`set_colour_property("grid", "grade")` gives array `"grade"`; `array=None` leaves the key `array` out of the input of `colour.set_active` (read `fake.tool_requests`) and the fake shows no array); scenario 2 (`set_colour_map("grid", "coolwarm")`; `colour_map_names()` is `("viridis", "grey", "coolwarm")`); scenario 3 (`set_colour_range("grid", 0.1, 0.9)` sends `min` `"0.1"` and `max` `"0.9"`; `clear_colour_range` clears it; `minimum >= maximum` gives `InvalidInputError` and no frame goes out); scenario 4 (`set_colour_scale("grid", "log")`; `"sqrt"` gives `InvalidInputError` and no frame; the input key is `object`, not `id`); `set_object_colour("grid", (255, 0, 0))` sends `[255, 0, 0]` and `None` sends null; `(256, 0, 0)` gives `InvalidInputError`; scenario 5 (after `set_colour_property("grid", "unit")`, each of `set_colour_map`, `set_colour_range` and `set_colour_scale` gives `InvalidInputError` that has the viewer text; the same for `flag`); scenario 6 (`set_colour_map("grid", "nope")` gives `InvalidInputError` whose text lists `viridis`, `grey` and `coolwarm`, and the fake received no `colour.set_map` (`fake.tool_requests`); a fake with `protocol_minor=17` (no names) gets the name and decides: it accepts `"nope"`; `colour_map_names()` with `protocol_minor=17` gives `MissingFeatureError`); `apply_colour(id, ColourSettings(array="grade", colour_map="viridis", range=(0.0, 1.0), scale="linear", plain=(1, 2, 3)))` sets all fields, and the 5 commands arrive in this order: `colour.set_active`, `colour.set_map`, `colour.set_range`, `colour.set_scale`, `object.set_colour`; a `ColourSettings` with only `scale` sends 1 command; a bad field of `ColourSettings` gives `InvalidInputError` before any frame; `FakeViewer(kind="headless")`: `set_colour_scale` gives `MissingFeatureError` and the other colour calls work; an unknown object gives `InvalidInputError`.

- [x] T074 [US1] [M] Add the filter members to `Viewer`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T072, T068, T055
  - Wave: 16
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_filter.py packages/loop_link/tests/test_selection.py packages/loop_link/tests/test_viewer.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add the 5 methods with the signatures of T072 (numpy-style docstrings with Parameters and Raises). They use `self._tool()` (as `select` does) and `commands.run_command(tool, "filter.edit", {"object": object_id, "edit": {...}}, object_id=object_id)`; they never open the data connection. Check the arguments before any frame: `object_id` and `array` are non-empty `str` (`check_text`); `minimum` and `maximum` by `check_number`, and `minimum <= maximum` else `InvalidInputError`; `values` by `check_codes`; `exclude`, `keep_nan` and `on` are `bool`. `filter_range` sends `{"kind": "add", "array": a}` (if `run_command` raises `InvalidInputError` for `add`, go on: the filter can exist), then `{"kind": "set_range", "array": a, "min": format_number(minimum), "max": format_number(maximum), "keep_nan": keep_nan}`, then `{"kind": "set_on", "array": a, "on": True}`. `filter_values` sends `add` (same rule), `{"kind": "set_values", "array": a, "values": codes}`, `{"kind": "set_exclude", "array": a, "exclude": exclude}`, then `set_on` true. `filter_enable`: one `set_on`. `filter_remove`: one `{"kind": "remove", "array": a}`. `cancel_filter`: one `{"kind": "cancel"}`. The first failure stops the call and goes up. `wait`, `cancel` and `timeout` are accepted and checked (`timeout > 0`) but not used yet (T081). `commands.py` does not import `viewer.py`, so a module-level import has no cycle.

- [x] T075 [US1] [M] Add the colour members and `apply_colour` to `Viewer`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T073, T074
  - Wave: 17
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_colour.py -k "not send"` passes and `ruff check packages/loop_link` passes.
  - Context: Add the methods of T073 (`colour_map_names`, `set_colour_property`, `set_colour_map`, `set_colour_range`, `clear_colour_range`, `set_colour_scale`, `set_object_colour`, `apply_colour`). Commands and inputs: `colour.set_active` `{"id", "array"}` (leave `array` out for `None`); `colour.set_map` `{"id", "name"}`; `colour.set_range` `{"id", "min": format_number(minimum), "max": format_number(maximum)}` with `minimum < maximum` else `InvalidInputError`; `colour.clear_range` `{"id"}`; `colour.set_scale` `{"object": id, "scale"}` with `scale` in `("linear", "log")`; `object.set_colour` `{"id", "colour": [r, g, b] or None}` (`check_rgb`). `colour_map_names()` calls `commands.list_colour_map_names(tool)` and keeps the tuple in `self._colour_map_names` for the life of the tool connection (reset it in `close`); a `None` result gives `MissingFeatureError("colour map names (geoviewer protocol 1.18)")`. `set_colour_map`: `name` by `check_text`; read the names with a private helper that returns `None` instead of raising when the viewer has none; when names exist and `name` is not one of them, raise `InvalidInputError(f"The viewer has no colour map '{name}'. The names are: {', '.join(names)}.")` before any `colour.set_map`. `apply_colour(object_id, settings)` sends, in this order and only for the fields that are not None: `colour.set_active`, `colour.set_map` (with the name check), `colour.set_range`, `colour.set_scale`, `object.set_colour`. A `settings` that is not a `ColourSettings` gives `InvalidInputError`. All use `self._tool()` and `run_command`; none opens the data connection.

- [x] T076 [US1] [M] Test: `send(..., colour=...)` in `test_colour.py` (FR-021).
  - Files: `packages/loop_link/tests/test_colour.py`
  - Needs: T075
  - Wave: 18
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_colour.py -k send` fails (the `colour` argument is passed to the adapter hints and refused), then passes after T077.
  - Context: Add tests named `test_send_colour_*` to the file of T073. New keyword-only argument `colour: Optional[ColourSettings] = None` of `Viewer.send(obj, *, id, name=None, parent=None, colour=None, **hints)`. Tests: `send(BlockModel(...with "grade"), id="grid", colour=ColourSettings(array="grade", colour_map="coolwarm", range=(0.0, 1.0), scale="log", plain=(9, 8, 7)))` returns revision 1, the object is in `fake.objects`, and `fake.colouring["grid"]` has all five settings; the object is in the fake before the first colour command (monkeypatch `FakeViewer._tool_reply` to record `"grid" in fake.objects` at the first `colour.set_active`); the commands arrive in the order of `apply_colour`; a failed command (`ColourSettings(array="missing")`) raises `InvalidInputError`, and the object stays in the viewer (`"grid" in fake.objects`, `viewer.revisions["grid"] == 1`); a refused object (`FakeViewer(reply_errors={"grid": "bad"})`) raises `ViewerRefusedError` and no colour command is sent (`flush` raises first); `colour=None` sends no tool frame (`fake.tool_requests` is empty and `viewer.tool_connected` is false); a `colour` that is not a `ColourSettings` gives `InvalidInputError` before the object goes out (`fake.objects` stays empty).

- [x] T077 [US1] [M] Add the `colour` argument to `Viewer.send`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T076
  - Wave: 19
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_colour.py packages/loop_link/tests/test_viewer.py packages/loop_link/tests/test_push.py` passes and `ruff check packages/loop_link` passes.
  - Context: Change the signature to `send(self, obj, *, id, name=None, parent=None, colour: Optional[ColourSettings] = None, **hints)`; the argument is before `**hints`, so it does not go to the adapter. First check `colour` (`None` or a `ColourSettings`, else `InvalidInputError("colour must be a ColourSettings object.")`) before `convert` and before any frame. Send the object as before. Then, when `colour is not None`: call `self.flush()` (this waits up to 0.5 s for a late refusal of the object and raises `ViewerRefusedError`), then `self.apply_colour(id, colour)`. If a colour command fails, the object stays in the viewer and the error goes up. The method returns the revision as before. Add the argument to the docstring (Parameters, Raises). Import `ColourSettings` from `loop_link.commands`. Do not change `push_objects`.

- [x] T078 [P] [Found] [S] Export `ColourSettings`.
  - Files: `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T068
  - Wave: 7
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "import loop_link as l; [getattr(l, n) for n in l.__all__]; assert 'ColourSettings' in l.__all__"` passes and `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_import.py` passes.
  - Context: Add `from .commands import ColourSettings` before the `from .errors import (...)` block (the module imports are in alphabetical order, `.commands` is before `.errors`) and add `"ColourSettings"` to `__all__` after `"Categorical"` (the list is sorted). Keep all other lines. Do not import pyvista, pandas, Qt or LoopStructural (SC-005, `test_import.py`).

- [x] T079 [US1] [M] Test: SC-007 speed of filter and colour commands in `test_colour.py`.
  - Files: `packages/loop_link/tests/test_colour.py`
  - Needs: T077
  - Wave: 20
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_colour.py -k speed` passes (the test checks finished code, so it has no failing first run).
  - Context: Mark `@pytest.mark.slow`. Name `test_sc007_speed_*`. Use `fake_viewer` and `viewer_factory`. Send `BlockModel(origin=(0, 0, 0), cell_sizes=(1, 1, 1), cells=(100, 100, 100), cell_properties={"grade": np.random.default_rng(0).random(1_000_000, dtype=np.float32)})` with id `grid`; call `viewer.flush()`. Warm up the tool connection with one `set_colour_property("grid", "grade")` (not timed). Then time each call with `time.perf_counter()` and assert that each is below 0.250 s: `set_colour_property`, `set_colour_map("grid", "viridis")`, `set_colour_range(0.0, 1.0)`, `clear_colour_range`, `set_colour_scale("linear")`, `set_object_colour((1, 2, 3))`, `filter_range("grid", "grade", 0.5, 1.0)` (3 commands: the whole call is below 0.250 s), `filter_enable("grid", "grade", False)`, `filter_remove("grid", "grade")`. The fake has `filter_delay=0.0`. The time of a filter build is not counted (spec SC-007).

### Stage S4d: wait for a filter build (FR-019; answered by geoviewer spec 072, protocol 1.20)

Spec 072 of geoviewer (merged 2026-11-09) makes `filter.edit` a long command. The viewer sends `command.progress` while it builds and `command.result` when the build ends. The result can have `passing`, `total` and `state` (`building` or `unbuilt`). See "Answers" of `plan.md` (OP15 and OP20). The wait is the wait for that result. A viewer of protocol 1.19 or less answers at once with no counts.

- [x] T088 [US1] [M] Test: the fake viewer gives a long `filter.edit` (protocol 1.20).
  - Files: `packages/loop_link/tests/test_fake_viewer.py`
  - Needs: T069, T071
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py -k filter_long` fails (the fake has no long `filter.edit`), then passes after T089.
  - Context: Read the helpers of the file first (how a test sends a `command.run` frame and reads frames). Tests are named `test_filter_long_*`. Use `FakeViewer(protocol_minor=20, filter_delay=0.4)`. Send a `BlockModel` of 1,000 cells with a float32 cell property `grade` = `np.arange(1000, dtype=np.float32)`, then these `filter.edit` requests on the tool connection path: `add` for `grade`, `set_range` min `"0.0"` max `"99.0"` keep_nan false, `set_on` true. Checks: (1) the `set_on` request gets at least one `command.progress` frame with the same `request` number, `fraction` from 0.0 to 1.0 and `text` `"filter.edit"`, before its `command.result`, and the result is `{"status": "done", "passing": 100, "total": 1000}`. (2) `set_on` false gives `passing` 1000 and `total` 1000. (3) `remove` and `add` give `{"status": "done"}` with no progress frame and no `passing`. (4) A `command.cancel` with the request number of a running `set_on` gives a `request.error` with code `cancelled`, and `fake.filters` shows the filter as it was before that edit (`on` false). (5) With `FakeViewer(protocol_minor=20, filter_delay=0.4, filter_stop_after=0.1)` the `set_on` result is `{"status": "done", "state": "building"}` after about 0.1 s. (6) The `commands` frame has `"long": true` for `filter.edit` and no `long` key for the other commands. (7) With `protocol_minor=18` the result of `set_on` is `{"status": "done"}` with no progress frame and no `long` key. The old tests must still pass.

- [x] T089 [US1] [M] Implement the long `filter.edit` in the fake viewer.
  - Files: `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T088, T071
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fake_viewer.py packages/loop_link/tests/test_filter.py packages/loop_link/tests/test_colour.py` passes and `ruff check packages/loop_link` passes.
  - Context: Change only the code for `protocol_minor >= 20`; the behaviour for a lower number stays as it is (busy time, result `{"status": "done"}`). Add the argument `filter_stop_after: float | None = None` (seconds; document it). In the `commands` frame add `"long": True` to the entry of `filter.edit` when `protocol_minor >= 20`. Route a `filter.edit` request to a new method `_stream_filter(conn, reader, frame)`, in the way `_stream_fetch` is routed (see the call of `_stream_fetch` in `_serve`; reuse `_cancelled`). Checks of input and of the object stay as in `_command_run`. A build runs when the edit is `set_range`, `set_values`, `set_exclude` or `set_on` and the filter is on after the edit and `filter_delay > 0`. During a build send `command.progress` frames (`{"type": "command.progress", "request": n, "fraction": f, "text": "filter.edit"}`) about every 0.1 s (at least one), poll `_cancelled` between them, and use the tool protocol version of `protocol.TOOL_PROTOCOL_VERSION`. On cancel: restore the filters of the object to the state before the edit, send `_request_error(request, "cancelled", ...)` and stop. If `filter_stop_after` is set and shorter than the build, send `{"status": "done", "state": "building"}` at that time (the edit stays applied). At the end of a build, send the result with `passing` and `total`: take the array of the edit from `StoredObject.buffers`, use its `association` from the `data` entry of the header, and count the elements that pass all filters of that association that are on (a range filter keeps `minimum <= v <= maximum`, and NaN only if `keep_nan`; a values filter keeps the listed codes, or all others if `exclude`); `total` is the size of the array. An edit that starts no build gives the same counts when the filter is on or off (a filter that is off passes all), `{"status": "done"}` for `add`, `remove` and `cancel`, and `{"status": "refused", ...}` as before. Do not set the busy time for `protocol_minor >= 20`. Docstrings in numpydoc style.

- [x] T090 [US1] [M] Test: `ToolConnection.post` sends a request and does not wait.
  - Files: `packages/loop_link/tests/test_tool.py`
  - Needs: T046, T089
  - Wave: 12
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_tool.py -k post` fails (no `post`), then passes after T091.
  - Context: Read the file for the fixtures that start a `FakeViewer` and a `ToolConnection`. Tests are named `test_tool_post_*`. Use `FakeViewer(protocol_minor=20, filter_delay=0.5)` with a `BlockModel` of 1,000 cells and the array `grade` (as in `test_fake_viewer.py`). `post(build)` takes the same `build` argument as `request`, returns the request number (an `int`) and returns in less than 0.2 s while the fake builds. A later `request` on the same connection (for example `commands.list`) gets its own reply and skips the late frames of the posted request. A `post` with no viewer raises `NoViewerError`; a closed connection raises `ViewerClosedError`.

- [x] T091 [US1] [M] Add `ToolConnection.post` to `tool.py`.
  - Files: `packages/loop_link/src/loop_link/tool.py`
  - Needs: T090, T047
  - Wave: 13
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_tool.py` passes and `ruff check packages/loop_link` passes.
  - Context: Add `post(self, build: Callable[[int], bytes]) -> int` next to `request`. It takes the lock, connects like `request` does (`_connect_locked`), adds 1 to the request number, sends `build(number)` with `send_binary`, and returns the number. It reads no frame. A link error (`_LINK_ERRORS`) drops the connection and raises `ViewerClosedError`, as `request` does. The frames of the posted request that come later are skipped by `_read_locked` (`header.get("request") != number`). Docstring in numpydoc style: say that a refusal or an error of the posted request is not seen.

- [x] T080 [US1] [M] Test: `FilterResult` and the wait of the filter calls.
  - Files: `packages/loop_link/tests/test_filter.py`, `packages/loop_link/tests/test_commands.py`
  - Needs: T074, T089, T091, T068
  - Wave: 17
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_filter.py packages/loop_link/tests/test_commands.py -k "wait or filter_result"` fails (no `FilterResult`), then passes after T081.
  - Context: Read OP15 and OP20 and their "Answers" in `plan.md` first. `loop_link.commands.FilterResult` is a frozen dataclass with `passing: int | None = None`, `total: int | None = None`, `state: str | None = None`, and the class method `from_result(result: Mapping)`. In `test_commands.py` (tests named `test_filter_result_*`): `{"status": "done", "passing": 100, "total": 1000}` gives `FilterResult(100, 1000, None)`; `{"status": "done"}` gives all `None`; `{"status": "done", "state": "unbuilt"}` and `"building"` give that state; `passing` that is a `bool`, a negative number or not an `int` raises `loop_link.protocol.ProtocolError`; `state` that is not a `str` raises `ProtocolError`. In `test_filter.py` (tests named `test_filter_wait_*`) use `FakeViewer(protocol_minor=20, filter_delay=0.4)` with `BlockModel` of 1,000 cells, `grade` = `np.arange(1000, dtype=np.float32)` and `viewer_factory`: `filter_range("grid", "grade", 0.0, 99.0)` returns `FilterResult(100, 1000, None)` after at least 0.3 s; `filter_values` on a category array returns the counts too; `filter_enable` and `filter_remove` return a `FilterResult`; `progress=callback` gets floats from 0.0 to 1.0 (at least one); `wait=False` returns `FilterResult()` in less than 0.2 s and a next call (`filter_enable`) works; a `cancel` event that a timer sets after 0.1 s raises `loop_link.errors.CancelledError` and `fake.filters` shows the filter off; with `filter_stop_after=0.1` the result has `state == "building"`; with `protocol_minor=18` the call returns `FilterResult()` at once. The old tests of both files must still pass (they ignore the return value).

- [x] T081 [US1] [M] Implement `FilterResult` and the wait of the filter calls.
  - Files: `packages/loop_link/src/loop_link/commands.py`, `packages/loop_link/src/loop_link/viewer.py`, `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T080, T077, T091
  - Wave: 20
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_filter.py packages/loop_link/tests/test_commands.py packages/loop_link/tests/test_colour.py packages/loop_link/tests/test_import.py` passes and `ruff check packages/loop_link` passes.
  - Context: Follow "Answers" of `plan.md` (OP15, 2026-11-09). In `commands.py` add `FilterResult` (add it to `__all__`) and give `run_command` the argument `progress: Callable[[float], None] | None = None` that it passes to `tool.request(progress=...)`. In `viewer.py`: `filter_range` and `filter_values` return `FilterResult` and use the arguments `wait`, `cancel`, `timeout` and the new `progress: Callable[[float], None] | None = None` (keyword only, after `timeout`). With `wait=True` every edit goes with `run_command(..., cancel=cancel, timeout=timeout, progress=progress)` and the call returns `FilterResult.from_result(result of the last set_on)`. With `wait=False` the edits before the last go as before, and the last `set_on` goes with `tool.post(...)` (build the frame with `encode.command_frame`); the call returns `FilterResult()`. `filter_enable` and `filter_remove` return `FilterResult.from_result(...)` of their result. `cancel_filter` stays as it is. `cancel` (a `threading.Event`) is passed to `tool.request`, which sends one `command.cancel`; the error `cancelled` becomes `CancelledError` already. A client `timeout` is the time to wait for the next frame, as in `ToolConnection.request` (a late frame raises `ViewerClosedError`); the viewer itself ends a long wait after 60 s with `state == "building"`. Remove the words "Not used yet" from the docstrings and describe `wait`, `cancel`, `timeout`, `progress` and the return value (`state` `building`: the build goes on; `unbuilt`: the viewer has no counts for the new filters; a viewer of protocol 1.19 or less gives no counts). Export `FilterResult` in `__init__.py` (keep `__all__` sorted).

- [x] T092 [US1] [M] Update M4 step 7 and the README for the filter wait.
  - Files: `packages/loop_link/tests/manual_checks.md`, `packages/loop_link/README.md`
  - Needs: T081, T084, T085
  - Wave: 21
  - Done when: `git grep -n -c "OP15" packages/loop_link/tests/manual_checks.md` gives 0 and `git grep -n -c "FilterResult" packages/loop_link/README.md` gives 1 or more.
  - Context: Read the section `## M4` of `manual_checks.md` and the section `## Filter and colour` of the README first. In M4: say in the intro that step 7 needs a viewer of protocol 1.20 (geoviewer spec 072) and that the other steps work with 1.18. Replace the part of step 7 about OP15 ("check what a second command does while the viewer builds (OP15)") with these checks: send an object of 1,000,000 cells; call `filter_range("grid", "grade", 0.5, 1.0)` with the default `wait=True`; check that the call returns after the build, that `FilterResult.passing` and `total` are numbers, and that `passing` fits the picture; call it again with `wait=False` and check that it returns at once; set a `threading.Event` during a build and check `CancelledError` and that the filter is as before; record the time of the build. In the closing paragraph change the list of open points to OP16, OP18 and OP19. In the README: show the result in the example (`result = viewer.filter_range("grid", "grade", 0.5, 1.0)` and `print(result.passing, result.total)`) and add 3 short sentences: the call waits for the end of the build by default (`wait=False` does not wait, `cancel` is a `threading.Event`, `progress` gets the fraction); the counts need a viewer of protocol 1.20, and an older viewer gives a `FilterResult` with no counts; `state` is `building` when the viewer stopped to wait (the build goes on) and `unbuilt` when the viewer has no counts for the new filters. ASD-STE100 English.
### Stage S4e: real codec, fixtures, manual check M4, docs

- [x] T082 [P] [Finish] [M] Test: command frames against the real geoviewer codec.
  - Files: `packages/loop_link/tests/test_codec_real.py`
  - Needs: T066, T058
  - Wave: 5
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_codec_real.py -rs` passes or skips with the reason that the variable is not set; with `LOOP_LINK_GEOVIEWER_REPO=C:\Users\lachl\dev\geoviewer` and `cargo` on the PATH it passes with no skip.
  - Context: In `_cases()` add these frames (the key is the case name): `commands_list` = `encode.commands_list_frame(1)`; `command_filter_add`, `command_filter_set_range`, `command_filter_set_values`, `command_filter_set_on`, `command_filter_set_exclude`, `command_filter_remove`, `command_filter_cancel` = `encode.command_frame(2, "filter.edit", {"object": "grid", "edit": {...}})` with the edits `{"kind": "add", "array": "grade"}`, `{"kind": "set_range", "array": "grade", "min": "0.5", "max": "1.0", "keep_nan": False}`, `{"kind": "set_values", "array": "unit", "values": [0, 2]}`, `{"kind": "set_on", "array": "grade", "on": True}`, `{"kind": "set_exclude", "array": "unit", "exclude": True}`, `{"kind": "remove", "array": "grade"}`, `{"kind": "cancel"}`; `command_colour_set_active` (`colour.set_active`, `{"id": "grid", "array": "grade"}`), `command_colour_set_map` (`colour.set_map`, `{"id": "grid", "name": "coolwarm"}`), `command_colour_set_range` (`colour.set_range`, `{"id": "grid", "min": "0.0", "max": "1.0"}`), `command_colour_clear_range` (`colour.clear_range`, `{"id": "grid"}`), `command_colour_set_scale` (`colour.set_scale`, `{"object": "grid", "scale": "log"}`), `command_object_set_colour` (`object.set_colour`, `{"id": "grid", "colour": [1, 2, 3]}`), `command_object_set_colour_none` (`object.set_colour`, `{"id": "grid", "colour": None}`). The helper, the validator program and the negative control stay. Do not change `codec_validator/`. If the validator reports FAIL for a case, make the task fail with the validator output; do not drop the case.

- [x] T083 [P] [Finish] [S] Test: command fixtures of the geoviewer repository.
  - Files: `packages/loop_link/tests/test_fixtures.py`
  - Needs: T066, T059
  - Wave: 5
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fixtures.py -rs` passes or skips with a reason; with `LOOP_LINK_GEOVIEWER_FIXTURES=C:\Users\lachl\dev\geoviewer\fixtures\protocol\python` it passes with no skip for the new tests.
  - Context: Use the helpers `_tool_header(name)` and `_check_header(made, fixture)` that the file already has (they skip when the folder `<variable>/../tool` is missing). Add two tests. `test_commands_list_matches_fixture`: `fixture = _tool_header("commands.list.json")`; `assert "request" in fixture`; build `commands_list_frame(fixture["request"])`; call `_check_header`. `test_command_run_matches_fixture`: `fixture = _tool_header("command.run.json")`; build `command_frame(fixture["request"], fixture["command"], fixture["input"])`; call `_check_header`. Add `commands_list_frame` and `command_frame` to the import from `loop_link.encode` (keep the list sorted). The fixtures are `{"type":"commands.list","version":{"major":1,"minor":11},"request":1,"buffers":[]}` and `{"type":"command.run","version":{"major":1,"minor":11},"request":2,"command":"object.set_visible","input":{"id":"p1","visible":false},"buffers":[]}`. Use assert messages in the style of the other tests of the file.

- [x] T084 [P] [Finish] [S] Manual check M4 in `manual_checks.md`.
  - Files: `packages/loop_link/tests/manual_checks.md`
  - Needs: T060
  - Wave: 2
  - Done when: `git grep -n -c "^## M4" packages/loop_link/tests/manual_checks.md` gives 1.
  - Context: Add a section `## M4: filter and colour with the real window viewer (stories 6 and 7)` before `## Record`, in ASD-STE100 English, in the style of M3. Intro: use a viewer of protocol 1.18, start with `Viewer.clear()`, and place the test objects near the origin of the scene. Steps: (1) send a `BlockModel` with a number property `grade`, call `filter_range("grid", "grade", 0.5, 1.0)`, check that only the cells in the range show; (2) send an object with a category property, call `filter_values` with two codes, then with `exclude=True`, check the picture each time; (3) call `filter_enable(False)`, then `filter_remove`, check that all cells show again; (4) call `filter_range` for a text property and for a property that does not exist, check the `InvalidInputError` and its text; (5) call `set_colour_property`, `set_colour_map("coolwarm")`, `set_colour_range`, `set_colour_scale("log")` and `set_object_colour`, check the picture after each call; (6) call the map, range and scale methods for a category property and check the error, then call `colour_map_names()` and check that it gives the names of the viewer; (7) send the same id again (a new revision) and check that the filter and colours stay; then send an object of 1,000,000 cells, call `filter_range` with `wait=False`, check what a second command does while the viewer builds (OP15), and call `cancel_filter` during the build; (8) send an object of 1,000,000 cells with `colour=ColourSettings(...)` and check that the settings apply (OP19: the 0.5 s of `flush` can be too short); (9) check the assumptions: `add` for an existing filter (OP16), and that a new filter is on after `add` and `set_range`. Closing paragraph: M4 checks OP15, OP16, OP18 and OP19, and records the time of one `set_range` on the object of 1,000,000 cells. Pass line: all 9 steps give the expected result.

- [x] T085 [P] [Finish] [S] README section "Filter and colour".
  - Files: `packages/loop_link/README.md`
  - Needs: T077, T078
  - Wave: 20
  - Done when: `git grep -n -c "^## Filter and colour" packages/loop_link/README.md` gives 1.
  - Context: Add `## Filter and colour` after `## Select and fetch` and before `## Test against the real codec`, with an example of about 10 lines: `from loop_link import ColourSettings, Viewer`; `viewer = Viewer()`; `viewer.send(grid, id="grid", colour=ColourSettings(array="grade", colour_map="viridis", range=(0.0, 1.0), scale="linear"))`; `viewer.filter_range("grid", "grade", 0.5, 1.0)`; `viewer.filter_values("grid", "unit", [0, 2], exclude=True)`; `viewer.filter_remove("grid", "grade")`; `viewer.set_colour_map("grid", "coolwarm")`; `viewer.colour_map_names()`. Add 3 short sentences: filter and colour use the tool connection, so they work while another program holds the data connection; the colour map names come from the viewer (protocol 1.18), and an older viewer gives `MissingFeatureError` for `colour_map_names()`; a value filter takes the integer ids of the categories (the codes of `Categorical`). Do not write about `wait` (it is not built yet). ASD-STE100 English.

## Phase 11: Finish S4

- [x] T086 [Finish] [S] Update `ROADMAP.md`: stage S4 done.
  - Files: `ROADMAP.md`
  - Needs: T081, T082, T083, T084, T085, T079, T092
  - Wave: 22
  - Done when: `git grep -n "stage S4" ROADMAP.md` shows 1 or more lines and `git grep -n -c "loop_link" ROADMAP.md` gives 3 or more.
  - Context: In the entry `- [ ] **Stage 7 — `loop_link`: Python link to geoviewer.**` replace the sentence "Filter and colour (stage S4) are not planned yet." with "Stage S4 (filter and colour) is done." (read the entry first; the sentence can differ a little). Add one status log entry at the end of "Status log", in the style of the neighbours, dated with the date of the run (YYYY-MM-DD): "Stage 7 stage S4 implemented: filter and colour calls, colour in `send` and the fake commands endpoint in `packages/loop_link`." Tick the Stage 7 box only if no other stage of the entry is open. If you are not sure, leave the box and say so in the result.

- [!] T087 [Finish] [M] Run the merge gates for S4 and fix the findings.  (root pytest finished: 7 failed in tests/unit/interpolator/test_2d_p1_p2_support.py, numpy 2D cross product, outside loop_link; ruff: 2 RUF100 in packages/loop_common)
  - Files: `packages/loop_link/`
  - Needs: T086, T062, T038, T063
  - Wave: 23
  - Done when: `ruff check LoopStructural packages && pytest` passes, `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests` passes, and `black --check packages/loop_link` passes.
  - Context: Constitution merge gate: `pytest` and `ruff check LoopStructural packages` must both pass. The root `pytest` does not run the package tests, so run the package command too. The root `pytest` was too slow on this machine for S5 (T063): run it, and if it cannot finish, say so in the result. Run `black packages/loop_link` to fix the format (line length 100, quotes kept). Fix only problems in `packages/loop_link/`. Do not change a gate or add a Ruff ignore for `packages/loop_link/src`. Also run `packages/loop_link/tests/test_import.py` alone (SC-005). If `LOOP_LINK_GEOVIEWER_REPO` and `LOOP_LINK_GEOVIEWER_FIXTURES` are set, run `test_codec_real.py` and `test_fixtures.py` with them. Report the result of each command in the commit message.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001, T002, T004, T005, T039, T040, T060 | 5 | 2 |
| 2 | T003, T006, T041, T051, T064, T084 | 3 | 3 |
| 3 | T007, T009, T011, T016, T042, T052, T065 | 1 | 6 |
| 4 | T008, T010, T017, T043, T058, T059, T066 | 1 | 6 |
| 5 | T012, T015, T019, T067, T069, T082, T083 | 1 | 6 |
| 6 | T013, T018, T044, T068 | 0 | 4 |
| 7 | T014, T020, T045, T078 | 1 | 3 |
| 8 | T021, T046, T070 | 0 | 3 |
| 9 | T023, T024, T025, T027, T034, T037, T047, T071 | 3 | 5 |
| 10 | T028, T030, T032, T035, T048, T072, T073, T088 | 0 | 8 |
| 11 | T029, T031, T033, T036, T049, T089 | 1 | 5 |
| 12 | T022, T050, T090 | 1 | 2 |
| 13 | T038, T053, T091 | 0 | 3 |
| 14 | T054 | 0 | 1 |
| 15 | T055 | 0 | 1 |
| 16 | T056, T057, T074 | 1 | 2 |
| 17 | T061, T075, T080 | 1 | 2 |
| 18 | T062, T076 | 1 | 1 |
| 19 | T063, T077 | 0 | 2 |
| 20 | T079, T081, T085 | 1 | 2 |
| 21 | T092 | 0 | 1 |
| 22 | T086 | 1 | 0 |
| 23 | T087 | 0 | 1 |

## Critical path

T001 -> T006 -> T007 -> T008 -> T019 -> T044 -> T045 -> T046 -> T047 -> T048 -> T049 -> T050 -> T053 -> T054 -> T055 -> T074 -> T075 -> T076 -> T077 -> T081 -> T092 -> T086 -> T087

## Dependencies

T001 to T038 (stages S1 to S3) and T039 to T063 (stage S5) are built, except the root `pytest` result of T063.
Stage S4 adds T064 to T092.

- S4a: `COMMAND_PROTOCOL_VERSION` (T064), then the frame builders as test and code (T065, T066), the
  commands module as test and code (T067, T068), and the fake commands endpoint as one test task and two
  code tasks that share `fake_viewer.py` (T069, T070 colour part, T071 filter part).
- S4b: the filter and colour tests run apart (T072, T073). The `Viewer` members follow in a chain,
  because they share `viewer.py` (T074 filter, T075 colour). The export (T078) needs only T068.
- S4c: tests and code for `send(..., colour=...)` (T076, T077), then the speed test of SC-007 (T079).
- S4d: geoviewer spec 072 (protocol 1.20) answers OP15 and OP20. The fake gets a long `filter.edit` with counts
  (T088, T089), `ToolConnection.post` is added (T090, T091), `FilterResult` and the wait follow (T080, T081),
  and the docs and M4 step 7 change (T092).
- S4e: real codec cases (T082), fixture checks (T083), manual check M4 (T084), README (T085), ROADMAP (T086)
  and the gate (T087). T087 also needs T063, the S5 gate, which failed on the root `pytest` only.
- Open points: OP15 and OP20 are answered (spec 072). OP16, OP18 and OP19 are checked by M4 and block no task.
