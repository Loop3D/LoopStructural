# Tasks: loop_link, a Python link to geoviewer

**Feature**: `006-loop-link` | **Plan**: [plan.md](plan.md)

<!--
Scope: stages S1 to S3 of the plan (stories 1 and 5). Stages S4 (stories 6, 7;
OP6) and S5 (stories 2 to 4; OP7) get tasks after their own planner runs.
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
- Waves: 13. Tiers: [S] = simple model, [M] = standard model.
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

- [ ] T007 [P] [Found] [M] Test: the codec in `test_protocol.py`.
  - Files: `packages/loop_link/tests/test_protocol.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_protocol.py` fails with ImportError of `loop_link.protocol`, then passes after T008.
  - Context: API under test: `encode_frame(header: dict, buffers: Sequence[tuple[str, np.ndarray]] = ()) -> bytes`, `decode_frame(data: bytes) -> Frame` (`Frame(header: dict, buffers: dict[str, np.ndarray])`), `decode_header(data: bytes) -> dict`, `ProtocolError`, constants `MAX_FRAME_BYTES=256*1024*1024`, `MAX_HEADER_BYTES=1024*1024`, `PROTOCOL_VERSION={"major":1,"minor":14}`. Frame layout: `[u32 LE header_len][JSON header][buffers]`, little-endian. `encode_frame` adds `version` and `buffers` (`[{"name","dtype","shape"}]`). dtype names `f32 f64 u32 i32 u8`. Tests: a round trip for each of the 5 dtypes with shape (2,3); a zero-element array is valid; a frame whose length differs from `4+header_len+sum(sizes)` gives `ProtocolError`; a header over 1 MiB gives `ProtocolError` on encode; a frame under 4 bytes gives `ProtocolError`; NaN in the header gives an error (`allow_nan=False`); an unsupported dtype (`complex128`) gives `ProtocolError`. A byte-literal test: `encode_frame({"type":"heartbeat","seq":1})` starts with `struct.pack('<I', len(header_json))` and the JSON has no spaces. If the environment variable `LOOP_LINK_GEOVIEWER_FIXTURES` is set, also decode each `*.bin` file in it and assert `decode_frame` does not raise (skip with a reason when not set).

- [ ] T008 [Found] [M] Implement the codec `protocol.py`.
  - Files: `packages/loop_link/src/loop_link/protocol.py`
  - Needs: T007
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_protocol.py` passes and `ruff check packages/loop_link` passes.
  - Context: Port the idea of `encode_frame` and `decode_header` from `C:\Users\lachl\dev\plugin_loopstructural-geoviewer\loopstructural\main\geoviewer\protocol.py` (read it). Public names: `PROTOCOL_VERSION`, `MAX_FRAME_BYTES`, `MAX_HEADER_BYTES`, `MAX_OBJECT_BYTES = 2*1024**3`, `F32_LIMIT = 4000.0`, `Frame` (frozen dataclass), `encode_frame`, `decode_frame`, `decode_header`, `ProtocolError(LoopLinkError)`. The codec raises `ProtocolError` for size problems. A frame over `MAX_FRAME_BYTES` raises `ProtocolError` in `encode_frame` (the object-level `ObjectTooLargeError` comes from `encode.py`, T014). JSON: `json.dumps(full, separators=(",", ":"), allow_nan=False)`. Buffers: `np.ascontiguousarray`, little-endian `tobytes`. `decode_frame` returns arrays with the shape of the header (`np.frombuffer(...).reshape(shape)`, a copy is not needed). Numpy and the standard library only.

- [ ] T009 [P] [Found] [M] Test: vector conversion in `test_orientation.py`.
  - Files: `packages/loop_link/tests/test_orientation.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_orientation.py` fails with ImportError, then passes after T010.
  - Context: Under test: `normals_to_orientation(normals) -> (dip_direction, dip)`, `directions_to_lineation(directions) -> (trend, plunge)`, `orientation_to_normals(dip_direction, dip) -> (n,3)`, `lineation_to_directions(trend, plunge) -> (n,3)`. Degrees. x east, y north, z up. Independent expected value in the test: a plane that dips to azimuth `a` with dip `d` has the upward normal `(sin d sin a, sin d cos a, cos d)` (write this formula in the test, not by calling the library). Tests: 12 planes (`a` in 0, 45, 90, 135, 180, 270, 359; `d` in 10, 45, 80) convert to `(a, d)` within 1e-9 degrees; a downward normal gives the same result as its upward twin; a vertical plane (n_z=0) gives dip 90; a horizontal plane gives dip 0 and dip direction 0.0; lineations: `(trend, plunge)` to direction `(cos p sin t, cos p cos t, -sin p)` and back; a direction with z>0 flips (axis); vertical line gives trend 0.0 and plunge 90; the inverse functions round trip; a zero vector and a NaN vector each raise `InvalidInputError`. Dip direction and trend are in [0, 360).

- [ ] T010 [Found] [M] Implement `orientation.py`.
  - Files: `packages/loop_link/src/loop_link/orientation.py`
  - Needs: T009
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_orientation.py` passes.
  - Context: Signatures take array-like (n,3) or (3,) input and return numpy arrays: `normals_to_orientation(normals) -> tuple[np.ndarray, np.ndarray]`; `directions_to_lineation(directions) -> tuple[np.ndarray, np.ndarray]`; `orientation_to_normals(dip_direction, dip) -> np.ndarray`; `lineation_to_directions(trend, plunge) -> np.ndarray`. Rules: `n = v/|v|`; for normals, if `n_z < 0` use `-n`; `dip = degrees(arccos(clip(n_z,-1,1)))`; `dip_direction = degrees(arctan2(n_x, n_y)) % 360`; if `n_x == n_y == 0` the dip direction is 0.0. For lineations, if `d_z > 0` use `-d`; `plunge = degrees(arcsin(clip(-d_z,-1,1)))`; `trend = degrees(arctan2(d_x, d_y)) % 360`; if `d_x == d_y == 0` the trend is 0.0. Inverse of normals: `(sin d sin a, sin d cos a, cos d)`. Inverse of lineations: `(cos p sin t, cos p cos t, -sin p)`. A zero-length or non-finite vector raises `InvalidInputError` (from `loop_link.errors`).

- [ ] T011 [P] [Found] [M] Test: the object model in `test_objects.py`.
  - Files: `packages/loop_link/tests/test_objects.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_objects.py` fails with ImportError of `loop_link.objects`, then passes after T012.
  - Context: Frozen dataclasses in `loop_link.objects`: `Points(vertices, properties={})`; `Vectors(vertices, orientation=None, lineation=None, properties={})` with `Vectors.from_directions(vertices, directions, *, kind, properties=None)` (`kind` is `"normal"` or `"lineation"`, else `InvalidInputError`); `Lines(vertices, segments, properties, cell_properties)`; `Surface(vertices, triangles, properties, cell_properties)`; `BlockModel(origin, cell_sizes, cells, properties, cell_properties)`; `TetraMesh(vertices, tetrahedra, properties, cell_properties)`; `Drillhole(collar, survey, intervals, end_of_hole=None, properties, interval_properties)`; `Categorical(codes, names, colours=None)`; `DEFAULT_CATEGORY_COLOURS` (10 RGBA tuples). Rules to test (each raises `InvalidInputError` that names the field): vertices not (n,3) or not finite; an index array with a negative value or a value >= n; a property of the wrong length (vertex: n; cell: m; `BlockModel` nodes: (nx+1)(ny+1)(nz+1), cells: nx*ny*nz, x fastest; Drillhole `properties`: k stations, `interval_properties`: j rows); a 2-D property; `cell_sizes` <= 0 or `cells` < 1 or non-integer; `Points` or `Vectors` with `cell_properties` (they have no such field, so the test checks `TypeError`); `Categorical` with a code >= `len(names)`. Valid input of all 7 types is accepted and arrays become numpy arrays.

- [ ] T012 [Found] [M] Implement `objects.py`.
  - Files: `packages/loop_link/src/loop_link/objects.py`
  - Needs: T011, T010
  - Wave: 5
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_objects.py` passes.
  - Context: Implement the dataclasses that T011 names, with the same fields and defaults (`field(default_factory=dict)` for the mappings). Each class does its checks in `__post_init__` and converts arrays with `object.__setattr__(self, name, np.asarray(...))` (frozen). Vertex and index arrays: `float64` and integer. `Vectors.from_directions(..., kind="normal")` calls `normals_to_orientation` and `kind="lineation"` calls `directions_to_lineation` (both in `loop_link.orientation`) and stores the pair in `orientation` or `lineation`. `Vectors` accepts both `orientation` and `lineation`. Add `SceneObject = Union[Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole]`, `Property = Union[np.ndarray, Sequence, Categorical]` and `ArrayLike`. `DEFAULT_CATEGORY_COLOURS` = the 10 tab10 colours as RGBA with alpha 255: (31,119,180), (255,127,14), (44,160,44), (214,39,40), (148,103,189), (140,86,75), (227,119,194), (127,127,127), (188,189,34), (23,190,207). A category without a given colour uses the colour at index `i % 10`. Do not import pyvista, pandas or Qt.

- [ ] T013 [Found] [M] Test: object encoding in `test_encode.py`.
  - Files: `packages/loop_link/tests/test_encode.py`
  - Needs: T012, T008
  - Wave: 6
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_encode.py` fails with ImportError of `loop_link.encode`, then passes after T014.
  - Context: Under test: `encode_object(obj, *, id, name, parent, revision, origin) -> bytes` (decode with `decode_frame`), `hello_frame(crs, origin, version)`, `heartbeat_frame(seq)`, `remove_frame(id)`, `clear_frame()`, `group_frame(id, name, parent, revision)`. Header fields of `scene.upsert`: `type id object_type name parent revision edit([]) replaces(null) geometry params data property_groups drape(null) style({})`. Type map: Points and Vectors -> `Points` (geometry `vertices`); Lines -> `Curve` (`vertices`, `segments` u32); Surface -> `Surface` (`vertices`, `triangles` u32); TetraMesh -> `TetraMesh` (`vertices`, `tetrahedra` u32); BlockModel -> `BlockModel` (no geometry; `params` `origin`, `cell_sizes`, `cells` as JSON lists, origin minus scene origin); Drillhole -> `Drillhole` (`collar` minus origin, `survey` f64 (k,3), `intervals` f64 (j,2); `params.end_of_hole=[v]` only when set). Vertices are f32 when max abs after the shift <= 4000, else f64; the shift is done in f64 (test with origin 500000, 7000000, 0 and vertices near it: f32 buffer, values equal to input minus origin within 1e-3). Data types: float32 -> `float`/f32; float64 -> `float`/f64; bool -> `boolean`/u8; int64 inside i32 -> `integer`/i32; int64 outside i32 -> `InvalidInputError`; `U` strings -> `text` (u8 bytes + u32 offsets named `data/<name>/offsets`, header field `offsets`); `Categorical` -> `referenced`/u32 with `value_map` entries `{"id","name","colour"}`; `object` dtype -> `InvalidInputError`; two properties that give the same buffer name -> `InvalidInputError`. Data array: `id`=`name`=property name, `buffer`=`data/<name>`, `association` `vertex` or `cell` (Drillhole `interval_properties`: `interval`). Vectors: `property_groups` `[{"id":"orientation","name":"orientation","role":"orientation","members":[{"role":"dip_direction","data":"orientation/dip_direction"},{"role":"dip","data":"orientation/dip"}]}]` and the f64 vertex arrays; lineation uses `lineation`, `trend`, `plunge`. Too large: build an object whose declared size exceeds 2 GiB with `np.broadcast_to` views (no real allocation) and expect `ObjectTooLargeError` with the id, size and limit.

- [ ] T014 [Found] [M] Implement `encode.py`.
  - Files: `packages/loop_link/src/loop_link/encode.py`
  - Needs: T013
  - Wave: 7
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_encode.py` passes.
  - Context: Public: `encode_object`, `hello_frame(crs: str, origin: Sequence[float], version: str) -> bytes` (header `type hello`, `application {"name":"loop_link","version":version}`, `crs`, `origin` as floats), `heartbeat_frame(seq: int)`, `remove_frame(id)` (`scene.remove`), `clear_frame()` (`scene.clear`), `group_frame(id, name, parent, revision)` (a `scene.upsert` of `object_type` `Group`). Private `_data_array(name, value, association) -> (desc: dict, buffers: list[tuple[str, np.ndarray]])` follows the type table in the test of T013. Use `loop_link.protocol.encode_frame`. Before building the frame, sum the declared bytes of all buffers plus `entries*8+name bytes` of each `value_map`; if the sum is over `MAX_OBJECT_BYTES`, or the frame is over `MAX_FRAME_BYTES` (catch `ProtocolError` from `encode_frame`), raise `ObjectTooLargeError(id, size, limit)`. Use broadcast views without copying when computing sizes (read `array.size * itemsize`, not `tobytes`). Coordinates: reject NaN and infinity with `InvalidInputError` (objects already check this). Principle IV: subtract the origin in float64 before any cast to f32.

- [ ] T015 [Found] [M] Port the stdlib WebSocket client to `_websocket.py`.
  - Files: `packages/loop_link/src/loop_link/_websocket.py`
  - Needs: T008
  - Wave: 5
  - Done when: `ruff check packages/loop_link` passes and `packages/loop_link/.venv/Scripts/python.exe -c "import loop_link._websocket"` passes. Behaviour is tested through `test_viewer.py` (T020).
  - Context: Read and port `C:\Users\lachl\dev\plugin_loopstructural-geoviewer\loopstructural\main\geoviewer\websocket_client.py` (245 lines, standard library only). Keep the class `WebSocketClient` with `connect(cls, host, port, *, token="", timeout=5.0, path="/")`, `send_binary(data: bytes)`, `receive(timeout: float) -> Optional[bytes]`, `close()`. Changes: (1) `WebSocketError(LoopLinkError)` and `WebSocketClosed(WebSocketError)`; (2) `MAX_INCOMING_BYTES = MAX_FRAME_BYTES + MAX_HEADER_BYTES` (import from `loop_link.protocol`; the old limit was 16 MiB); (3) in `connect`, an HTTP 401 answer raises `TokenError`, and a refused or timed-out TCP connection raises `NoViewerError` (both from `loop_link.errors`); (4) add numpy-style docstrings and type hints to every public method (Ruff `D` and `ANN` apply to new files); (5) the handshake header stays `Authorization: Bearer <token>`. Keep the thread rule: not thread safe.

- [ ] T016 [P] [Found] [M] Test: discovery in `test_discovery.py`.
  - Files: `packages/loop_link/tests/test_discovery.py`
  - Needs: T006
  - Wave: 3
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_discovery.py` fails with ImportError, then passes after T017.
  - Context: Under test: `find_viewer(*, port=None, token=None, host="127.0.0.1", settings=None, timeout=0.5) -> ViewerAddress` (frozen dataclass `host, port, token`), `settings_dir() -> Path` (`%APPDATA%\geoviewer` on Windows, else `~/.geoviewer`), `DEFAULT_PORT = 47615`. Use `tmp_path` as `settings` and `monkeypatch` for `GEOVIEWER_TOKEN`. Token order: the argument; then `GEOVIEWER_TOKEN` (an empty value is not a token); then the field `token` of `settings/link.json`; none gives `TokenError`; a damaged `link.json` (not JSON) gives `TokenError`. Instance files: `settings/instances/<pid>.json` with fields `version (1), pid, port, kind ("window" or "headless"), started, focused` (ms). To make a live port, open a `socket` listener on port 0 in the test. Rules: a file whose port has no listener is dropped; `window` before `headless`; then the largest `focused`; a damaged file is ignored; an explicit `port` skips the files; no file and no listener on the default port gives `NoViewerError`. The function never starts a program. `settings_dir()` test: monkeypatch `sys.platform`, `APPDATA` and `HOME`.

- [ ] T017 [Found] [M] Implement `discovery.py`.
  - Files: `packages/loop_link/src/loop_link/discovery.py`
  - Needs: T016
  - Wave: 4
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_discovery.py` passes.
  - Context: Public: `DEFAULT_PORT = 47615`, `ViewerAddress`, `settings_dir()`, `find_viewer(...)`. Token text on failure: "Set GEOVIEWER_TOKEN, or save a token in the link popup of the viewer." (`TokenError`). Alive check: `socket.create_connection((host, port), timeout)` then close. Sort key of live instances: `(kind != "window", -focused)`. Read JSON with `json.loads`; catch `OSError` and `ValueError` per file. The `NoViewerError` text says: "No viewer found. Start geoviewer and turn the link on. This library never starts a viewer." Use only the standard library.

- [ ] T019 [Found] [M] Implement the fake viewer.
  - Files: `packages/loop_link/src/loop_link/testing/__init__.py`, `packages/loop_link/src/loop_link/testing/fake_viewer.py`
  - Needs: T008, T006
  - Wave: 5
  - Done when: `ruff check packages/loop_link` passes and `packages/loop_link/.venv/Scripts/python.exe -c "from loop_link.testing import FakeViewer"` passes. Behaviour is tested by T020 and T025.
  - Context: Base: `C:\Users\lachl\dev\plugin_loopstructural-geoviewer\tests\unit\geoviewer_fake_server.py` (179 lines, read it). New class `FakeViewer(token="secret", crs="EPSG:28350", origin=(0.0, 0.0, 0.0), kind="window", heartbeat_timeout=15.0, reply_errors=None, read_delay=0.0)`. Members: `.port`, `.objects: dict[str, StoredObject]` (`StoredObject(header: dict, buffers: dict[str, np.ndarray])`), `.received: list[Frame]`, `.inventory() -> list[dict]` (`id`, `object_type`, `revision`), `.write_instance_file(settings: Path) -> Path` (writes `instances/<pid>.json` and `link.json` is not written), `.accept_frame(data: bytes) -> list[Frame]` (apply one data-client frame and return the reply frames), `.close()`, context manager. Rules: listen on `127.0.0.1` port 0 in a thread; many connections; bearer token check (HTTP 401 when wrong); answer `hello` with `hello` (`application {"name":"fake-viewer"}`, `crs`, `origin`) then `scene.inventory` (`entries`); the second data connection (path not `/tool`) gets an `error` frame code `busy` then is closed; close a data connection with no frame for `heartbeat_timeout` seconds; the path `/tool` is accepted and every message gets `error` code `unknown_message`; `scene.upsert` rules (geometry roles present for the type, index < vertex count, data array length from association, finite coordinates, buffer sizes match, object size <= 2 GiB) else `error` code `invalid_object` with field `object`; `reply_errors` maps id to a forced error text; `scene.remove` and `scene.clear` change `objects`; `read_delay` waits that many seconds after each frame. The object keeps the data until a `scene.remove`. Frames use `loop_link.protocol`. Add docstrings and type hints.

- [ ] T018 [Found] [M] Test support: `conftest.py`.
  - Files: `packages/loop_link/tests/conftest.py`
  - Needs: T019, T017
  - Wave: 6
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests --collect-only -q` passes.
  - Context: Fixture `fake_viewer` (function scope): start `FakeViewer()`, write its instance file into `tmp_path / "geoviewer"`, `monkeypatch.setenv("GEOVIEWER_TOKEN", "secret")`, yield an object with `.viewer` (the `FakeViewer`) and `.settings` (the settings `Path`), and close at the end. Fixture `viewer_factory`: returns a function that makes `Viewer(settings_dir=<settings>, heartbeat_interval=0.1)` objects and closes all of them at teardown. Helper `geoviewer_fixtures() -> Optional[Path]`: the path in `LOOP_LINK_GEOVIEWER_FIXTURES` when it is a folder, else `None`. No test code in this file.

- [ ] T020 [Found] [M] Test: `Viewer` in `test_viewer.py` (story 1).
  - Files: `packages/loop_link/tests/test_viewer.py`
  - Needs: T018, T012
  - Wave: 7
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_viewer.py` fails with ImportError of `loop_link.viewer`, then passes after T021.
  - Context: Use fixtures `fake_viewer` and `viewer_factory` of `conftest.py`. API: `Viewer(*, port=None, token=None, host="127.0.0.1", crs="", origin=(0,0,0), settings_dir=None, timeout=5.0, heartbeat_interval=5.0)`; `connect()`, `send(obj, *, id, name=None, parent=None) -> int`, `send_group`, `remove`, `clear`, `inventory() -> list[InventoryEntry]`, `flush(settle=0.5, *, raise_on_refusal=True) -> list[Refusal]`, `release_data()`, `close()`, properties `revisions`, `scene_origin`, `scene_crs`, `connected`. Tests (one per scenario of story 1): 1 discovery with no port or token argument and the object in `fake_viewer.viewer.objects`; 2 same id twice gives revisions 1 then 2 and the stored header has revision 2; 3 no viewer gives `NoViewerError` and no process starts; 4 properties are the same names and values in the stored buffers; 5 one object of each of the 7 types gives 7 objects with the right `object_type`; 6 a `Vectors.from_directions` object gives stored `dip_direction` and `dip` equal to the expected values; 7 a second `Viewer` gets `BusyError` and the objects of the first do not change, then after `release_data()` of the first the second sends and the objects of the first stay. Also: a wrong token gives `TokenError`; `reply_errors` makes `flush()` raise `ViewerRefusedError` with a `Refusal` that has the id; the scene origin of the viewer differs from the proposed origin: coordinates use the scene origin; a viewer that closes in the middle of a send (`FakeViewer.close()` in a thread) gives `ViewerClosedError`; keepalive: `FakeViewer(heartbeat_timeout=0.6)` with `heartbeat_interval=0.1` stays connected for 2 s and a later `send` works; an `error` code `unsupported_version` gives `MissingFeatureError`.

- [ ] T021 [Found] [M] Implement `Viewer` in `viewer.py`.
  - Files: `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T020, T014, T015, T017
  - Wave: 8
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_viewer.py` passes.
  - Context: Class and method signatures as in T020 (also `@dataclass(frozen=True) class InventoryEntry: id: str; object_type: str; revision: int`). `connect()`: `find_viewer(...)`, `WebSocketClient.connect(host, port, token=token, timeout=timeout, path="/")`, send `hello_frame(crs, origin, __version__)`, read until both `hello` and `scene.inventory` arrive (else `ViewerClosedError`); an `error` frame with code `busy` raises `BusyError` and closes the socket; `unsupported_version` raises `MissingFeatureError`; other `error` raises `ViewerRefusedError`. Scene origin: the viewer `hello` origin when not null, else the proposed origin; keep `scene_crs` from the viewer `hello`. All frames use the scene origin. Revision of an id: `max(local, inventory) + 1`, 1 for a new id. `send(obj, ...)` calls `adapters.convert` only when `obj` is not a `SceneObject` (T027 adds the hook; before that, a non-`SceneObject` raises `InvalidInputError` that names the accepted inputs). After each frame, drain the socket with a 1 ms wait and keep each `error` frame that has an `object` as a `Refusal`. `flush` waits up to `settle` for late frames. Keepalive: daemon thread `loop-link-heartbeat`; every `heartbeat_interval` it tries a non-blocking `lock.acquire`; if the last send was older than the interval, it sends `heartbeat_frame(seq)`; an `OSError` stops the thread and marks the connection closed. `release_data` and `close` set a `threading.Event`, join with timeout `2*heartbeat_interval`, close the socket, and keep `revisions`. OP2: default `crs=""`; first read the hello handling of the geoviewer repository (`C:\Users\lachl\dev\geoviewer\crates`, search `Hello`) and, if the viewer refuses an empty CRS, raise `InvalidInputError("give crs=...")` for an empty scene. Convert `OSError` and `WebSocketClosed` on send to `ViewerClosedError`.

- [ ] T023 [US1] [S] Export the public names in `__init__.py`.
  - Files: `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "import loop_link as l; [getattr(l, n) for n in l.__all__]"` passes and `ruff check packages/loop_link` passes.
  - Context: Keep `__version__ = "0.1.0"`. Import and list in `__all__`: from `.viewer`: `Viewer`, `InventoryEntry`; from `.objects`: `Points`, `Vectors`, `Lines`, `Surface`, `BlockModel`, `TetraMesh`, `Drillhole`, `Categorical`; from `.errors`: `LoopLinkError`, `NoViewerError`, `TokenError`, `BusyError`, `MissingFeatureError`, `InvalidInputError`, `ObjectTooLargeError`, `ViewerClosedError`, `ViewerRefusedError`, `MissingExtraError`, `Refusal`. Do not export `PushResult` or `push_objects` yet (T036). Do not import pyvista, pandas, Qt or LoopStructural here.

- [ ] T024 [P] [US1] [S] Test: SC-001 speed in `test_performance.py`.
  - Files: `packages/loop_link/tests/test_performance.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_performance.py` passes.
  - Context: Mark the test `@pytest.mark.slow`. Use fixtures `fake_viewer` and `viewer_factory`. Build `BlockModel(origin=(0,0,0), cell_sizes=(1,1,1), cells=(100,100,100), cell_properties={"grade": np.random.default_rng(0).random(1_000_000, dtype=np.float32)})`. Time with `time.perf_counter()` from the start of `viewer.send(obj, id="grid")` to the moment that `fake_viewer.viewer.inventory()` has an entry with `id == "grid"` (poll every 10 ms, give up after 10 s). Assert the time is below 5.0 s (SC-001). The stored buffer must equal the input (`np.array_equal`).

- [ ] T025 [P] [US1] [M] Test: frames of the geoviewer repository (SC-004).
  - Files: `packages/loop_link/tests/test_fixtures.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_fixtures.py -rs` passes, or skips with the reason "set LOOP_LINK_GEOVIEWER_FIXTURES"; with the variable set to `C:\Users\lachl\dev\geoviewer\fixtures\protocol\python` it passes with no skip.
  - Context: Use the helper `geoviewer_fixtures()` of `conftest.py` (OP8). `pytest.skip("set LOOP_LINK_GEOVIEWER_FIXTURES to the folder fixtures/protocol/python of the geoviewer repository")` when it is `None`. Test 1: for each `*.bin` file in the folder, `FakeViewer.accept_frame(bytes)` returns no `error` reply (use a fresh `FakeViewer` and send a `hello_frame` first). Test 2: for the file `scene_upsert_drillhole_end_of_hole.bin`, decode it with `decode_frame`, rebuild the same object with `Drillhole` and `encode_object` (read the header for `id`, `name`, collar, survey, intervals, `end_of_hole`), decode the result, and assert that the headers are equal except `version`, and the buffers are equal with `np.array_equal`. Other plugin frames: compare the header keys `type`, `object_type`, `geometry` only.

- [ ] T027 [US1] [M] Add the adapter dispatch `adapters/__init__.py` and hook it into `Viewer.send`.
  - Files: `packages/loop_link/src/loop_link/adapters/__init__.py`, `packages/loop_link/src/loop_link/viewer.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_viewer.py packages/loop_link/tests/test_import.py` passes (no regression) and `ruff check packages/loop_link` passes.
  - Context: `adapters/__init__.py`: `convert(obj: Any, **hints: Any) -> SceneObject`. Return `obj` when it is a `SceneObject`. Else choose by the first part of `type(obj).__module__`: `pyvista` -> `adapters.pyvista.from_pyvista(obj, **hints)`; `loop_common` -> `adapters.loopstructural.from_loopstructural(obj)`; import the adapter module lazily inside the branch. Any other type raises `InvalidInputError("Accepted inputs: loop_link objects (Points, Vectors, Lines, Surface, BlockModel, TetraMesh, Drillhole); pyvista meshes and grids; LoopStructural Surface and StructuredGrid; pandas tables for drillholes (use loop_link.adapters.pandas.drillholes_from_pandas).")`. The adapter modules are created by T029 and T031; `convert` imports them lazily, so a missing module raises `ImportError` only when used. Do not create stubs for them. In `viewer.py`, `Viewer.send` calls `convert(obj)` first. Add `**hints` to `Viewer.send` (`send(obj, *, id, name=None, parent=None, **hints)`) and pass them on.

- [ ] T034 [P] [US5] [M] Test: plugin push API in `test_push.py` (story 5).
  - Files: `packages/loop_link/tests/test_push.py`
  - Needs: T021
  - Wave: 9
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_push.py` fails with ImportError of `loop_link.push`, then passes after T035.
  - Context: Under test: `push_objects(viewer, objects, *, replace_prefix=None, progress=None, cancel=None) -> PushResult` and `PushResult(sent=0, removed=0, refused=[], skipped=[], cancelled=False)` with `summary() -> str`. `objects` is an iterable of `(id, object)`. Use `fake_viewer` and `viewer_factory`. Tests: (1) frames equal: push 3 objects and send the same 3 with `Viewer.send` to a second `FakeViewer`; the stored headers (without `version`) and buffers are equal (story 5, scenario 1); (2) `replace_prefix="model/"`: the fake starts with ids `model/old` and `other/keep` (send them first); after the push of `model/new` the ids `model/old` is removed, `other/keep` stays, `removed == 1`; (3) an object that raises `InvalidInputError` (for example an unsupported type `42`) goes to `skipped` as `(id, reason)` and the push goes on; (4) `cancel` set after the first object: `cancelled is True`, `sent == 1`; (5) `progress` gets the text `Sending <id>...` before each object, on the calling thread; (6) a slow viewer: `FakeViewer(read_delay=0.05)`; run `push_objects` in a `threading.Thread` for 5 objects; in the main thread count 20 ticks of 0.01 s with `time.sleep`; all 20 ticks finish while the push still runs or within 0.5 s of its start (the main thread is not blocked; story 5, scenario 2); (7) `summary()` text: "3 objects sent to GeoViewer", plus ", 1 old objects removed" and ", 1 objects skipped" when set, and ". The viewer refused N messages" when `refused` is not empty, ending with ".".

- [ ] T037 [P] [Finish] [S] Manual checks document.
  - Files: `packages/loop_link/tests/manual_checks.md`
  - Needs: T021
  - Wave: 9
  - Done when: `git grep -n -c "^## M" packages/loop_link/tests/manual_checks.md` gives 2.
  - Context: Two sections. `## M1: speed with the real headless viewer (SC-001)`: start `geoviewer-headless --port 0` with `GEOVIEWER_TOKEN=secret` (it prints the port); run a script that sends a `BlockModel` of 100 x 100 x 100 cells with one `float32` cell property to it and prints the time until the object is in `Viewer.inventory()` after a reconnect; pass if under 5 s; record the time. `## M2: seven types in the window viewer`: start the viewer with a window, turn the link on, save a token, send one object of each of the 7 types with one property each, check by eye that 7 objects with the right types show, and that the property names show in the property table. Keep each section to 8 steps or fewer, in ASD-STE100 English.

## Phase 3: Stage S2 adapters (FR-010)

- [ ] T028 [P] [US1] [M] Test: the pyvista adapter.
  - Files: `packages/loop_link/tests/test_adapter_pyvista.py`
  - Needs: T027
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pyvista.py` fails with ImportError of `loop_link.adapters.pyvista`, then passes after T029.
  - Context: `pv = pytest.importorskip("pyvista")`. Under test: `from_pyvista(mesh, *, kind=None, vectors=None, vectors_kind="normal") -> SceneObject` and `Viewer.send(mesh, id=...)`. Cases: `pv.PolyData(points)` with only vertices -> `Points`, `point_data` -> `properties`; with a point array `"n"` of shape (n,3) and `vectors="n"` -> `Vectors` (check `dip` against the formula arccos(n_z) for upward normals); `pv.lines_from_points` PolyData -> `Lines` with segments (m,2); `pv.Sphere().triangulate()` -> `Surface` (vertices, triangles (m,3), cell data -> `cell_properties`); a quad `pv.Plane()` face set is triangulated -> `Surface`; `pv.ImageData(dimensions=(4,5,6), spacing=(1,2,3), origin=(10,20,30))` -> `BlockModel(origin=(10,20,30), cell_sizes=(1,2,3), cells=(3,4,5))` with point data on nodes and cell data on cells; `pv.RectilinearGrid` with uneven spacing -> `InvalidInputError` (OP9); `pv.StructuredGrid` -> `InvalidInputError`; an `UnstructuredGrid` of all tetra cells (`pv.UnstructuredGrid(cells, [pv.CellType.TETRA], points)`) -> `TetraMesh` with tetrahedra (m,4); a mixed-cell `UnstructuredGrid` -> `InvalidInputError`; `kind="Drillhole"` -> `InvalidInputError` whose text names `Drillhole` and `drillholes_from_pandas` (OP3); a 2-D array in `point_data` of shape (n,2) -> `InvalidInputError`. Analytic grid test: `f = x + 10 y + 100 z` on the nodes of the ImageData grid; send it and read the stored buffer from the `FakeViewer`; it equals `f` at the node order x fastest.

- [ ] T029 [US1] [M] Implement the pyvista adapter.
  - Files: `packages/loop_link/src/loop_link/adapters/pyvista.py`
  - Needs: T028
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pyvista.py` passes.
  - Context: `from_pyvista(mesh, *, kind: Optional[str] = None, vectors: Optional[str] = None, vectors_kind: str = "normal") -> SceneObject`. Import `pyvista` inside the module; if it is missing raise `MissingExtraError("pyvista")` at the first call of `from_pyvista`, not at `import loop_link`. Rules: `PolyData` with `n_lines > 0` -> `Lines` (segments from the line cell array: each line of k points gives k-1 segments); `PolyData` with faces -> `Surface` (call `mesh.triangulate()` when `not mesh.is_all_triangles`; triangles `faces.reshape(-1, 4)[:, 1:]`); `PolyData` of vertices only -> `Points`, or `Vectors.from_directions(points, mesh.point_data[vectors], kind=vectors_kind)` when `vectors` is set; `ImageData` -> `BlockModel(origin=mesh.origin, cell_sizes=mesh.spacing, cells=np.array(mesh.dimensions) - 1)`; `RectilinearGrid` with equal spacing on each axis (`np.allclose(np.diff(axis), step)`) -> `BlockModel`, else `InvalidInputError` (OP9); `UnstructuredGrid` with all cells `pyvista.CellType.TETRA` -> `TetraMesh` (`cells_dict[pyvista.CellType.TETRA]`); any other type or `kind="Drillhole"` -> `InvalidInputError` that names the accepted inputs (OP3 text: "A pyvista mesh cannot hold a drillhole. Use loop_link.Drillhole or loop_link.adapters.pandas.drillholes_from_pandas."). `point_data` -> `properties`, `cell_data` -> `cell_properties`, as numpy arrays; a 2-D array raises `InvalidInputError` that names the array. Order of pyvista image data is x fastest: no reorder. Add docstrings and type hints.

- [ ] T030 [P] [US1] [M] Test: the LoopStructural adapter.
  - Files: `packages/loop_link/tests/test_adapter_loopstructural.py`
  - Needs: T027
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_loopstructural.py` fails with ImportError of `loop_link.adapters.loopstructural`, then passes after T031.
  - Context: `geometry = pytest.importorskip("loop_common.geometry")`. Classes: `loop_common.geometry.Surface` (fields `vertices`, `triangles`, `properties: dict|None`, `cell_properties: dict|None`, `name`) and `loop_common.geometry.StructuredGrid` (fields `origin`, `step_vector`, `nsteps` (node counts), `properties` (per node, x fastest), `cell_properties` (per cell, x fastest), property `nsteps_cells`). Under test: `from_loopstructural(obj)` and `Viewer.send(obj, id=...)`. Cases: `Surface` with 4 vertices, 2 triangles, `properties={"a": ...}`, `cell_properties={"b": ...}` -> `Surface` with the same arrays; `properties=None` -> empty mapping; `StructuredGrid(origin=[1,2,3], step_vector=[1,2,3], nsteps=[4,5,6])` -> `BlockModel(origin=(1,2,3), cell_sizes=(1,2,3), cells=(3,4,5))`; analytic check: set `properties["f"]` to `grid.nodes @ [1, 10, 100]` and send it; the stored buffer from the `FakeViewer` equals it (so no reorder is needed); any other class from `loop_common.geometry` (for example `ValuePoints`) -> `InvalidInputError` that names the accepted inputs.

- [ ] T031 [US1] [M] Implement the LoopStructural adapter.
  - Files: `packages/loop_link/src/loop_link/adapters/loopstructural.py`
  - Needs: T030
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_loopstructural.py` passes.
  - Context: `from_loopstructural(obj: Any) -> SceneObject`. Do not import `loop_common` at module level. Find the classes with `module = sys.modules.get("loop_common.geometry")`; if `module is None` or `obj` is not an instance of `module.Surface` or `module.StructuredGrid`, raise `InvalidInputError` that names the accepted inputs. If `loop_common` cannot be imported when the extra is needed, raise `MissingExtraError("loopstructural")`. `Surface` -> `Surface(vertices=obj.vertices, triangles=obj.triangles, properties=obj.properties or {}, cell_properties=obj.cell_properties or {})`. `StructuredGrid` -> `BlockModel(origin=obj.origin, cell_sizes=obj.step_vector, cells=obj.nsteps_cells, properties=obj.properties, cell_properties=obj.cell_properties)`. Coordinates are the public attributes (world coordinates, spec 005). No other transform. Add docstrings and type hints.

- [ ] T032 [P] [US1] [M] Test: the pandas adapter.
  - Files: `packages/loop_link/tests/test_adapter_pandas.py`
  - Needs: T027
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pandas.py` fails with ImportError of `loop_link.adapters.pandas`, then passes after T033.
  - Context: `pd = pytest.importorskip("pandas")`. Under test: `drillholes_from_pandas(collar, survey, logs=None, *, hole_id="hole_id", x="x", y="y", z="z", depth="depth", azimuth="azimuth", dip="dip", from_="from", to="to", dip_positive_down=True) -> dict[str, Drillhole]`. Data: collar `hole_id,x,y,z` for `DH1`, `DH2`; survey `hole_id,depth,azimuth,dip` (2 rows each); logs `hole_id,from,to,grade,lith` (grade float, lith strings). Cases: result keys `DH1`, `DH2`; `Drillhole.collar` equals the row; `survey` shape (2,3) in the order depth, azimuth, dip; `intervals` shape (j,2); `interval_properties["grade"]` float array and `["lith"]` a numpy `U` array; `dip_positive_down=False` changes the sign of `dip`; a custom column name through `x="east"` works; a hole in `survey` that is not in `collar` -> `InvalidInputError`; an `object` column with a non-string value in `logs` -> `InvalidInputError`; `logs=None` gives holes with 0 intervals (`intervals` shape (0,2)); the result sends with `Viewer.send(hole, id=hole_id)` to the `FakeViewer` and the stored `object_type` is `Drillhole`.

- [ ] T033 [US1] [M] Implement the pandas adapter.
  - Files: `packages/loop_link/src/loop_link/adapters/pandas.py`
  - Needs: T032
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_adapter_pandas.py` passes.
  - Context: Signature as in T032. Import `pandas` inside the function body of the first use; if it is missing raise `MissingExtraError("pandas")`. Group the three tables by `hole_id` (`groupby`, `sort=False`); sort the survey rows by `depth`; `Drillhole(collar=[x,y,z], survey=array of [depth, azimuth, dip], intervals=array of [from, to], properties={}, interval_properties={extra columns of logs})`. `dip_positive_down=False` multiplies `dip` by -1 (the `geoh5` sign is the opposite of the viewer, which has positive down). A column of strings (all `str`) becomes `np.asarray(values, dtype="U")`; a numeric or bool column keeps its dtype; any other `object` column raises `InvalidInputError` that names the column. A hole id in `survey` or `logs` that is not in `collar` raises `InvalidInputError` that names the id. Return a `dict` that keeps the order of `collar`. Add docstrings and type hints.

## Phase 4: Stage S3 plugin API (story 5)

- [ ] T035 [US5] [M] Implement `push.py`.
  - Files: `packages/loop_link/src/loop_link/push.py`
  - Needs: T034
  - Wave: 10
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_push.py` passes.
  - Context: `@dataclass class PushResult: sent: int = 0; removed: int = 0; refused: list[Refusal] = field(default_factory=list); skipped: list[tuple[str, str]] = field(default_factory=list); cancelled: bool = False` and `summary(self) -> str` (text as in the test of T034: "<sent> objects sent to GeoViewer", ", <n> old objects removed", ", <n> objects skipped", ". The viewer refused <n> messages", final "."). `push_objects(viewer: Viewer, objects: Iterable[tuple[str, Any]], *, replace_prefix: Optional[str] = None, progress: Optional[Callable[[str], None]] = None, cancel: Optional[threading.Event] = None) -> PushResult`. Steps: `viewer.connect()` if `not viewer.connected`; for each `(id, obj)`: if `cancel` is set, set `cancelled=True` and stop; call `progress(f"Sending {id}...")`; `viewer.send(obj, id=id)` and `sent += 1`; `InvalidInputError` or `ObjectTooLargeError` -> `skipped.append((id, str(error)))` and go on; `ViewerClosedError` and `BusyError` are not caught. After the loop, if `replace_prefix`, for each id in `viewer.inventory()` that starts with the prefix and was not sent in this call: `viewer.remove(id)` and `removed += 1`. Then `refused = viewer.flush(raise_on_refusal=False)`. No Qt, no module state; `progress` runs on the calling thread. Docstrings and type hints.

- [ ] T036 [US5] [S] Export `PushResult` and `push_objects`.
  - Files: `packages/loop_link/src/loop_link/__init__.py`
  - Needs: T035, T023
  - Wave: 11
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -c "from loop_link import PushResult, push_objects"` passes and `ruff check packages/loop_link` passes.
  - Context: Add `from .push import PushResult, push_objects` to `__init__.py` and add both names to `__all__`. Keep all other lines.

- [ ] T022 [Finish] [S] Test: the core import needs no Qt or optional extras (SC-005).
  - Files: `packages/loop_link/tests/test_import.py`
  - Needs: T036
  - Wave: 12
  - Done when: `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests/test_import.py` passes.
  - Context: Run in a subprocess: `subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)` where `code` imports `loop_link`, then prints `sorted(m for m in sys.modules if m.split(".")[0] in {"PyQt5","PyQt6","PySide2","PySide6","qgis","pyvista","pandas","loop_common","LoopStructural"})`. The test asserts the printed list is `[]` and the return code is 0. Add a second test: `loop_link.__all__` names all exist (`getattr`).

## Phase 5: Finish

- [ ] T038 [Finish] [M] Run the merge gates and fix the findings.
  - Files: `packages/loop_link/`
  - Needs: T002, T003, T004, T005, T022, T024, T025, T029, T031, T033, T037
  - Wave: 13
  - Done when: `ruff check LoopStructural packages && pytest` passes, `packages/loop_link/.venv/Scripts/python.exe -m pytest packages/loop_link/tests` passes, and `black --check packages/loop_link` passes.
  - Context: Constitution merge gate: `pytest` and `ruff check LoopStructural packages` must both pass. The root `pytest` does not run the package tests, so run the package command too. Run `black packages/loop_link` to fix format (line length 100, quotes kept). Fix only problems in `packages/loop_link/`. Do not change a gate or add a Ruff ignore for `packages/loop_link/src`. Run `git grep -n loop_link ROADMAP.md .release-please-config.json .release-please-manifest.json` and check 2 files or more show the name (SC-006). Report the result of each command in the commit message.

## Wave table

| Wave | Tasks | Tier [S] | Tier [M] |
|---|---|---|---|
| 1 | T001, T002, T004, T005 | 2 | 2 |
| 2 | T003, T006 | 1 | 1 |
| 3 | T007, T009, T011, T016 | 0 | 4 |
| 4 | T008, T010, T017 | 0 | 3 |
| 5 | T012, T015, T019 | 0 | 3 |
| 6 | T013, T018 | 0 | 2 |
| 7 | T014, T020 | 0 | 2 |
| 8 | T021 | 0 | 1 |
| 9 | T023, T024, T025, T027, T034, T037 | 3 | 3 |
| 10 | T028, T030, T032, T035 | 0 | 4 |
| 11 | T029, T031, T033, T036 | 1 | 3 |
| 12 | T022 | 1 | 0 |
| 13 | T038 | 0 | 1 |

## Critical path

T001 -> T006 -> T009 -> T010 -> T012 -> T013 -> T014 -> T021 -> T034 -> T035 -> T036 -> T022 -> T038

## Dependencies

Setup (T001 to T005) first. T006 (errors) opens the foundation. Then three chains in parallel:
codec (T007, T008), vector math (T009, T010), objects (T011, T012; needs T010), and discovery
(T016, T017). The chain objects, then encode (T013, T014), joins the WebSocket port (T015) and
discovery in `Viewer` (T020, T021; needs the fake viewer T019 and `conftest.py` T018). After T021:
exports (T023), speed (T024), fixtures (T025), adapter dispatch (T027), push tests (T034) and
manual checks (T037) run in parallel. S2 adapters follow T027 as test-then-code pairs
(T028 to T033). S3 is T034, T035, T036. T022 and T038 close the work.
