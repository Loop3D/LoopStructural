"""Give the frames of loop_link to the real geoviewer codec.

The codec is the Rust crate ``geoviewer-protocol``. The geoviewer repository is private, so
these tests are optional. Set ``LOOP_LINK_GEOVIEWER_REPO`` to a geoviewer checkout. The tests
skip if the variable is not set or if ``cargo`` is not on the PATH.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from loop_link import encode
from loop_link.objects import (
    BlockModel,
    Categorical,
    Drillhole,
    Lines,
    Points,
    Surface,
    TetraMesh,
    Vectors,
)

_PACKAGE_DIR = Path(__file__).resolve().parent.parent
_SOURCE_DIR = Path(__file__).resolve().parent / "codec_validator"
_BUILD_DIR = _PACKAGE_DIR / ".codec-validator-build"
_BUILD_TIMEOUT_SECONDS = 15 * 60
_ENV_NAME = "LOOP_LINK_GEOVIEWER_REPO"


@pytest.fixture(scope="session")
def validator_exe() -> Path:
    """Build the codec validator and return the path of the program."""
    repo = os.environ.get(_ENV_NAME)
    if not repo:
        pytest.skip(f"The variable {_ENV_NAME} is not set. It must name a geoviewer checkout.")
    if shutil.which("cargo") is None:
        pytest.skip("The program 'cargo' is not on the PATH.")
    protocol = Path(repo).resolve() / "crates" / "geoviewer-protocol"
    if not protocol.is_dir():
        pytest.fail(f"The folder {protocol} does not exist. Check {_ENV_NAME}.")
    (_BUILD_DIR / "src").mkdir(parents=True, exist_ok=True)
    template = (_SOURCE_DIR / "Cargo.toml.template").read_text(encoding="utf-8")
    manifest = template.replace("@GEOVIEWER_PROTOCOL_PATH@", protocol.as_posix())
    (_BUILD_DIR / "Cargo.toml").write_text(manifest, encoding="utf-8")
    shutil.copyfile(_SOURCE_DIR / "main.rs", _BUILD_DIR / "src" / "main.rs")
    result = subprocess.run(
        # The flag keeps the target folder in the build folder, also if the user config of
        # cargo names another target folder.
        ["cargo", "build", "--release", "--target-dir", str(_BUILD_DIR / "target")],
        cwd=_BUILD_DIR,
        check=False,
        capture_output=True,
        text=True,
        timeout=_BUILD_TIMEOUT_SECONDS,
    )
    if result.returncode != 0:
        pytest.fail(f"cargo build failed:\n{result.stdout}\n{result.stderr}")
    suffix = ".exe" if sys.platform == "win32" else ""
    exe = _BUILD_DIR / "target" / "release" / f"codec-validator{suffix}"
    if not exe.is_file():
        pytest.fail(f"The program {exe} does not exist after the build.")
    return exe


def _encode(obj: object, **kwargs: object) -> bytes:
    """Encode an object with default arguments."""
    arguments: dict[str, object] = {
        "id": "obj",
        "name": "Object",
        "parent": None,
        "revision": 1,
        "origin": (0.0, 0.0, 0.0),
    }
    arguments.update(kwargs)
    return encode.encode_object(obj, **arguments)  # type: ignore[arg-type]


def _cases() -> dict[str, bytes]:
    """Return the frames to check. The key is the name of the case."""
    rng = np.random.default_rng(1)
    pts = rng.uniform(0.0, 10.0, size=(6, 3))
    tri = [[0, 1, 2], [1, 2, 3], [3, 4, 5]]
    angles = (np.linspace(0, 300, 6), np.linspace(10, 80, 6))
    cases: dict[str, bytes] = {}
    cases["points"] = _encode(Points(pts))
    cases["points_named_parent"] = _encode(Points(pts), parent="group", revision=7)
    cases["points_properties"] = _encode(
        Points(
            pts,
            properties={
                "float": rng.uniform(size=6),
                "float32": rng.uniform(size=6).astype(np.float32),
                "integer": np.arange(6),
                "int32": np.arange(6, dtype=np.int32),
                "unsigned": np.arange(6, dtype=np.uint8),
                "text": np.array(["a", "b", "a", "c", "b", "a"]),
                "bool": np.array([True, False, True, True, False, False]),
            },
        )
    )
    cases["points_categorical"] = _encode(
        Points(
            pts,
            properties={
                "unit": Categorical(
                    [0, 1, 2, 1, 0, 2],
                    ["sand", "clay", "rock"],
                    colours=[(255, 0, 0, 255), (0, 255, 0, 255), (0, 0, 255, 128)],
                ),
                "default_colours": Categorical([0, 1, 0, 1, 0, 1], ["a", "b"]),
            },
        )
    )
    cases["points_f64_far_from_origin"] = _encode(
        Points(np.array([[5000.0, 6000.0, 10.0], [5001.5, 6002.5, 12.0]])),
        origin=(1000.0, 2000.0, 0.0),
    )
    cases["points_f64_far_with_property"] = _encode(
        Points(np.array([[4500.5, 0.25, 1.0], [4600.5, 1.25, 2.0]]), properties={"v": [1.0, 2.0]}),
        origin=(0.5, 0.25, 0.0),
    )
    cases["vectors_normal"] = _encode(
        Vectors.from_directions(pts, rng.normal(size=(6, 3)), kind="normal")
    )
    cases["vectors_lineation"] = _encode(
        Vectors.from_directions(pts, rng.normal(size=(6, 3)), kind="lineation")
    )
    cases["vectors_orientation_angles"] = _encode(Vectors(pts, orientation=angles))
    cases["vectors_lineation_angles"] = _encode(Vectors(pts, lineation=angles))
    cases["vectors_both"] = _encode(
        Vectors(
            pts,
            orientation=angles,
            lineation=(np.linspace(5, 305, 6), np.linspace(1, 50, 6)),
            properties={"w": np.arange(6.0)},
        )
    )
    cases["lines"] = _encode(Lines(pts, [[0, 1], [1, 2], [3, 4]]))
    cases["lines_properties"] = _encode(
        Lines(
            pts,
            [[0, 1], [1, 2], [3, 4]],
            properties={"v": np.arange(6.0)},
            cell_properties={"c": np.arange(3), "t": ["x", "y", "x"]},
        )
    )
    cases["surface"] = _encode(Surface(pts, tri))
    cases["surface_properties"] = _encode(
        Surface(
            pts,
            tri,
            properties={"v": np.arange(6.0)},
            cell_properties={
                "c": np.arange(3.0),
                "unit": Categorical([0, 1, 0], ["a", "b"]),
            },
        )
    )
    cases["block_model"] = _encode(BlockModel([0, 0, 0], [1, 1, 1], [2, 3, 4]))
    cases["block_model_properties"] = _encode(
        BlockModel(
            [100.0, 200.0, 0.0],
            [10.0, 10.0, 5.0],
            [2, 3, 4],
            properties={"node": np.arange(3 * 4 * 5, dtype=float)},
            cell_properties={
                "cell": np.arange(24.0),
                "label": Categorical(np.arange(24) % 3, ["a", "b", "c"]),
            },
        ),
        origin=(100.0, 200.0, 0.0),
    )
    cases["tetra_mesh"] = _encode(TetraMesh(pts[:5], [[0, 1, 2, 3], [1, 2, 3, 4]]))
    cases["tetra_mesh_properties"] = _encode(
        TetraMesh(
            pts[:5],
            [[0, 1, 2, 3], [1, 2, 3, 4]],
            properties={"v": np.arange(5.0)},
            cell_properties={"c": np.arange(2.0), "id": np.arange(2)},
        )
    )
    survey = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, -10.0], [1.0, 0.0, -20.0], [2.0, 0.0, -30.0]])
    cases["drillhole_vertex_properties"] = _encode(
        Drillhole(
            [100.0, 200.0, 50.0],
            survey,
            [[0.0, 10.0], [10.0, 30.0]],
            properties={"azimuth": np.arange(4.0)},
        ),
        origin=(100.0, 200.0, 0.0),
    )
    cases["drillhole_interval_properties"] = _encode(
        Drillhole(
            [100.0, 200.0, 50.0],
            survey,
            [[0.0, 10.0], [10.0, 20.0], [20.0, 30.0]],
            end_of_hole=30.0,
            properties={"azimuth": np.arange(4.0)},
            interval_properties={
                "grade": np.array([0.1, 0.5, 0.2]),
                "rock": Categorical([0, 1, 0], ["shale", "chert"]),
                "note": ["a", "b", "c"],
                "count": np.arange(3),
            },
        ),
        origin=(100.0, 200.0, 0.0),
    )
    cases["drillhole_no_interval_properties"] = _encode(
        Drillhole([0.0, 0.0, 0.0], survey, [[0.0, 30.0]])
    )
    cases["group"] = encode.group_frame("grp", "Group", None, 1)
    cases["group_with_parent"] = encode.group_frame("grp2", "Group 2", "grp", 2)
    cases["remove"] = encode.remove_frame("obj")
    cases["clear"] = encode.clear_frame()
    cells = np.array([1, 2, 3])
    nodes = np.array([4])
    cases["selection_set"] = encode.selection_set_frame(3, "grid", cells, nodes, 2)
    cases["selection_set_no_revision"] = encode.selection_set_frame(3, "grid", cells, nodes, None)
    cases["selection_set_cells_only"] = encode.selection_set_frame(3, "grid", cells, None, 2)
    cases["selection_set_clear"] = encode.selection_set_frame(3, "grid", None, None, None)
    cases["selection_get"] = encode.selection_get_frame(4)
    cases["object_fetch"] = encode.fetch_frame(5, "grid")
    cases["command_cancel"] = encode.cancel_frame(6)
    cases["commands_list"] = encode.commands_list_frame(1)
    filter_edits = {
        "add": {"kind": "add", "array": "grade"},
        "set_range": {
            "kind": "set_range",
            "array": "grade",
            "min": "0.5",
            "max": "1.0",
            "keep_nan": False,
        },
        "set_values": {"kind": "set_values", "array": "unit", "values": [0, 2]},
        "set_on": {"kind": "set_on", "array": "grade", "on": True},
        "set_exclude": {"kind": "set_exclude", "array": "unit", "exclude": True},
        "remove": {"kind": "remove", "array": "grade"},
        "cancel": {"kind": "cancel"},
    }
    for name, edit in filter_edits.items():
        cases[f"command_filter_{name}"] = encode.command_frame(
            2, "filter.edit", {"object": "grid", "edit": edit}
        )
    colour_commands = {
        "command_colour_set_active": ("colour.set_active", {"id": "grid", "array": "grade"}),
        "command_colour_set_map": ("colour.set_map", {"id": "grid", "name": "coolwarm"}),
        "command_colour_set_range": (
            "colour.set_range",
            {"id": "grid", "min": "0.0", "max": "1.0"},
        ),
        "command_colour_clear_range": ("colour.clear_range", {"id": "grid"}),
        "command_colour_set_scale": ("colour.set_scale", {"object": "grid", "scale": "log"}),
        "command_object_set_colour": ("object.set_colour", {"id": "grid", "colour": [1, 2, 3]}),
        "command_object_set_colour_none": ("object.set_colour", {"id": "grid", "colour": None}),
    }
    for case, (command, arguments) in colour_commands.items():
        cases[case] = encode.command_frame(2, command, arguments)
    return cases


_CASES = _cases()


def _run(exe: Path, paths: list[Path]) -> subprocess.CompletedProcess[str]:
    """Run the validator over the files."""
    return subprocess.run(
        [str(exe), *map(str, paths)], check=False, capture_output=True, text=True, timeout=120
    )


def test_real_codec_accepts_all_frames(validator_exe: Path, tmp_path: Path) -> None:
    """The real codec accepts each frame that loop_link makes."""
    paths: dict[str, Path] = {}
    for case, frame in _CASES.items():
        paths[case] = tmp_path / f"{case}.frame"
        paths[case].write_bytes(frame)
    result = _run(validator_exe, list(paths.values()))
    report = f"exit code {result.returncode}\n{result.stdout}{result.stderr}"
    lines = result.stdout.splitlines()
    assert result.returncode == 0, report
    assert len(lines) == len(paths), report
    for case, path in paths.items():
        assert f"OK {path}" in lines, f"{case} is not accepted.\n{report}"


@pytest.mark.parametrize("case", sorted(_CASES))
def test_real_codec_accepts_frame(validator_exe: Path, tmp_path: Path, case: str) -> None:
    """The real codec accepts one frame. A failure names the case."""
    path = tmp_path / f"{case}.frame"
    path.write_bytes(_CASES[case])
    result = _run(validator_exe, [path])
    assert result.returncode == 0, f"{case}: {result.stdout}{result.stderr}"
    assert result.stdout.strip() == f"OK {path}"


def test_real_codec_rejects_truncated_frame(validator_exe: Path, tmp_path: Path) -> None:
    """Negative control: the validator reports FAIL for a frame with a short buffer."""
    frame = _CASES["points_properties"]
    good = tmp_path / "good.frame"
    bad = tmp_path / "truncated.frame"
    good.write_bytes(frame)
    bad.write_bytes(frame[:-5])
    result = _run(validator_exe, [good, bad])
    lines = result.stdout.splitlines()
    assert result.returncode == 1, f"{result.stdout}{result.stderr}"
    assert f"OK {good}" in lines
    assert any(line.startswith(f"FAIL {bad}:") for line in lines), result.stdout


_NULL_NAME_FRAMES = {
    "object": lambda: _encode(Points(np.zeros((1, 3))), name=None),
    "group": lambda: encode.group_frame("grp", None, None, 1),
}


@pytest.mark.parametrize("kind", sorted(_NULL_NAME_FRAMES))
def test_real_codec_null_name(validator_exe: Path, tmp_path: Path, kind: str) -> None:
    """The real codec accepts a frame that was made with name=None."""
    path = tmp_path / "null_name.frame"
    path.write_bytes(_NULL_NAME_FRAMES[kind]())
    assert _run(validator_exe, [path]).returncode == 0


def test_validator_usage_error(validator_exe: Path) -> None:
    """The validator exits with code 2 if it gets no file."""
    result = _run(validator_exe, [])
    assert result.returncode == 2
