"""Frame keys, world-value round trips and repeated save/load cycles (U3.1, U3.2, U3.8; SC-006)."""

import copy

import numpy as np
import pytest
from loop_interpolation import GeologicalInterpolator, InterpolatorFactory, InterpolatorFormatError

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    PERSIST_ENGINES,
    QUERY_OFFSETS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    r1_interpolator,
    r1_value,
    rot_z,
    strip_frame_keys,
)

CELLS = [
    pytest.param(engine, anchor_index, rotation, id=f"{engine}-a{anchor_index}-r{int(rotation)}")
    for engine in PERSIST_ENGINES
    for anchor_index in range(len(ANCHORS))
    for rotation in ROTATIONS_DEG
]


def _build(engine, anchor_index, rotation):
    """Build one exact R1 interpolator for a cell."""
    anchor = ANCHORS[anchor_index]
    interp, _ = r1_interpolator(engine, anchor, rotation, exact=True)
    return interp, anchor


def _assert_same(left, right):
    """Compare two payloads that can hold arrays, dicts and lists."""
    if isinstance(left, dict):
        assert isinstance(right, dict)
        assert list(left) == list(right)
        for key in left:
            _assert_same(left[key], right[key])
    elif isinstance(left, (list, tuple)) and not np.isscalar(left):
        assert len(left) == len(right)
        for a, b in zip(left, right):
            _assert_same(a, b)
    else:
        np.testing.assert_array_equal(left, right)


@pytest.mark.parametrize("engine, anchor_index, rotation", CELLS)
def test_payload_declares_frame(engine, anchor_index, rotation):
    """The payload carries the format version, the space, the transform and the box."""
    interp, _ = _build(engine, anchor_index, rotation)
    payload = interp.to_dict()
    box = interp.bounding_box
    assert payload["format_version"] == 2
    assert payload["constraint_space"] == "local"
    transform = payload["transform"]
    assert transform["dimensions"] == 3
    np.testing.assert_array_equal(transform["local_origin"], box.local_origin)
    np.testing.assert_array_equal(transform["local_rotation"], box.local_rotation)
    box_payload = payload["bounding_box"]
    assert {"origin", "maximum", "nsteps"} <= set(box_payload)
    assert not any(key.startswith("local_") for key in box_payload)
    keys = list(payload)
    assert keys.index("format_version") < keys.index("constraint_space")
    assert keys.index("format_version") < keys.index("transform")
    assert keys.index("format_version") < keys.index("bounding_box")


@pytest.mark.parametrize("engine, anchor_index, rotation", CELLS)
def test_round_trip_world_values(engine, anchor_index, rotation):
    """Check that a loaded interpolator matches its source within the SC-006 limits."""
    interp, anchor = _build(engine, anchor_index, rotation)
    loaded = GeologicalInterpolator.from_dict(interp.to_dict())
    _assert_world_values(loaded, interp, anchor)


@pytest.mark.parametrize("engine, anchor_index, rotation", CELLS)
def test_local_arrays_unchanged_after_two_cycles(engine, anchor_index, rotation):
    """Two save/load cycles keep the local arrays, the payload and the frame."""
    interp, _ = _build(engine, anchor_index, rotation)
    first_payload = interp.to_dict()
    first = GeologicalInterpolator.from_dict(first_payload)
    second_payload = first.to_dict()
    second = GeologicalInterpolator.from_dict(second_payload)
    for key in interp.data:
        np.testing.assert_array_equal(second.data[key], interp.data[key])
    assert second_payload["format_version"] == first_payload["format_version"]
    _assert_same(second_payload["transform"], first_payload["transform"])
    assert second_payload["constraint_space"] == first_payload["constraint_space"]
    assert first_payload["format_version"] == 2
    _assert_same(second_payload, first_payload)


@pytest.mark.parametrize("engine", PERSIST_ENGINES)
def test_json_and_yaml_carry_frame(engine):
    """JSON and YAML text keep the transform of the source."""
    interp, _ = _build(engine, 3, 30.0)
    source = interp.to_dict()["transform"]
    assert source["dimensions"] == 3
    for loaded in (
        GeologicalInterpolator.from_json(interp.to_json()),
        GeologicalInterpolator.from_yaml(interp.to_yaml()),
    ):
        _assert_same(loaded.to_dict()["transform"], source)


ENGINE_CELLS = [pytest.param(engine, id=engine) for engine in PERSIST_ENGINES]


def _with_frame_keys(payload, interp):
    """Return a copy of a payload with the frame keys that a v2 `to_dict` writes."""
    result = copy.deepcopy(payload)
    box = interp.bounding_box
    result["format_version"] = 2
    result["constraint_space"] = "local"
    result["transform"] = {
        "dimensions": 3,
        "local_origin": np.asarray(box.local_origin).tolist(),
        "local_rotation": np.asarray(box.local_rotation).tolist(),
    }
    result["bounding_box"] = {
        "origin": np.asarray(box.origin).tolist(),
        "maximum": np.asarray(box.maximum).tolist(),
        "nsteps": np.asarray(box.nsteps).tolist(),
    }
    return result


def _assert_format_error(excinfo, reason, prefix):
    """Check the reason and the start of the message of a format error."""
    assert excinfo.value.reason == reason
    assert str(excinfo.value).startswith(prefix)


def _assert_world_values(loaded, source, anchor):
    """Check that a loaded interpolator matches its source within the SC-006 limits."""
    points = anchor + QUERY_OFFSETS
    scale = max(1.0, np.abs(r1_value(points, anchor)).max())
    value_error = np.abs(loaded.evaluate_value(points) - source.evaluate_value(points)).max()
    gradient_error = np.abs(
        loaded.evaluate_gradient(points) - source.evaluate_gradient(points)
    ).max()
    assert value_error <= 1e-12 * scale
    assert gradient_error <= 1e-12 * max(1.0, np.abs(R1_GRADIENT).max())


@pytest.mark.parametrize("engine", ENGINE_CELLS)
def test_legacy_payload_requires_choice(engine):
    """A legacy payload with no choice raises `legacy_no_transform`."""
    interp, _ = _build(engine, 1, 30.0)
    legacy = strip_frame_keys(interp.to_dict())
    before = copy.deepcopy(legacy)
    with pytest.raises(InterpolatorFormatError) as excinfo:
        GeologicalInterpolator.from_dict(legacy)
    _assert_format_error(
        excinfo, "legacy_no_transform", "legacy interpolator payload has no transform"
    )
    assert "local_origin" in str(excinfo.value)
    assert "declare_identity" in str(excinfo.value)
    _assert_same(legacy, before)


@pytest.mark.parametrize("engine, anchor_index, rotation", CELLS)
def test_legacy_with_supplied_transform(engine, anchor_index, rotation):
    """A legacy payload loads with a supplied transform and gives the source values."""
    interp, anchor = _build(engine, anchor_index, rotation)
    box = interp.bounding_box
    legacy = strip_frame_keys(interp.to_dict())
    before = copy.deepcopy(legacy)
    loaded = GeologicalInterpolator.from_dict(
        legacy, local_origin=box.local_origin, local_rotation=box.local_rotation
    )
    _assert_world_values(loaded, interp, anchor)
    _assert_same(legacy, before)


@pytest.mark.parametrize("engine", ENGINE_CELLS)
def test_legacy_with_identity_declaration(engine):
    """A legacy payload loads with `declare_identity=True` for an identity frame."""
    interp, anchor = _build(engine, 0, 0.0)
    legacy = strip_frame_keys(interp.to_dict())
    before = copy.deepcopy(legacy)
    loaded = GeologicalInterpolator.from_dict(legacy, declare_identity=True)
    _assert_world_values(loaded, interp, anchor)
    _assert_same(legacy, before)


CONFLICT_CASES = [
    pytest.param(
        0.0, {"declare_identity": True, "local_origin": np.zeros(3)}, id="identity-and-origin"
    ),
    pytest.param(0.0, {"local_origin": np.array([1.0, 2.0, 3.0])}, id="other-origin"),
    pytest.param(0.0, {"local_rotation": rot_z(10.0)}, id="other-rotation"),
    pytest.param(30.0, {"declare_identity": True}, id="identity-on-rotated-payload"),
]


@pytest.mark.parametrize("rotation, kwargs", CONFLICT_CASES)
def test_conflicting_inputs(rotation, kwargs):
    """Inputs that conflict with the transform raise `conflict`."""
    interp, _ = _build("fdi", 0, rotation)
    payload = _with_frame_keys(interp.to_dict(), interp)
    if "local_origin" in kwargs and "declare_identity" in kwargs:
        payload = strip_frame_keys(payload)
    before = copy.deepcopy(payload)
    with pytest.raises(InterpolatorFormatError) as excinfo:
        GeologicalInterpolator.from_dict(payload, **kwargs)
    _assert_format_error(excinfo, "conflict", "conflicting transform")
    _assert_same(payload, before)


def _drop_transform(payload):
    """Remove the transform."""
    del payload["transform"]


def _world_space(payload):
    """Set the constraint space to world."""
    payload["constraint_space"] = "world"


def _scale_rotation(payload):
    """Scale the rotation matrix by 2."""
    rotation = np.asarray(payload["transform"]["local_rotation"])
    payload["transform"]["local_rotation"] = (2.0 * rotation).tolist()


def _nan_origin(payload):
    """Put NaN in the origin."""
    payload["transform"]["local_origin"][0] = float("nan")


def _short_origin(payload):
    """Cut the origin to length 2."""
    payload["transform"]["local_origin"] = payload["transform"]["local_origin"][:2]


def _no_box(payload):
    """Remove the bounding box."""
    payload["bounding_box"] = None


def _legacy_keeps_transform(payload):
    """Make a legacy payload that keeps its transform."""
    del payload["format_version"]


def _string_version(payload):
    """Write the format version as a string."""
    payload["format_version"] = "2"


MALFORMED_CASES = [
    pytest.param(_drop_transform, id="missing-transform"),
    pytest.param(_world_space, id="constraint-space-world"),
    pytest.param(_scale_rotation, id="rotation-scaled-by-2"),
    pytest.param(_nan_origin, id="nan-origin"),
    pytest.param(_short_origin, id="origin-length-2"),
    pytest.param(_no_box, id="no-box-with-rotated-transform"),
    pytest.param(_legacy_keeps_transform, id="legacy-payload-keeps-transform"),
    pytest.param(_string_version, id="format-version-string"),
]


@pytest.mark.parametrize("mutate", MALFORMED_CASES)
def test_malformed_metadata(mutate):
    """A payload with bad frame metadata raises `malformed`."""
    interp, _ = _build("fdi", 1, 30.0)
    payload = _with_frame_keys(interp.to_dict(), interp)
    mutate(payload)
    before = copy.deepcopy(payload)
    with pytest.raises(InterpolatorFormatError) as excinfo:
        GeologicalInterpolator.from_dict(payload)
    assert excinfo.value.reason == "malformed"
    _assert_same(payload, before)


@pytest.mark.parametrize("version", [3, 99])
def test_unsupported_future_version(version):
    """A payload from a future format version raises `unsupported_version`."""
    interp, _ = _build("fdi", 1, 30.0)
    payload = _with_frame_keys(interp.to_dict(), interp)
    payload["format_version"] = version
    with pytest.raises(InterpolatorFormatError) as excinfo:
        GeologicalInterpolator.from_dict(payload)
    assert excinfo.value.reason == "unsupported_version"


def test_creation_call_unchanged():
    """The creation call builds an empty interpolator and refuses frame keywords."""
    interp, _ = _build("fdi", 1, 30.0)
    box = interp.bounding_box
    empty = InterpolatorFactory.from_dict({"type": "FDI", "boundingbox": box, "nelements": 1000})
    assert empty.n_g == 0
    assert empty.n_i == 0
    with pytest.raises(InterpolatorFormatError) as excinfo:
        InterpolatorFactory.from_dict(
            {"type": "FDI", "boundingbox": box, "nelements": 1000}, declare_identity=True
        )
    assert excinfo.value.reason == "malformed"


def test_old_reader_would_name_format_version():
    """An old reader that pops the known keys meets `format_version` first."""
    interp, _ = _build("fdi", 1, 30.0)
    payload = interp.to_dict()
    rest = copy.deepcopy(payload)
    for key in ("type", "support", "data", "c", "up_to_date", "valid"):
        rest.pop(key)
    assert next(iter(rest)) == "format_version"
