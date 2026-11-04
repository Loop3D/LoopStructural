"""Frame keys, world-value round trips and repeated save/load cycles (U3.1, U3.2, U3.8; SC-006)."""

import numpy as np
import pytest
from loop_interpolation import GeologicalInterpolator

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    PERSIST_ENGINES,
    QUERY_OFFSETS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    r1_interpolator,
    r1_value,
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
    """A loaded interpolator gives the same world values and gradients."""
    interp, anchor = _build(engine, anchor_index, rotation)
    loaded = GeologicalInterpolator.from_dict(interp.to_dict())
    points = anchor + QUERY_OFFSETS
    scale = max(1.0, np.abs(r1_value(points, anchor)).max())
    value_error = np.abs(loaded.evaluate_value(points) - r1_value(points, anchor)).max()
    gradient_error = np.abs(loaded.evaluate_gradient(points) - R1_GRADIENT).max()
    assert value_error <= 1e-12 * scale
    assert gradient_error <= 1e-12 * max(1.0, np.abs(R1_GRADIENT).max())


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
