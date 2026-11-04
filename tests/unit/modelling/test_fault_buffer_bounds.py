"""Check buffered fault supports in rotated model coordinate frames."""

from itertools import product

import numpy as np
import pytest

from tests.fixtures.coordinate_frames import ROTATIONS_DEG, model_with_frame, r5_fault_data, tp

ORIGIN = np.array([500000.0, 6500000.0, 100.0])
HALF_SIZE = np.full(3, 10.0)
BUFFERS = (0.0, 0.2, 0.5)
ENGINES = (("FDI", 20_000), ("PLI", 1_000))


def _make_fault(buffer, rotation_deg, explicit_center, engine, nelements):
    model = model_with_frame(ORIGIN - HALF_SIZE, ORIGIN + HALF_SIZE, rotation_deg)
    model.data = r5_fault_data(ORIGIN)
    fault = model.create_and_add_fault(
        "fault",
        1.0,
        data=model.data,
        interpolatortype=engine,
        nelements=nelements,
        fault_normal_vector=np.array([1.0, 0.0, 0.0]),
        fault_slip_vector=np.array([0.0, 0.0, -1.0]),
        fault_center=ORIGIN if explicit_center else None,
        major_axis=12,
        minor_axis=4,
        intermediate_axis=12,
        force_mesh_geometry=True,
        fault_buffer=buffer,
    )
    return model, fault.builder


def _expected_bounds(model, buffer):
    radius = 10.0 + 20.0 * buffer
    world_corners = np.array(
        list(product(*[(ORIGIN[i] - radius, ORIGIN[i] + radius) for i in range(3)]))
    )
    local_corners = model.bounding_box.project(world_corners)
    return np.min(local_corners, axis=0), np.max(local_corners, axis=0), world_corners


def _supports(builder):
    return [coordinate_builder.interpolator.support for coordinate_builder in builder.builders]


def _assert_support_matches(model, builder, buffer):
    expected_origin, expected_maximum, world_corners = _expected_bounds(model, buffer)
    length = np.max(expected_maximum - expected_origin)
    world_magnitude = np.max(np.abs(world_corners))
    for coordinate_builder, support in zip(builder.builders, _supports(builder)):
        spacing = np.asarray(support.step_vector)
        tolerance = tp(length, world_magnitude)
        actual_origin = np.asarray(support.origin)
        actual_maximum = np.asarray(support.maximum)
        assert np.all(expected_origin - actual_origin >= -tolerance)
        assert np.all(expected_origin - actual_origin <= spacing + tolerance)
        assert np.all(actual_maximum - expected_maximum >= -tolerance)
        assert np.all(actual_maximum - expected_maximum <= spacing + tolerance)

        rows = coordinate_builder.data[["X", "Y", "Z"]].to_numpy(dtype=float)
        rows = rows[np.isfinite(rows).all(axis=1)]
        assert rows.size > 0
        local_rows = coordinate_builder.interpolator.bounding_box.project(rows)
        assert np.all(local_rows >= actual_origin - tolerance)
        assert np.all(local_rows <= actual_maximum + tolerance)


@pytest.mark.parametrize("engine,nelements", ENGINES)
@pytest.mark.parametrize("rotation_deg", ROTATIONS_DEG)
@pytest.mark.parametrize("explicit_center", [True, False])
@pytest.mark.parametrize("buffer", BUFFERS)
def test_buffer_bounds_match_request(buffer, explicit_center, rotation_deg, engine, nelements):
    model, builder = _make_fault(buffer, rotation_deg, explicit_center, engine, nelements)
    _assert_support_matches(model, builder, buffer)


@pytest.mark.parametrize("engine,nelements", ENGINES)
@pytest.mark.parametrize("rotation_deg", ROTATIONS_DEG)
@pytest.mark.parametrize("explicit_center", [True, False])
def test_buffers_give_distinct_bounds(explicit_center, rotation_deg, engine, nelements):
    bounds = []
    for buffer in BUFFERS:
        model, builder = _make_fault(buffer, rotation_deg, explicit_center, engine, nelements)
        supports = _supports(builder)
        bounds.append([(support.origin.copy(), support.maximum.copy()) for support in supports])
        _assert_support_matches(model, builder, buffer)

    for index in range(3):
        assert not np.allclose(bounds[0][index][0], bounds[1][index][0])
        assert not np.allclose(bounds[1][index][0], bounds[2][index][0])
        assert not np.allclose(bounds[0][index][1], bounds[1][index][1])
        assert not np.allclose(bounds[1][index][1], bounds[2][index][1])


@pytest.mark.parametrize("engine,nelements", ENGINES)
@pytest.mark.parametrize("rotation_deg", ROTATIONS_DEG)
@pytest.mark.parametrize("explicit_center", [True, False])
@pytest.mark.parametrize("buffer", BUFFERS)
def test_property_change_keeps_buffer_and_forced_geometry(
    buffer, explicit_center, rotation_deg, engine, nelements
):
    model, builder = _make_fault(buffer, rotation_deg, explicit_center, engine, nelements)
    before = [(support.origin.copy(), support.maximum.copy()) for support in _supports(builder)]

    builder.fault_normal_vector = np.array([1.0, 0.0, 0.0])

    after = _supports(builder)
    for (origin_before, maximum_before), support in zip(before, after):
        np.testing.assert_array_equal(support.origin, origin_before)
        np.testing.assert_array_equal(support.maximum, maximum_before)
    _assert_support_matches(model, builder, buffer)
