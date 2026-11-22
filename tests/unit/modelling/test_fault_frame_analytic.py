"""Check analytically defined fault frames in model coordinate frames."""

import numpy as np
import pytest

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    ROTATIONS_DEG,
    model_with_frame,
    r2_fault_data,
    r2_frame,
    r2_probes,
    r4_frame,
)

HALF_SIZE = np.full(3, 10.0)
NORMAL = np.array([1.0, 0.0, 0.0])
SLIP = np.array([0.0, 0.0, -1.0])


def _fault(anchor, rotation_deg, interpolatortype, explicit_center=True, slip=SLIP):
    model = model_with_frame(anchor - HALF_SIZE, anchor + HALF_SIZE, rotation_deg)
    data = r2_fault_data(anchor)
    model.data = data
    fault = model.create_and_add_fault(
        "fault",
        0,
        data=data,
        interpolatortype=interpolatortype,
        fault_center=anchor if explicit_center else None,
        fault_normal_vector=NORMAL,
        fault_slip_vector=slip,
        major_axis=20,
        minor_axis=10,
        intermediate_axis=20,
        fault_buffer=0.2,
        nelements=20_000 if interpolatortype == "FDI" else 1_000,
        fault_trace_anisotropy=0,
        fault_dip_anisotropy=0,
    )
    # These are the CG solve settings of SC-003. The builder ignores
    # top-level rtol, atol and maxiter, so they go in solver_kwargs.
    fault.builder.build(
        solver="cg",
        solver_kwargs={"rtol": 1e-10, "atol": 1e-12, "maxiter": 20_000},
        w1=1,
        w2=1,
        w3=1,
        step=10,
    )
    return fault


def test_explicit_centre_rows_not_projected():
    anchor = ANCHORS[1]
    fault = _fault(anchor, 30, "FDI")
    for builder in fault.builder.builders:
        generated = builder.data
        rows = generated.loc[generated["val"] == 0, ["X", "Y", "Z"]].to_numpy()
        assert np.any(np.all(rows == anchor, axis=1))


def test_inferred_centre_keeps_elevation():
    anchor = ANCHORS[1]
    fault = _fault(anchor, 30, "FDI", explicit_center=False)
    np.testing.assert_allclose(fault.builder.fault_centre, anchor, rtol=0, atol=1e-12)
    assert fault.builder.fault_centre[2] == anchor[2]


def test_slip_vector_made_orthogonal():
    fault = _fault(ANCHORS[0], 0, "FDI", slip=np.array([1.0, 0.0, -1.0]))
    normal = fault.builder.fault_normal_vector
    slip = fault.builder.fault_slip_vector
    np.testing.assert_allclose(slip, [0, 0, -1], rtol=0, atol=1e-12)
    assert abs(np.dot(normal, slip)) <= 1e-12
    coordinate1 = fault.builder.builders[1].data
    gradients = coordinate1[["nx", "ny", "nz"]].to_numpy()
    # The coordinate 1 gradient is the unit in-plane slip divided by the intermediate axis.
    expected = np.array([0.0, 0.0, -1.0]) / fault.builder.fault_intermediate_axis
    assert np.any(np.all(np.isclose(gradients, expected, rtol=0, atol=1e-12), axis=1))


def test_parallel_slip_raises():
    model = model_with_frame(ANCHORS[0] - HALF_SIZE, ANCHORS[0] + HALF_SIZE, 0)
    data = r2_fault_data(ANCHORS[0])
    with pytest.raises(ValueError):
        model.create_and_add_fault(
            "fault",
            0,
            data=data,
            fault_center=ANCHORS[0],
            fault_normal_vector=NORMAL,
            fault_slip_vector=NORMAL,
            major_axis=20,
            minor_axis=10,
            intermediate_axis=20,
        )


def test_extent_direction_row_at_centre():
    anchor = ANCHORS[0]
    fault = _fault(anchor, 0, "FDI", explicit_center=False)
    coordinate2 = fault.builder.builders[2].data
    directions = coordinate2[["nx", "ny", "nz"]].to_numpy()
    finite = np.all(np.isfinite(directions), axis=1)
    assert np.count_nonzero(finite) == 1
    assert np.all(coordinate2.loc[finite, ["X", "Y", "Z"]].to_numpy()[0] == anchor)
    row = coordinate2.loc[finite]
    np.testing.assert_allclose(
        row[["nx", "ny", "nz"]].to_numpy()[0], [0, 1 / 20, 0], rtol=0, atol=1e-10
    )
    _, gradients = r4_frame(anchor[None, :], anchor)
    np.testing.assert_allclose(
        row[["nx", "ny", "nz"]].to_numpy()[0], gradients[2][0], rtol=0, atol=1e-10
    )
    points = anchor + np.array([[0, -6, 0], [0, 6, 0]])
    directions = r4_frame(points, anchor)[1][2]
    cosine = (
        np.dot(directions[0], directions[1])
        / np.linalg.norm(directions[0])
        / np.linalg.norm(directions[1])
    )
    assert np.rad2deg(np.arccos(np.clip(cosine, -1, 1))) >= 15


@pytest.mark.parametrize("anchor", ANCHORS.tolist())
@pytest.mark.parametrize("rotation_deg", ROTATIONS_DEG)
@pytest.mark.parametrize("explicit_center", [True, False])
@pytest.mark.parametrize("interpolatortype", ["FDI", "PLI"])
def test_planar_fault_frame_matches_analytic(
    anchor, rotation_deg, explicit_center, interpolatortype
):
    anchor = np.asarray(anchor)
    fault = _fault(anchor, rotation_deg, interpolatortype, explicit_center)
    points = r2_probes(anchor)
    expected, expected_gradients = r2_frame(points, anchor)
    actual = np.column_stack([fault[i].evaluate_value(points) for i in range(3)])
    actual_gradients = np.stack([fault[i].evaluate_gradient(points) for i in range(3)], axis=1)
    assert np.all(np.isfinite(actual))
    assert np.all(np.isfinite(actual_gradients))
    assert np.max(np.abs(actual - expected)) <= 1e-6
    expected_magnitude = np.linalg.norm(expected_gradients, axis=2)
    actual_magnitude = np.linalg.norm(actual_gradients, axis=2)
    assert (
        np.max(np.abs(actual_magnitude - expected_magnitude)) / np.max(expected_magnitude) <= 1e-4
    )
    dot = np.sum(actual_gradients * expected_gradients, axis=2)
    denom = np.linalg.norm(actual_gradients, axis=2) * expected_magnitude
    angles = np.rad2deg(np.arccos(np.clip(dot / denom, -1, 1)))
    assert np.all(np.isfinite(angles))
    assert np.max(angles) <= 0.01


def test_caller_normal_not_mutated():
    # The builder must copy the normal vector of the caller.
    anchor = ANCHORS[1]
    model = model_with_frame(anchor - HALF_SIZE, anchor + HALF_SIZE, 30)
    data = r2_fault_data(anchor)
    model.data = data
    normal = np.array([2.0, 0.0, 0.0])
    slip = np.array([0.0, 0.0, -2.0])
    model.create_and_add_fault(
        "fault",
        0,
        data=data,
        interpolatortype="FDI",
        fault_center=anchor,
        fault_normal_vector=normal,
        fault_slip_vector=slip,
        major_axis=20,
        minor_axis=10,
        intermediate_axis=20,
        fault_buffer=0.2,
        nelements=20_000,
        fault_trace_anisotropy=0,
        fault_dip_anisotropy=0,
    )
    np.testing.assert_array_equal(normal, [2.0, 0.0, 0.0])
    np.testing.assert_array_equal(slip, [0.0, 0.0, -2.0])


def test_fault_ellipsoid_centre_is_world(monkeypatch):
    # fault_centre is a world point, so the ellipsoid must not convert it again.
    import sys
    import types

    received = {}

    class FakePolyData:
        def __init__(self, pts):
            received["pts"] = np.array(pts, dtype=float)
            self._data = {}

        def __setitem__(self, key, value):
            self._data[key] = value

        def __getitem__(self, key):
            return self._data[key]

        def glyph(self, **kwargs):
            return self

    fake = types.ModuleType("pyvista")
    fake.PolyData = FakePolyData
    fake.ParametricEllipsoid = lambda *args: object()
    monkeypatch.setitem(sys.modules, "pyvista", fake)

    anchor = ANCHORS[1]
    model = model_with_frame(anchor - HALF_SIZE, anchor + HALF_SIZE, 30)
    data = r2_fault_data(anchor)
    model.data = data
    fault = model.create_and_add_fault(
        "fault",
        0,
        data=data,
        interpolatortype="FDI",
        fault_center=anchor,
        fault_normal_vector=NORMAL,
        fault_slip_vector=SLIP,
        major_axis=20,
        minor_axis=10,
        intermediate_axis=20,
        fault_buffer=0.2,
        nelements=20_000,
        fault_trace_anisotropy=0,
        fault_dip_anisotropy=0,
    )
    fault.fault_ellipsoid()
    np.testing.assert_allclose(received["pts"][0], anchor, atol=1e-9)
