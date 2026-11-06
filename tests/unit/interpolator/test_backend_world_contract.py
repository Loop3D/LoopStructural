"""World contract of the native and Surfe adapters (FR-030, SC-018).

The adapters must send LOCAL data to the backend and return WORLD results.
These tests use recording stubs for the backend modules, so they test the Python
side only. The stubs return the exact R1 field. No test here needs a real
extension. All tests run.
"""

from __future__ import annotations

import copy
import importlib
import sys
import types

import numpy as np
import pytest
from loop_interpolation import P1Interpolator, _native

from tests.fixtures.coordinate_frames import (
    ANCHORS,
    QUERY_OFFSETS,
    R1_GRADIENT,
    ROTATIONS_DEG,
    r1_gradient,
    r1_interpolator,
    r1_value,
    rot_z,
    tp,
)

CELLS = [(a, r) for a in range(len(ANCHORS)) for r in ROTATIONS_DEG]
NATIVE_IDS = [f"native-boundary-anchor{a}-rot{int(r)}" for a, r in CELLS]
SURFE_IDS = [f"surfe-boundary-anchor{a}-rot{int(r)}" for a, r in CELLS]
UNIT_GRADIENT = R1_GRADIENT / np.linalg.norm(R1_GRADIENT)
GRADIENT_OFFSETS = np.array([(2, 3, 4), (5, 6, 7), (3, 8, 5), (6, 4, 9)], dtype=float)
LENGTH = 12.0

# Largest world error of the Surfe cells (filled by the Surfe test, reported at the end).
SURFE_ERRORS: list[float] = []


def _width(anchor: np.ndarray) -> float:
    """Return the largest absolute world coordinate of a cell."""
    return float(np.max(np.abs(anchor))) + LENGTH


def _scale(anchor: np.ndarray, query: np.ndarray) -> float:
    """Return F: the larger of 1 and the largest absolute field value at the query points."""
    return max(1.0, float(np.max(np.abs(r1_value(query, anchor)))))


def _gradient_rows(anchor: np.ndarray) -> np.ndarray:
    """Return 4 gradient rows (X, Y, Z, gx, gy, gz, w) in world coordinates."""
    points = anchor + GRADIENT_OFFSETS
    vectors = np.tile(UNIT_GRADIENT, (len(points), 1))
    return np.column_stack((points, vectors, np.ones(len(points))))


def _local_vectors_by_hand(box, anchor: np.ndarray, vectors: np.ndarray) -> np.ndarray:
    """Return local vectors as the difference of two projected points (no vector helper)."""
    start = anchor + 1.0
    return box.project(start + vectors) - box.project(start)


@pytest.fixture(scope="module", autouse=True)
def report_surfe_error():
    """Print the largest world error of the Surfe cells at the end of the module."""
    yield
    worst = max(SURFE_ERRORS) if SURFE_ERRORS else float("nan")
    print(f"\nMEASURED maximum world error of the Surfe cells: {worst:.3e}")


class _RecordingProblem:
    """Stand-in for the core class ``Problem``. It records every call."""

    instances: list = []
    field: np.ndarray | None = None

    def __init__(self, origin, nsteps, step_vector):
        self.args = (origin, nsteps, step_vector)
        self.calls: list[tuple[str, tuple]] = []
        type(self).instances.append(self)

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)

        def _record(*args):
            self.calls.append((name, args))

        return _record

    def solve(self, settings):
        self.calls.append(("solve", (settings,)))
        return type(self).field.copy(), {"converged": True, "iterations": 1, "residual": 0.0}

    def named(self, name: str) -> list[tuple]:
        return [args for call, args in self.calls if call == name]


@pytest.fixture
def stub_core(monkeypatch):
    """Put a recording module ``loop_interpolation_core`` in ``sys.modules``."""
    problem_class = type("Problem", (_RecordingProblem,), {"instances": [], "field": None})
    module = types.ModuleType("loop_interpolation_core")
    module.Problem = problem_class
    monkeypatch.setitem(sys.modules, "loop_interpolation_core", module)
    return problem_class


_ABSENT = object()


class _RecordingSurfeApi:
    """Stand-in for ``surfepy.Surfe_API``. It records every call."""

    instances: list = []
    value_function = None
    gradient_function = None

    def __init__(self, method=1):
        self.method = method
        self.calls: list[tuple[str, tuple]] = []
        type(self).instances.append(self)

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)

        def _record(*args):
            self.calls.append((name, tuple(np.array(a, copy=True) for a in args)))

        return _record

    def EvaluateInterpolantAtPoints(self, points):
        self.calls.append(("EvaluateInterpolantAtPoints", (np.array(points, copy=True),)))
        return type(self).value_function(np.asarray(points))

    def EvaluateVectorInterpolantAtPoints(self, points):
        self.calls.append(("EvaluateVectorInterpolantAtPoints", (np.array(points, copy=True),)))
        return type(self).gradient_function(np.asarray(points))

    def named(self, name: str) -> list[tuple]:
        return [args for call, args in self.calls if call == name]


@pytest.fixture
def stub_surfe():
    """Put a recording module ``surfepy`` in ``sys.modules`` and import a new Surfe wrapper.

    The fixture restores the three global states at teardown and asserts the restore.
    """
    package = importlib.import_module("loop_interpolation")
    wrapper_name = "loop_interpolation._surfe_wrapper"
    old_surfepy = sys.modules.get("surfepy", _ABSENT)
    old_wrapper = sys.modules.get(wrapper_name, _ABSENT)
    old_attribute = package.__dict__.get("_surfe_wrapper", _ABSENT)

    api_class = type(
        "Surfe_API", (_RecordingSurfeApi,), {"instances": [], "value_function": None}
    )
    module = types.ModuleType("surfepy")
    module.Surfe_API = api_class
    sys.modules["surfepy"] = module
    sys.modules.pop(wrapper_name, None)
    try:
        wrapper = importlib.import_module(wrapper_name)
        yield wrapper, api_class
    finally:
        for key, old in (("surfepy", old_surfepy), (wrapper_name, old_wrapper)):
            if old is _ABSENT:
                sys.modules.pop(key, None)
            else:
                sys.modules[key] = old
        if old_attribute is _ABSENT:
            package.__dict__.pop("_surfe_wrapper", None)
        else:
            package._surfe_wrapper = old_attribute
    assert sys.modules.get("surfepy", _ABSENT) is old_surfepy
    assert sys.modules.get(wrapper_name, _ABSENT) is old_wrapper
    assert package.__dict__.get("_surfe_wrapper", _ABSENT) is old_attribute


def _nodal_field(interp, anchor: np.ndarray) -> np.ndarray:
    return r1_value(interp.local_to_world_points(interp.support.nodes), anchor)


@pytest.mark.parametrize(("anchor_index", "rotation"), CELLS, ids=NATIVE_IDS)
def test_native_adapter_sends_local_rows_and_grid(stub_core, anchor_index, rotation):
    anchor = ANCHORS[anchor_index]
    interp, observations = r1_interpolator("fdi", anchor, rotation, exact=True)
    interp.set_gradient_constraints(_gradient_rows(anchor))
    stub_core.field = _nodal_field(interp, anchor)

    assert _native.solve(interp, 1e-8, {})

    assert len(stub_core.instances) == 1
    problem = stub_core.instances[0]
    support = interp.support
    np.testing.assert_array_equal(problem.args[0], support.origin)
    np.testing.assert_array_equal(problem.args[1], support.nsteps)
    np.testing.assert_array_equal(problem.args[2], support.step_vector)

    world_points = observations[["X", "Y", "Z"]].to_numpy()
    expected_points = interp.bounding_box.project(world_points)
    (sent_points, sent_values, _), = problem.named("add_value")
    sent_points = sent_points.reshape(-1, 3)
    assert sent_points.shape == (27, 3)
    limit = tp(LENGTH, _width(anchor))
    assert np.max(np.abs(sent_points - expected_points)) <= limit
    np.testing.assert_allclose(sent_values, observations["val"].to_numpy(), rtol=0, atol=limit)

    (grad_points, grad_vectors, _), = problem.named("add_gradient")
    grad_points = grad_points.reshape(-1, 3)
    grad_vectors = grad_vectors.reshape(-1, 3)
    assert grad_vectors.shape == (4, 3)
    expected_vectors = interp.world_to_local_vectors(np.tile(UNIT_GRADIENT, (4, 1)))
    assert np.max(np.abs(grad_vectors - expected_vectors)) <= 1e-12
    by_hand = _local_vectors_by_hand(
        interp.bounding_box, anchor, np.tile(UNIT_GRADIENT, (4, 1))
    )
    assert np.max(np.abs(grad_vectors - by_hand)) <= limit
    expected_grad_points = interp.bounding_box.project(anchor + GRADIENT_OFFSETS)
    assert np.max(np.abs(grad_points - expected_grad_points)) <= limit

    # Wrong-answer guard: world data would differ from the local data.
    if rotation in (30.0, 90.0):
        assert np.max(np.abs(sent_points - world_points)) > 1
        assert np.max(np.abs(grad_vectors - UNIT_GRADIENT)) > 1e-3


@pytest.mark.parametrize(("anchor_index", "rotation"), CELLS, ids=NATIVE_IDS)
def test_native_solved_field_evaluates_in_world(stub_core, anchor_index, rotation):
    anchor = ANCHORS[anchor_index]
    interp, _ = r1_interpolator("fdi", anchor, rotation, exact=True)
    interp.set_gradient_constraints(_gradient_rows(anchor))
    field = _nodal_field(interp, anchor)
    stub_core.field = field
    interp.c = np.zeros_like(field)

    converged = _native.solve(interp, 1e-8, {})
    assert converged
    # ``solve_system`` sets this flag after the solve. Without it, evaluation solves again.
    interp.up_to_date = converged

    np.testing.assert_array_equal(interp.c, field)
    query = anchor + QUERY_OFFSETS
    scale = _scale(anchor, query)
    value_error = np.max(np.abs(interp.evaluate_value(query) - r1_value(query, anchor)))
    assert value_error <= 1e-8 * scale
    expected = r1_gradient(query, anchor)
    gradient_error = np.abs(interp.evaluate_gradient(query) - expected)
    assert np.all(gradient_error <= 1e-8 * np.maximum(1.0, np.abs(expected)))


def _rotated_support(interp):
    interp.support.rotation_xy = rot_z(30.0)
    return interp, "rotation_xy"


def _region(interp):
    interp._setup_kwargs = {"region": lambda points: np.ones(len(points), dtype=bool)}
    return interp, "region"


def _p1(_interp):
    interp, _ = r1_interpolator("pli", ANCHORS[1], 30.0, exact=True)
    assert isinstance(interp, P1Interpolator)
    return interp, "interpolator class"


@pytest.mark.parametrize("make", [_rotated_support, _region, _p1], ids=["rotation", "region", "p1"])
def test_native_unsupported_settings_keep_rejection_and_state(stub_core, make):
    interp, _ = r1_interpolator("fdi", ANCHORS[1], 30.0, exact=True)
    interp, word = make(interp)
    stub_core.field = np.zeros(interp.support.n_nodes)
    before_c = np.array(interp.c, copy=True)
    before_data = copy.deepcopy(interp.data)
    before_flag = interp.up_to_date

    with pytest.raises(NotImplementedError) as caught:
        _native.solve(interp, 1e-8, {})

    assert word in str(caught.value)
    np.testing.assert_array_equal(interp.c, before_c)
    assert interp.up_to_date == before_flag
    assert interp.data.keys() == before_data.keys()
    for key, rows in before_data.items():
        np.testing.assert_array_equal(interp.data[key], rows)
    # The rejection must come before any core call.
    assert stub_core.instances == []


def _surfe_interpolator(wrapper, anchor, rotation):
    """Return a Surfe wrapper with R1 world data, and the R1 FDI interpolator for the frame."""
    reference, observations = r1_interpolator("fdi", anchor, rotation, exact=True)
    interp = wrapper.SurfeRBFInterpolator()
    interp.bounding_box = reference.bounding_box
    interp.set_value_constraints(observations[["X", "Y", "Z", "val", "w"]].to_numpy())
    normal_rows = _gradient_rows(anchor)
    interp.set_normal_constraints(normal_rows)
    return interp, observations


@pytest.mark.parametrize(("anchor_index", "rotation"), CELLS, ids=SURFE_IDS)
def test_surfe_adapter_sends_local_points(stub_surfe, anchor_index, rotation):
    wrapper, api_class = stub_surfe
    anchor = ANCHORS[anchor_index]
    interp, observations = _surfe_interpolator(wrapper, anchor, rotation)

    interp.setup_interpolator()

    assert len(api_class.instances) == 1
    api = api_class.instances[0]
    sent = np.array([args[:4] for args in api.named("AddInterfaceConstraint")])
    sent = sent.reshape(-1, 4)
    assert sent.shape == (27, 4)
    world_points = observations[["X", "Y", "Z"]].to_numpy()
    expected = interp.bounding_box.project(world_points)
    limit = tp(LENGTH, _width(anchor))
    assert np.max(np.abs(sent[:, :3] - expected)) <= limit
    np.testing.assert_allclose(sent[:, 3], observations["val"].to_numpy(), rtol=0, atol=limit)

    (planar,) = [args[0] for args in api.named("SetPlanarConstraints")]
    assert planar.shape == (4, 6)
    expected_vectors = _local_vectors_by_hand(
        interp.bounding_box, anchor, np.tile(UNIT_GRADIENT, (4, 1))
    )
    assert np.max(np.abs(planar[:, 3:6] - expected_vectors)) <= limit
    expected_planar_points = interp.bounding_box.project(anchor + GRADIENT_OFFSETS)
    assert np.max(np.abs(planar[:, :3] - expected_planar_points)) <= limit

    # Wrong-answer guard: world data would differ from the local data.
    if rotation in (30.0, 90.0):
        assert np.max(np.abs(sent[:, :3] - world_points)) > 1
        assert np.max(np.abs(planar[:, 3:6] - UNIT_GRADIENT)) > 1e-3


@pytest.mark.parametrize(("anchor_index", "rotation"), CELLS, ids=SURFE_IDS)
def test_surfe_adapter_returns_world_values_and_gradients(stub_surfe, anchor_index, rotation):
    wrapper, api_class = stub_surfe
    anchor = ANCHORS[anchor_index]
    interp, _ = _surfe_interpolator(wrapper, anchor, rotation)
    api_class.value_function = lambda local: r1_value(interp.local_to_world_points(local), anchor)
    api_class.gradient_function = lambda local: interp.world_to_local_vectors(
        r1_gradient(interp.local_to_world_points(local), anchor)
    )
    interp.setup_interpolator()

    query = anchor + QUERY_OFFSETS
    values = interp.evaluate_value(query)
    gradients = interp.evaluate_gradient(query)

    api = api_class.instances[0]
    (sent_local,) = [args[0] for args in api.named("EvaluateInterpolantAtPoints")]
    expected_local = interp.bounding_box.project(query)
    assert np.max(np.abs(sent_local - expected_local)) <= tp(LENGTH, _width(anchor))
    if anchor.any() or rotation:
        assert not np.allclose(sent_local, query, atol=1e-3)

    scale = _scale(anchor, query)
    value_error = float(np.max(np.abs(values - r1_value(query, anchor))))
    expected = r1_gradient(query, anchor)
    gradient_error = np.abs(gradients - expected)
    SURFE_ERRORS.append(max(value_error, float(np.max(gradient_error))))
    assert value_error <= 1e-8 * scale
    assert np.all(gradient_error <= 1e-8 * np.maximum(1.0, np.abs(expected)))
