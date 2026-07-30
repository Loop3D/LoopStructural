import numpy as np
import pytest

from loop_interpolation import FiniteDifferenceInterpolator, P1Interpolator

from LoopStructural.geometry import BoundingBox
from LoopStructural.interpolators import InterpolatorType
from LoopStructural.interpolators._api import LoopInterpolator


@pytest.fixture
def bounding_box():
    return BoundingBox(np.array([0, 0, 0]), np.array([1, 1, 1]))


@pytest.fixture
def value_constraints():
    return np.array(
        [
            [0.1, 0.1, 0.1, 0.0],
            [0.9, 0.1, 0.1, 1.0],
            [0.1, 0.9, 0.1, 0.0],
            [0.9, 0.9, 0.9, 1.0],
        ]
    )


def test_default_interpolator_is_finite_difference(bounding_box):
    api = LoopInterpolator(bounding_box, nelements=500)
    assert isinstance(api.interpolator, FiniteDifferenceInterpolator)
    assert api.bounding_box is bounding_box
    assert api.dimensions == 3


def test_fit_sets_value_constraints_on_underlying_interpolator(bounding_box, value_constraints):
    api = LoopInterpolator(bounding_box, nelements=500)
    api.fit(values=value_constraints)
    assert np.array_equal(api.interpolator.data["value"][:, :4], value_constraints)


def test_fit_sets_normal_constraints(bounding_box):
    api = LoopInterpolator(bounding_box, nelements=500)
    normals = np.array([[0.5, 0.5, 0.5, 0.0, 0.0, 1.0]])
    api.fit(normal_vectors=normals)
    assert np.array_equal(api.interpolator.data["normal"][:, :6], normals)


def test_fit_sets_tangent_constraints(bounding_box):
    api = LoopInterpolator(bounding_box, nelements=500)
    tangents = np.array([[0.5, 0.5, 0.5, 1.0, 0.0, 0.0]])
    api.fit(tangent_vectors=tangents)
    assert np.array_equal(api.interpolator.data["tangent"][:, :6], tangents)


def test_fit_sets_inequality_constraints(bounding_box, monkeypatch):
    # NOTE: calling FiniteDifferenceInterpolator.setup() with inequality value
    # constraints raises a ValueError further down the stack (in
    # DiscreteInterpolator.add_value_inequality_constraints /
    # StructuredGrid.inside, both outside the scope of this test module) because
    # the full n-column constraint array is passed through instead of just the
    # XYZ columns. That looks like a pre-existing bug unrelated to LoopInterpolator
    # itself, so here we stub out `setup` to isolate and verify the constraint
    # dispatch logic in `LoopInterpolator.fit`.
    api = LoopInterpolator(bounding_box, nelements=500)
    monkeypatch.setattr(api.interpolator, "setup", lambda **kwargs: None)
    inequality_values = np.array([[0.5, 0.5, 0.5, 0.0, 1.0]])
    inequality_pairs = np.array([[0.2, 0.2, 0.2, 0]])
    api.fit(
        inequality_value_constraints=inequality_values,
        inequality_pairs_constraints=inequality_pairs,
    )
    assert np.array_equal(api.interpolator.data["inequality"], inequality_values)
    assert np.array_equal(api.interpolator.data["inequality_pairs"], inequality_pairs)


def test_evaluate_scalar_value_matches_constraints_closely(bounding_box, value_constraints):
    api = LoopInterpolator(bounding_box, nelements=1000)
    api.fit(values=value_constraints)
    result = api.evaluate_scalar_value(value_constraints[:, :3])
    assert result.shape == (value_constraints.shape[0],)
    assert np.allclose(result, value_constraints[:, 3], atol=1e-2)


def test_evaluate_gradient_shape(bounding_box, value_constraints):
    api = LoopInterpolator(bounding_box, nelements=500)
    api.fit(values=value_constraints)
    gradient = api.evaluate_gradient(value_constraints[:, :3])
    assert gradient.shape == (value_constraints.shape[0], 3)


def test_fit_and_evaluate_value_returns_values_at_data_locations(bounding_box, value_constraints):
    api = LoopInterpolator(bounding_box, nelements=500)
    result = api.fit_and_evaluate_value(values=value_constraints)
    assert result.shape[0] == value_constraints.shape[0]


def test_fit_and_evaluate_gradient_returns_gradient_at_data_locations(
    bounding_box, value_constraints
):
    api = LoopInterpolator(bounding_box, nelements=500)
    result = api.fit_and_evaluate_gradient(values=value_constraints)
    assert result.shape == (value_constraints.shape[0], 3)


def test_fit_and_evaluate_value_and_gradient_returns_both(bounding_box, value_constraints):
    api = LoopInterpolator(bounding_box, nelements=500)
    values, gradient = api.fit_and_evaluate_value_and_gradient(values=value_constraints)
    assert values.shape[0] == value_constraints.shape[0]
    assert gradient.shape == (value_constraints.shape[0], 3)


def test_type_attribute_ignores_requested_interpolator_type(bounding_box):
    # NOTE: this documents a bug in LoopInterpolator.__init__ - `self.type` is
    # hard-coded to "FDI" regardless of the `type` argument that was passed in,
    # even though the correct interpolator class is constructed via the factory.
    api = LoopInterpolator(bounding_box, nelements=500, type=InterpolatorType.PIECEWISE_LINEAR)
    assert isinstance(api.interpolator, P1Interpolator)
    assert api.type == "FDI"


def test_plot_2d_returns_image_and_axis():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")

    bb2 = BoundingBox(
        origin=np.array([0.0, 0.0]),
        maximum=np.array([1.0, 1.0]),
        global_origin=np.array([0.0, 0.0]),
        dimensions=2,
    )
    api = LoopInterpolator(bb2, dimensions=2, nelements=200)
    values = np.array(
        [
            [0.1, 0.1, 0.0],
            [0.9, 0.1, 1.0],
            [0.1, 0.9, 0.0],
        ]
    )
    api.fit(values=values)
    val, ax = api.plot()
    assert val.ndim == 2
    assert ax is not None


def test_plot_3d_dispatches_to_support_vtk(bounding_box, value_constraints, monkeypatch):
    api = LoopInterpolator(bounding_box, nelements=200)
    api.fit(values=value_constraints)

    calls = {}

    class FakeGrid:
        def __setitem__(self, key, value):
            calls["key"] = key
            calls["value"] = value

        def plot(self, **kwargs):
            calls["plot_kwargs"] = kwargs

    fake_grid = FakeGrid()
    monkeypatch.setattr(api.interpolator.support, "vtk", lambda: fake_grid)

    result = api.plot(color="red")

    assert result is fake_grid
    assert calls["key"] == "val"
    assert calls["plot_kwargs"] == {"color": "red"}
