import numpy as np
import pytest

from LoopStructural.interpolators import P2Interpolator


class FakeP2Support:
    """Minimal stand-in for a P2UnstructuredTetMesh, providing just enough
    surface area to exercise P2Interpolator's constraint-building methods
    without needing a geometrically valid quadratic tetrahedral mesh.
    """

    dimension = 3

    def __init__(self, n_elements=2, dof_per_element=10):
        self.n_elements = n_elements
        self.n_nodes = n_elements * dof_per_element
        self.nodes = np.zeros((self.n_nodes, 3))
        self.elements = np.arange(self.n_nodes).reshape(n_elements, dof_per_element)
        self.element_size = np.ones(n_elements)

    def evaluate_shape_derivatives(self, points):
        n = points.shape[0]
        grad = np.ones((n, self.elements.shape[1], 3))
        elements = np.zeros(n, dtype=int)
        return grad, elements

    def evaluate_shape(self, points):
        n = points.shape[0]
        N = np.ones((n, self.elements.shape[1]))
        elements = np.zeros(n, dtype=int)
        mask = np.ones(n, dtype=bool)
        return N, elements, mask

    def evaluate_d2(self, points, c):
        assert not np.any(np.isnan(points)), "nan rows leaked through to support.evaluate_d2"
        return np.full(points.shape[0], 42.0)


@pytest.fixture
def interpolator():
    return P2Interpolator(FakeP2Support())


def test_add_gradient_constraints_uses_correct_support_indexing(interpolator):
    """Regression test: add_gradient_constraints used to index
    `self.support[elements[inside]]` instead of `self.support.elements[...]`.
    No support class implements __getitem__, so this raised a TypeError as
    soon as any gradient constraint was added.
    """
    interpolator.set_gradient_constraints(np.array([[0.1, 0.1, 0.1, 1.0, 0.0, 0.0, 1.0]]))
    interpolator.add_gradient_constraints(w=1.0)
    assert "gradient" in interpolator.constraints
    matrix = interpolator.constraints["gradient"]["matrix"]
    assert matrix.shape == (1, interpolator.dof)


def test_evaluate_d2_masks_nan_rows(interpolator):
    """Regression test: evaluate_d2's nan mask was computed as
    `evaluation_points == np.nan`, which is always False (nan != nan), so
    rows containing nan were never filtered out before being passed to
    support.evaluate_d2.
    """
    interpolator.c = np.zeros(interpolator.support.n_nodes)
    points = np.array(
        [
            [0.1, 0.1, 0.1],
            [np.nan, 0.2, 0.2],
            [0.3, 0.3, 0.3],
        ]
    )
    result = interpolator.evaluate_d2(points)
    assert result.shape == (3,)
    assert result[1] == 0.0
    assert result[0] == 42.0
    assert result[2] == 42.0


def test_add_value_constraints_single_point_not_dropped():
    """Regression test: add_value_constraints required
    `points.shape[0] > 1`, silently discarding a single value constraint
    (inconsistent with the finite-difference interpolator's `> 0` check).
    """
    interp = P2Interpolator(FakeP2Support(n_elements=1))
    interp.set_value_constraints(np.array([[0.1, 0.1, 0.1, 5.0, 1.0]]))
    interp.add_value_constraints(w=1.0)
    assert "value" in interp.constraints
    matrix = interp.constraints["value"]["matrix"]
    assert matrix.shape[0] == 1
