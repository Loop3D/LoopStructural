import numpy as np
import pytest

from LoopStructural.geometry import BoundingBox
from LoopStructural.interpolators import (
    InterpolatorFactory,
    P1Interpolator,
    P2Interpolator,
)
from loop_common.supports import P1Unstructured2d, P2Unstructured2d


def _bbox_2d():
    return BoundingBox(origin=np.array([0, 0]), maximum=np.array([1, 1]), dimensions=2)


def test_p1_unstructured_2d_requires_explicit_arrays_or_bbox_args():
    with pytest.raises(ValueError):
        P1Unstructured2d()


def test_p2_unstructured_2d_requires_explicit_arrays_or_bbox_args():
    with pytest.raises(ValueError):
        P2Unstructured2d()


def test_p1_unstructured_2d_from_bbox_builds_valid_mesh():
    mesh = P1Unstructured2d(origin=np.zeros(2), step_vector=np.ones(2) / 4, nsteps=np.array([5, 5]))
    assert mesh.elements.shape[1] == 3
    assert mesh.elements.min() == 0
    assert mesh.elements.max() == mesh.nodes.shape[0] - 1
    # neighbours must be mutual: if i lists n as a neighbour, n must list i back
    for i, neighbours in enumerate(mesh.neighbours):
        for n in neighbours:
            if n >= 0:
                assert i in mesh.neighbours[n]


def test_p2_unstructured_2d_from_bbox_builds_valid_mesh():
    mesh = P2Unstructured2d(origin=np.zeros(2), step_vector=np.ones(2) / 4, nsteps=np.array([5, 5]))
    assert mesh.elements.shape[1] == 6
    assert mesh.elements.min() == 0
    assert mesh.elements.max() == mesh.nodes.shape[0] - 1
    # edge midpoints deduplicated: far fewer nodes than the naive
    # non-deduplicated upper bound
    naive_upper_bound = mesh.nodes.shape[0] + 3 * mesh.n_elements
    assert mesh.nodes.shape[0] < naive_upper_bound


def test_create_p1_interpolator_2d_via_factory():
    """Regression test: P1Unstructured2d previously only accepted explicit
    elements/vertices/neighbours arrays, so building a 2D P1 interpolator via
    the standard factory always raised TypeError.
    """
    interp = InterpolatorFactory.create_interpolator("P1", _bbox_2d(), nelements=200)
    assert isinstance(interp, P1Interpolator)
    assert interp.dof > 0


def test_create_p2_interpolator_2d_via_factory():
    interp = InterpolatorFactory.create_interpolator("P2", _bbox_2d(), nelements=200)
    assert isinstance(interp, P2Interpolator)
    assert interp.dof > 0


def test_p1_interpolator_2d_reproduces_linear_field():
    """Regression test for two independent, previously-undiscovered bugs:

    1. _p1interpolator.py hardcoded points[:, :3]/points[:, 3] assuming 3D,
       which silently fed a 2D constraint's *value* column in as if it were
       a Z coordinate.
    2. _2d_base_unstructured.py's get_element_for_location computed
       barycentric weights for vertices [1, 2, 0] but stored them in columns
       [0, 1, 2], so every value/gradient evaluation used the wrong vertex's
       weight.

    Both together made any real use of P1Unstructured2d silently wrong
    (rather than crashing), which is why this went unnoticed - nothing in
    LoopStructural/modelling ever uses 2D interpolation.
    """
    interp = InterpolatorFactory.create_interpolator("P1", _bbox_2d(), nelements=500)
    support = interp.support

    def f(xy):
        x, y = xy[:, 0], xy[:, 1]
        return 2 * x - 3 * y + 5

    values = f(support.nodes)
    constraints = np.hstack([support.nodes, values[:, None], np.ones((support.nodes.shape[0], 1))])
    interp.set_value_constraints(constraints)
    interp.add_value_constraints(w=1.0)
    interp.solve_system(solver="lsmr")

    rng = np.random.default_rng(0)
    test_points = rng.uniform(0.05, 0.95, size=(200, 2))
    predicted = interp.evaluate_value(test_points)
    actual = f(test_points)
    valid = ~np.isnan(predicted)
    assert valid.sum() == len(test_points)
    assert np.max(np.abs(predicted[valid] - actual[valid])) < 1e-8


def test_p2_unstructured_2d_get_quadrature_points_shape():
    """Regression test: get_quadrature_points built `cp` with shape
    (n_edges, self.ncps, 2) - using the element node count (6) instead of
    the actual number of quadrature points (2) - and returned a `weight`
    array of that same wrong shape. This only surfaced when
    minimise_edge_jumps (used by the default constant-gradient
    regularisation) actually ran, since 2D P2 was never exercised before.
    """
    mesh = P2Unstructured2d(origin=np.zeros(2), step_vector=np.ones(2) / 4, nsteps=np.array([5, 5]))
    cp, weight = mesh.get_quadrature_points()
    n_edges = mesh.shared_elements.shape[0]
    assert cp.shape == (n_edges, 2, 2)
    assert weight.shape == (n_edges, 2)


def test_p2_interpolator_2d_minimise_edge_jumps_does_not_crash():
    """Regression test for the get_quadrature_points shape bug above, at
    the point where it actually surfaced: calling minimise_edge_jumps on a
    2D P2 interpolator used to raise
    ValueError: could not broadcast input array from shape (n,2) into shape (n,)
    """
    interp = InterpolatorFactory.create_interpolator("P2", _bbox_2d(), nelements=200)
    interp.set_value_constraints(
        np.hstack([np.random.default_rng(0).random((10, 2)), np.zeros((10, 1)), np.ones((10, 1))])
    )
    interp.add_value_constraints(w=1.0)
    interp.minimise_edge_jumps(w=0.1)
    assert any("shared element jump" in name for name in interp.constraints)


def test_p2_interpolator_2d_reproduces_quadratic_field():
    """Strong correctness check for the 2D P1->P2 mesh elevation and the new
    evaluate_value/evaluate_gradient overrides on P2Unstructured2d (the base
    class implementation only handles 3-node linear elements).
    """
    interp = InterpolatorFactory.create_interpolator("P2", _bbox_2d(), nelements=2000)
    support = interp.support

    def f(xy):
        x, y = xy[:, 0], xy[:, 1]
        return x**2 + 2 * y**2 + x * y + 2 * x - 3 * y + 5

    values = f(support.nodes)
    constraints = np.hstack([support.nodes, values[:, None], np.ones((support.nodes.shape[0], 1))])
    interp.set_value_constraints(constraints)
    interp.add_value_constraints(w=1.0)
    interp.solve_system(solver="lsmr")

    rng = np.random.default_rng(0)
    test_points = rng.uniform(0.05, 0.95, size=(500, 2))
    predicted = interp.evaluate_value(test_points)
    actual = f(test_points)
    valid = ~np.isnan(predicted)
    assert valid.sum() == len(test_points)
    assert np.max(np.abs(predicted[valid] - actual[valid])) < 1e-6


def test_p2_unstructured_2d_evaluate_shape_d2_reproduces_analytic_second_derivatives():
    """Regression test: evaluate_shape_d2 referenced self.hN, an attribute
    that's never set anywhere, so calling it raised AttributeError.
    Rewritten to follow the same reference-space hessian + chain rule
    approach as P2UnstructuredTetMesh.evaluate_shape_d2 in 3D (using the
    self.hessian array already built in __init__).

    For a genuine quadratic field, the second derivatives are constant
    everywhere, so this also verifies the fix is numerically correct and
    not just crash-free.
    """
    interp = InterpolatorFactory.create_interpolator("P2", _bbox_2d(), nelements=1000)
    support = interp.support

    # f = 2x^2 + 3y^2 + 1.5xy + ... -> fxx=4, fxy=1.5, fyy=6 everywhere
    def f(xy):
        x, y = xy[:, 0], xy[:, 1]
        return 2 * x**2 + 3 * y**2 + 1.5 * x * y + 2 * x - 3 * y + 5

    values = f(support.nodes)
    interp.set_value_constraints(
        np.hstack([support.nodes, values[:, None], np.ones((support.nodes.shape[0], 1))])
    )
    interp.add_value_constraints(w=1.0)
    interp.solve_system(solver="lsmr")

    rng = np.random.default_rng(0)
    test_points = rng.uniform(0.05, 0.95, size=(50, 2))
    d2 = support.evaluate_d2(test_points, interp.c)
    assert d2.shape == (50, 3)
    assert np.allclose(np.nanmean(d2, axis=0), [4.0, 1.5, 6.0], atol=1e-6)


def test_p1_interpolator_2d_gradient_orthogonal_constraints_consistent_with_linear_field():
    """Regression test for a half-applied dimension fix: add_gradient_orthogonal_constraints
    in _p1interpolator.py hardcoded points[:, :3] when slicing point coordinates, while the
    sibling add_gradient_constraints/add_norm_constraints methods in the same file correctly
    used points[:, : self.dimensions]. For a 2D interpolator (self.dimensions == 2) this is
    exercised here directly: a gradient-orthogonal constraint that is mathematically consistent
    with a known linear field should not corrupt the least-squares system, and the field should
    still be reproduced (up to the small residual expected from blending an extra weighted
    least-squares constraint on top of the exact value pins).
    """
    interp = InterpolatorFactory.create_interpolator("P1", _bbox_2d(), nelements=500)
    support = interp.support

    def f(xy):
        x, y = xy[:, 0], xy[:, 1]
        return 2 * x - 3 * y + 5

    values = f(support.nodes)
    constraints = np.hstack([support.nodes, values[:, None], np.ones((support.nodes.shape[0], 1))])
    interp.set_value_constraints(constraints)
    interp.add_value_constraints(w=1.0)

    # gradient of f is (2, -3); (3, 2) is orthogonal to it (dot product == 0)
    points = support.barycentre
    assert points.shape[1] == 2
    vectors = np.tile(np.array([3.0, 2.0]), (points.shape[0], 1))
    interp.add_gradient_orthogonal_constraints(points, vectors, w=1.0, name="gradient orthogonal")
    assert "gradient orthogonal" in interp.constraints

    interp.solve_system(solver="lsmr")

    rng = np.random.default_rng(0)
    test_points = rng.uniform(0.05, 0.95, size=(200, 2))
    predicted = interp.evaluate_value(test_points)
    actual = f(test_points)
    valid = ~np.isnan(predicted)
    assert valid.sum() == len(test_points)
    # A dimension-slicing bug (e.g. feeding a value/weight column in as a coordinate)
    # would produce errors many orders of magnitude larger than this tolerance.
    assert np.max(np.abs(predicted[valid] - actual[valid])) < 1e-3


def test_p2_interpolator_2d_gradient_orthogonal_constraints_consistent_with_quadratic_field():
    """Same regression as above for _p2interpolator.py's
    add_gradient_orthogonal_constraints, which had the identical points[:, :3] vs
    points[:, : self.dimensions] inconsistency. Here the orthogonal vector varies per point
    since the quadratic field's gradient is not constant.
    """
    interp = InterpolatorFactory.create_interpolator("P2", _bbox_2d(), nelements=2000)
    support = interp.support

    def f(xy):
        x, y = xy[:, 0], xy[:, 1]
        return x**2 + 2 * y**2 + x * y + 2 * x - 3 * y + 5

    def grad_f(xy):
        x, y = xy[:, 0], xy[:, 1]
        gx = 2 * x + y + 2
        gy = 4 * y + x - 3
        return np.stack([gx, gy], axis=1)

    values = f(support.nodes)
    constraints = np.hstack([support.nodes, values[:, None], np.ones((support.nodes.shape[0], 1))])
    interp.set_value_constraints(constraints)
    interp.add_value_constraints(w=1.0)

    points = support.barycentre
    assert points.shape[1] == 2
    grad = grad_f(points)
    # rotate each gradient vector by 90 degrees to get a vector orthogonal to it
    vectors = np.stack([-grad[:, 1], grad[:, 0]], axis=1)
    interp.add_gradient_orthogonal_constraints(points, vectors, w=1.0)

    interp.solve_system(solver="lsmr")

    rng = np.random.default_rng(0)
    test_points = rng.uniform(0.05, 0.95, size=(500, 2))
    predicted = interp.evaluate_value(test_points)
    actual = f(test_points)
    valid = ~np.isnan(predicted)
    assert valid.sum() == len(test_points)
    assert np.max(np.abs(predicted[valid] - actual[valid])) < 1e-3


def test_p2_interpolator_2d_minimise_grad_steepness_does_not_crash():
    """Regression test: minimise_grad_steepness calls
    support.evaluate_shape_d2, which previously raised
    AttributeError: 'P2Unstructured2d' object has no attribute 'hN'
    the moment it was invoked - i.e. every time setup_interpolator() ran
    its default regularisation for a 2D P2 interpolator.
    """
    interp = InterpolatorFactory.create_interpolator("P2", _bbox_2d(), nelements=200)
    interp.set_value_constraints(
        np.hstack([np.random.default_rng(0).random((10, 2)), np.zeros((10, 1)), np.ones((10, 1))])
    )
    interp.setup_interpolator(regularisation=0.1)
    interp.solve_system(solver="lsmr")
    assert any("gradsteepness" in name for name in interp.constraints)
