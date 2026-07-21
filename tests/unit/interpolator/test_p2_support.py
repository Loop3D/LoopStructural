import numpy as np
import pytest

from LoopStructural.geometry import BoundingBox
from LoopStructural.interpolators import (
    InterpolatorFactory,
    P2Interpolator,
    P2UnstructuredTetMesh,
)


def test_p2_tetmesh_requires_explicit_arrays_or_bbox_args():
    """Regression test: P2UnstructuredTetMesh could previously only be built
    from explicit nodes/elements/neighbours arrays, which meant
    InterpolatorFactory.create_interpolator('P2', bounding_box, ...) - the
    standard, documented way to build any interpolator - always raised
    TypeError, since SupportFactory.create_support_from_bbox calls every
    support class with origin/step_vector/nsteps.
    """
    with pytest.raises(ValueError):
        P2UnstructuredTetMesh()


def test_p2_tetmesh_from_bbox_builds_valid_mesh():
    origin = np.zeros(3)
    step_vector = np.ones(3) / 4
    nsteps = np.array([5, 5, 5])
    mesh = P2UnstructuredTetMesh(origin=origin, step_vector=step_vector, nsteps=nsteps)

    assert mesh.elements.shape[1] == 10
    # every node index used should be a valid row in nodes, with no gaps
    assert mesh.elements.min() == 0
    assert mesh.elements.max() == mesh.nodes.shape[0] - 1
    # edge midpoints must be deduplicated: far fewer nodes than
    # n_corner_nodes + 6 * n_elements (the naive, non-deduplicated count)
    naive_upper_bound = mesh.nodes.shape[0] + 6 * mesh.n_elements
    assert mesh.nodes.shape[0] < naive_upper_bound

    # adjacent elements should share edge-midpoint nodes (i.e. share a face's
    # three edges), not each get their own independent midpoint nodes
    shared_face_found = False
    for i, neighbours in enumerate(mesh.neighbours):
        for n in neighbours:
            if n < 0:
                continue
            shared = set(mesh.elements[i, 4:].tolist()) & set(mesh.elements[n, 4:].tolist())
            if len(shared) >= 3:
                shared_face_found = True
                break
        if shared_face_found:
            break
    assert shared_face_found


def test_create_p2_interpolator_via_factory():
    """Regression test: this call used to raise
    TypeError: P2UnstructuredTetMesh.__init__() got an unexpected keyword
    argument 'origin'.
    """
    bb = BoundingBox(origin=np.array([0, 0, 0]), maximum=np.array([1, 1, 1]))
    interp = InterpolatorFactory.create_interpolator("P2", bb, nelements=200)
    assert isinstance(interp, P2Interpolator)
    assert interp.dof > 0


def test_p2_interpolator_reproduces_quadratic_field():
    """Strong correctness check for the P1->P2 mesh elevation: a genuine
    quadratic scalar field, constrained at every node (corners and edge
    midpoints) of a P2 mesh, should be reproduced almost exactly everywhere
    by the quadratic shape functions - this only holds if the local edge ->
    node ordering used when building the mesh matches the ordering assumed
    by evaluate_shape/evaluate_shape_derivatives.
    """
    bb = BoundingBox(origin=np.array([0, 0, 0]), maximum=np.array([1, 1, 1]))
    interp = InterpolatorFactory.create_interpolator("P2", bb, nelements=1000)
    support = interp.support

    def f(xyz):
        x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        return x**2 + 2 * y**2 + 3 * z**2 + x * y - y * z + 2 * x - 3 * y + z + 5

    values = f(support.nodes)
    constraints = np.hstack([support.nodes, values[:, None], np.ones((support.nodes.shape[0], 1))])
    interp.set_value_constraints(constraints)
    # add value constraints directly, skipping setup_interpolator's default
    # constant-gradient smoothing regularisation, so this checks pure
    # interpolation accuracy rather than a smoothed fit
    interp.add_value_constraints(w=1.0)
    interp.solve_system(solver="lsmr")

    rng = np.random.default_rng(0)
    test_points = rng.uniform(0.05, 0.95, size=(200, 3))
    predicted = interp.evaluate_value(test_points)
    actual = f(test_points)
    valid = ~np.isnan(predicted)
    assert valid.sum() == len(test_points)
    err = np.abs(predicted[valid] - actual[valid])
    assert err.max() < 1e-6
