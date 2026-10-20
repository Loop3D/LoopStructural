"""Tests for P2FoldInterpolator."""
import numpy as np
import pytest
from loop_interpolation import P2Interpolator, P2UnstructuredTetMesh
from loop_interpolation._p2_fold_interpolator import P2FoldInterpolator
from loop_interpolation._quadrature import tet_quadrature_points


class SimpleFold:
    """Fold with constant orientation (z), axis (y) and dgz (x)."""

    def get_deformed_orientation(self, points):
        n = points.shape[0]
        ori = np.tile([0.0, 0.0, 1.0], (n, 1))
        axis = np.tile([0.0, 1.0, 0.0], (n, 1))
        dgz = np.tile([1.0, 0.0, 0.0], (n, 1))
        return ori, axis, dgz


@pytest.fixture
def mesh():
    return P2UnstructuredTetMesh(
        origin=np.zeros(3), step_vector=np.ones(3), nsteps_cells=np.array([2, 2, 2])
    )


def _interp(mesh):
    return P2FoldInterpolator(mesh, fold=SimpleFold())


def _setup(mesh, **fw):
    i = _interp(mesh)
    i.set_value_constraints(np.array([[0.5, 0.5, 0.5, 0.0], [1.5, 1.5, 1.5, 1.0]]))
    i.setup_interpolator(fold_weights=fw)
    return i


@pytest.mark.parametrize("rule,n", [("gauss4", 4), ("centroid", 1)])
def test_counts_and_weights(mesh, rule, n):
    i = _setup(mesh, sample_rule=rule)
    names = {"fold orientation": 10.0, "fold axis": 10.0, "fold normalisation": 1.0}
    _, _, q = tet_quadrature_points(mesh.nodes, mesh.elements, rule)
    for name, w in names.items():
        c = i.constraints[name]
        assert c["matrix"].shape[0] == mesh.elements.shape[0] * n
        np.testing.assert_allclose(c["w"], w * np.sqrt(q))


def test_rows_have_no_volume_scaling(mesh):
    i = _setup(mesh, sample_rule="centroid")
    c = i.constraints["fold normalisation"]
    A = c["matrix"].toarray()
    np.testing.assert_allclose(np.linalg.norm(A, axis=1), 1.0)
    pts, tet, _ = tet_quadrature_points(mesh.nodes, mesh.elements, "centroid")
    dN = mesh.evaluate_shape_derivatives(pts, elements=tet)[0]
    raw = dN[:, 0, :]  # dgz = x
    # B = target / |row|, with no volume factor
    np.testing.assert_allclose(c["b"], -1.0 / np.linalg.norm(raw, axis=1))


def test_invalid_rule_raises(mesh):
    with pytest.raises(ValueError):
        _interp(mesh).setup_interpolator(fold_weights={"sample_rule": "bad"})


def test_copy_type(mesh):
    c = _interp(mesh).copy()
    assert type(c) is P2FoldInterpolator
    assert c.fold is not None
    assert isinstance(c, P2Interpolator)


def test_update_keeps_fold_weights(mesh):
    i = _setup(mesh, sample_rule="centroid", fold_axis_w=3.0)
    assert i._setup_kwargs["fold_weights"]["fold_axis_w"] == 3.0
    i.setup_interpolator(**i._setup_kwargs)
    np.testing.assert_allclose(i.constraints["fold axis"]["w"], 3.0)
