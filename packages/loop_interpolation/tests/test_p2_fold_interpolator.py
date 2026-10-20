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


def test_update_replays_setup_kwargs(mesh):
    i = _interp(mesh)
    i.set_value_constraints(np.array([[0.5, 0.5, 0.5, 0.0], [1.5, 1.5, 1.5, 1.0]]))
    fw = {"sample_rule": "centroid", "fold_orientation": 7.0, "fold_axis_w": 3.0}
    i.setup_interpolator(fold_weights=fw, regularisation=0.3)
    assert i.solve_system("cg")
    i.up_to_date = False
    # drop the constraints, then update() must rebuild them from the stored kwargs
    i.reset()
    assert i.update()
    for name, w in (("fold orientation", 7.0), ("fold axis", 3.0)):
        c = i.constraints[name]
        assert c["matrix"].shape[0] == mesh.elements.shape[0]  # centroid: 1 per tet
        np.testing.assert_allclose(c["w"], w * np.sqrt(0.0 + _q_centroid(mesh)))
    assert i._setup_kwargs["regularisation"] == 0.3
    assert i._setup_kwargs["fold_weights"] == fw
    assert i.interpolation_weights["cgw"] == i.resolve_regularisation_config(
        regularisation=0.3, directional_regularisation=None
    ).isotropic


def _q_centroid(mesh):
    return tet_quadrature_points(mesh.nodes, mesh.elements, "centroid")[2]


def test_copy_is_independent_fold_interpolator(mesh):
    i = _setup(mesh, sample_rule="centroid")
    c = i.copy()
    assert c is not i and c.fold is i.fold and c.support is i.support
    assert c.constraints is not i.constraints
    # copy starts fresh: it can be set up with its own kwargs
    c.set_value_constraints(np.array([[0.5, 0.5, 0.5, 0.0]]))
    c.setup_interpolator(fold_weights={"sample_rule": "gauss4"})
    assert c.constraints["fold axis"]["matrix"].shape[0] == 4 * mesh.elements.shape[0]
    assert i.constraints["fold axis"]["matrix"].shape[0] == mesh.elements.shape[0]


def test_weights_stay_aligned_with_rows_unequal_weights(mesh):
    # a weight that differs per sample point, so a row/weight mix-up shows up
    pts, tet, q = tet_quadrature_points(mesh.nodes, mesh.elements, "gauss4")

    def mask_fn(p):
        return np.asarray(p)[:, 0] > 1.2  # zero weight on part of the mesh

    i = _interp(mesh)
    i.set_value_constraints(np.array([[0.5, 0.5, 0.5, 0.0]]))
    # a non constant dgz so rows differ and can be matched back to their sample
    i.fold.get_deformed_orientation = lambda p: (
        np.c_[p[:, 1] + 1, p[:, 2] + 2, p[:, 0] + 3],
        np.tile([0.0, 1.0, 0.0], (len(p), 1)),
        np.c_[p[:, 0] + 1, np.ones(len(p)), p[:, 2] + 1],
    )
    i.setup_interpolator(
        fold_weights={"sample_rule": "gauss4", "fold_axis_w": None, "fold_normalisation": None}
    )
    base = i.constraints["fold orientation"]
    # rebuild with the mask, from scratch, and compare to the expected per-sample weights
    i.reset()
    i.add_fold_constraints(fold_axis_w=None, fold_normalisation=None, mask_fn=mask_fn)
    c = i.constraints["fold orientation"]
    dN = mesh.evaluate_shape_derivatives(pts, elements=tet)[0]
    vec = i.fold.get_deformed_orientation(pts)[0]
    raw = np.einsum("ij,ijk->ik", vec, dN)
    raw /= np.linalg.norm(raw, axis=1)[:, None]
    idc = mesh.elements[tet]
    expected_w = 10.0 * np.sqrt(q) * (~mask_fn(pts))
    # unequal weights: zero on masked samples, in an irregular pattern
    assert len(np.unique(np.round(expected_w, 12))) == 2
    A = c["matrix"].toarray()
    assert A.shape[0] == len(pts)
    exp_A = np.zeros_like(A)
    for r in range(len(pts)):
        exp_A[r, idc[r]] = raw[r]
    np.testing.assert_allclose(A, exp_A, atol=1e-12)
    np.testing.assert_allclose(c["w"], expected_w)
    assert base["matrix"].shape == c["matrix"].shape


def test_weights_follow_rows_when_nan_rows_dropped(mesh):
    from loop_interpolation._discrete_interpolator import DiscreteInterpolator  # noqa: F401

    i = _interp(mesh)
    A = np.array([[1.0, 0.0], [np.nan, 1.0], [0.0, 2.0], [3.0, 4.0]])
    idc = np.array([[0, 1], [1, 2], [2, 3], [3, 4]])
    B = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.array([10.0, 20.0, 30.0, 40.0])
    i.add_constraints_to_least_squares(A, B, idc, w=w, name="t")
    c = i.constraints["t"]
    np.testing.assert_allclose(c["w"], [10.0, 30.0, 40.0])
    np.testing.assert_allclose(c["b"], [1.0, 1.5, 0.8])
    M = c["matrix"].toarray()
    assert M[1, 3] == 1.0 and M[2, 4] == 0.8


def test_per_row_weight_with_3d_input_is_repeated_not_tiled(mesh):
    # (n, k, m) input is flattened row-major: row (i, j) -> i*k + j
    i = _interp(mesh)
    A = np.ones((2, 3, 2))
    idc = np.tile(np.array([[0, 1]]), (2, 3, 1))
    i.add_constraints_to_least_squares(A, np.ones((2, 3)), idc, w=np.array([1.0, 2.0]), name="t3")
    np.testing.assert_allclose(i.constraints["t3"]["w"], [1, 1, 1, 2, 2, 2])


@pytest.mark.parametrize("family", ["fold orientation", "fold axis", "fold normalisation"])
def test_quadrature_preserves_total_family_energy(mesh, family):
    # Rows are normalised before weighting. The squared Frobenius norm is
    # E[||WAx||^2] for an isotropic unit-variance nodal field x, so this tests
    # the family energy without assuming identical residuals at different points.
    energies = []
    for rule in ("centroid", "gauss4"):
        c = _setup(mesh, sample_rule=rule).constraints[family]
        WA = c["matrix"].multiply(c["w"][:, None])
        energies.append(float(WA.multiply(WA).sum()))
    np.testing.assert_allclose(energies[0], energies[1], rtol=1e-12)


def test_nan_fold_rows_build_and_log_count(mesh, caplog):
    class PartialFold(SimpleFold):
        def get_deformed_orientation(self, points):
            ori, axis, dgz = super().get_deformed_orientation(points)
            ori[points[:, 0] < 1.0] = np.nan
            return ori, axis, dgz

    i = P2FoldInterpolator(mesh, fold=PartialFold())
    pts, _, _ = tet_quadrature_points(mesh.nodes, mesh.elements, "gauss4")
    dropped = int(np.sum(pts[:, 0] < 1.0))
    i.setup_interpolator()
    assert i.constraints["fold orientation"]["matrix"].shape[0] == len(pts) - dropped
    assert i.constraints["fold axis"]["matrix"].shape[0] == len(pts)
    assert f"Dropped {dropped} non-finite rows from fold orientation" in caplog.text
    assert i.solve_system("cg")
    assert np.isfinite(i.evaluate_value(np.array([[0.5, 0.5, 0.5]]))).all()


@pytest.mark.parametrize("alignment,expected_sign", [("none", -1), ("warn", -1), ("correct", 1)])
def test_fold_norm_sign_matches_dfi(mesh, alignment, expected_sign):
    from loop_interpolation import DiscreteFoldInterpolator, TetMesh
    from loop_interpolation._fold_norm_alignment import resolve_fold_norm_target
    from unittest.mock import Mock

    normals = np.array([[0.5, 0.5, 0.5, -1., 0., 0.]])
    target = resolve_fold_norm_target(
        fold=SimpleFold(), normal_constraints=normals, fold_norm=-1.,
        dgz_alignment=alignment, logger=Mock(),
    )
    assert np.sign(target) == expected_sign
    dfi = DiscreteFoldInterpolator(
        TetMesh(origin=np.zeros(3), step_vector=np.ones(3), nsteps=np.full(3, 3)),
        fold=SimpleFold(),
    )
    for i in (_interp(mesh), dfi):
        i.set_normal_constraints(normals)
        i.setup_interpolator(fold_weights={"dgz_alignment": alignment})
        assert np.all(np.sign(i.constraints["fold normalisation"]["b"]) == np.sign(target))


@pytest.mark.parametrize("angle", [0.0, 30.0])
def test_constant_fold_without_orientation_data(mesh, angle):
    theta = np.deg2rad(angle)
    orientation = np.array([np.sin(theta), 0., np.cos(theta)])
    dgz = np.array([np.cos(theta), 0., -np.sin(theta)])
    expected = -dgz  # normalisation target -1

    class ConstantAngleFold:
        def get_deformed_orientation(self, points):
            return tuple(np.tile(v, (len(points), 1)) for v in
                         (orientation, np.array([0., 1., 0.]), dgz))

    i = P2FoldInterpolator(mesh, fold=ConstantAngleFold())
    contacts = np.array([[0.5, 0.5, 0.5], [1.5, 1.5, 1.5]])
    i.set_value_constraints(np.c_[contacts, contacts @ expected])
    i.setup_interpolator()
    assert i.solve_system("cg")
    points = np.random.default_rng(83).uniform(0.1, 1.9, (100, 3))
    gradient = i.evaluate_gradient(points)
    assert np.isfinite(gradient).all()
    angles = np.degrees(np.arccos(np.clip(
        gradient @ expected / np.linalg.norm(gradient, axis=1), -1., 1.
    )))
    assert np.median(angles) < 1.0
