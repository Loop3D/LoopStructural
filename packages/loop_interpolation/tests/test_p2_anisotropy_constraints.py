"""Tests that the element-based anisotropy constraints (and their underlying
generic constraint methods) also work for P2Interpolator, not just P1.

Coverage
--------
1. Regression tests for P2Interpolator.add_gradient_orthogonal_constraints:
   the ``b``/target parameter used to be silently ignored (always forced to
   zero), and ``name`` wasn't a parameter at all.
2. Regression test for P2Interpolator directional regularisation: previously
   unimplemented (NotImplementedError from the DiscreteInterpolator base),
   now backed by minimise_edge_jumps via a face-centroid sample-point grid.
3. add_element_anisotropy_constraints works end-to-end against a
   P2Interpolator for all four constraint families.
"""

import numpy as np
import pytest
from loop_common.supports import P2UnstructuredTetMesh
from loop_interpolation import P2Interpolator
from loop_interpolation._anisotropy import add_element_anisotropy_constraints
from loop_interpolation._regularisation import DirectionalRegularisation


def _make_p2_mesh(nsteps_cells=(3, 3, 3)):
    return P2UnstructuredTetMesh(
        origin=np.zeros(3),
        step_vector=np.ones(3) / 4,
        nsteps_cells=np.array(nsteps_cells),
    )


def _constant_anisotropy(normal_dir=(0, 0, 1), primary_dir=(1, 0, 0), secondary_dir=(0, 1, 0)):
    class _AnisotropyStub:
        def get_directions(self, points):
            n = points.shape[0]
            primary = np.tile(np.array(primary_dir, dtype=float), (n, 1))
            secondary = np.tile(np.array(secondary_dir, dtype=float), (n, 1))
            normal = np.tile(np.array(normal_dir, dtype=float), (n, 1))
            return primary, secondary, normal

    return _AnisotropyStub()


class TestP2GradientOrthogonalConstraintsRegression:
    """P2Interpolator.add_gradient_orthogonal_constraints used to silently
    force the target ('B') to zero and had no 'name' parameter at all."""

    def test_target_b_is_actually_applied(self):
        """b used to be hardcoded to zero regardless of the argument passed in.

        add_constraints_to_least_squares row-normalises A and B together
        (divides both by the row norm of A), so the stored 'b' is not the
        literal target -- but it must be strictly positive whenever a
        positive target is requested (it was always exactly zero before
        this fix).
        """
        mesh = _make_p2_mesh()
        interp = P2Interpolator(mesh)
        points = mesh.barycentre
        vectors = np.tile(np.array([0.0, 0.0, 1.0]), (points.shape[0], 1))

        interp.add_gradient_orthogonal_constraints(points, vectors, w=1.0, b=2.0, name="target_b")
        assert np.all(interp.constraints["target_b"]["b"] > 0.0)

        interp.reset()
        interp.add_gradient_orthogonal_constraints(points, vectors, w=1.0, b=-2.0, name="target_b")
        assert np.all(interp.constraints["target_b"]["b"] < 0.0)

    def test_name_parameter_is_used(self):
        mesh = _make_p2_mesh()
        interp = P2Interpolator(mesh)
        points = mesh.barycentre
        vectors = np.tile(np.array([0.0, 0.0, 1.0]), (points.shape[0], 1))

        interp.add_gradient_orthogonal_constraints(points, vectors, w=1.0, name="custom name")
        assert "custom name" in interp.constraints

    def test_default_b_is_zero(self):
        mesh = _make_p2_mesh()
        interp = P2Interpolator(mesh)
        points = mesh.barycentre
        vectors = np.tile(np.array([0.0, 0.0, 1.0]), (points.shape[0], 1))

        interp.add_gradient_orthogonal_constraints(points, vectors, w=1.0, name="default_b")
        assert np.all(interp.constraints["default_b"]["b"] == 0.0)


class TestP2DirectionalRegularisationRegression:
    """P2Interpolator previously had no get_regularisation_sample_points/
    _add_directional_regularisation, so add_directional_regularisation with
    a non-empty config raised NotImplementedError."""

    def test_directional_regularisation_no_longer_raises(self):
        mesh = _make_p2_mesh()
        interp = P2Interpolator(mesh)
        vector = np.array([0.0, 0.0, 1.0])

        interp.add_directional_regularisation(
            (
                DirectionalRegularisation(
                    weight=0.1,
                    direction=lambda points: np.tile(vector, (points.shape[0], 1)),
                    name="p2 directional reg",
                ),
            )
        )
        matching = [k for k in interp.constraints if "p2 directional reg" in k]
        assert matching
        assert all(interp.constraints[k]["matrix"].shape[0] > 0 for k in matching)

    def test_isotropic_cgw_regularisation_still_works(self):
        """cgw-driven minimise_edge_jumps (unaffected by the name/b changes)."""
        mesh = _make_p2_mesh()
        interp = P2Interpolator(mesh)
        interp.setup_interpolator(cgw=0.1, gpw=0.0, npw=0.0, tpw=0.0, cpw=0.0)
        assert any("shared element jump" in k for k in interp.constraints)


class TestAddElementAnisotropyConstraintsOnP2:
    @pytest.fixture
    def setup_interp(self):
        mesh = _make_p2_mesh()
        anisotropy = _constant_anisotropy()
        interp = P2Interpolator(mesh)
        interp.setup_interpolator(cgw=0.0, gpw=0.0, npw=0.0, tpw=0.0, cpw=0.0)
        return interp, anisotropy

    def test_primary_constraints_added(self, setup_interp):
        interp, anisotropy = setup_interp
        add_element_anisotropy_constraints(
            interp,
            anisotropy,
            primary_w=5.0,
            secondary_w=None,
            regularisation_weights=None,
            normal_w=None,
        )
        assert any("anisotropy primary" in k for k in interp.constraints)

    def test_secondary_constraints_added(self, setup_interp):
        interp, anisotropy = setup_interp
        add_element_anisotropy_constraints(
            interp,
            anisotropy,
            primary_w=None,
            secondary_w=5.0,
            regularisation_weights=None,
            normal_w=None,
        )
        assert any("anisotropy secondary" in k for k in interp.constraints)

    def test_normal_constraints_use_target(self, setup_interp):
        interp, anisotropy = setup_interp
        add_element_anisotropy_constraints(
            interp,
            anisotropy,
            primary_w=None,
            secondary_w=None,
            regularisation_weights=None,
            normal_w=1.0,
            norm_target=2.0,
            norm_alignment="none",
        )
        keys = [k for k in interp.constraints if "anisotropy normal" in k]
        assert keys
        # b is row-normalised alongside A (see test_target_b_is_actually_applied),
        # so it won't equal norm_target exactly -- but it must be strictly
        # positive for a positive target (it was always exactly zero before
        # the add_gradient_orthogonal_constraints fix).
        assert np.all(interp.constraints[keys[0]]["b"] > 0.0)

    def test_regularisation_constraints_added(self, setup_interp):
        interp, anisotropy = setup_interp
        add_element_anisotropy_constraints(
            interp,
            anisotropy,
            primary_w=None,
            secondary_w=None,
            regularisation_weights=[0.1, 0.01, 0.01],
            normal_w=None,
        )
        assert any("anisotropy regularisation" in k for k in interp.constraints)

    def test_all_families_together(self, setup_interp):
        interp, anisotropy = setup_interp
        add_element_anisotropy_constraints(
            interp,
            anisotropy,
            primary_w=5.0,
            secondary_w=5.0,
            regularisation_weights=[0.1, 0.01, 0.01],
            normal_w=1.0,
        )
        assert any("anisotropy primary" in k for k in interp.constraints)
        assert any("anisotropy secondary" in k for k in interp.constraints)
        assert any("anisotropy normal" in k for k in interp.constraints)
        assert any("anisotropy regularisation" in k for k in interp.constraints)
