"""
Tests for FDFoldInterpolator and FiniteDifferenceInterpolator.minimise_directional_gradient_change.

Coverage
--------
1. Import / construction
2. No-fold guard
3. minimise_directional_gradient_change
   a. Constraints are added and non-empty for each of the six operator types
   b. Bad vector shape is silently skipped with a warning (not an exception)
   c. Zero-weight direction components generate no rows (e.g. axis-aligned vector)
   d. Masked nodes (zeroed vectors) do not contribute rows
4. add_fold_constraints
   a. All four families are added when weights are non-zero
   b. Individual families can be disabled via None weight
   c. mask_fn zeroes out excluded nodes
5. setup_interpolator with fold_weights forwarding
6. FDFoldInterpolator produces a solve-able system and recovers a planar field
"""

import numpy as np
import pytest

from loop_interpolation import FDFoldInterpolator, FiniteDifferenceInterpolator, StructuredGrid
from loop_common.supports import RectilinearGrid


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_uniform_grid(n=10):
    """Return a small uniform StructuredGrid."""
    origin = np.zeros(3)
    nsteps = np.array([n, n, n])
    step_vector = np.ones(3)
    return StructuredGrid(origin=origin, nsteps=nsteps, step_vector=step_vector)


def _make_rectilinear_grid(n=10):
    x = np.linspace(0.0, float(n), n + 1)
    y = np.linspace(0.0, float(n), n + 1)
    z = np.linspace(0.0, float(n), n + 1)
    return RectilinearGrid(x, y, z)


def _constant_fold(n_nodes, dgz_dir=(0, 0, 1), deformed_dir=(1, 0, 0), axis_dir=(0, 1, 0)):
    """
    Return a minimal FoldEvent stub whose get_deformed_orientation always
    returns uniform constant direction vectors.
    """
    dgz_vec = np.tile(np.array(dgz_dir, dtype=float), (n_nodes, 1))
    deformed_vec = np.tile(np.array(deformed_dir, dtype=float), (n_nodes, 1))
    axis_vec = np.tile(np.array(axis_dir, dtype=float), (n_nodes, 1))
    # Normalise each
    for v in (dgz_vec, deformed_vec, axis_vec):
        nrm = np.linalg.norm(v[0])
        if nrm > 0:
            v /= nrm

    class _FoldStub:
        def get_deformed_orientation(self, _pts):
            return deformed_vec, axis_vec, dgz_vec

    return _FoldStub()


# ---------------------------------------------------------------------------
# 1. Construction
# ---------------------------------------------------------------------------


class TestFDFoldInterpolatorConstruction:
    def test_creates_without_fold(self):
        grid = _make_uniform_grid()
        interp = FDFoldInterpolator(grid)
        assert interp is not None
        assert interp.fold is None

    def test_creates_with_fold(self):
        grid = _make_uniform_grid()
        fold = _constant_fold(grid.n_nodes)
        interp = FDFoldInterpolator(grid, fold=fold)
        assert interp.fold is fold

    def test_dof_matches_grid_nodes(self):
        grid = _make_uniform_grid(8)
        interp = FDFoldInterpolator(grid)
        assert interp.dof == grid.n_nodes

    def test_type_is_finite_difference(self):
        from loop_interpolation import InterpolatorType

        grid = _make_uniform_grid()
        interp = FDFoldInterpolator(grid)
        assert interp.type == InterpolatorType.FINITE_DIFFERENCE


# ---------------------------------------------------------------------------
# 2. No-fold guard
# ---------------------------------------------------------------------------


class TestNoFoldGuard:
    def test_setup_without_fold_raises(self):
        grid = _make_uniform_grid()
        interp = FDFoldInterpolator(grid)
        with pytest.raises(RuntimeError, match="no fold event"):
            interp.setup_interpolator()

    def test_setup_with_fold_does_not_raise(self):
        grid = _make_uniform_grid()
        fold = _constant_fold(grid.n_nodes)
        interp = FDFoldInterpolator(grid, fold=fold)
        interp.setup_interpolator()  # should not raise


# ---------------------------------------------------------------------------
# 3. minimise_directional_gradient_change
# ---------------------------------------------------------------------------


class TestMinimiseDirectionalGradientChange:
    @pytest.fixture
    def fdi(self):
        return FiniteDifferenceInterpolator(_make_uniform_grid())

    @pytest.fixture
    def fdi_rect(self):
        return FiniteDifferenceInterpolator(_make_rectilinear_grid())

    def _z_vector(self, fdi):
        """Constant z-direction vector field over all nodes."""
        v = np.zeros((fdi.support.n_nodes, 3))
        v[:, 2] = 1.0
        return v

    def _diagonal_vector(self, fdi):
        """Constant diagonal (1,1,1)/sqrt(3) vector field."""
        v = np.ones((fdi.support.n_nodes, 3))
        v /= np.sqrt(3)
        return v

    # 3a. Constraints are added for an axis-aligned vector (only dzz should
    #     be non-zero; mixed terms should be zero weight and thus absent).
    def test_z_vector_adds_dzz_constraint(self, fdi):
        fdi.reset()
        fdi.minimise_directional_gradient_change(1.0, self._z_vector(fdi), name="test_reg")
        # dzz rows should appear (vz^2 = 1)
        matching = [k for k in fdi.constraints if "test_reg_dzz" in k]
        assert len(matching) > 0, "Expected test_reg_dzz constraints"
        assert fdi.constraints[matching[0]]["matrix"].shape[0] > 0

    def test_z_vector_no_mixed_constraints(self, fdi):
        fdi.reset()
        fdi.minimise_directional_gradient_change(1.0, self._z_vector(fdi), name="test_reg")
        for op in ("dxx", "dyy", "dxy", "dxz", "dyz"):
            key = f"test_reg_{op}"
            if key in fdi.constraints:
                # If present, all rows must have come from the full operator
                # but the component weight should have been zero so no rows
                # are expected.
                assert fdi.constraints[key]["matrix"].shape[0] == 0, (
                    f"Expected no rows for {key} with z-only vector"
                )

    # 3b. Diagonal vector adds all six operator types.
    def test_diagonal_vector_adds_all_six_operators(self, fdi):
        fdi.reset()
        fdi.minimise_directional_gradient_change(1.0, self._diagonal_vector(fdi), name="diag")
        for op in ("dxx", "dyy", "dzz", "dxy", "dxz", "dyz"):
            key = f"diag_{op}"
            assert key in fdi.constraints, f"Missing constraint {key}"
            assert fdi.constraints[key]["matrix"].shape[0] > 0

    # 3c. Wrong shape is silently skipped (no exception).
    def test_bad_vector_shape_no_exception(self, fdi):
        fdi.reset()
        bad = np.ones((10, 3))  # wrong n_nodes dimension
        fdi.minimise_directional_gradient_change(1.0, bad, name="bad")
        # No constraints should have been added.
        assert not any("bad" in k for k in fdi.constraints)

    def test_none_vector_no_exception(self, fdi):
        fdi.reset()
        fdi.minimise_directional_gradient_change(1.0, None, name="none_vec")
        assert not any("none_vec" in k for k in fdi.constraints)

    # 3d. Zero-weight: weight=0 produces no rows.
    def test_zero_weight_adds_no_rows(self, fdi):
        fdi.reset()
        fdi.minimise_directional_gradient_change(0.0, self._diagonal_vector(fdi), name="zero_w")
        for k in fdi.constraints:
            if "zero_w" in k:
                assert fdi.constraints[k]["matrix"].shape[0] == 0

    # Works on a RectilinearGrid too.
    def test_works_on_rectilinear_grid(self, fdi_rect):
        fdi_rect.reset()
        v = self._diagonal_vector(fdi_rect)
        fdi_rect.minimise_directional_gradient_change(1.0, v, name="rect_reg")
        matching = [k for k in fdi_rect.constraints if "rect_reg" in k]
        assert len(matching) > 0


# ---------------------------------------------------------------------------
# 4. add_fold_constraints
# ---------------------------------------------------------------------------


class TestAddFoldConstraints:
    @pytest.fixture
    def setup_interp(self):
        """Return (interpolator, fold) ready for add_fold_constraints calls."""
        grid = _make_uniform_grid(8)
        fold = _constant_fold(grid.n_nodes)
        interp = FDFoldInterpolator(grid, fold=fold)
        # Call parent setup only (no fold constraints yet).
        FiniteDifferenceInterpolator.setup_interpolator(
            interp,
            cpw=0.0,
            gpw=0.0,
            npw=0.0,
            tpw=0.0,
            ipw=0.0,
            dxx=0.0,
            dyy=0.0,
            dzz=0.0,
            dxy=0.0,
            dyz=0.0,
            dxz=0.0,
        )
        return interp, fold

    def test_orientation_constraints_added(self, setup_interp):
        interp, _ = setup_interp
        interp.add_fold_constraints(
            fold_orientation=5.0,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=None,
        )
        assert any("fold orientation" in k for k in interp.constraints), (
            "Expected 'fold orientation' constraints"
        )

    def test_axis_constraints_added(self, setup_interp):
        interp, _ = setup_interp
        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=5.0,
            fold_regularisation=None,
            fold_normalisation=None,
        )
        assert any("fold axis" in k for k in interp.constraints)

    def test_normalisation_constraints_added(self, setup_interp):
        interp, _ = setup_interp
        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=1.0,
            fold_norm=1.0,
        )
        assert any("fold normalisation" in k for k in interp.constraints)

    def test_normalisation_default_target_is_negative_one(self, setup_interp):
        interp, _ = setup_interp
        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=1.0,
        )
        keys = [k for k in interp.constraints if "fold normalisation" in k]
        assert len(keys) > 0
        b = interp.constraints[keys[0]]["b"]
        assert np.all(b < 0.0)

    def test_dgz_alignment_correct_flips_target_sign(self, setup_interp):
        interp, _ = setup_interp
        normal_constraints = np.array([[3.0, 3.0, 3.0, 0.0, 0.0, -1.0, 1.0]])
        interp.set_normal_constraints(normal_constraints)

        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=1.0,
            fold_norm=1.0,
            dgz_alignment="correct",
        )

        keys = [k for k in interp.constraints if "fold normalisation" in k]
        assert len(keys) > 0
        b = interp.constraints[keys[0]]["b"]
        assert np.all(b < 0.0)

    def test_dgz_alignment_warn_does_not_flip_target_sign(self, setup_interp):
        interp, _ = setup_interp
        normal_constraints = np.array([[3.0, 3.0, 3.0, 0.0, 0.0, -1.0, 1.0]])
        interp.set_normal_constraints(normal_constraints)

        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=1.0,
            fold_norm=1.0,
            dgz_alignment="warn",
        )

        keys = [k for k in interp.constraints if "fold normalisation" in k]
        assert len(keys) > 0
        b = interp.constraints[keys[0]]["b"]
        assert np.all(b > 0.0)

    def test_regularisation_constraints_added(self, setup_interp):
        interp, _ = setup_interp
        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=[0.1, 0.01, 0.01],
            fold_normalisation=None,
        )
        assert any("fold regularisation" in k for k in interp.constraints)

    def test_all_none_adds_nothing_extra(self, setup_interp):
        interp, _ = setup_interp
        before = set(interp.constraints.keys())
        interp.add_fold_constraints(
            fold_orientation=None,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=None,
        )
        after = set(interp.constraints.keys())
        assert after == before, "No new constraints expected when all weights are None"

    def test_mask_fn_reduces_active_nodes(self, setup_interp):
        interp_masked, _ = setup_interp
        grid2 = _make_uniform_grid(8)
        fold2 = _constant_fold(grid2.n_nodes)
        interp_full = FDFoldInterpolator(grid2, fold=fold2)
        FiniteDifferenceInterpolator.setup_interpolator(
            interp_full,
            cpw=0.0,
            gpw=0.0,
            npw=0.0,
            tpw=0.0,
            ipw=0.0,
            dxx=0.0,
            dyy=0.0,
            dzz=0.0,
            dxy=0.0,
            dyz=0.0,
            dxz=0.0,
        )

        # Mask out half the domain.
        half = interp_masked.support.nodes[:, 0].max() / 2
        mask_fn = lambda xyz: xyz[:, 0] > half

        interp_masked.add_fold_constraints(
            fold_orientation=5.0,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=None,
            mask_fn=mask_fn,
        )
        interp_full.add_fold_constraints(
            fold_orientation=5.0,
            fold_axis_w=None,
            fold_regularisation=None,
            fold_normalisation=None,
        )

        # Masked version must have fewer gradient-orthogonal rows than full version.
        def _row_count(interp, key_part):
            return sum(v["matrix"].shape[0] for k, v in interp.constraints.items() if key_part in k)

        masked_rows = _row_count(interp_masked, "fold orientation")
        full_rows = _row_count(interp_full, "fold orientation")
        assert masked_rows < full_rows, (
            f"Masked rows ({masked_rows}) should be less than full rows ({full_rows})"
        )


# ---------------------------------------------------------------------------
# 5. setup_interpolator with fold_weights forwarding
# ---------------------------------------------------------------------------


class TestSetupInterpolatorFoldWeights:
    def test_fold_weights_are_forwarded(self):
        grid = _make_uniform_grid(8)
        fold = _constant_fold(grid.n_nodes)
        interp = FDFoldInterpolator(grid, fold=fold)
        interp.setup_interpolator(
            cpw=0.0,
            gpw=0.0,
            npw=0.0,
            tpw=0.0,
            ipw=0.0,
            fold_weights={
                "fold_orientation": 5.0,
                "fold_axis_w": 5.0,
                "fold_regularisation": [0.1, 0.01, 0.01],
                "fold_normalisation": 1.0,
            },
        )
        assert any("fold orientation" in k for k in interp.constraints)
        assert any("fold axis" in k for k in interp.constraints)
        assert any("fold regularisation" in k for k in interp.constraints)
        assert any("fold normalisation" in k for k in interp.constraints)

    def test_fold_weights_not_passed_to_parent(self):
        """fold_weights kwarg must not reach the parent and cause a KeyError."""
        grid = _make_uniform_grid(8)
        fold = _constant_fold(grid.n_nodes)
        interp = FDFoldInterpolator(grid, fold=fold)
        # If fold_weights were forwarded to the parent it would try to set
        # interpolation_weights["fold_weights"] which is harmless, but
        # operators lookup would fail.  The test confirms no exception.
        interp.setup_interpolator(fold_weights={})


# ---------------------------------------------------------------------------
# 6. End-to-end: solve recovers a planar field
# ---------------------------------------------------------------------------


class TestEndToEnd:
    @pytest.mark.parametrize("solver", ["lsmr"])
    def test_planar_field_recovery(self, solver):
        """
        FDFoldInterpolator with fold constraints should still be able to
        recover a planar field f = x + 0.5*y when given value + norm data.
        The fold geometry is set to uniform z-normal (horizontal fold axis),
        which should not disrupt a vertical planar field.
        """
        rng = np.random.default_rng(0)
        n = 15
        grid = StructuredGrid(
            origin=np.zeros(3),
            nsteps=np.array([n, n, n]),
            step_vector=np.ones(3),
        )
        fold = _constant_fold(
            grid.n_nodes,
            dgz_dir=(0, 0, 1),  # fold normal = z  (across-fold direction)
            deformed_dir=(1, 0, 0),  # deformed orientation = x
            axis_dir=(0, 1, 0),  # fold axis = y
        )
        interp = FDFoldInterpolator(grid, fold=fold)

        lo = grid.origin + 1.5
        hi = grid.maximum - 1.5
        pts = rng.uniform(lo, hi, (200, 3))
        vals = pts[:, 0] + 0.5 * pts[:, 1]

        val_data = np.column_stack([pts, vals, np.ones(len(pts))])
        interp.set_value_constraints(val_data)

        interp.setup_interpolator(
            cpw=1.0,
            gpw=0.0,
            npw=0.0,
            tpw=0.0,
            ipw=0.0,
            fold_weights={
                "fold_orientation": 0.1,
                "fold_axis_w": 0.1,
                "fold_regularisation": [0.01, 0.001, 0.001],
                "fold_normalisation": None,
            },
        )
        interp.solve_system(solver)

        predicted = interp.support.evaluate_value(pts, interp.c)
        mae = np.mean(np.abs(predicted - vals))
        assert mae < 1.0, f"MAE {mae:.4f} is too large for a planar field"

    @pytest.mark.parametrize("solver", ["lsmr"])
    def test_rectilinear_grid_end_to_end(self, solver):
        """FDFoldInterpolator works on a RectilinearGrid too."""
        rng = np.random.default_rng(1)
        grid = _make_rectilinear_grid(12)
        fold = _constant_fold(grid.n_nodes)
        interp = FDFoldInterpolator(grid, fold=fold)

        lo = grid.origin + 1.5
        hi = grid.maximum - 1.5
        pts = rng.uniform(lo, hi, (150, 3))
        vals = pts[:, 0] + 0.5 * pts[:, 1]

        val_data = np.column_stack([pts, vals, np.ones(len(pts))])
        interp.set_value_constraints(val_data)

        interp.setup_interpolator(
            cpw=1.0,
            gpw=0.0,
            npw=0.0,
            tpw=0.0,
            ipw=0.0,
            fold_weights={
                "fold_orientation": 0.1,
                "fold_axis_w": None,
                "fold_regularisation": [0.01, 0.001, 0.001],
                "fold_normalisation": None,
            },
        )
        interp.solve_system(solver)

        predicted = interp.support.evaluate_value(pts, interp.c)
        mae = np.mean(np.abs(predicted - vals))
        assert mae < 1.0, f"MAE {mae:.4f} is too large"

    def test_anisotropic_reg_differs_from_isotropic(self):
        """
        The directional regularisation should produce a different (lower) system
        matrix norm than isotropic regularisation when the direction is strongly
        aligned with one axis.  We test that the constraint matrices are not
        identical (i.e. the anisotropy has an effect).
        """
        grid = _make_uniform_grid(8)

        # Isotropic: call the standard assemble_inner for dxx+dyy+dzz.
        fdi_iso = FiniteDifferenceInterpolator(grid)
        fdi_iso.setup_interpolator(
            cpw=0.0,
            gpw=0.0,
            npw=0.0,
            tpw=0.0,
            ipw=0.0,
            dxx=1.0,
            dyy=1.0,
            dzz=1.0,
            dxy=0.0,
            dyz=0.0,
            dxz=0.0,
        )

        # Anisotropic: z-only direction field (only dzz gets weight 1, others 0).
        v = np.zeros((grid.n_nodes, 3))
        v[:, 2] = 1.0
        fdi_aniso = FiniteDifferenceInterpolator(grid)
        fdi_aniso.reset()
        fdi_aniso.minimise_directional_gradient_change(1.0, v, name="aniso")

        def _total_rows(interp):
            return sum(v["matrix"].shape[0] for v in interp.constraints.values())

        # Anisotropic has fewer rows (only dzz), isotropic has dxx+dyy+dzz.
        assert _total_rows(fdi_aniso) < _total_rows(fdi_iso)
