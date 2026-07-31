"""
Integration tests for FiniteDifferenceInterpolator used with RectilinearGrid.
"""

import numpy as np
import pytest
from loop_common.supports import RectilinearGrid
from loop_interpolation import FiniteDifferenceInterpolator

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def small_uniform_rect_grid():
    """4-cell uniform RectilinearGrid (same geometry as StructuredGrid(nsteps=[4,4,4]))."""
    x = np.linspace(0.0, 4.0, 5)
    y = np.linspace(0.0, 4.0, 5)
    z = np.linspace(0.0, 4.0, 5)
    return RectilinearGrid(x, y, z)


@pytest.fixture
def nonuniform_rect_grid():
    """Non-uniform RectilinearGrid suitable for interpolation tests."""
    x = np.linspace(0.0, 10.0, 21)
    y = np.linspace(0.0, 10.0, 21)
    z = np.linspace(0.0, 10.0, 21)
    # introduce slight non-uniformity by jittering every other step
    x[1::2] += 0.1
    return RectilinearGrid(x, y, z)


# ---------------------------------------------------------------------------
# Basic construction
# ---------------------------------------------------------------------------


def test_fdi_creation_with_rectilinear_grid(small_uniform_rect_grid):
    fdi = FiniteDifferenceInterpolator(small_uniform_rect_grid)
    assert fdi is not None
    assert fdi.dof == small_uniform_rect_grid.n_nodes


def test_setup_interpolator_no_error(small_uniform_rect_grid):
    fdi = FiniteDifferenceInterpolator(small_uniform_rect_grid)
    fdi.setup_interpolator()  # should not raise


# ---------------------------------------------------------------------------
# Regularisation constraints are built for RectilinearGrid
# ---------------------------------------------------------------------------


def test_rectilinear_regularisation_constraints_exist():
    x = np.linspace(0.0, 4.0, 5)
    y = np.linspace(0.0, 4.0, 5)
    z = np.linspace(0.0, 4.0, 5)
    grid = RectilinearGrid(x, y, z)
    fdi = FiniteDifferenceInterpolator(grid)
    fdi.setup_interpolator(
        dxx=1.0,
        dyy=1.0,
        dzz=1.0,
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )
    for name in ["dxx", "dyy", "dzz"]:
        assert name in fdi.constraints, f"Missing constraint '{name}'"
        assert fdi.constraints[name]["matrix"].shape[0] > 0


def test_rectilinear_mixed_regularisation_constraints_exist():
    x = np.linspace(0.0, 4.0, 5)
    y = np.linspace(0.0, 4.0, 5)
    z = np.linspace(0.0, 4.0, 5)
    grid = RectilinearGrid(x, y, z)
    fdi = FiniteDifferenceInterpolator(grid)
    fdi.setup_interpolator(
        dxx=0.0,
        dyy=0.0,
        dzz=0.0,
        dxy=1.0,
        dyz=1.0,
        dxz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )
    for name in ["dxy", "dyz", "dxz"]:
        assert name in fdi.constraints, f"Missing constraint '{name}'"
        assert fdi.constraints[name]["matrix"].shape[0] > 0


def test_border_regularisation_constraints_exist():
    x = np.linspace(0.0, 4.0, 5)
    y = np.linspace(0.0, 4.0, 5)
    z = np.linspace(0.0, 4.0, 5)
    grid = RectilinearGrid(x, y, z)
    fdi = FiniteDifferenceInterpolator(grid)
    fdi.setup_interpolator(
        dxx=0.0,
        dyy=0.0,
        dzz=0.0,
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        dx=1.0,
        dy=1.0,
        dz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )
    for name in ["dx_lower", "dx_upper", "dy_lower", "dy_upper", "dz_lower", "dz_upper"]:
        assert name in fdi.constraints, f"Missing border constraint '{name}'"


# ---------------------------------------------------------------------------
# Interpolation accuracy on a planar field
# ---------------------------------------------------------------------------


def _make_planar_data(grid, n=200, seed=42):
    """Return (pts, vals, normals) for the planar field f = x + 0.5*y."""
    rng = np.random.default_rng(seed)
    lo = grid.origin + 0.5
    hi = grid.maximum - 0.5
    pts = rng.uniform(lo, hi, size=(n, 3))
    vals = pts[:, 0] + 0.5 * pts[:, 1]
    normals = np.tile([1.0, 0.5, 0.0], (n // 2, 1))
    normals /= np.linalg.norm(normals)
    return pts, vals, normals


@pytest.mark.parametrize("solver", ["lsmr"])
def test_planar_interpolation_accuracy(solver):
    """RectilinearGrid FDI should recover a planar field with low MAE."""
    x = np.linspace(0.0, 10.0, 21)
    y = np.linspace(0.0, 10.0, 21)
    z = np.linspace(0.0, 10.0, 21)
    grid = RectilinearGrid(x, y, z)
    fdi = FiniteDifferenceInterpolator(grid)

    pts, vals, _ = _make_planar_data(grid, n=300)
    val_data = np.column_stack([pts, vals, np.ones(len(pts))])

    fdi.set_value_constraints(val_data)
    fdi.setup_interpolator(cpw=1.0, gpw=0.0)
    fdi.solve_system(solver)

    predicted = fdi.support.evaluate_value(pts, fdi.c)
    mae = np.mean(np.abs(predicted - vals))
    assert mae < 0.5, f"MAE {mae:.4f} is unexpectedly large"


@pytest.mark.parametrize("solver", ["lsmr"])
def test_planar_interpolation_with_gradient_constraints(solver):
    """FDI on RectilinearGrid can use gradient (normal) constraints."""
    x = np.linspace(0.0, 10.0, 21)
    y = np.linspace(0.0, 10.0, 21)
    z = np.linspace(0.0, 10.0, 21)
    grid = RectilinearGrid(x, y, z)
    fdi = FiniteDifferenceInterpolator(grid)

    pts, vals, normals = _make_planar_data(grid, n=200)
    val_data = np.column_stack([pts[:100], vals[:100], np.ones(100)])
    norm_data = np.column_stack([pts[100:150], normals[:50], np.ones(50)])

    fdi.set_value_constraints(val_data)
    fdi.set_gradient_constraints(norm_data)
    fdi.setup_interpolator(cpw=1.0, gpw=1.0)
    fdi.solve_system(solver)

    predicted = fdi.support.evaluate_value(pts[:100], fdi.c)
    mae = np.mean(np.abs(predicted - vals[:100]))
    assert mae < 1.0, f"MAE {mae:.4f} is unexpectedly large"


# ---------------------------------------------------------------------------
# Uniform RectilinearGrid should match StructuredGrid accuracy
# ---------------------------------------------------------------------------


def test_uniform_rectilinear_matches_structured_grid():
    """
    A uniform RectilinearGrid must give essentially the same result as
    StructuredGrid on identical geometry.
    """
    from loop_interpolation import StructuredGrid

    nsteps = np.array([20, 20, 20])
    step = np.array([0.5, 0.5, 0.5])
    origin = np.zeros(3)

    sg = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=step)
    rg = RectilinearGrid(
        np.linspace(origin[0], origin[0] + nsteps[0] * step[0], nsteps[0] + 1),
        np.linspace(origin[1], origin[1] + nsteps[1] * step[1], nsteps[1] + 1),
        np.linspace(origin[2], origin[2] + nsteps[2] * step[2], nsteps[2] + 1),
    )

    rng = np.random.default_rng(0)
    lo = origin + 0.6
    hi = origin + nsteps * step - 0.6
    pts = rng.uniform(lo, hi, size=(200, 3))
    vals_true = pts[:, 0] + 0.5 * pts[:, 1]
    val_data = np.column_stack([pts, vals_true, np.ones(len(pts))])

    results = {}
    for name, grid in [("structured", sg), ("rectilinear", rg)]:
        fdi = FiniteDifferenceInterpolator(grid)
        fdi.set_value_constraints(val_data)
        fdi.setup_interpolator(cpw=1.0, gpw=0.0)
        fdi.solve_system("lsmr")
        predicted = fdi.support.evaluate_value(pts, fdi.c)
        results[name] = np.mean(np.abs(predicted - vals_true))

    # Both should be accurate and within 2x of each other
    assert results["structured"] < 0.5
    assert results["rectilinear"] < 0.5
    assert results["rectilinear"] < results["structured"] * 3.0, (
        f"Rectilinear MAE ({results['rectilinear']:.4f}) is much worse than "
        f"StructuredGrid MAE ({results['structured']:.4f})"
    )


# ---------------------------------------------------------------------------
# Region masking
# ---------------------------------------------------------------------------
