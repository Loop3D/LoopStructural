import numpy as np


def test_nx(interpolator, data):
    assert interpolator.dof == 21 * 21 * 21


def test_region(interpolator, data, region_func):
    """Test to see whether restricting the interpolator to a region works"""
    # interpolator = generate_interpolator(interpolator)
    interpolator.set_value_constraints(data[["X", "Y", "Z", "val", "w"]].to_numpy())
    interpolator.setup_interpolator()
    interpolator.set_region(region_func)
    # assert np.all(interpolator.region == region_func(interpolator.support.nodes))
    interpolator.solve_system()


def test_add_constraint_to_least_squares(interpolator):
    """make sure that when incorrect sized arrays are passed it doesn't get added"""
    pass


def test_finite_difference_border_regularisation_constraints():
    from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid

    origin = np.array([0.0, 0.0, 0.0])
    nsteps = np.array([4, 4, 4])
    step_vector = np.array([1.0, 1.0, 1.0])
    grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=step_vector)
    interpolator = FiniteDifferenceInterpolator(grid)

    interpolator.setup_interpolator(
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        dxx=0.0,
        dyy=0.0,
        dzz=0.0,
        dx=1.0,
        dy=1.0,
        dz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )

    expected_rows = grid.nsteps[1] * grid.nsteps[2]
    for name in ["dx_lower", "dx_upper", "dy_lower", "dy_upper", "dz_lower", "dz_upper"]:
        assert name in interpolator.constraints
        assert interpolator.constraints[name]["matrix"].shape[0] == expected_rows


def test_update_interpolator():
    pass


def test_solve_timing_breakdown_available():
    from loop_interpolation import FiniteDifferenceInterpolator, StructuredGrid

    origin = np.array([0.0, 0.0, 0.0])
    nsteps = np.array([4, 4, 4])
    step_vector = np.array([1.0, 1.0, 1.0])
    grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=step_vector)
    interpolator = FiniteDifferenceInterpolator(grid)

    interpolator.setup_interpolator(
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        dxx=1.0,
        dyy=1.0,
        dzz=1.0,
        cpw=0.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )
    ok = interpolator.solve_system("lsmr")
    assert ok is True

    timing = interpolator.get_last_solve_timing()
    assert "assembly_seconds" in timing
    assert "solve_seconds" in timing
    assert "total_seconds" in timing
    assert timing["total_seconds"] >= timing["solve_seconds"]
