import numpy as np
from loop_interpolation import (
    ConstraintDiagnosticsReport,
    FiniteDifferenceInterpolator,
    PiecewiseLinearInterpolator,
    StructuredGrid,
    TetMesh,
)


def test_fdi_setup_returns_constraint_diagnostics_report():
    grid = StructuredGrid(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([4, 4, 4]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )
    interpolator = FiniteDifferenceInterpolator(grid)

    value_constraints = np.array(
        [
            [0.5, 0.5, 0.5, 1.0, 2.0],
            [20.0, 20.0, 20.0, 2.0, 1.0],
        ]
    )
    interpolator.set_value_constraints(value_constraints)

    report = interpolator.setup_interpolator(
        dxy=0.0,
        dyz=0.0,
        dxz=0.0,
        dxx=0.0,
        dyy=0.0,
        dzz=0.0,
        dx=0.0,
        dy=0.0,
        dz=0.0,
        cpw=1.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )

    assert isinstance(report, ConstraintDiagnosticsReport)
    assert "value" in report.families
    assert report.families["value"].row_count == 1
    assert report.families["value"].dropped_rows == 1
    assert report.outside_model_points["value"] == 1
    assert isinstance(report.summary(), str)
    assert "Active families" in report.summary()


def test_pli_setup_and_setup_wrapper_return_diagnostics_report():
    mesh = TetMesh(
        origin=np.array([0.0, 0.0, 0.0]),
        nsteps=np.array([4, 4, 4]),
        step_vector=np.array([1.0, 1.0, 1.0]),
    )
    interpolator = PiecewiseLinearInterpolator(mesh)

    value_constraints = np.array(
        [
            [0.5, 0.5, 0.5, 1.0, 1.5],
            [10.0, 10.0, 10.0, 0.0, 1.0],
        ]
    )
    interpolator.set_value_constraints(value_constraints)

    report = interpolator.setup(
        cgw=0.0,
        cpw=1.0,
        gpw=0.0,
        npw=0.0,
        tpw=0.0,
        ipw=0.0,
    )

    assert isinstance(report, ConstraintDiagnosticsReport)
    assert report is interpolator.latest_diagnostics_report
    assert report.families["value"].source_point_count == 2
    assert report.families["value"].row_count == 1
    assert report.families["value"].dropped_rows == 1
    assert report.outside_model_points["value"] == 1
    assert "Region coverage" in report.summary()
