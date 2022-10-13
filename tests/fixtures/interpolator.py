from LoopStructural.interpolators import (
    FiniteDifferenceInterpolator as FDI,
    PiecewiseLinearInterpolator as PLI,
)
from LoopStructural.interpolators import StructuredGrid, TetMesh

import pytest
import numpy as np


@pytest.fixture(params=["FDI", "PLI"])
def interpolator(request):
    interpolator = request.param
    origin = np.array([-0.1, -0.1, -0.1])
    maximum = np.array([1.1, 1.1, 1.1])
    nsteps = np.array([20, 20, 20])
    step_vector = (maximum - origin) / nsteps
    if interpolator == "FDI":
        grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=step_vector)
        interpolator = FDI(grid)
        return interpolator
    elif interpolator == "PLI":
        grid = TetMesh(origin=origin, nsteps=nsteps, step_vector=step_vector)
        interpolator = PLI(grid)
        return interpolator
    else:
        raise ValueError(f"Invalid interpolator: {interpolator}")
