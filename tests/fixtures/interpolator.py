from LoopStructural.interpolators import (
    FiniteDifferenceInterpolator as FDI,
    PiecewiseLinearInterpolator as PLI,
)
from LoopStructural.interpolators import StructuredGrid, TetMesh
from LoopStructural.utils import BoundingBox
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


@pytest.fixture(params=["PLI", "FDI"])
def interpolatortype(request):
    return request.param


@pytest.fixture(params=[1e3, 1e4, 4e4])
def nelements(request):
    nelements = request.param
    return nelements


@pytest.fixture()
def bounding_box():
    return BoundingBox(np.array([0, 0, 0]), np.array([1, 1, 1]))


@pytest.fixture(params=["PLI", "FDI"])
def interpolator_type(request):
    interpolator_type = request.param
    return interpolator_type


@pytest.fixture(params=["grid", "tetra"])
def support(request):
    support_type = request.param
    if support_type == "grid":
        return StructuredGrid()
    if support_type == "tetra":
        return TetMesh()


@pytest.fixture(params=["grid", "tetra"])
def support_class(request):
    support_type = request.param
    if support_type == "grid":
        return StructuredGrid
    if support_type == "tetra":
        return TetMesh


@pytest.fixture(params=["everywhere", "restricted"])
def region_func(request):
    region_type = request.param

    if region_type == "restricted":
        return lambda xyz: xyz[:, 0] > 0.5
    if region_type == "everywhere":
        return lambda xyz: np.ones(xyz.shape[0], dtype=bool)
