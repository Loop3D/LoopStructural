import pytest

from loop_common.supports import StructuredGrid, TetMesh


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
