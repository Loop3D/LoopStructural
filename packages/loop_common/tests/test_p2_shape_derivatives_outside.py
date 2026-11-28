"""evaluate_shape_derivatives must handle points outside the mesh."""

import numpy as np
import pytest
from loop_common.supports._3d_p2_tetra import P2UnstructuredTetMesh
from loop_common.supports._p2_structured_tetra import P2TetMesh


def _unstructured():
    return P2UnstructuredTetMesh(
        origin=np.zeros(3), step_vector=np.ones(3) / 4, nsteps_cells=np.array([4, 4, 4])
    )


def _structured():
    return P2TetMesh(origin=np.zeros(3), nsteps=np.array([3, 3, 3]), step_vector=np.ones(3))


@pytest.mark.parametrize(
    "make,inside_pt,outside_pt",
    [
        (_unstructured, [0.41, 0.52, 0.33], [5.0, 5.0, 5.0]),
        (_structured, [1.1, 1.2, 1.3], [50.0, 50.0, 50.0]),
    ],
)
def test_outside_point_gives_nan_and_does_not_poison_batch(make, inside_pt, outside_pt):
    mesh = make()
    ref = mesh.evaluate_shape_derivatives(np.array([inside_pt]))
    pts = np.array([inside_pt, outside_pt, inside_pt])
    res = mesh.evaluate_shape_derivatives(pts)
    assert len(res) == 3
    dN, elements, inside = res
    assert dN.shape == (3, 3, 10)
    assert elements[1] == -1
    assert not inside[1]
    assert np.all(np.isnan(dN[1]))
    assert inside[0] and inside[2]
    assert np.all(np.isfinite(dN[0])) and np.all(np.isfinite(dN[2]))
    np.testing.assert_allclose(dN[0], ref[0][0])
