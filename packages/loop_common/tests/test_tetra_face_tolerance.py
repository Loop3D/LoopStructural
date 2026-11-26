"""A point on a tetrahedron face must be found despite rounding.

On a face the barycentric coordinate across the face is exactly 0. Rounding can
give a value a few 1e-17 below 0, and the point was lost: its value was NaN. This
happened on some platforms, for points on a lattice plane of the support.
"""

import numpy as np
import pytest
from loop_common.supports import P2TetMesh, TetMesh

FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def _meshes():
    origin = np.array([0.1, 0.2, 0.3])
    step = np.array([0.7, 0.9, 1.1])
    nsteps = np.array([6, 6, 6])
    return {
        "TetMesh": TetMesh(origin=origin, nsteps_cells=nsteps, step_vector=step),
        "P2TetMesh": P2TetMesh(origin=origin, nsteps=nsteps + 1, step_vector=step),
    }, origin, step * nsteps


def _interior_face_centres(mesh, origin, extent):
    centres = []
    for element in mesh.elements[:, :4]:
        for face in FACES:
            centres.append(mesh.nodes[element[list(face)]].mean(axis=0))
    centres = np.unique(np.round(np.array(centres), 9), axis=0)
    keep = np.all((centres > origin + 1e-3) & (centres < origin + extent - 1e-3), axis=1)
    return centres[keep]


@pytest.mark.parametrize("name", ["TetMesh", "P2TetMesh"])
def test_points_on_interior_faces_are_found_with_rounding_noise(name):
    meshes, origin, extent = _meshes()
    mesh = meshes[name]
    centres = _interior_face_centres(mesh, origin, extent)
    assert len(centres) > 100
    rng = np.random.default_rng(1)
    scale = np.abs(centres).max()
    for _ in range(5):
        # The noise is about 1 ulp of the coordinates: far below the tolerance.
        noisy = centres + rng.uniform(-1e-16, 1e-16, size=centres.shape) * scale
        inside = mesh.get_element_for_location(noisy)[3]
        assert inside.all(), f"{(~inside).sum()} of {len(inside)} face points lost"


def test_points_outside_the_mesh_stay_outside():
    meshes, origin, extent = _meshes()
    mesh = meshes["TetMesh"]
    outside = np.array(
        [
            origin + extent * [0.5, 0.5, 0.5] + [-1.0, 0.0, 0.0] * extent,
            origin + extent * [0.5, 0.5, 0.5] + [0.0, 1.0, 0.0] * extent,
            origin + extent * [0.5, 0.5, 0.5] + [0.0, 0.0, 1.0] * extent,
        ]
    )
    inside = mesh.get_element_for_location(outside)[3]
    assert not inside.any()
