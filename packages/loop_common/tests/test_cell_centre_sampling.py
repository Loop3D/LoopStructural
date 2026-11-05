"""Cell centres of a BoundingBox agree with the cells of its structured grids, also for rotated boxes."""

import itertools

import numpy as np
import pytest
from loop_common.geometry import BoundingBox

BOXES = {
    "unit": dict(origin=[0.0, 0.0, 0.0], maximum=[1.0, 1.0, 1.0], nsteps=[10, 10, 10]),
    "offset_anisotropic": dict(
        origin=[100.0, -20.0, 5.0], maximum=[112.0, 4.0, 35.0], nsteps=[7, 5, 4]
    ),
}
ANGLES = [0.0, 30.0, 90.0]
LOCAL_ORIGINS = {"origin0": (0.0, 0.0, 0.0), "origin_shifted": (7.0, -3.0, 11.0)}


def _rotation_about_z(angle):
    radians = np.deg2rad(angle)
    rotation = np.eye(3)
    rotation[:2, :2] = [[np.cos(radians), -np.sin(radians)], [np.sin(radians), np.cos(radians)]]
    return rotation


@pytest.fixture(
    params=list(itertools.product(BOXES, ANGLES, LOCAL_ORIGINS)),
    ids=lambda p: "-".join(map(str, p)),
)
def box(request):
    name, angle, origin_name = request.param
    box = BoundingBox(**BOXES[name])
    box.set_local_transform(
        local_origin=np.array(LOCAL_ORIGINS[origin_name]), rotation_matrix=_rotation_about_z(angle)
    )
    return box


def _tolerance(box):
    """Point tolerance, the same as in test_rigid_transform_contract."""
    length = np.linalg.norm(box.maximum - box.origin)
    world_magnitude = max(np.linalg.norm(box.origin), np.linalg.norm(box.maximum))
    return 1e-10 * length + 32 * np.finfo(float).eps * max(world_magnitude, length)


def _world_centres(box):
    return box.reproject(box.cell_centres(order="F"))


def test_each_cell_has_exactly_one_centre(box):
    """Each cell of the box local grid contains exactly one cell centre.

    This test is evidence for FR-015 and SC-011.

    The test works in the local frame of the box. The enclosure grid of
    structured_grid(local_coordinates=True) has axis-aligned cells. These
    cells differ from the cells of a rotated box. A one-to-one count against
    that grid does not hold (measured: 8 of 12 cells fail with the first
    version of this test).
    """
    n_cells = int(np.prod(box.nsteps - 1))
    centres = box.cell_centres(order="F")
    # World round trip must give the same local centres.
    local = box.project(box.reproject(centres))
    np.testing.assert_allclose(local, centres, atol=_tolerance(box), rtol=0.0)
    # The grid nodes are linspace(origin, maximum) in the frame of the box, and
    # project() maps them to the local frame. So go back with reproject(), then
    # find the cell index on that grid: x changes fastest (order F).
    origin = box.regular_grid(order="F", local=False)[0]
    step = box.node_step_vector
    cell = np.floor((box.reproject(local) - origin) / step).astype(int)
    nx, ny, nz = (int(n) - 1 for n in box.nsteps)
    assert len(cell) == n_cells
    assert np.all(cell >= 0)
    assert np.all(cell < np.array([nx, ny, nz]))
    index = cell[:, 0] + nx * (cell[:, 1] + ny * cell[:, 2])
    np.testing.assert_array_equal(
        np.bincount(index, minlength=n_cells), np.ones(n_cells, dtype=int)
    )
    np.testing.assert_array_equal(index, np.arange(n_cells))


def test_local_coordinates_grid_associates_own_centres(box):
    """Each cell centre falls in a cell of the local grid.

    This test is evidence for FR-015 and SC-011.
    """
    grid = box.structured_grid(local_coordinates=True)
    assert np.all(grid.cell_index_at(box.project(_world_centres(box))) != -1)


def test_centres_are_cell_midpoints_in_world(box):
    """Each world centre is the mean of the 8 corner nodes of its cell."""
    grid = box.structured_grid()
    nx, ny, nz = (int(n) for n in box.nsteps)
    k, j, i = np.meshgrid(np.arange(nz - 1), np.arange(ny - 1), np.arange(nx - 1), indexing="ij")
    i, j, k = i.ravel(), j.ravel(), k.ravel()  # x varies fastest, as the flat cell index
    corners = np.zeros((len(i), 3))
    for di, dj, dk in itertools.product((0, 1), repeat=3):
        corners += grid.nodes[(i + di) + nx * ((j + dj) + ny * (k + dk))]
    np.testing.assert_allclose(_world_centres(box), corners / 8.0, atol=_tolerance(box), rtol=0.0)


def test_centres_match_exported_cells(box):
    """The world centres equal the cell centres of the exported grid."""
    grid = box.structured_grid()
    np.testing.assert_allclose(
        _world_centres(box), grid.cell_centres, atol=_tolerance(box), rtol=0.0
    )


def test_centres_inside_box(box):
    """Each world centre is inside the box limits."""
    world = _world_centres(box)
    tolerance = _tolerance(box)
    assert np.all(world >= box.origin - tolerance)
    assert np.all(world <= box.maximum + tolerance)
