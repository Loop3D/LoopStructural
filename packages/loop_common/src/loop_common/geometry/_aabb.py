import numpy as np
from scipy import sparse


def _initialise_aabb(grid):
    minx = np.min(grid.nodes[grid.elements[:, :4], 0], axis=1)
    maxx = np.max(grid.nodes[grid.elements[:, :4], 0], axis=1)
    miny = np.min(grid.nodes[grid.elements[:, :4], 1], axis=1)
    maxy = np.max(grid.nodes[grid.elements[:, :4], 1], axis=1)

    cell_indexes = grid.aabb_grid.global_index_to_cell_index(np.arange(grid.aabb_grid.n_elements))
    corners = grid.aabb_grid.cell_corner_indexes(cell_indexes)
    positions = grid.aabb_grid.node_indexes_to_position(corners)
    x_boundary = positions[:, [0, 1], 0]
    y_boundary = positions[:, [0, 2], 1]
    a = np.logical_and(minx[None, :] > x_boundary[:, None, 0], minx[None, :] < x_boundary[:, None, 1])
    b = np.logical_and(maxx[None, :] < x_boundary[:, None, 1], maxx[None, :] > x_boundary[:, None, 0])
    c = np.logical_and(minx[None, :] < x_boundary[:, None, 0], maxx[None, :] > x_boundary[:, None, 0])
    x_logic = np.logical_or(np.logical_or(a, b), c)

    a = np.logical_and(miny[None, :] > y_boundary[:, None, 0], miny[None, :] < y_boundary[:, None, 1])
    b = np.logical_and(maxy[None, :] < y_boundary[:, None, 1], maxy[None, :] > y_boundary[:, None, 0])
    c = np.logical_and(miny[None, :] < y_boundary[:, None, 0], maxy[None, :] > y_boundary[:, None, 0])
    y_logic = np.logical_or(np.logical_or(a, b), c)
    logic = np.logical_and(x_logic, y_logic)

    if grid.dimension == 3:
        z_boundary = positions[:, [0, 6], 2]
        minz = np.min(grid.nodes[grid.elements[:, :4], 2], axis=1)
        maxz = np.max(grid.nodes[grid.elements[:, :4], 2], axis=1)
        a = np.logical_and(minz[None, :] > z_boundary[:, None, 0], minz[None, :] < z_boundary[:, None, 1])
        b = np.logical_and(maxz[None, :] < z_boundary[:, None, 1], maxz[None, :] > z_boundary[:, None, 0])
        c = np.logical_and(minz[None, :] < z_boundary[:, None, 0], maxz[None, :] > z_boundary[:, None, 0])
        z_logic = np.logical_or(np.logical_or(a, b), c)
        logic = np.logical_and(logic, z_logic)

    grid._aabb_table = sparse.csr_matrix(logic)
