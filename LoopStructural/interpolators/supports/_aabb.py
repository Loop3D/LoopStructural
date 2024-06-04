import numpy as np
from scipy import sparse


def _initialise_aabb(grid):
    """assigns the tetras to the grid cells where the bounding box
    of the tetra element overlaps the grid cell.
    It could be changed to use the separating axis theorem, however this would require
    significantly more calculations. (12 more I think).. #TODO test timing
    """
    # calculate the bounding box for all tetraherdon in the mesh
    # find the min/max extents for xyz
    # tetra_bb = np.zeros((grid.elements.shape[0], 19, 3))
    minx = np.min(grid.nodes[grid.elements[:, :4], 0], axis=1)
    maxx = np.max(grid.nodes[grid.elements[:, :4], 0], axis=1)
    miny = np.min(grid.nodes[grid.elements[:, :4], 1], axis=1)
    maxy = np.max(grid.nodes[grid.elements[:, :4], 1], axis=1)

    cell_indexes = grid.aabb_grid.global_index_to_cell_index(np.arange(grid.aabb_grid.n_elements))
    corners = grid.aabb_grid.cell_corner_indexes(cell_indexes)
    positions = grid.aabb_grid.node_indexes_to_position(corners)
    ## Because we known the node orders just select min/max from each
    # coordinate. Use these to check whether the tetra is in the cell
    x_boundary = positions[:, [0, 1], 0]
    y_boundary = positions[:, [0, 2], 1]
    a = np.logical_and(
        minx[None, :] > x_boundary[:, None, 0],
        minx[None, :] < x_boundary[:, None, 1],
    )  # min point between cell
    b = np.logical_and(
        maxx[None, :] < x_boundary[:, None, 1],
        maxx[None, :] > x_boundary[:, None, 0],
    )  # max point between cell
    c = np.logical_and(
        minx[None, :] < x_boundary[:, None, 0],
        maxx[None, :] > x_boundary[:, None, 0],
    )  # min point < than cell & max point > cell

    x_logic = np.logical_or(np.logical_or(a, b), c)

    a = np.logical_and(
        miny[None, :] > y_boundary[:, None, 0],
        miny[None, :] < y_boundary[:, None, 1],
    )  # min point between cell
    b = np.logical_and(
        maxy[None, :] < y_boundary[:, None, 1],
        maxy[None, :] > y_boundary[:, None, 0],
    )  # max point between cell
    c = np.logical_and(
        miny[None, :] < y_boundary[:, None, 0],
        maxy[None, :] > y_boundary[:, None, 0],
    )  # min point < than cell & max point > cell

    y_logic = np.logical_or(np.logical_or(a, b), c)
    logic = np.logical_and(x_logic, y_logic)

    if grid.dimension == 3:
        z_boundary = positions[:, [0, 6], 2]
        minz = np.min(grid.nodes[grid.elements[:, :4], 2], axis=1)
        maxz = np.max(grid.nodes[grid.elements[:, :4], 2], axis=1)
        a = np.logical_and(
            minz[None, :] > z_boundary[:, None, 0],
            minz[None, :] < z_boundary[:, None, 1],
        )  # min point between cell
        b = np.logical_and(
            maxz[None, :] < z_boundary[:, None, 1],
            maxz[None, :] > z_boundary[:, None, 0],
        )  # max point between cell
        c = np.logical_and(
            minz[None, :] < z_boundary[:, None, 0],
            maxz[None, :] > z_boundary[:, None, 0],
        )  # min point < than cell & max point > cell

        z_logic = np.logical_or(np.logical_or(a, b), c)
        logic = np.logical_and(logic, z_logic)

    grid._aabb_table = sparse.csr_matrix(logic)
