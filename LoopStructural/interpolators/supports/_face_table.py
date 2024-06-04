import numpy as np
from scipy import sparse


def _init_face_table(grid):
    """
    Fill table containing elements that share a face, and another
    table that contains the nodes for a face.
    """
    # need to identify the shared nodes for pairs of elements
    # we do this by creating a sparse matrix that has N rows (number of elements)
    # and M columns (number of nodes).
    # We then fill the location where a node is in an element with true
    # Then we create a table for the pairs of elements in the mesh
    # we have the neighbour relationships, which are the 4 neighbours for each element
    # create a new table that shows the element index repeated four times
    # flatten both of these arrays so we effectively have a table with pairs of neighbours
    # disgard the negative neighbours because these are border neighbours
    rows = np.tile(np.arange(grid.n_elements)[:, None], (1, grid.dimension + 1))
    elements = grid.elements
    neighbours = grid.neighbours
    # add array of bool to the location where there are elements for each node

    # use this to determine shared faces

    element_nodes = sparse.coo_matrix(
        (
            np.ones(elements.shape[0] * (grid.dimension + 1)),
            (rows.ravel(), elements[:, : grid.dimension + 1].ravel()),
        ),
        shape=(grid.n_elements, grid.n_nodes),
        dtype=bool,
    ).tocsr()
    n1 = np.tile(np.arange(neighbours.shape[0], dtype=int)[:, None], (1, grid.dimension + 1))
    n1 = n1.flatten()
    n2 = neighbours.flatten()
    n1 = n1[n2 >= 0]
    n2 = n2[n2 >= 0]
    el_rel = np.zeros((grid.neighbours.flatten().shape[0], 2), dtype=int)
    el_rel[:] = -1
    el_rel[np.arange(n1.shape[0]), 0] = n1
    el_rel[np.arange(n1.shape[0]), 1] = n2
    el_rel = el_rel[el_rel[:, 0] >= 0, :]

    # el_rel2 = np.zeros((grid.neighbours.flatten().shape[0], 2), dtype=int)
    grid._shared_element_relationships[:] = -1
    el_pairs = sparse.coo_matrix((np.ones(el_rel.shape[0]), (el_rel[:, 0], el_rel[:, 1]))).tocsr()
    i, j = sparse.tril(el_pairs).nonzero()
    grid._shared_element_relationships[: len(i), 0] = i
    grid._shared_element_relationships[: len(i), 1] = j

    grid._shared_element_relationships = grid.shared_element_relationships[
        grid.shared_element_relationships[:, 0] >= 0, :
    ]

    faces = element_nodes[grid.shared_element_relationships[:, 0], :].multiply(
        element_nodes[grid.shared_element_relationships[:, 1], :]
    )
    shared_faces = faces[np.array(np.sum(faces, axis=1) == grid.dimension).flatten(), :]
    row, col = shared_faces.nonzero()
    row = row[row.argsort()]
    col = col[row.argsort()]
    shared_face_index = np.zeros((shared_faces.shape[0], grid.dimension), dtype=int)
    shared_face_index[:] = -1
    shared_face_index[row.reshape(-1, grid.dimension)[:, 0], :] = col.reshape(-1, grid.dimension)
    grid._shared_elements[np.arange(grid.shared_element_relationships.shape[0]), :] = (
        shared_face_index
    )
    # resize
    grid._shared_elements = grid.shared_elements[: len(grid.shared_element_relationships), :]
