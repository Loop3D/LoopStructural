"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""

from typing import Tuple

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, tril

from loop_common.logging import get_logger as getLogger

from . import StructuredGrid, SupportType
from ._base_support import BaseSupport

logger = getLogger(__name__)


def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cross product for stacks of 3-vectors, ~10x faster than np.cross
    for this shape because it skips np.cross's generic axis handling.
    """
    return np.stack(
        (
            a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
            a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
            a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0],
        ),
        axis=1,
    )


class UnStructuredTetMesh(BaseSupport):
    """ """

    dimension = 3

    def __init__(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        neighbours: np.ndarray,
        aabb_nsteps=None,
    ):
        """An unstructured mesh defined by nodes, elements and neighbours
        An axis aligned bounding box (AABB) grid is used to speed up finding
        which tetra a point is in: each grid cell is sized from the tetra's
        bounding box extent (see _initialise_aabb) so it only needs to test
        a handful of candidate tetra per query point.

        Parameters
        ----------
        nodes : array or array like
            container of vertex locations
        elements : array or array like, dtype cast to long
            container of tetra indicies
        neighbours : array or array like, dtype cast to long
            array containing element neighbours
        aabb_nsteps : list, optional
            force nsteps for aabb, by default None
        """
        self.type = SupportType.UnStructuredTetMesh
        self._nodes = np.array(nodes)
        if self._nodes.shape[1] != 3:
            raise ValueError("Nodes must be 3D")
        self.neighbours = np.array(neighbours, dtype=np.int64)
        if self.neighbours.shape[1] != 4:
            raise ValueError("Neighbours array is too big")
        self._elements = np.array(elements, dtype=np.int64)
        if self.elements.shape[0] != self.neighbours.shape[0]:
            raise ValueError("Number of elements and neighbours do not match")
        self._barycentre = np.sum(self.nodes[self.elements[:, :4]][:, :, :], axis=1) / 4.0
        self.minimum = np.min(self.nodes, axis=0)
        self.maximum = np.max(self.nodes, axis=0)
        length = self.maximum - self.minimum
        self.minimum -= length * 0.1
        self.maximum += length * 0.1
        # cache each tetra's bounding box; reused here to size the aabb grid
        # and again in _initialise_aabb to bucket tetra into grid cells.
        tetra_nodes = self.nodes[self.elements[:, :4]]
        self._tetra_bbox_min = np.min(tetra_nodes, axis=1)
        self._tetra_bbox_max = np.max(tetra_nodes, axis=1)
        if aabb_nsteps is None:
            # Size grid cells from the tetrahedra's *bounding box* extent, not
            # their volume: a thin/skewed tetrahedron's AABB can be many times
            # larger than its volume would suggest, so a volume-based estimate
            # undersizes cells and each tetra ends up straddling many cells.
            # Cells about half the typical (median) tetra AABB edge keep the
            # per-cell candidate count low (a handful of tetra per cell)
            # without letting the sparse table blow up in size.
            bbox_extent = self._tetra_bbox_max - self._tetra_bbox_min
            median_extent = np.median(bbox_extent, axis=0)
            median_extent = np.maximum(median_extent, np.max(length) * 1e-6)
            step_vector = median_extent * 0.5
            # number of steps is the length of the box / step vector
            aabb_nsteps = np.ceil((self.maximum - self.minimum) / step_vector).astype(int)
            # make sure there is at least one cell in every dimension
            aabb_nsteps[aabb_nsteps < 2] = 2
        aabb_nsteps = np.array(aabb_nsteps, dtype=int)
        step_vector = (self.maximum - self.minimum) / (aabb_nsteps - 1)
        self.aabb_grid = StructuredGrid(self.minimum, nsteps=aabb_nsteps, step_vector=step_vector)
        # make a big table to store which tetra are in which element.
        # if this takes up too much memory it could be simplified by using sparse matrices or dict but
        # at the expense of speed
        self.aabb_table = csr_matrix((self.aabb_grid.n_elements, len(self.elements)), dtype=bool)
        self.shared_element_relationships = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 2), dtype=int
        )
        self.shared_elements = np.zeros(
            (self.neighbours[self.neighbours >= 0].flatten().shape[0], 3), dtype=int
        )
        self._init_face_table()
        self._initialise_aabb()

    def set_nelements(self, nelements):
        raise NotImplementedError("Cannot set number of elements for unstructured mesh")

    @property
    def nodes(self):
        return self._nodes

    @property
    def elements(self):
        return self._elements

    @property
    def barycentre(self):
        return self._barycentre

    @property
    def n_nodes(self):
        return self.nodes.shape[0]

    def onGeometryChange(self):
        pass

    def _init_face_table(self):
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
        rows = np.tile(np.arange(self.n_elements)[:, None], (1, 4))
        elements = self.get_elements()
        neighbours = self.get_neighbours()
        # add array of bool to the location where there are elements for each node

        # use this to determine shared faces

        element_nodes = coo_matrix(
            (np.ones(elements.shape[0] * 4), (rows.ravel(), elements[:, :4].ravel())),
            shape=(self.n_elements, self.n_nodes),
            dtype=bool,
        ).tocsr()
        n1 = np.tile(np.arange(neighbours.shape[0], dtype=int)[:, None], (1, 4))
        n1 = n1.flatten()
        n2 = neighbours.flatten()
        n1 = n1[n2 >= 0]
        n2 = n2[n2 >= 0]
        el_rel = np.zeros((self.neighbours.flatten().shape[0], 2), dtype=int)
        el_rel[:] = -1
        el_rel[np.arange(n1.shape[0]), 0] = n1
        el_rel[np.arange(n1.shape[0]), 1] = n2
        el_rel = el_rel[el_rel[:, 0] >= 0, :]

        # el_rel2 = np.zeros((self.neighbours.flatten().shape[0], 2), dtype=int)
        self.shared_element_relationships[:] = -1
        el_pairs = coo_matrix((np.ones(el_rel.shape[0]), (el_rel[:, 0], el_rel[:, 1]))).tocsr()
        i, j = tril(el_pairs).nonzero()
        self.shared_element_relationships[: len(i), 0] = i
        self.shared_element_relationships[: len(i), 1] = j

        self.shared_element_relationships = self.shared_element_relationships[
            self.shared_element_relationships[:, 0] >= 0, :
        ]

        faces = element_nodes[self.shared_element_relationships[:, 0], :].multiply(
            element_nodes[self.shared_element_relationships[:, 1], :]
        )
        shared_faces = faces[np.array(np.sum(faces, axis=1) == 3).flatten(), :]
        row, col = shared_faces.nonzero()
        row = row[row.argsort()]
        col = col[row.argsort()]
        shared_face_index = np.zeros((shared_faces.shape[0], 3), dtype=int)
        shared_face_index[:] = -1
        shared_face_index[row.reshape(-1, 3)[:, 0], :] = col.reshape(-1, 3)

        self.shared_elements[np.arange(self.shared_element_relationships.shape[0]), :] = (
            shared_face_index
        )
        # resize
        self.shared_elements = self.shared_elements[: len(self.shared_element_relationships), :]
        # flag = np.zeros(self.elements.shape[0])
        # face_index = 0
        # for i, t in enumerate(self.elements):
        #     flag[i] = True
        #     for n in self.neighbours[i]:
        #         if n < 0:
        #             continue
        #         if flag[n]:
        #             continue
        #         face_node_index = 0
        #         self.shared_element_relationships[face_index, 0] = i
        #         self.shared_element_relationships[face_index, 1] = n
        #         for v in t:
        #             if v in self.elements[n, :4]:
        #                 self.shared_elements[face_index, face_node_index] = v
        #                 face_node_index += 1

        #         face_index += 1
        # self.shared_elements = self.shared_elements[:face_index, :]
        # self.shared_element_relationships = self.shared_element_relationships[
        #     :face_index, :
        # ]

    def _initialise_aabb(self):
        """Builds a sparse mapping from AABB grid cells to the tetrahedra whose
        axis-aligned bounding box overlaps that cell.

        Rather than testing every (cell, tetra) pair -- which costs
        O(n_cells * n_elements) time and memory and does not scale to large
        meshes -- each tetra's bounding box is converted directly into the
        (small) range of grid cells it spans, and only those (cell, tetra)
        pairs are recorded. Cells are sized (in __init__) at about half the
        median tetra bounding-box edge, so most tetra touch only a handful of
        cells and the resulting table stays close to O(n_elements) in size,
        rather than O(n_cells * n_elements).
        """
        bbox_min = self._tetra_bbox_min
        bbox_max = self._tetra_bbox_max

        origin = self.aabb_grid.origin
        step = self.aabb_grid.step_vector
        nsteps_cells = self.aabb_grid.nsteps_cells

        cell_min = np.floor((bbox_min - origin[None, :]) / step[None, :]).astype(np.int64)
        cell_max = np.floor((bbox_max - origin[None, :]) / step[None, :]).astype(np.int64)
        for d in range(3):
            cell_min[:, d] = np.clip(cell_min[:, d], 0, nsteps_cells[d] - 1)
            cell_max[:, d] = np.clip(cell_max[:, d], 0, nsteps_cells[d] - 1)

        # number of grid cells each tetra's bounding box spans per axis (>=1)
        span = cell_max - cell_min + 1
        counts = span[:, 0] * span[:, 1] * span[:, 2]
        total = int(counts.sum())

        tetra_id = np.repeat(np.arange(self.n_elements), counts)
        block_start = np.repeat(np.concatenate(([0], np.cumsum(counts)[:-1])), counts)
        local_idx = np.arange(total) - block_start

        ny_rep = np.repeat(span[:, 1], counts)
        nz_rep = np.repeat(span[:, 2], counts)
        k_off = local_idx % nz_rep
        j_off = (local_idx // nz_rep) % ny_rep
        i_off = local_idx // (nz_rep * ny_rep)

        i = np.repeat(cell_min[:, 0], counts) + i_off
        j = np.repeat(cell_min[:, 1], counts) + j_off
        k = np.repeat(cell_min[:, 2], counts) + k_off

        global_cell = i + nsteps_cells[0] * j + nsteps_cells[0] * nsteps_cells[1] * k

        self.aabb_table = csr_matrix(
            (np.ones(total, dtype=bool), (global_cell, tetra_id)),
            shape=(self.aabb_grid.n_elements, self.n_elements),
            dtype=bool,
        )

    @property
    def ntetra(self):
        return self.elements.shape[0]

    @property
    def n_elements(self):
        return self.ntetra

    @property
    def n_cells(self):
        return None

    @property
    def shared_element_norm(self):
        """
        Get the normal to all of the shared elements
        """
        elements = self.shared_elements
        v1 = self.nodes[elements[:, 1], :] - self.nodes[elements[:, 0], :]
        v2 = self.nodes[elements[:, 2], :] - self.nodes[elements[:, 0], :]
        return np.cross(v1, v2, axisa=1, axisb=1)

    @property
    def shared_element_size(self):
        """
        Get the area of the share triangle
        """
        norm = self.shared_element_norm
        return 0.5 * np.linalg.norm(norm, axis=1)

    @property
    def element_size(self):
        """Calculate the volume of a tetrahedron using the 4 corners
        volume = abs(det(A))/6 where A is the jacobian of the corners

        Returns
        -------
        _type_
            _description_
        """
        vecs = (
            self.nodes[self.elements[:, :4], :][:, 1:, :]
            - self.nodes[self.elements[:, :4], :][:, 0, None, :]
        )
        return np.abs(np.linalg.det(vecs)) / 6

    def evaluate_shape_derivatives(self, locations, elements=None):
        """
        Get the gradients of all tetras

        Parameters
        ----------
        elements

        Returns
        -------

        """
        inside = None
        if elements is not None:
            inside = np.zeros(self.n_elements, dtype=bool)
            inside[elements] = True
        if elements is None:
            verts, c, elements, inside = self.get_element_for_location(locations)
            # elements = np.arange(0, self.n_elements, dtype=int)
        ps = self.nodes[self.elements, :]
        m = np.array(
            [
                [
                    (ps[:, 1, 0] - ps[:, 0, 0]),
                    (ps[:, 1, 1] - ps[:, 0, 1]),
                    (ps[:, 1, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 2, 0] - ps[:, 0, 0]),
                    (ps[:, 2, 1] - ps[:, 0, 1]),
                    (ps[:, 2, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 3, 0] - ps[:, 0, 0]),
                    (ps[:, 3, 1] - ps[:, 0, 1]),
                    (ps[:, 3, 2] - ps[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I

        return element_gradients[elements, :, :], elements, inside

    def evaluate_shape(self, locations):
        """
        Convenience function returning barycentric coords

        """
        locations = np.array(locations)
        verts, c, elements, inside = self.get_element_for_location(locations)
        return c, elements, inside

    def evaluate_value(self, pos, property_array):
        """
        Evaluate value of interpolant

        Parameters
        ----------
        pos - numpy array
            locations
        prop - string
            property name

        Returns
        -------

        """
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        vertices, c, tetras, inside = self.get_element_for_location(pos)
        values[inside] = np.sum(
            c[inside, :] * property_array[self.elements[tetras[inside], :]], axis=1
        )
        return values

    def evaluate_gradient(self, pos, property_array):
        """
        Evaluate the gradient of an interpolant at the locations

        Parameters
        ----------
        pos - numpy array
            locations
        prop - string
            property to evaluate


        Returns
        -------

        """
        values = np.zeros(pos.shape)
        values[:] = np.nan
        (
            vertices,
            element_gradients,
            tetras,
            inside,
        ) = self.get_element_gradient_for_location(pos)
        # grads = np.zeros(tetras.shape)
        values[inside, :] = (
            element_gradients[inside, :, :] * property_array[self.elements[tetras][inside, None, :]]
        ).sum(2)
        # length = np.sum(values[inside, :], axis=1)
        # values[inside,:] /= length[:,None]
        return values

    def inside(self, pos):
        if pos.shape[1] > 3:
            logger.warning(f"Converting {pos.shape[1]} to 3d using first 3 columns")
            pos = pos[:, :3]

        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i] + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
        return inside

    def get_elements(self):
        return self.elements

    def get_element_for_location(self, points: np.ndarray) -> Tuple:
        """
        Determine the tetrahedron from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        points = np.asarray(points)
        npoints = points.shape[0]
        verts = np.zeros((npoints, 4, 3))
        bc = np.zeros((npoints, 4))
        tetras = np.full(npoints, -1, dtype="int64")
        inside = np.zeros(npoints, dtype=bool)
        npts_step = int(1e4)
        # process points in blocks to bound the size of the candidate table
        for start in range(0, npoints, npts_step):
            end = min(start + npts_step, npoints)
            block = points[start:end, :]

            cell_index, block_inside = self.aabb_grid.position_to_cell_index(block)
            # block_inside indicates which points fall within the aabb grid bounds
            block_inside_idx = np.flatnonzero(block_inside)
            if block_inside_idx.size == 0:
                continue
            global_index = (
                cell_index[block_inside_idx, 0]
                + self.aabb_grid.nsteps_cells[0] * cell_index[block_inside_idx, 1]
                + self.aabb_grid.nsteps_cells[0]
                * self.aabb_grid.nsteps_cells[1]
                * cell_index[block_inside_idx, 2]
            )

            tetra_indices = self.aabb_table[global_index, :].tocoo()
            # tetra_indices.row indexes into block_inside_idx (the compacted
            # set of in-bounds points), map it back to indices within block
            point_idx = block_inside_idx[tetra_indices.row]
            col = tetra_indices.col
            # using returned indexes calculate barycentric coords to determine which tetra the points are in
            vertices = self.nodes[self.elements[col, :4]]
            pos = block[point_idx, :]
            vap = pos[:, :] - vertices[:, 0, :]
            vbp = pos[:, :] - vertices[:, 1, :]
            #         # vcp = p - points[:, 2, :]
            #         # vdp = p - points[:, 3, :]
            vab = vertices[:, 1, :] - vertices[:, 0, :]
            vac = vertices[:, 2, :] - vertices[:, 0, :]
            vad = vertices[:, 3, :] - vertices[:, 0, :]
            vbc = vertices[:, 2, :] - vertices[:, 1, :]
            vbd = vertices[:, 3, :] - vertices[:, 1, :]

            va = np.einsum("ij, ij->i", vbp, _cross(vbd, vbc)) / 6.0
            vb = np.einsum("ij, ij->i", vap, _cross(vac, vad)) / 6.0
            vc = np.einsum("ij, ij->i", vap, _cross(vad, vab)) / 6.0
            vd = np.einsum("ij, ij->i", vap, _cross(vab, vac)) / 6.0
            v = np.einsum("ij, ij->i", vab, _cross(vac, vad)) / 6.0
            c = np.zeros((va.shape[0], 4))
            c[:, 0] = va / v
            c[:, 1] = vb / v
            c[:, 2] = vc / v
            c[:, 3] = vd / v
            mask = np.all(c >= 0, axis=1)

            found_point_idx = point_idx[mask]
            verts[start + found_point_idx, :, :] = vertices[mask, :, :]
            bc[start + found_point_idx, :] = c[mask, :]
            tetras[start + found_point_idx] = col[mask]
            inside[start + found_point_idx] = True

        return verts, bc, tetras, inside

    def get_element_gradients(self, elements=None):
        """
        Get the gradients of all tetras

        Parameters
        ----------
        elements

        Returns
        -------

        """
        # points = np.zeros((5, 4, self.n_cells, 3))
        # points[:, :, even_mask, :] = nodes[:, even_mask, :][self.tetra_mask_even, :, :]
        # points[:, :, ~even_mask, :] = nodes[:, ~even_mask, :][self.tetra_mask, :, :]

        # # changing order to points, tetra, nodes, coord
        # points = points.swapaxes(0, 2)
        # points = points.swapaxes(1, 2)
        if elements is None:
            elements = np.arange(0, self.n_elements, dtype=int)
        ps = self.nodes[
            self.elements, :
        ]  # points.reshape(points.shape[0] * points.shape[1], points.shape[2], points.shape[3])
        # vertices = self.nodes[self.elements[col,:]]
        m = np.array(
            [
                [
                    (ps[:, 1, 0] - ps[:, 0, 0]),
                    (ps[:, 1, 1] - ps[:, 0, 1]),
                    (ps[:, 1, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 2, 0] - ps[:, 0, 0]),
                    (ps[:, 2, 1] - ps[:, 0, 1]),
                    (ps[:, 2, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 3, 0] - ps[:, 0, 0]),
                    (ps[:, 3, 1] - ps[:, 0, 1]),
                    (ps[:, 3, 2] - ps[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I

        return element_gradients[elements, :, :]

    def get_element_gradient_for_location(self, pos):
        """
        Get the gradient of the tetra for a location

        Parameters
        ----------
        pos

        Returns
        -------

        """
        vertices, bc, tetras, inside = self.get_element_for_location(pos)
        ps = vertices
        m = np.array(
            [
                [
                    (ps[:, 1, 0] - ps[:, 0, 0]),
                    (ps[:, 1, 1] - ps[:, 0, 1]),
                    (ps[:, 1, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 2, 0] - ps[:, 0, 0]),
                    (ps[:, 2, 1] - ps[:, 0, 1]),
                    (ps[:, 2, 2] - ps[:, 0, 2]),
                ],
                [
                    (ps[:, 3, 0] - ps[:, 0, 0]),
                    (ps[:, 3, 1] - ps[:, 0, 1]),
                    (ps[:, 3, 2] - ps[:, 0, 2]),
                ],
            ]
        )
        I = np.array([[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I
        return vertices, element_gradients, tetras, inside

    def get_neighbours(self):
        """
        This function goes through all of the elements in the mesh and assembles a numpy array
        with the neighbours for each element

        Returns
        -------

        """
        return self.neighbours

    def vtk(self, node_properties=None, cell_properties=None):
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for vtk support")

        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}
        from pyvista import CellType

        celltype = np.full(self.elements.shape[0], CellType.TETRA, dtype=np.uint8)
        elements = np.hstack(
            [np.zeros(self.elements.shape[0], dtype=int)[:, None] + 4, self.elements]
        )
        elements = elements.flatten()
        grid = pv.UnstructuredGrid(elements, celltype, self.nodes)
        for key, value in node_properties.items():
            grid[key] = value
        for key, value in cell_properties.items():
            grid.cell_arrays[key] = value

        return grid
