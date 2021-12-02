"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
import logging

import numpy as np
from LoopStructural.interpolators.cython.dsi_helper import cg, constant_norm, fold_cg
from .base_structured_3d_support import BaseStructuredSupport

from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class UnStructuredTetMesh:
    """ """

    def __init__(self, nodes, elements, neighbours, aabb_nsteps=None):
        """An unstructured mesh defined by nodes, elements and neighbours
        An axis aligned bounding box (AABB) is used to speed up finding
        which tetra a point is in.
        The aabb grid is calculated so that there are approximately 10 tetra per
        element.

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
        self.nodes = np.array(nodes)
        self.n_nodes = self.nodes.shape[0]
        self.neighbours = np.array(neighbours, dtype=np.int64)
        self.elements = np.array(elements, dtype=np.int64)
        self.barycentre = np.sum(self.nodes[self.elements][:, :, :], axis=1) / 4.0
        self.minimum = np.min(self.nodes, axis=0)
        self.maximum = np.max(self.nodes, axis=0)
        length = self.maximum - self.minimum
        self.minimum -= length * 0.1
        self.maximum += length * 0.1
        if aabb_nsteps == None:
            box_vol = np.product(self.maximum - self.minimum)
            element_volume = box_vol / (len(self.elements) / 20)
            # calculate the step vector of a regular cube
            step_vector = np.zeros(3)
            step_vector[:] = element_volume ** (1.0 / 3.0)
            # number of steps is the length of the box / step vector
            aabb_nsteps = np.ceil((self.maximum - self.minimum) / step_vector).astype(
                int
            )
            # make sure there is at least one cell in every dimension
            aabb_nsteps[aabb_nsteps < 2] = 2
        aabb_nsteps = np.array(aabb_nsteps, dtype=int)
        step_vector = (self.maximum - self.minimum) / (aabb_nsteps - 1)
        self.aabb_grid = BaseStructuredSupport(
            self.minimum, nsteps=aabb_nsteps, step_vector=step_vector
        )
        # make a big table to store which tetra are in which element.
        # if this takes up too much memory it could be simplified by using sparse matrices or dict but
        # at the expense of speed
        self.aabb_table = np.zeros(
            (self.aabb_grid.n_elements, len(self.elements)), dtype=bool
        )
        self.aabb_table[:] = False
        self._initialise_aabb()

    def _initialise_aabb(self):
        """assigns the tetras to the grid cells where the bounding box
        of the tetra element overlaps the grid cell.
        It could be changed to use the separating axis theorem, however this would require
        significantly more calculations. (12 more I think).. #TODO test timing
        """
        # calculate the bounding box for all tetraherdon in the mesh
        # find the min/max extents for xyz
        tetra_bb = np.zeros((self.elements.shape[0], 19, 3))
        minx = np.min(self.nodes[self.elements, 0], axis=1)
        maxx = np.max(self.nodes[self.elements, 0], axis=1)
        miny = np.min(self.nodes[self.elements, 1], axis=1)
        maxy = np.max(self.nodes[self.elements, 1], axis=1)
        minz = np.min(self.nodes[self.elements, 2], axis=1)
        maxz = np.max(self.nodes[self.elements, 2], axis=1)
        ix, iy, iz = self.aabb_grid.global_index_to_cell_index(
            np.arange(self.aabb_grid.n_elements)
        )
        cix, ciy, ciz = self.aabb_grid.cell_corner_indexes(ix, iy, iz)
        px, py, pz = self.aabb_grid.node_indexes_to_position(cix, ciy, ciz)
        x_boundary = px[:, [0, 1]]
        y_boundary = py[:, [0, 2]]
        z_boundary = pz[:, [0, 3]]
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
        logic = np.logical_and(x_logic, y_logic)
        logic = np.logical_and(logic, z_logic)
        # inside_x =
        # # add the corners of the cube
        # tetra_bb[:,0,:] = np.array([minx,miny,minz]).T
        # tetra_bb[:,1,:] = np.array([maxx,miny,minz]).T
        # tetra_bb[:,2,:] = np.array([maxx,maxy,minz]).T
        # tetra_bb[:,3,:] = np.array([minx,maxy,minz]).T
        # tetra_bb[:,4,:] = np.array([minx,miny,maxz]).T
        # tetra_bb[:,5,:] = np.array([maxx,miny,minz]).T
        # tetra_bb[:,6,:] = np.array([maxx,maxy,maxz]).T
        # tetra_bb[:,7,:] = np.array([minx,maxy,maxz]).T
        # # add centre
        # tetra_bb[:,8,:] = np.array([(maxx-minx)/2,(maxy-miny)/2,(maxy-miny)/2]).T

        # # find which
        # cell_index_ijk = np.array(self.aabb_grid.position_to_cell_index(tetra_bb.reshape((-1,3)))).swapaxes(0,1)
        # cell_index_global = cell_index_ijk[:, 0] + self.aabb_grid.nsteps_cells[ None, 0] \
        #        * cell_index_ijk[:, 1] + self.aabb_grid.nsteps_cells[ None, 0] * \
        #        self.aabb_grid.nsteps_cells[ None, 1] * cell_index_ijk[:, 2]
        # bbcorners_grid_cell = cell_index_global.reshape((tetra_bb.shape[0],tetra_bb.shape[1]))
        # tetra_index = np.arange(self.elements.shape[0],dtype=int)
        # tetra_index = np.tile(tetra_index,(9,1)).T
        self.aabb_table = logic
        # [bbcorners_grid_cell,tetra_index] = True

    @property
    def ntetra(self):
        return self.elements.shape[0]

    @property
    def n_elements(self):
        return self.ntetra

    @property
    def n_cells(self):
        return None

    def barycentre(self, elements=None):
        """
        Return the barycentres of all tetrahedrons or of specified tetras using
        global index

        Parameters
        ----------
        elements - numpy array
            global index

        Returns
        -------

        """
        np.sum(self.nodes[self.elements][:, :, :], axis=1) / 4.0

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
            c[inside, :] * property_array[tetras[inside, :]], axis=1
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
            element_gradients[inside, :, :] * property_array[tetras[inside, None, :]]
        ).sum(2)
        length = np.sum(values[inside, :], axis=1)
        # values[inside,:] /= length[:,None]
        return values

    def inside(self, pos):
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= (
                pos[:, i]
                < self.origin[None, i]
                + self.step_vector[None, i] * self.nsteps_cells[None, i]
            )
        return inside

    def get_elements(self):
        return self.elements

    def get_element_for_location(self, points):
        """
        Determine the tetrahedron from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        cell_index = np.array(self.aabb_grid.position_to_cell_index(points)).swapaxes(
            0, 1
        )
        global_index = (
            cell_index[:, 0]
            + self.aabb_grid.nsteps_cells[None, 0] * cell_index[:, 1]
            + self.aabb_grid.nsteps_cells[None, 0]
            * self.aabb_grid.nsteps_cells[None, 1]
            * cell_index[:, 2]
        )
        tetra_indices = self.aabb_table[global_index, :]
        row, col = np.where(tetra_indices)
        vertices = self.nodes[self.elements[col, :]]
        pos = points[row, :]
        vap = pos[:, :] - vertices[:, 0, :]
        vbp = pos[:, :] - vertices[:, 1, :]
        #         # vcp = p - points[:, 2, :]
        #         # vdp = p - points[:, 3, :]
        vab = vertices[:, 1, :] - vertices[:, 0, :]
        vac = vertices[:, 2, :] - vertices[:, 0, :]
        vad = vertices[:, 3, :] - vertices[:, 0, :]
        vbc = vertices[:, 2, :] - vertices[:, 1, :]
        vbd = vertices[:, 3, :] - vertices[:, 1, :]

        va = np.einsum("ij, ij->i", vbp, np.cross(vbd, vbc, axisa=1, axisb=1)) / 6.0
        vb = np.einsum("ij, ij->i", vap, np.cross(vac, vad, axisa=1, axisb=1)) / 6.0
        vc = np.einsum("ij, ij->i", vap, np.cross(vad, vab, axisa=1, axisb=1)) / 6.0
        vd = np.einsum("ij, ij->i", vap, np.cross(vab, vac, axisa=1, axisb=1)) / 6.0
        v = np.einsum("ij, ij->i", vab, np.cross(vac, vad, axisa=1, axisb=1)) / 6.0
        c = np.zeros((va.shape[0], 4))
        c[:, 0] = va / v
        c[:, 1] = vb / v
        c[:, 2] = vc / v
        c[:, 3] = vd / v
        # inside = np.ones(c.shape[0],dtype=bool)
        mask = np.all(c >= 0, axis=1)
        verts = np.zeros((points.shape[0], 4, 3))
        bc = np.zeros((points.shape[0], 4))
        tetras = np.zeros(points.shape[0], dtype="int64")
        inside = np.zeros(points.shape[0], dtype=bool)

        verts[row[mask], :, :] = vertices[mask, :, :]
        bc[row[mask], :] = c[mask, :]
        tetras[row[mask]] = col[mask]
        inside[row[mask]] = True
        return verts, bc, self.elements[tetras], inside

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
        I = np.array(
            [[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]]
        )
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
        I = np.array(
            [[-1.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]]
        )
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
