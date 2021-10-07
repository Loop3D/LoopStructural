"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
import logging

import numpy as np
from LoopStructural.interpolators.cython.dsi_helper import cg, constant_norm, fold_cg
from .base_structured_3d_support import BaseStructuredSupport

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class TetMesh:
    """

    """
    def __init__(self, nodes, elements, neighbours):
        self.nodes = np.array(nodes)
        self.neighbours = np.array(neighbours,dtype=np.int64)
        self.elements = np.array(elements,dtype=np.int64)
        
        self.barycentre = np.sum(self.nodes[self.elements][:, :, :],
                                 axis=1) / 4.
    @property
    def ntetra(self):
        return np.product(self.nsteps_cells) * 5
    
    @property
    def n_elements(self):
        return self.ntetra

    @property
    def n_cells(self):
        return np.product(self.nsteps_cells)

    def barycentre(self, elements = None):
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
        np.sum(self.nodes[self.elements][:, :, :],
                                 axis=1) / 4.
        # if elements is None:
        #     elements = np.arange(0,self.ntetra)
        # tetra = self.get_elements()[elements]
        # barycentre = np.sum(self.nodes[tetra][:, :, :],
        #                          axis=1) / 4.
        # return barycentre

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
        vertices, c, tetras, inside = self.get_tetra_for_location(pos)
        values[inside] = np.sum(c[inside,:]*property_array[tetras[inside,:]],axis=1)
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
        vertices, element_gradients, tetras, inside = self.get_tetra_gradient_for_location(pos)
        #grads = np.zeros(tetras.shape)
        values[inside,:] = (element_gradients[inside,:,:]*property_array[tetras[inside,None,:]]).sum(2)
        length = np.sum(values[inside,:],axis=1)
        # values[inside,:] /= length[:,None]
        return values

    def inside(self, pos):
        inside = np.ones(pos.shape[0]).astype(bool)
        for i in range(3):
            inside *= pos[:, i] > self.origin[None, i]
            inside *= pos[:, i] < self.origin[None, i] + \
                      self.step_vector[None, i] * self.nsteps_cells[None, i]
        return inside

    def get_tetra_for_location(self, pos):
        """
        Determine the tetrahedron from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        

    def get_constant_gradient(self, region, direction=None):
        """
        Get the constant gradient for the specified nodes

        Parameters
        ----------
        region : np.array(dtype=bool)
            mask of nodes to calculate cg for

        Returns
        -------

        """
        """
        Add the constant gradient regularisation to the system

        Parameters
        ----------
        w (double) - weighting of the cg parameter

        Returns
        -------

        """
        if direction is not None:
            print('using cg direction')
            logger.info("Running constant gradient")
            elements_gradients = self.get_element_gradients(np.arange(self.ntetra))
            if elements_gradients.shape[0] != direction.shape[0]:
                logger.error('Cannot add directional CG, vector field is not the correct length')
                return
            region = region.astype('int64')

            neighbours = self.get_neighbours()
            elements = self.get_elements()
            idc, c, ncons = fold_cg(elements_gradients, direction, neighbours.astype('int64'), elements.astype('int64'), self.nodes)

            idc = np.array(idc[:ncons, :])
            c = np.array(c[:ncons, :])
            B = np.zeros(c.shape[0])
            return c, idc, B
        if self.cg is None:
            logger.info("Running constant gradient")
            elements_gradients = self.get_element_gradients(np.arange(self.ntetra))
            region = region.astype('int64')

            neighbours = self.get_neighbours()
            elements = self.get_elements()
            idc, c, ncons = cg(elements_gradients, neighbours.astype('int64'), elements.astype('int64'), self.nodes,
                               region.astype('int64'))

            idc = np.array(idc[:ncons, :])
            c = np.array(c[:ncons, :])
            B = np.zeros(c.shape[0])
            self.cg = (c,idc,B)
        return self.cg[0], self.cg[1], self.cg[2]

    def get_elements(self):
        """
        Get a numpy array of all of the elements in the mesh

        Returns
        -------
        numpy array elements

        """

        

    def get_element_gradients(self, elements = None):
        """
        Get the gradients of all tetras

        Parameters
        ----------
        elements

        Returns
        -------

        """
        points = np.zeros((5, 4, self.n_cells, 3))
        points[:, :, even_mask, :] = nodes[:, even_mask, :][self.tetra_mask_even, :, :]
        points[:, :, ~even_mask, :] = nodes[:, ~even_mask, :][self.tetra_mask, :, :]

        # changing order to points, tetra, nodes, coord
        points = points.swapaxes(0, 2)
        points = points.swapaxes(1, 2)

        ps = points.reshape(points.shape[0] * points.shape[1], points.shape[2], points.shape[3])

        m = np.array(
            [[(ps[:, 1, 0] - ps[:, 0, 0]), (ps[:, 1, 1] - ps[:, 0, 1]),
              (ps[:, 1, 2] - ps[:, 0, 2])],
             [(ps[:, 2, 0] - ps[:, 0, 0]), (ps[:, 2, 1] - ps[:, 0, 1]),
              (ps[:, 2, 2] - ps[:, 0, 2])],
             [(ps[:, 3, 0] - ps[:, 0, 0]), (ps[:, 3, 1] - ps[:, 0, 1]),
              (ps[:, 3, 2] - ps[:, 0, 2])]])
        I = np.array(
            [[-1., 1., 0., 0.],
             [-1., 0., 1., 0.],
             [-1., 0., 0., 1.]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I

        return element_gradients[elements,:,:]

    def get_tetra_gradient_for_location(self, pos):
        """
        Get the gradient of the tetra for a location

        Parameters
        ----------
        pos

        Returns
        -------

        """
        vertices, bc, tetras, inside = self.get_tetra_for_location(pos)
        ps = vertices
        m = np.array(
            [[(ps[:, 1, 0] - ps[:, 0, 0]), (ps[:, 1, 1] - ps[:, 0, 1]),(ps[:, 1, 2] - ps[:, 0, 2])],
             [(ps[:, 2, 0] - ps[:, 0, 0]), (ps[:, 2, 1] - ps[:, 0, 1]),(ps[:, 2, 2] - ps[:, 0, 2])],
             [(ps[:, 3, 0] - ps[:, 0, 0]), (ps[:, 3, 1] - ps[:, 0, 1]),(ps[:, 3, 2] - ps[:, 0, 2])]])
        I = np.array(
            [[-1., 1., 0., 0.],
             [-1., 0., 1., 0.],
             [-1., 0., 0., 1.]])
        m = np.swapaxes(m, 0, 2)
        element_gradients = np.linalg.inv(m)

        element_gradients = element_gradients.swapaxes(1, 2)
        element_gradients = element_gradients @ I
        return vertices, element_gradients, tetras, inside



    def global_node_indicies(self, indexes):
        """
        Convert from node indexes to global node index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        indexes = np.array(indexes).swapaxes(0, 2)
        return indexes[:, :, 0] + self.nsteps[None, None, 0] \
               * indexes[:, :, 1] + self.nsteps[None, None, 0] * \
               self.nsteps[None, None, 1] * indexes[:, :, 2]

    def global_cell_indicies(self, indexes):
        """
        Convert from cell indexes to global cell index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        indexes = np.array(indexes).swapaxes(0, 2)
        return indexes[:, :, 0] + self.nsteps_cells[None, None, 0] \
               * indexes[:, :, 1] + self.nsteps_cells[None, None, 0] * \
               self.nsteps_cells[None, None, 1] * indexes[:, :, 2]

    def cell_corner_indexes(self, x_cell_index, y_cell_index, z_cell_index):
        """
        Returns the indexes of the corners of a cell given its location xi,
        yi, zi

        Parameters
        ----------
        x_cell_index
        y_cell_index
        z_cell_index

        Returns
        -------

        """
        x_cell_index = np.array(x_cell_index)
        y_cell_index = np.array(y_cell_index)
        z_cell_index = np.array(z_cell_index)

        xcorner = np.array([0, 1, 0, 1, 0, 1, 0, 1])
        ycorner = np.array([0, 0, 1, 1, 0, 0, 1, 1])
        zcorner = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        xcorners = x_cell_index[:, None] + xcorner[None, :]
        ycorners = y_cell_index[:, None] + ycorner[None, :]
        zcorners = z_cell_index[:, None] + zcorner[None, :]
        return xcorners, ycorners, zcorners

    def position_to_cell_corners(self, pos):
        """
        Find the nodes that belong to a cell which contains a point

        Parameters
        ----------
        pos

        Returns
        -------

        """
        inside = self.inside(pos)
        ix, iy, iz = self.position_to_cell_index(pos)
        cornersx, cornersy, cornersz = self.cell_corner_indexes(ix, iy, iz)
        globalidx = self.global_cell_indicies(
            np.dstack([cornersx, cornersy, cornersz]).T)
        return globalidx, inside


    def node_indexes_to_position(self, xindex, yindex, zindex):
        """
        Get the xyz position from the node coordinates

        Parameters
        ----------
        xindex
        yindex
        zindex

        Returns
        -------

        """
        x = self.origin[0] + self.step_vector[0] * xindex
        y = self.origin[1] + self.step_vector[1] * yindex
        z = self.origin[2] + self.step_vector[2] * zindex

        return np.array([x, y, z])

    def global_index_to_node_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        # determine the ijk indices for the global index.
        # remainder when dividing by nx = i
        # remained when dividing modulus of nx by ny is j
        x_index = global_index % self.nsteps[0, None]
        y_index = global_index // self.nsteps[0, None] % \
                  self.nsteps[1, None]
        z_index = global_index // self.nsteps[0, None] // \
                  self.nsteps[1, None]
        return x_index, y_index, z_index

    def global_index_to_cell_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        # determine the ijk indices for the global index.
        # remainder when dividing by nx = i
        # remained when dividing modulus of nx by ny is j

        x_index = global_index % self.nsteps_cells[0, None]
        y_index = global_index // self.nsteps_cells[0, None] % \
                  self.nsteps_cells[1, None]
        z_index = global_index // self.nsteps_cells[0, None] // \
                  self.nsteps_cells[1, None]
        return x_index, y_index, z_index

    def get_neighbours(self):
        """
        This function goes through all of the elements in the mesh and assembles a numpy array
        with the neighbours for each element

        Returns
        -------

        """
        



