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
    def __init__(self, nodes, elements, neighbours,aabb_nsteps=[10,10,10]):
        self.nodes = np.array(nodes)
        self.neighbours = np.array(neighbours,dtype=np.int64)
        self.elements = np.array(elements,dtype=np.int64)
        self.barycentre = np.sum(self.nodes[self.elements][:, :, :],
                                 axis=1) / 4.
        self.minimum = np.min(self.nodes,axis=0)
        self.maximum = np.max(self.nodes,axis=0)
        length = (self.maximum-self.minimum)
        self.minimum -= length*0.2
        self.maximum += length*0.2

        aabb_nsteps = np.array(aabb_nsteps,dtype=int)
        step_vector = (self.maximum-self.minimum)/aabb_nsteps
        self.aabb_grid = BaseStructuredSupport(self.minimum,nsteps=aabb_nsteps,step_vector=step_vector)
        # make a big table to store which tetra are in which element.
        # if this takes up too much memory it could be simplified by using sparse matrices or dict but 
        # at the expense of speed 
        self.aabb_table = np.zeros((self.aabb_grid.n_elements,len(self.elements)),dtype=bool)
        self.aabb_table[:] = False
    
    def initialise_aabb(self):
        # calculate the bounding box for all tetraherdon in the mesh
        # find the min/max extents for xyz
        tetra_bb = np.zeros((self.elements.shape[0],8,3))
        minx = np.min(self.nodes[self.elements,0],axis=1)
        maxx = np.max(self.nodes[self.elements,0],axis=1)
        miny = np.min(self.nodes[self.elements,1],axis=1)
        maxy = np.max(self.nodes[self.elements,1],axis=1)
        minz = np.min(self.nodes[self.elements,2],axis=1)
        maxz = np.max(self.nodes[self.elements,2],axis=1)
        
        # add the corners of the cube
        tetra_bb[:,0,:] = np.array([minx,miny,minz]).T
        tetra_bb[:,1,:] = np.array([maxx,miny,minz]).T
        tetra_bb[:,2,:] = np.array([maxx,maxy,minz]).T
        tetra_bb[:,3,:] = np.array([minx,maxy,minz]).T
        tetra_bb[:,4,:] = np.array([minx,miny,maxz]).T
        tetra_bb[:,5,:] = np.array([maxx,miny,minz]).T
        tetra_bb[:,6,:] = np.array([maxx,maxy,maxz]).T
        tetra_bb[:,7,:] = np.array([minx,maxy,maxz]).T

        # find which 
        cell_index_ijk = np.array(self.aabb_grid.position_to_cell_index(tetra_bb.reshape((-1,3)))).swapaxes(0,1)
        cell_index_global = cell_index_ijk[:, 0] + self.aabb_grid.nsteps_cells[ None, 0] \
               * cell_index_ijk[:, 1] + self.aabb_grid.nsteps_cells[ None, 0] * \
               self.aabb_grid.nsteps_cells[ None, 1] * cell_index_ijk[:, 2]
        bbcorners_grid_cell = cell_index_global.reshape((tetra_bb.shape[0],tetra_bb.shape[1]))
        tetra_index = np.arange(self.elements.shape[0],dtype=int)
        tetra_index = np.tile(tetra_index,(8,1)).T
        self.aabb_table[bbcorners_grid_cell,tetra_index] = True

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

    def get_tetra_for_location(self, points):
        """
        Determine the tetrahedron from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        cell_index = np.array(grid.position_to_cell_index(points)).swapaxes(0,1)
        global_index = cell_index[:, 0] + grid.nsteps_cells[ None, 0] \
                    * cell_index[:, 1] + grid.nsteps_cells[ None, 0] * \
                    grid.nsteps_cells[ None, 1] * cell_index[:, 2]
        tetra_indices = grid_to_tetra[global_index,:]
        row, col = np.where(tetra_indices)
        vertices = nodes[elements[col,:]]
        pos = points[row,:]
        vap = pos[:, :] - vertices[:, 0, :]
        vbp = pos[:, :] - vertices[:, 1, :]
        #         # vcp = p - points[:, 2, :]
        #         # vdp = p - points[:, 3, :]
        vab = vertices[:, 1, :] - vertices[:, 0, :]
        vac = vertices[:, 2, :] - vertices[:, 0, :]
        vad = vertices[:, 3, :] - vertices[:, 0, :]
        vbc = vertices[:, 2, :] - vertices[:, 1, :]
        vbd = vertices[:, 3, :] - vertices[:, 1, :]

        va = np.einsum('ij, ij->i', vbp, np.cross(vbd, vbc, axisa=1, axisb=1)) / 6.
        vb = np.einsum('ij, ij->i', vap, np.cross(vac, vad, axisa=1, axisb=1)) / 6.
        vc = np.einsum('ij, ij->i', vap, np.cross(vad, vab, axisa=1, axisb=1)) / 6.
        vd = np.einsum('ij, ij->i', vap, np.cross(vab, vac, axisa=1, axisb=1)) / 6.
        v = np.einsum('ij, ij->i', vab, np.cross(vac, vad, axisa=1, axisb=1)) / 6.
        c = np.zeros((va.shape[0], 4))
        c[:,  0] = va / v
        c[:,  1] = vb / v
        c[:,  2] = vc / v
        c[:,  3] = vd / v

        mask = np.all(c>0,axis=1)
        return vertices[mask,:,:], c[mask], col[mask], np.ones(col[inside],dtype=bool)

        

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
        



