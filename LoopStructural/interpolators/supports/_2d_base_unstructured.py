"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


class BaseUnstructured2d:
    """ """

    def __init__(self, elements, vertices, neighbours):
        self.elements = elements
        self.vertices = vertices
        if self.elements.shape[1] == 3:
            self.order = 1
        elif self.elements.shape[1] == 6:
            self.order = 2
        self.nelements = self.elements.shape[0]
        self.nx = self.vertices.shape[0]
        self.n_nodes = self.nx
        self.neighbours = neighbours

        self.properties = {}
        # build an array of edges and edge relationships
        self.edges = np.zeros((self.nelements * 3, 2), dtype=int)
        self.edge_relationships = np.zeros((self.nelements * 3, 2), dtype=int)
        edge_index = 0
        flag = np.zeros(self.elements.shape[0], dtype=bool)
        for i, t in enumerate(self.elements):
            flag[i] = True
            for n in self.neighbours[i]:
                if n < 0:
                    continue
                if flag[n]:
                    continue
                edge_node_index = 0
                self.edge_relationships[edge_index, 0] = i
                self.edge_relationships[edge_index, 1] = n
                for v in t:
                    if v in self.elements[n, :3]:
                        self.edges[edge_index, edge_node_index] = v
                        edge_node_index += 1

                edge_index += 1
        self.edges = self.edges[:edge_index, :]
        self.edge_relationships = self.edge_relationships[:edge_index, :]

    @property
    def ncps(self):
        """
        Returns the number of nodes for an element in the mesh
        """
        return self.elements.shape[1]

    @property
    def nodes(self):
        """
        Gets the nodes of the mesh as a property rather than using a function, accessible as a property! Python magic!

        Returns
        -------
        nodes : np.array((N,3))
            Fortran ordered
        """
        return self.vertices

    @property
    def barycentre(self):
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
        element_idx = np.arange(0, self.nelements)
        elements = self.elements[element_idx]
        barycentre = np.sum(self.nodes[elements][:, :3, :], axis=1) / 3.0
        return barycentre

    def element_area(self, elements):
        tri_points = self.nodes[self.elements[elements, :], :]
        M_t = np.ones((tri_points.shape[0], 3, 3))
        M_t[:, :, 1:] = tri_points[:, :3, :]
        area = np.abs(np.linalg.det(M_t)) * 0.5
        return area

    def evaluate_value(self, pos, prop):
        """
        Evaluate value of interpolant

        Parameters
        ----------
        pos - numpy array
            locations
        prop - numpy array
            property values at nodes

        Returns
        -------

        """
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        c, tri = self.evaluate_shape(pos[:, :2])
        inside = tri >= 0
        # vertices, c, elements, inside = self.get_elements_for_location(pos)
        values[inside] = np.sum(
            c[inside, :] * self.properties[prop][self.elements[tri[inside], :]], axis=1
        )
        return values

    def evaluate_gradient(self, pos, prop):
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
        element_gradients, tri = self.evaluate_shape_derivatives(pos[:, :2])
        inside = tri >= 0
        # ?vertices, element_gradients, elements, inside = self.get_element_gradient_for_location(pos[:,:2])
        # vertex_vals = self.properties[prop][elements]
        # grads = np.zeros(tetras.shape)
        # v = (element_gradients[inside,:,:]*tmesh.properties['defaultproperty'][tmesh.elements[tri[inside],:,None]]).sum(1)
        values[inside, :] = (
            element_gradients[inside, :, :]
            * self.properties[prop][self.elements[tri[inside], :, None]]
        ).sum(1)
        length = np.sum(values[inside, :], axis=1)
        # values[inside,:] /= length[:,None]
        return values

    def get_element_for_location(self, pos):
        """
        Determine the elements from a numpy array of points

        Parameters
        ----------
        pos : np.array



        Returns
        -------

        """
        M = np.ones((self.elements.shape[0], 3, 3))
        M[:, :, 1:] = self.vertices[self.elements, :][:, :3, :]
        points_ = np.ones((pos.shape[0], 3))
        points_[:, 1:] = pos
        minv = np.linalg.inv(M)
        c = np.einsum("kij,li->lkj", minv, points_)
        isinside = np.all(c >= 0, axis=2)
        ix, iy = np.where(isinside == True)
        element_idx = np.zeros(pos.shape[0], dtype=int)
        element_idx[:] = -1
        element_idx[ix] = iy
        c_return = np.zeros((pos.shape[0], 3))
        c_return[:] = np.nan
        c_return[ix, :] = c[isinside, :]
        return c_return, element_idx
