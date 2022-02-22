"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
import logging

import numpy as np
from ._2d_base_unstructured import BaseUnstructured2d

logger = logging.getLogger(__name__)


class P1Unstructured2d(BaseUnstructured2d):
    """ """

    def __init__(self, elements, vertices, neighbours):
        BaseUnstructured2d.__init__(self, elements, vertices, neighbours)

    def evaluate_shape_derivatives(self, locations, elements=None):
        """
        compute dN/ds (1st row), dN/dt(2nd row)
        """
        locations = np.array(locations)
        if elements is None:
            c, tri = self.get_element_for_location(locations)
        else:
            tri = elements
            M = np.ones((elements.shape[0], 3, 3))
            M[:, :, 1:] = self.vertices[self.elements[elements], :][:, :3, :]
            points_ = np.ones((locations.shape[0], 3))
            points_[:, 1:] = locations
            minv = np.linalg.inv(M)
            c = np.einsum("lij,li->lj", minv, points_)

        vertices = self.nodes[self.elements[tri][:, :3]]
        jac = np.zeros((tri.shape[0], 2, 2))
        jac[:, 0, 0] = vertices[:, 1, 0] - vertices[:, 0, 0]
        jac[:, 0, 1] = vertices[:, 1, 1] - vertices[:, 0, 1]
        jac[:, 1, 0] = vertices[:, 2, 0] - vertices[:, 0, 0]
        jac[:, 1, 1] = vertices[:, 2, 1] - vertices[:, 0, 1]
        N = np.zeros((tri.shape[0], 6))

        # dN containts the derivatives of the shape functions
        dN = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])

        # find the derivatives in x and y by calculating the dot product between the jacobian^-1 and the
        # derivative matrix
        #         d_n = np.einsum('ijk,ijl->ilk',np.linalg.inv(jac),dN)
        d_n = np.linalg.inv(jac)
        #         d_n = d_n.swapaxes(1,2)
        d_n = d_n @ dN
        d_n = d_n.swapaxes(2, 1)
        # d_n = np.dot(np.linalg.inv(jac),dN)
        return d_n, tri

    def evaluate_shape(self, locations):
        locations = np.array(locations)
        c, tri = self.get_element_for_location(locations)
        # c = np.dot(np.array([1,x,y]),np.linalg.inv(M)) # convert to barycentric coordinates
        # order of bary coord is (1-s-t,s,t)
        N = c  # np.zeros((c.shape[0],3)) #evaluate shape functions at barycentric coordinates
        return N, tri
