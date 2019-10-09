from LoopStructural.supports.base_grid import BaseGrid
from scipy.spatial import cKDTree
import numpy as np

import logging
logger = logging.getLogger(__name__)


class TetMesh(BaseGrid):
    """
    An unstructured grid
    """
    def __init__(self, corners):
        super().__init__(self, corners)
        self.node_properties = {}
        self.element_properties = {}
        self.node_regions = {}
        self.element_regions = {}
        self.nodes = None
        self.n_nodes = 0
        self.elements = None
        self.n_elements = 0
        self.neighbours = None
        self.barycentre = None
        self.tree = None

    def set_nodes(self, nodes):
        """
        Set the node values of the tetrahedral mesh. This should be called from a builder
        :param nodes:
        :return:
        """
        self.nodes = nodes
        self.n_nodes = nodes.shape[0]

    def set_elements(self, elements):
        """
        set the array for the elements
        :param elements:
        :return:
        """
        self.elements = elements
        self.n_elements = elements.shape[0]

    def set_neighbours(self, neighbours):
        """
        set the neighour array
        :param neighbours:
        :return:
        """
        self.neighbours = neighbours

    def set_barycentre(self, barycentre):
        """
        set the barycentre arrays
        :param barycentre:
        :return:
        """
        self.barycentre = barycentre
        self.tree = cKDTree(self.barycentre)

    def update_property(self, propertyname, property):
        """
        updates or adds a property to the property db
        :param propertyname:
        :param property:
        :return:
        """
        self.node_properties[propertyname] = property
        grads = self.get_elements_gradients(np.arange(self.n_elements))
        props = self.node_properties[propertyname][self.elements[np.arange(self.n_elements)]]
        grad = np.einsum('ikj,ij->ik', grads, props)
        self.element_properties[propertyname] = grad

    def update_region(self, regionname, region):
        """
        updates or adds a region to the region db
        :param regionname:
        :param region:
        :return:
        """
        self.node_regions[regionname] = region
        self.node_properties['REGION_' + regionname] = region.astype(float)

    def elements_for_array(self, evaluation_points):
        """
        Get the elements for an array of points
        :param evaluation_points:
        :param k:
        :return:
        """
        # find the nearest k elements for the arrays
        # reducing k could speed up calculation but migh add errors

        d, ee = self.tree.query(evaluation_points)
        inside = self.is_inside(evaluation_points)
        return ee, inside

    def evaluate_gradient(self, evaluation_points, property_name):
        e, inside = self.elements_for_array(evaluation_points)
        ps = self.nodes[self.elements[e[inside]]]
        m = np.array(
            [[(ps[:, 1, 0] - ps[:, 0, 0]), (ps[:, 1, 1] - ps[:, 0, 1]), (ps[:, 1, 2] - ps[:, 0, 2])],
             [(ps[:, 2, 0] - ps[:, 0, 0]), (ps[:, 2, 1] - ps[:, 0, 1]), (ps[:, 2, 2] - ps[:, 0, 2])],
             [(ps[:, 3, 0] - ps[:, 0, 0]), (ps[:, 3, 1] - ps[:, 0, 1]), (ps[:, 3, 2] - ps[:, 0, 2])]])
        I = np.array(
            [[-1., 1., 0., 0.],
             [-1., 0., 1., 0.],
             [-1., 0., 0., 1.]])
        m = np.swapaxes(m, 0, 2)
        grads = np.linalg.inv(m)

        grads = grads.swapaxes(1, 2)
        grads = grads @ I
        vals = self.node_properties[property_name][self.elements[e[inside]]]
        a = np.zeros(evaluation_points.shape)
        a[inside] = (grads * vals[:, None, :]).sum(2) / 4.  # @ vals.T/
        a[~inside, :] = np.nan
        a /= np.sum(a, axis=1)[:, None]
        return a

    def evaluate_value(self, evaluation_points, property_name):
        e, inside = self.elements_for_array(evaluation_points)
        bc = self.calculate_barycentric_coordinates(e[inside], evaluation_points[inside])
        prop_int = np.zeros(e.shape)

        props = self.node_properties[property_name][self.elements[e[inside]]]
        prop_int[inside] = np.sum((bc.T * props), axis=1)
        prop_int[~inside] = np.nan
        return prop_int

    def calculate_barycentric_coordinates(self, e, p):
        """
        Calculate the barycentre coordinates for an array of n elements
        and n points
        """
        points = self.nodes[self.elements[e]]
        vap = p - points[:, 0, :]
        vbp = p - points[:, 1, :]
        vab = points[:, 1, :] - points[:, 0, :]
        vac = points[:, 2, :] - points[:, 0, :]
        vad = points[:, 3, :] - points[:, 0, :]
        vbc = points[:, 2, :] - points[:, 1, :]
        vbd = points[:, 3, :] - points[:, 1, :]

        va = np.sum(vbp * np.cross(vbd, vbc, axisa=1, axisb=1), axis=1) / 6.
        vb = np.sum(vap * np.cross(vac, vad, axisa=1, axisb=1), axis=1) / 6.
        vc = np.sum(vap * np.cross(vad, vab, axisa=1, axisb=1), axis=1) / 6.
        vd = np.sum(vap * np.cross(vab, vac, axisa=1, axisb=1), axis=1) / 6.
        v = np.sum(vab * np.cross(vac, vad, axisa=1, axisb=1), axis=1) / 6.
        c = np.zeros((4, p.shape[0]))
        c[0, :] = va / v
        c[1, :] = vb / v
        c[2, :] = vc / v
        c[3, :] = vd / v
        return c