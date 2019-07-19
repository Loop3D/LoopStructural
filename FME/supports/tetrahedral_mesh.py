from .base_support import BaseGrid
from scipy.spatial import cKDTree


class TetMeshSupport(BaseGrid):
    """
    An unstructured grid
    """
    def __init__(self,corners):
        super().__init__(self,corners)
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

    def set_nodes(self,nodes):
        """
        Set the node values of the tetrahedral mesh. This should be called from a builder
        :param nodes:
        :return:
        """
        self.nodes = nodes
        self.n_nodes = nodes.shape[0]

    def set_elements(self,elements):
        """
        set the array for the elements
        :param elements:
        :return:
        """
        self.elements = elements
        self.n_elements = elements.shape[0]

    def set_neighbours(self,neighbours):
        """
        set the neighour array
        :param neighbours:
        :return:
        """
        self.neighbours = neighbours

    def set_barycentre(self,barycentre):
        """
        set the barycentre arrays
        :param barycentre:
        :return:
        """
        self.barycentre = barycentre
        self.tree =cKDTree(self.barycentre)

    def update_property(self,propertyname,property):
        """
        updates or adds a property to the property db
        :param propertyname:
        :param property:
        :return:
        """
        pass

    def update_region(self, regionname, region):
        """
        updates or adds a region to the region db
        :param regionname:
        :param region:
        :return:
        """
        pass

    def evaluate_gradient(self,evaluation_points):

        pass

    def evaluate_value(self,evaluation_points):

        pass