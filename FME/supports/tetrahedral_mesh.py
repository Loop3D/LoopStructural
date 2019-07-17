from .base_support import BaseGrid
from scipy.spatial import cKDTree


class TetMeshSupport(BaseGrid):

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
        self.nodes = nodes
        self.n_nodes = nodes.shape[0]

    def set_elements(self,elements):
        self.elements = elements
        self.n_elements = elements.shape[0]

    def set_neighbours(self,neighbours):

        self.neighbours = neighbours

    def set_barycentre(self,barycentre):

        self.barycentre = barycentre
        self.tree =cKDTree(self.barycentre)

    def update_property(self,propertyname,property):

        pass

    def update_region(self,regionname,region):

        pass

