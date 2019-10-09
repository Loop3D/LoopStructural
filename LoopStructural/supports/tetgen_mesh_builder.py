import meshpy.tet
import numpy as np
from sklearn.decomposition import PCA


import logging
logger = logging.getLogger(__name__)


class TetMeshBuilder:

    def __init__(self, mesh, **kwargs):
        self.mesh = mesh
        self.pca = PCA(n_components=3)
        self.pca.fit(corners)
        self.transformed_corners = self.pca.transform(corners)
        self.minpc0 = np.min(self.transformed_corners[:, 0])
        self.maxpc0 = np.max(self.transformed_corners[:, 0])
        self.minpc1 = np.min(self.transformed_corners[:, 1])
        self.maxpc1 = np.max(self.transformed_corners[:, 1])
        self.minpc2 = np.min(self.transformed_corners[:, 2])
        self.maxpc2 = np.max(self.transformed_corners[:, 2])
        self.properties = {}
        self.regions = {}
        self.cell_properties = {}
        self.n = None
        self.n_cell = None

    def build_mesh(self, **kwargs):
        self.n_tetra = 4000  # maxvol=0.001
        if 'n_tetra' in kwargs:
            self.n_tetra = kwargs['n_tetra']
        # calculate the volume of the bounding box
        self.maxvol = 0.001
        if 'maxvol' not in kwargs:
            lengthU = maxx - minx
            lengthV = maxy - miny
            lengthW = maxz - minz
            boxVol = lengthU * lengthW * lengthV
            correction_factor = 1.91
            self.maxvol = correction_factor * boxVol / n_tetra
        if 'maxvol' in kwargs:
            self.maxvol = kwargs['maxvol']
        facets = [
            [0, 1, 2, 3],
            [4, 5, 6, 7],
            [0, 4, 5, 1],
            [1, 5, 6, 2],
            [2, 6, 7, 3],
            [3, 7, 4, 0]
        ]
        # create the mesh
        info = meshpy.tet.MeshInfo()
        # use the projected points to build the mesh
        info.set_points(self.mesh.transformed_corners)
        info.set_facets(facets)
        meshpy_mesh = meshpy.tet.build(info, max_volume=self.maxvol, options=meshpy.tet.Options('pqn'))
        self.mesh.set_nodes(self.mesh.pca.inverse_transform(np.array(meshpy_mesh.points)))
        self.mesh.set_elements(np.array(meshpy_mesh.elements))
        self.mesh.set_neighbours(np.array(meshpy_mesh.neighbors))

        #set faces?
        self.support.set_barycentre(np.sum(self.nodes[self.elements][:, :, :], axis=1) / 4.)
        return
