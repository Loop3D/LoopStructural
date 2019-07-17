import meshpy.tet
import numpy as np

class TetMeshBuilder:

    def __init__(self,support,**kwargs):
        self.support = support
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

    def build_mesh(self):

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
        info.set_points(self.support.transformed_corners)
        info.set_facets(facets)
        meshpy_mesh = meshpy.tet.build(info, max_volume=self.maxvol, options=meshpy.tet.Options('pqn'))
        self.support.set_nodes(self.support.pca.inverse_transform(np.array(meshpy_mesh.points)))
        self.support.set_elements(np.array(meshpy_mesh.elements))
        self.support.set_neighbours(np.array(meshpy_mesh.neighbors))

        #set faces?
        self.support.set_barycentre(np.sum(self.nodes[self.elements][:, :, :], axis=1) / 4.)

        return
