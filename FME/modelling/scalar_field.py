import numpy as np
from ..cython.marching_tetra import marching_tetra

class ScalarField:
    """
    Generic support to represent a scalar field
    """
    def __init__(self):
        pass

    def evaluate_value(self,evaluation_points):

        pass


    def evaluate_gradient(self,evaluation_points):

        pass


class TetrahedralMeshScalarField:
    """
    Support to represent a scalar field using a tetrahedral mesh
    """
    def __init__(self, mesh, property_name):
        self.mesh = mesh
        self.property_name = property_name
        self.region = self.mesh.regions["everywhere"]
        self.interpolator = None
    @classmethod
    def from_node_values(cls, mesh, property_name, node_values):
        mesh.update_property(property_name, node_values)
        return cls(mesh, property_name)

    @classmethod
    def from_interpolator(cls, interpolator):
        interpolator.update()
        interpolator.mesh.update_property(interpolator.propertyname, interpolator.c)
        scalar_field = cls(
            interpolator.mesh,
            interpolator.propertyname)
        scalar_field.interpolator = interpolator
        return scalar_field

    def evaluate_value(self, evaluation_points):
        return self.mesh.evaluate_value(evaluation_points, self.property_name)

    def evaluate_gradient(self, evaluation_points):
        if evaluation_points.shape[0]>0:
            return self.mesh.evaluate_gradient(evaluation_points, self.property_name)
        return np.zeros((0,3))
    def mean_property_value(self):
        return np.nanmean(self.mesh.properties[self.property_name])

    def min_property_value(self):
        return np.nanmin(self.mesh.properties[self.property_name])

    def max_property_value(self):
        return np.nanmax(self.mesh.properties[self.property_name])

    def get_connected_nodes_for_mask(self, mask):
        return self.mesh.get_connected_nodes_for_mask(mask)

    def get_node_values(self):
        return self.mesh.properties[self.property_name]

    def number_of_nodes(self):
        return self.mesh.n_nodes

    def slice(self, isovalue):
        tri, ntri = marching_tetra(isovalue,
                                   self.mesh.elements,
                                   self.mesh.nodes,
                                   self.region,
                                   self.mesh.properties[self.property_name])

        ##convert from memoryview to np array
        tri = np.array(tri)
        ntri = np.array(ntri)[0]
        ##create a triangle indices array and initialise to -1
        tris = np.zeros((ntri, 3)).astype(int)
        tris[:, :] = -1
        ##create a dict for storing nodes index where the key is the node as as a tuple.
        # A dict is preferable because it is very quick to check if a key exists
        # assemble arrays for unique vertex and triangles defined by vertex indices
        nodes = {}
        n = 0  # counter
        for i in range(ntri):
            for j in range(3):
                if tuple(tri[i, j, :]) in nodes:
                    tris[i, j] = nodes[tuple(tri[i, j, :])]
                else:
                    nodes[tuple(tri[i, j, :])] = n
                    tris[i, j] = n
                    n += 1
        nodes_np = np.zeros((n, 3))
        for v in nodes.keys():
            nodes_np[nodes[v], :] = np.array(v)

        # find the normal vector to the faces using the vertex order
        a = nodes_np[tris[:, 0], :] - nodes_np[tris[:, 1], :]
        b = nodes_np[tris[:, 0], :] - nodes_np[tris[:, 2], :]

        crosses = np.cross(a, b)
        crosses = crosses / (np.sum(crosses ** 2, axis=1) ** (0.5))[:, np.newaxis]
        tribc = np.mean(nodes_np[tris, :], axis=1)
        mask = np.any(~np.isnan(tribc), axis=1)
        tribc = tribc[mask, :]

        propertygrad = self.evaluate_gradient(tribc)
        propertygrad /= np.linalg.norm(propertygrad, axis=1)[:, None]

        # dot product between gradient and normal indicates if faces are incorrectly ordered
        dotproducts = np.zeros(tris.shape[0])
        dotproducts[mask] = (propertygrad * crosses[mask]).sum(axis=1)
        # if dot product > 0 then adjust triangle indexing
        dotproducts[np.isnan(dotproducts)] = -1
        # todo need to check if nan
        indices = (dotproducts > 0)
        tris[indices] = tris[indices, ::-1]
        #
        #
        return tris, nodes_np
