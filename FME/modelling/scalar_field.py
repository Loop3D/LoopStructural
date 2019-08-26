import numpy as np
from ..cython.marching_tetra import marching_tetra


class StructuredGridScalarField:
    def __init__(self,grid,property_name):
        self.grid = grid
        self.property_name = property_name
        self.interpolator = None

    @classmethod
    def from_node_values(cls, grid, property_name, node_values):
        pass

    @classmethod
    def from_interpolator(cls, interpolator):
        pass

    def evaluate_value(self, evaluation_points):
        """
        Evaluate the scalar value fo the scalar field at the locations
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(np.isnan(evaluation_points),axis=1)

        if evaluation_points[~mask,:].shape[0]>0:
            evaluated[~mask] = self.grid.evaluate_value(
                evaluation_points[~mask], self.property_name)
        return evaluated

    def evaluate_gradient(self, evaluation_points):
        """
        Evaluate the gradient of the scalar field at the evaluation points
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        if evaluation_points.shape[0]>0:
            return self.grid.evaluate_gradient(evaluation_points, self.property_name)
        return np.zeros((0,3))

    def mean(self):
        """
        average value of the scalar field
        Returns
        -------

        """
        return np.nanmean(self.grid.properties[self.property_name])

    def min(self):
        """
        Min value of the scalar field
        Returns
        -------

        """
        return np.nanmin(self.grid.properties[self.property_name])

    def max(self):
        return np.nanmax(self.grid.properties[self.property_name])

    def get_node_values(self):
        """
        Node values from the mesh object
        Returns
        -------

        """
        return self.grid.properties[self.property_name]

    def number_of_nodes(self):
        """
        Number of nodes in the mesh
        Returns
        -------

        """
        return self.grid.n

    def update_property(self, values):
        """
        Updates the values of the scalar field on the mesh
        Parameters
        ----------
        values

        Returns
        -------

        """
        self.grid.properties[self.property_name] = values

    def slice(self, isovalue):
        pass
class TetrahedralMeshScalarField:
    """
    Support to represent a scalar field using a tetrahedral mesh
    """
    def __init__(self, mesh, property_name):
        """
        Constructor for an existing property
        Parameters
        ----------
        mesh - the mesh object to store this scalar field
        property_name - name of the property on the mesh
        """
        self.mesh = mesh
        self.property_name = property_name
        self.region = self.mesh.regions["everywhere"]
        self.interpolator = None

    @classmethod
    def from_node_values(cls, mesh, property_name, node_values):
        """
        Constructor for a new property to add to mesh
        Parameters
        ----------
        mesh -
        property_name
        node_values

        Returns
        -------

        """
        mesh.update_property(property_name, node_values)
        return cls(mesh, property_name)

    @classmethod
    def from_interpolator(cls, interpolator):
        """
        Constructor for scalar field from interpolator
        Parameters
        ----------
        interpolator

        Returns
        -------

        """
        interpolator.update()
        interpolator.mesh.update_property(interpolator.propertyname, interpolator.c)
        scalar_field = cls(
            interpolator.mesh,
            interpolator.propertyname)
        scalar_field.interpolator = interpolator
        return scalar_field

    def evaluate_value(self, evaluation_points):
        """
        Evaluate the scalar value fo the scalar field at the locations
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(np.isnan(evaluation_points),axis=1)

        if evaluation_points[~mask,:].shape[0]>0:
            evaluated[~mask] = self.mesh.evaluate_value(
                evaluation_points[~mask], self.property_name)
        return evaluated

    def evaluate_gradient(self, evaluation_points):
        """
        Evaluate the gradient of the scalar field at the evaluation points
        Parameters
        ----------
        evaluation_points

        Returns
        -------

        """
        if evaluation_points.shape[0]>0:
            return self.mesh.evaluate_gradient(evaluation_points, self.property_name)
        return np.zeros((0,3))

    def mean(self):
        """
        average value of the scalar field
        Returns
        -------

        """
        return np.nanmean(self.mesh.properties[self.property_name])

    def min(self):
        """
        Min value of the scalar field
        Returns
        -------

        """
        return np.nanmin(self.mesh.properties[self.property_name])

    def max(self):
        return np.nanmax(self.mesh.properties[self.property_name])

    def get_node_values(self):
        """
        Node values from the mesh object
        Returns
        -------

        """
        return self.mesh.properties[self.property_name]

    def number_of_nodes(self):
        """
        Number of nodes in the mesh
        Returns
        -------

        """
        return self.mesh.n_nodes

    def update_property(self, values):
        """
        Updates the values of the scalar field on the mesh
        Parameters
        ----------
        values

        Returns
        -------

        """
        self.mesh.properties[self.property_name] = values

    def slice(self, isovalue):
        """
        Create a triangular surface from an isovalue of the scalar field
        Parameters
        ----------
        isovalue - value for creating surface for

        Returns
        -------

        """
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
