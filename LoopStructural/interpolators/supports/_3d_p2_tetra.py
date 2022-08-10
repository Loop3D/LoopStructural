from ._3d_unstructured_tetra import UnStructuredTetMesh

import numpy as np


class P2UnstructuredTetMesh(UnStructuredTetMesh):
    def __init__(
        self,
        nodes: np.ndarray,
        elements: np.ndarray,
        neighbours: np.ndarray,
        aabb_nsteps=None,
    ):
        UnStructuredTetMesh.__init__(self, nodes, elements, neighbours, aabb_nsteps)
        if self.elements.shape[1] != 10:
            raise ValueError(
                f"P2 tetrahedron must have 8 nodes, has {self.elements.shape[1]}"
            )
        self.hessian = np.array(
            [
                [
                    [4, 4, 0, 0, 0, 0, -8, 0, 0, 0],
                    [4, 0, 0, 0, 0, 0, -4, 4, 0, -4],
                    [4, 0, 0, 0, 0, -4, -4, 0, 4, 0],
                ],
                [
                    [4, 0, 0, 0, 0, 0, -4, 4, 0, -4],
                    [4, 0, 4, 0, 0, 0, 0, 0, 0, -8],
                    [4, 0, 0, 0, 4, -4, 0, 0, 0, -4],
                ],
                [
                    [4, 0, 0, 0, 0, -4, -4, 0, 4, 0],
                    [4, 0, 0, 0, 4, -4, 0, 0, 0, -4],
                    [4, 0, 0, 4, 0, -8, 0, 0, 0, 0],
                ],
            ]
        )

    def get_quadrature_points(self, npts: int = 3):
        """Calculate the quadrature points for the triangle using 3 points
        these points are at the barycentric coordinates of (1/6,1/6), (1/6,2/3), (2/3,1/6)
        All points are weighted equally at 1/6

        Parameters
        ----------
        npts : int, optional
            _description_, by default 3

        Returns
        -------
        _type_
            _description_
        """
        if npts == 3:
            vertices = self.nodes[self.shared_elements]
            cp = np.zeros((vertices.shape[0], 3, 3))
            reference_points = np.array([[1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]])

            cp[:, 0, :] = (
                vertices[:, 0, :]
                * (1 - reference_points[0, 0] - reference_points[1, 0])
                + vertices[:, 1, :] * (reference_points[0, 0])
                + vertices[:, 2, :] * (reference_points[1, 0])
            )
            cp[:, 1, :] = (
                vertices[:, 0, :]
                * (1 - reference_points[0, 1] - reference_points[1, 1])
                + vertices[:, 1, :] * (reference_points[0, 1])
                + vertices[:, 2, :] * (reference_points[1, 1])
            )
            cp[:, 2, :] = (
                vertices[:, 0, :]
                * (1 - reference_points[0, 2] - reference_points[1, 2])
                + vertices[:, 1, :] * (reference_points[0, 2])
                + vertices[:, 2, :] * (reference_points[1, 2])
            )
            weights = np.zeros((vertices.shape[0], 3))
            weights[:, :] = 1 / 6
            return cp, weights
        if npts == 1:
            vertices = self.nodes[self.shared_elements]
            cp = np.zeros((vertices.shape[0], 1, 3))
            reference_points = np.array([[1 / 3], [1 / 3]])

            cp[:, 0, :] = (
                vertices[:, 0, :]
                * (1 - reference_points[0, 0] - reference_points[1, 0])
                + vertices[:, 1, :] * (reference_points[0, 0])
                + vertices[:, 2, :] * (reference_points[1, 0])
            )
            weights = np.zeros((vertices.shape[0], 1))
            weights[:, :] = 1 / 2
            return cp, weights

    def evaluate_shape_d2(self, indexes: np.ndarray) -> np.ndarray:
        """evaluate second derivatives of shape functions in s and t

        Parameters
        ----------
        indexes : np.ndarray
            array of indexes

        Returns
        -------
        np.array
            array of second derivative shape function
        """

        vertices = self.nodes[self.elements[indexes], :]

        jac = np.array(
            [
                [
                    (vertices[:, 1, 0] - vertices[:, 0, 0]),
                    (vertices[:, 1, 1] - vertices[:, 0, 1]),
                    (vertices[:, 1, 2] - vertices[:, 0, 2]),
                ],
                [
                    (vertices[:, 2, 0] - vertices[:, 0, 0]),
                    (vertices[:, 2, 1] - vertices[:, 0, 1]),
                    (vertices[:, 2, 2] - vertices[:, 0, 2]),
                ],
                [
                    (vertices[:, 3, 0] - vertices[:, 0, 0]),
                    (vertices[:, 3, 1] - vertices[:, 0, 1]),
                    (vertices[:, 3, 2] - vertices[:, 0, 2]),
                ],
            ]
        )
        jac = jac.swapaxes(0, 2)
        jac = jac.swapaxes(1, 2)
        jac = np.linalg.inv(jac)
        # calculate derivative by summation
        d2 = np.zeros((vertices.shape[0], 6, self.elements.shape[1]))
        ii = 0
        for i in range(3):
            for j in range(i, 3):
                for k in range(3):
                    for l in range(3):
                        d2[:, ii, :] += (
                            jac[:, i, k, None]
                            * jac[:, j, l, None]
                            * self.hessian[None, k, l, :]
                        )
                ii += 1
        return d2

    def evaluate_shape_derivatives(
        self, locations: np.ndarray, elements: np.ndarray = None
    ) -> np.ndarray:
        """
        compute dN/ds (1st row), dN/dt(2nd row)


        Parameters
        ----------
        locations : np.array
            location (n,3) array
        elements : np.array, optional
            indexes to calculate shape function for. Used when evaluating quad points
            on faces as two tetra hold the point by default None
            When it is none, the index is calculated from the location

        Returns
        -------
        np.array
            array of shape paramters
        """

        locations = np.array(locations)
        if elements is None:
            verts, c, elements, inside = self.get_element_for_location(locations)
        else:
            M = np.ones((elements.shape[0], 4, 4))
            M[:, :, 1:] = self.nodes[self.elements[elements], :][:, :4, :]
            points_ = np.ones((locations.shape[0], 4))
            points_[:, 1:] = locations
            minv = np.linalg.inv(M)
            c = np.einsum("lij,li->lj", minv, points_)
            verts = self.nodes[self.elements[elements][:, :4]]
        jac = np.array(
            [
                [
                    (verts[:, 1, 0] - verts[:, 0, 0]),
                    (verts[:, 1, 1] - verts[:, 0, 1]),
                    (verts[:, 1, 2] - verts[:, 0, 2]),
                ],
                [
                    (verts[:, 2, 0] - verts[:, 0, 0]),
                    (verts[:, 2, 1] - verts[:, 0, 1]),
                    (verts[:, 2, 2] - verts[:, 0, 2]),
                ],
                [
                    (verts[:, 3, 0] - verts[:, 0, 0]),
                    (verts[:, 3, 1] - verts[:, 0, 1]),
                    (verts[:, 3, 2] - verts[:, 0, 2]),
                ],
            ]
        )
        r = c[:, 1]
        s = c[:, 2]
        t = c[:, 3]
        jac = np.swapaxes(jac, 0, 2)
        # dN containts the derivatives of the shape functions
        dN = np.zeros((elements.shape[0], 3, 10))

        dN[:, 0, 0] = 4 * r + 4 * s + 4 * t - 3
        dN[:, 0, 1] = 4 * r - 1
        dN[:, 0, 2] = 0
        dN[:, 0, 3] = 0
        dN[:, 0, 4] = 0
        dN[:, 0, 5] = -4 * t
        dN[:, 0, 6] = -8 * r - 4 * s - 4 * t + 4
        dN[:, 0, 7] = 4 * s
        dN[:, 0, 8] = 4 * t
        dN[:, 0, 9] = -4 * s

        dN[:, 1, 0] = 4 * r + 4 * s + 4 * t - 3
        dN[:, 1, 1] = 0
        dN[:, 1, 2] = 4 * s - 1
        dN[:, 1, 3] = 0
        dN[:, 1, 4] = 4 * t
        dN[:, 1, 5] = -4 * t
        dN[:, 1, 6] = -4 * r
        dN[:, 1, 7] = 4 * r
        dN[:, 1, 8] = -0
        dN[:, 1, 9] = -4 * r - 8 * s - 4 * t + 4

        dN[:, 2, 0] = 4 * r + 4 * s + 4 * t - 3
        dN[:, 2, 1] = 0
        dN[:, 2, 2] = 0
        dN[:, 2, 3] = 4 * t - 1
        dN[:, 2, 4] = 4 * s
        dN[:, 2, 5] = -4 * r - 4 * s - 8 * t + 4
        dN[:, 2, 6] = -4 * r
        dN[:, 2, 7] = 0
        dN[:, 2, 8] = 4 * r
        dN[:, 2, 9] = -4 * s

        # find the derivatives in x and y by calculating the dot product between the jacobian^-1 and the
        # derivative matrix
        #         d_n = np.einsum('ijk,ijl->ilk',np.linalg.inv(jac),dN)
        d_n = np.linalg.inv(jac)
        #         d_n = d_n.swapaxes(1,2)
        d_n = d_n.swapaxes(1, 2)
        d_n = d_n @ dN
        # d_n = np.dot(np.linalg.inv(jac),dN)
        return d_n, elements

    def evaluate_shape(self, locations: np.ndarray):
        locations = np.array(locations)
        verts, c, elements, inside = self.get_element_for_location(locations)
        # order of bary coord is (1-s-t,s,t)
        N = np.zeros(
            (c.shape[0], 10)
        )  # evaluate shape functions at barycentric coordinates

        for i in range(c.shape[1]):
            N[:, i] = (2 * c[:, i] - 1) * c[:, i]

        N[:, 4] = 4 * c[:, 3] * c[:, 2]
        N[:, 5] = 4 * c[:, 0] * c[:, 3]
        N[:, 6] = 4 * c[:, 0] * c[:, 1]
        N[:, 7] = 4 * c[:, 1] * c[:, 2]
        N[:, 8] = 4 * c[:, 1] * c[:, 3]
        N[:, 9] = 4 * c[:, 0] * c[:, 2]
        # inside = np.all(c>0,axis=1)
        return N, elements, inside

    def evaluate_d2(self, pos: np.ndarray, prop: np.ndarray) -> np.ndarray:
        """
        Evaluate the second derivative of the interpolant
        d2x, dxdy, d2y, dxdz dydz d2dz

        Parameters
        ----------
        pos - numpy array
            locations
        prop - numpy array
            property values at nodes

        Returns
        -------

        """
        c, tri, inside = self.evaluate_shape(pos)
        d2 = self.evaluate_shape_d2(tri)
        values = np.zeros((pos.shape[0], d2.shape[1]))
        values[:] = np.nan

        for i in range(d2.shape[1]):
            values[inside, i] = np.sum(
                d2[inside, i, :] * prop[self.elements[tri[inside], :]],
                axis=1,
            )

        return values

    def evaluate_value(self, pos: np.ndarray, property_array: np.ndarray) -> np.ndarray:
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
        if len(pos.shape) != 2 or pos.shape[1] != 3:
            raise ValueError(
                f"pos must be a numpy array of shape (n,3), shape is {pos.shape}"
            )
        if property_array.shape[0] != self.n_nodes:
            raise ValueError("property array must have same length as nodes")
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        N, tetras, inside = self.evaluate_shape(pos)
        values[inside] = np.sum(
            N[inside, :] * property_array[self.elements[tetras][inside, :]], axis=1
        )
        return values

    def evaluate_gradient(
        self, pos: np.ndarray, property_array: np.ndarray
    ) -> np.ndarray:
        if len(pos.shape) != 2 or pos.shape[1] != 3:
            raise ValueError(
                f"pos must be a numpy array of shape (n,3), shape is {pos.shape}"
            )
        if property_array.shape[0] != self.n_nodes:
            raise ValueError("property array must have same length as nodes")
        values = np.zeros(pos.shape)
        values[:] = np.nan
        element_gradients, tetra = self.evaluate_shape_derivatives(pos[:, :3])
        inside = tetra >= 0
        values[inside, :] = (
            element_gradients[:, :, :]
            * property_array[self.elements[tetra[inside], None, :]]
        ).sum(2)

        return values
