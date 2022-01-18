from ._3d_unstructured_tetra import UnStructuredTetMesh

import numpy as np

class P2UnstructuredTetMesh(UnStructuredTetMesh):

    def __init__(self, nodes, elements, neighbours, aabb_nsteps=None):
        UnStructuredTetMesh.__init__(self, nodes, elements, neighbours, aabb_nsteps)

    # def evaluate_mixed_derivative(self, indexes):
    #     """
    #     evaluate partial of N with respect to st (to set u_xy=0)
    #     """

    #     vertices = self.nodes[self.elements[indexes], :]
    #     jac = np.array(
    #         [
    #             [
    #                 (vertices[:, 1, 0] - vertices[:, 0, 0]),
    #                 (vertices[:, 1, 1] - vertices[:, 0, 1]),
    #             ],
    #             [
    #                 vertices[:, 2, 0] - vertices[:, 0, 0],
    #                 vertices[:, 2, 1] - vertices[:, 0, 1],
    #             ],
    #         ]
    #     ).T
    #     Nst_coeff = jac[:, 0, 0] * jac[:, 1, 1] + jac[:, 0, 1] * jac[:, 1, 0]

    #     Nst = self.Nst[None, :] * Nst_coeff[:, None]
    #     return (
    #         Nst
    #         + self.hN[None, 0, :] * (jac[:, 0, 0] * jac[:, 1, 0])[:, None]
    #         + self.hN[None, 1, :] * (jac[:, 1, 0] * jac[:, 1, 1])[:, None]
    #     )

    # def evaluate_shape_d2(self, indexes):
    #     """evaluate second derivatives of shape functions in s and t

    #     Parameters
    #     ----------
    #     M : [type]
    #         [description]

    #     Returns
    #     -------
    #     [type]
    #         [description]
    #     """

    #     vertices = self.nodes[self.elements[indexes], :]

    #     jac = np.array(
    #         [
    #             [
    #                 (vertices[:, 1, 0] - vertices[:, 0, 0]),
    #                 (vertices[:, 1, 1] - vertices[:, 0, 1]),
    #                 (vertices[:, 1, 2] - vertices[:, 0, 2]),
    #             ],
    #             [
    #                 (vertices[:, 2, 0] - vertices[:, 0, 0]),
    #                 (vertices[:, 2, 1] - vertices[:, 0, 1]),
    #                 (vertices[:, 2, 2] - vertices[:, 0, 2]),
    #             ],
    #             [
    #                 (vertices[:, 3, 0] - verticess[:, 0, 0]),
    #                 (vertices[:, 3, 1] - vertices[:, 0, 1]),
    #                 (vertices[:, 3, 2] - vertices[:, 0, 2]),
    #             ],
    #         ]
    #     ).T
    #     jac = np.linalg.inv(jac)
    #     jac = jac * jac

    #     d2_prod = np.einsum("lij,ik->lik", jac, self.hN)
    #     d2Const = d2_prod[:, 0, :] + d2_prod[:, 1, :]
    #     xxConst = d2_prod[:, 0, :]
    #     yyConst = d2_prod[:, 1, :]

    #     return xxConst, yyConst

    def evaluate_shape_derivatives(self, locations, elements=None):
        """
        compute dN/ds (1st row), dN/dt(2nd row)
        """
        locations = np.array(locations)
        if elements is None:
            verts, c, elements, inside = self.get_element_for_location(locations)
        else:
            M = np.ones((elements.shape[0], 3, 3))
            M[:, :, 1:] = self.vertices[self.elements[elements], :][:, :3, :]
            points_ = np.ones((locations.shape[0], 3))
            points_[:, 1:] = locations
            minv = np.linalg.inv(M)
            c = np.einsum("lij,li->lj", minv, points_)
        vertices = self.nodes[self.elements[elements][:, :3]]
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
        ).T
        N = np.zeros((elements.shape[0], 6))
 
        # dN containts the derivatives of the shape functions
        dN = np.zeros((elements.shape[0], 3, 10))
        dN[:, 0, 0] = 4*c[0] - 1
        dN[:, 0, 1] =  0
        dN[:, 0, 2] = 0
        dN[:, 0, 3] = 4*c[0] + 4*c[1] + 4*c[2] - 3
        dN[:, 0, 4] = 0
        dN[:, 0, 5] = 4*c[2]
        dN[:, 0, 6] = 4*c[1]
        dN[:, 0, 7] = -8*c[0] - 4*c[1] - 4*c[2] + 4
        dN[:, 0, 8] = -4*c[1]
        dN[:, 0, 9] = -4*c[2]

        dN[:, 1, 0] = 0
        dN[:, 1, 1] = 4*c[1] - 1
        dN[:, 1, 2] = 0
        dN[:, 1, 3] = 4*c[0] + 4*c[1] + 4*c[2] - 3
        dN[:, 1, 4] = 4*c[2]
        dN[:, 1, 5] = 0
        dN[:, 1, 6] = 4*c[0]
        dN[:, 1, 7] = -4*c[0]
        dN[:, 1, 8] = -4*c[0] - 8*c[1] - 4*c[2] + 4
        dN[:, 1, 9] = -4*c[2]

        dN[:,2,0] = 0
        dN[:,2,1] = 0
        dN[:,2,2] = 4*c[2] - 1
        dN[:,2,3] = 4*c[0] + 4*c[1] + 4*c[2] - 3
        dN[:,2,4] = 4*c[1]
        dN[:,2,5] = 4*c[0]
        dN[:,2,6] = 0
        dN[:,2,7] = -4*c[0]
        dN[:,2,8] = -4*c[1]
        dN[:,2,9] = -4*c[0] - 4*c[1] - 8*c[2] + 4

        # find the derivatives in x and y by calculating the dot product between the jacobian^-1 and the
        # derivative matrix
        #         d_n = np.einsum('ijk,ijl->ilk',np.linalg.inv(jac),dN)
        d_n = np.linalg.inv(jac)
        #         d_n = d_n.swapaxes(1,2)
        d_n = d_n @ dN
        d_n = d_n.swapaxes(2, 1)
        # d_n = np.dot(np.linalg.inv(jac),dN)
        return d_n, elements

    def evaluate_shape(self, locations):
        locations = np.array(locations)
        verts, c, elements, inside = self.get_element_for_location(locations)
        # order of bary coord is (1-s-t,s,t)
        N = np.zeros(
            (c.shape[0], 10)
        )  # evaluate shape functions at barycentric coordinates
        for i in range(c.shape[1]):
            N[:, i] = (2*c[:,i]-1)*c[:,i]
        
        N[:,4] = 4*c[:,1]*c[:,2]
        N[:,5] = 4*c[:,0]*c[:,2]
        N[:,6] = 4*c[:,0]*c[:,1]
        N[:,7] = 4*c[:,0]*c[:,3]
        N[:,8] = 4*c[:,1]*c[:,3]
        N[:,9] = 4*c[:,2]*c[:,3]
        inside = np.all(c>0,axis=1)
        return N, elements, inside

    # def evaluate_d2(self, pos, prop):
    #     """
    #     Evaluate value of interpolant

    #     Parameters
    #     ----------
    #     pos - numpy array
    #         locations
    #     prop - numpy array
    #         property values at nodes

    #     Returns
    #     -------

    #     """
    #     values = np.zeros(pos.shape[0])
    #     values[:] = np.nan
    #     c, tri = self.evaluate_shape(pos[:, :2])
    #     xxConst, yyConst = self.evaluate_shape_d2(tri)
    #     xyConst = self.evaluate_mixed_derivative(tri)
    #     inside = tri > 0
    #     # vertices, c, elements, inside = self.get_elements_for_location(pos)
    #     values[inside] = np.sum(
    #         xxConst[inside, :] * self.properties[prop][self.elements[tri[inside], :]],
    #         axis=1,
    #     )
    #     values[inside] += np.sum(
    #         yyConst[inside, :] * self.properties[prop][self.elements[tri[inside], :]],
    #         axis=1,
    #     )
    #     values[inside] += np.sum(
    #         xyConst[inside, :] * self.properties[prop][self.elements[tri[inside], :]],
    #         axis=1,
    #     )

    #     return values

    def evaluate_value(self, pos, property_array):
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
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        N, tetras, inside = self.evaluate_shape(pos)
        values[inside] = np.sum(
            N[inside, :] * property_array[self.elements[tetras][inside, :]], axis=1
        )
        return values

    def evaluate_gradient(self, pos, property_array):
        values = np.zeros(pos.shape)
        values[:] = np.nan
        element_gradients, tetra = self.evaluate_shape_derivatives(pos[:, :3])
        inside = tri >= 0
        # ?vertices, element_gradients, elements, inside = self.get_element_gradient_for_location(pos[:,:2])
        # vertex_vals = self.properties[prop][elements]
        # grads = np.zeros(tetras.shape)
        # v = (element_gradients[inside,:,:]*tmesh.properties['defaultproperty'][tmesh.elements[tri[inside],:,None]]).sum(1)
        values[inside, :] = (
            element_gradients[inside, :, :]
            * property_array[self.elements[tetra[inside], :, None]]
        ).sum(1)
        length = np.sum(values[inside, :], axis=1)
        # values[inside,:] /= length[:,None]
        return values