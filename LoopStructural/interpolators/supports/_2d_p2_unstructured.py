"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
import logging

import numpy as np
from ._2d_base_unstructured import BaseUnstructured2d

logger = logging.getLogger(__name__)


class P2Unstructured2d(BaseUnstructured2d):
    """ """

    def __init__(self, elements, vertices, neighbours):
        BaseUnstructured2d.__init__(self, elements, vertices, neighbours)
        # hessian of shape functions
        self.H = np.array(
            [
                [[4, 4, 0, 0, 0, -8], [4, 0, 0, 4, -4, -4]],
                [[4, 0, 0, 4, -4, -4], [4, 0, 4, 0, -8, 0]],
            ]
        )

    def evaluate_d2_shape(self, indexes):
        vertices = self.nodes[self.elements[indexes], :]
        jac = np.array(
            [
                [
                    (vertices[:, 1, 0] - vertices[:, 0, 0]),
                    (vertices[:, 1, 1] - vertices[:, 0, 1]),
                ],
                [
                    vertices[:, 2, 0] - vertices[:, 0, 0],
                    vertices[:, 2, 1] - vertices[:, 0, 1],
                ],
            ]
        )
        jac = np.linalg.inv(jac)
        dxy = (
            self.H[None, 0, 1, :] * jac[:, 0, 0] * jac[:, 1, 1]
            + self.H[None, 0, 1, :] * jac[:, 1, 0] * jac[:, 0, 1]
            + self.H[None, 0, 0, :] * jac[:, 0, 0] * jac[:, 0, 1]
            + self.H[None, 1, 1, :] * jac[:, 1, 0] * jac[:, 1, 1]
        )
        dxx = (
            self.H[None, 0, 0, :] * jac[:, 0, 0] * jac[:, 0, 0]
            + jac[:, 0, 0] * jac[:, 1, 0] * self.H[None, 0, 1, :]
            + jac[:, 1, 0] * jac[:, 1, 0] * self.H[None, 1, 1]
        )
        dyy = (
            self.H[None, 0, 0, :] * jac[:, 1, 0] * jac[:, 1, 0]
            + jac[:, 1, 0] * jac[:, 1, 1] * self.H[None, 0, 1, :]
            + jac[:, 1, 1] * jac[:, 1, 1] * self.H[None, 1, 1]
        )
        return dxx, dyy, dxy

    #     vertices = np.zeros((3,2))
    #     vertices[0,:] = [M[0,1],M[0,2]]
    #     vertices[1,:] = [M[1,1],M[1,2]]
    #     vertices[2,:] = [M[2,1],M[2,2]]
    #     jac  = np.array([[(vertices[1,0]-vertices[0,0]),(vertices[1,1]-vertices[0,1])],
    #         [vertices[2,0]-vertices[0,0],vertices[2,1]-vertices[0,1]]])
    #     Nst_coeff = jac[0,0]*jac[1,1]+jac[0,1]*jac[1,0]

    #     #N_st
    #     Nst = np.zeros(6)
    #     Nst[0] = 4
    #     Nst[1] = 0
    #     Nst[2] = 0
    #     Nst[3] = 4
    #     Nst[4] = -4
    #     Nst[5] = -4

    #     hN = np.zeros((2,6))

    #     #N_ss
    #     hN[0,0] = 4
    #     hN[0,1] = 4
    #     hN[0,2] = 0
    #     hN[0,3] = 0
    #     hN[0,4] = 0
    #     hN[0,5] = -8

    #     #N_tt
    #     hN[1,0] = 4
    #     hN[1,1] = 0
    #     hN[1,2] = 4
    #     hN[1,3] = 0
    #     hN[1,4] = -8
    #     hN[1,5] = 0

    #     xyConst = Nst*Nst_coeff + hN[0] * jac[0,0]*jac[1,0] + hN[1] * jac[1,0]*jac[1,1]
    #     jac = np.linalg.inv(jac)
    #     jac = jac*jac

    #     d2_prod = np.dot(jac,hN)
    #     d2Const = d2_prod[0] + d2_prod[1]
    #     xxConst = d2_prod[0]
    #     yyConst = d2_prod[1]

    #     return xxConst,yyConst,xyConstz
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
    #             ],
    #             [
    #                 vertices[:, 2, 0] - vertices[:, 0, 0],
    #                 vertices[:, 2, 1] - vertices[:, 0, 1],
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
        dN = np.zeros((tri.shape[0], 2, 6))
        dN[:, 0, 0] = 4 * c[:, 1] + 4 * c[:, 2] - 3  # diff(N1,s).evalf(subs=vmap)
        dN[:, 0, 1] = 4 * c[:, 1] - 1  # diff(N2,s).evalf(subs=vmap)
        dN[:, 0, 2] = 0  # diff(N3,s).evalf(subs=vmap)
        dN[:, 0, 3] = 4 * c[:, 2]  # diff(N4,s).evalf(subs=vmap)
        dN[:, 0, 4] = -4 * c[:, 2]  # diff(N5,s).evalf(subs=vmap)
        dN[:, 0, 5] = -8 * c[:, 1] - 4 * c[:, 2] + 4  # diff(N6,s).evalf(subs=vmap)

        dN[:, 1, 0] = 4 * c[:, 1] + 4 * c[:, 2] - 3  # diff(N1,t).evalf(subs=vmap)
        dN[:, 1, 1] = 0  # diff(N2,t).evalf(subs=vmap)
        dN[:, 1, 2] = 4 * c[:, 2] - 1  # diff(N3,t).evalf(subs=vmap)
        dN[:, 1, 3] = 4 * c[:, 1]  # diff(N4,t).evalf(subs=vmap)
        dN[:, 1, 4] = -4 * c[:, 1] - 8 * c[:, 2] + 4  # diff(N5,t).evalf(subs=vmap)
        dN[:, 1, 5] = -4 * c[:, 1]  # diff(N6,t).evalf(subs=vmap)

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
        N = np.zeros(
            (c.shape[0], 6)
        )  # evaluate shape functions at barycentric coordinates
        N[:, 0] = c[:, 0] * (2 * c[:, 0] - 1)  # (1-s-t)(1-2s-2t)
        N[:, 1] = c[:, 1] * (2 * c[:, 1] - 1)  # s(2s-1)
        N[:, 2] = c[:, 2] * (2 * c[:, 2] - 1)  # t(2t-1)
        N[:, 3] = 4 * c[:, 1] * c[:, 2]  # 4st
        N[:, 4] = 4 * c[:, 2] * c[:, 0]  # 4t(1-s-t)
        N[:, 5] = 4 * c[:, 1] * c[:, 0]  # 4s(1-s-t)

        return N, tri

    def evaluate_d2(self, pos, prop):
        """
        Evaluate value of interpolant

        Parameters
        ----------
        pos - numpy array
            locations
        prop - numpy array
            property values at nodes

        Returns
        -------

        """
        values = np.zeros(pos.shape[0])
        values[:] = np.nan
        c, tri = self.evaluate_shape(pos[:, :2])
        xxConst, yyConst, xyConst = self.evaluate_shape_d2(tri)
        # xyConst = self.evaluate_mixed_derivative(tri)
        inside = tri > 0
        # vertices, c, elements, inside = self.get_elements_for_location(pos)
        values[inside] = np.sum(
            xxConst[inside, :] * self.properties[prop][self.elements[tri[inside], :]],
            axis=1,
        )
        values[inside] += np.sum(
            yyConst[inside, :] * self.properties[prop][self.elements[tri[inside], :]],
            axis=1,
        )
        values[inside] += np.sum(
            xyConst[inside, :] * self.properties[prop][self.elements[tri[inside], :]],
            axis=1,
        )

        return values

    def get_quadrature_points(self, npts=2):
        if npts == 2:

            v1 = self.nodes[self.edges][:, 0, :]
            v2 = self.nodes[self.edges][:, 1, :]
            cp = np.zeros((v1.shape[0], self.ncps, 2))
            cp[:, 0] = 0.25 * v1 + 0.75 * v2
            cp[:, 1] = 0.75 * v1 + 0.25 * v2
            return cp
        raise NotImplementedError("Only 2 point quadrature is implemented")

    def get_edge_normal(self, e):
        v = self.nodes[self.edges][:, 0, :] - self.nodes[self.edges][:, 1, :]
        # e_len = np.linalg.norm(v, axis=1)
        normal = np.array([v[:, 1], -v[:, 0]]).T
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        return normal
