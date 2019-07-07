import sys

import numpy.linalg as la
import scipy.sparse.linalg as sla
from scipy.sparse import coo_matrix, spdiags

from .dsi_helper import compute_cg_regularisation_constraint
from .geological_interpolator import GeologicalInterpolator
from .geological_points import *


class PiecewiseLinearInterpolator(DiscreteInterpolator):
    """
    Piecewise Linear Interpolator
    Approximates scalar field by finding coefficients to a piecewise linear
    equation on a tetrahedral mesh

    """

    def __init__(self, mesh, **kwargs):
        """

        :param mesh: the mesh to apply PLI on
        :param kwargs: possible kwargs are 'region' being the subset of the mesh to approximate
        the linear equations on
        'propertyname' the name of the property that is interpolated on the mesh
        """
        if 'region' in kwargs:
            region = kwargs['region']
        if 'region' not in kwargs:
            region = 'everywhere'
        if 'propertyname' in kwargs:
            self.propertyname = kwargs['propertyname']
        # whether to assemble a rectangular matrix or a square matrix
        self.shape = 'rectangular'
        if 'shape' in kwargs:
            self.shape = kwargs['shape']
        DiscreteInterpolator.__init__(self)
        self.interpolator_type = 'PLI'
        self.mesh = mesh
        self.region = self.mesh.regions[region]
        self.region_map = np.zeros(mesh.n_nodes).astype(int)
        self.region_map[self.region] = np.array(range(0, len(self.region_map[self.region])))

        self.nx = len(self.mesh.nodes[self.region])



    def _setup_interpolator(self, **kwargs):
        """
        adds all of the constraints to the interpolation matrix
        :param kwargs: 'cgw' is the constant gradient weight
        'cpw' control point weight
        'gpw' gradient control point weight
        'tpw' tangent control point weight
        'cg' boolean is cg being used
        :return:
        """
        cgw = 0.1
        cpw = 1.0
        gpw = 1.0
        tpw = 1.0
        cg = True
        if 'cgw' in kwargs:
            cgw = kwargs['cgw']
        if 'cpw' in kwargs:
            cpw = kwargs['cpw']
        if 'gpw' in kwargs:
            gpw = kwargs['gpw']
        if 'tpw' in kwargs:
            tpw = kwargs['tpw']
        if 'cg' in kwargs:
            cg = kwargs['cg']
        if cg:
            self.add_constant_gradient(cgw)
        print("Setting up interpolator with %i value control points \n\
        %i gradient control points and %i tangent control points and \n\
        constant gradient regularization with a weight of %f" % (self.n_i, self.n_g, self.n_t, cgw))
        self.add_gradient_ctr_pts(gpw)
        self.add_ctr_pts(cpw)
        self.add_tangent_ctr_pts(tpw)

    def add_constant_gradient(self, w=0.1):
        """
        adds constant gradient regularisation to the PLI interpolator
        :param w: weighting (per constraint) to give the constant gradient interpolation
        :return:
        """
        # iterate over all elements
        A, B, row, col, c_ = self.mesh.get_constant_gradient(region=self.region, shape=self.shape)
        self.add_constraints_to_least_squares(A,B,col)
        return

    def add_gradient_ctr_pts(self, w=1.0):  # for now weight all gradient points the same
        """
        add gradient norm constraints to the interpolator
        :param w: weight per constraint
        :return:
        """
        e,inside = self.mesh.elements_for_array(pos)
        d_t = self.mesh.get_element_gradients(e)
        norm/ = np.linalg.norm(norm,axis=1)
        #add in the element gradient matrix into the inte
        self.add_constraint_to_least_squares(d_t,norm,e)

    def add_tangent_ctr_pts(self, w=1.0):
        """
        Add tangent constraint to the scalar field
        :param w: weight per constraint
        :return:
        """


    def add_ctr_pts(self, w=1.0):  # for now weight all value points the same
        """
        add value data to the interpolator
        :param w: weight per constraint
        :return:
        """
        #get elements for points
        e = self.mesh.elements_for_array(pos)
        #get barycentric coordinates for points

        A = self.mesh.calc_bary_c(e,pos)
        A*=w
        idc = self.mesh.elements[e]
        self.add_constraint_to_least_squares(A,B,idc)

    def add_elements_gradient_orthogonal_constraint(self, elements, normals, w=1.0, B=0):
        """
        constraints scalar field to be orthogonal to a given vector
        :param elements: index of elements to apply constraint to
        :param normals: list of normals for elements
        :param w: global weighting per constraint
        :param B: norm value
        :return:
        """
        d_t = self.mesh.get_elements_gradients(elements)
        normals = normals / np.einsum('ij,ij->i', normals, normals)[:, None]
        c = np.einsum('ij,ijk->ik', normals, d_t)
        a = self.mesh.elements[elements]

        c *= w
        if self.shape == 'rectangular':
            rows = np.arange(self.c_, self.c_ + len(elements))
            rows = np.tile(rows, (4, 1)).T
            b = np.zeros(len(elements))
            ##add the existing number of rows to the row count
            # rows+=self.c_
            # update dsi row counter
            self.c_ += len(elements)
            self.A.extend(c.flatten().tolist())
            self.col.extend(a.flatten().tolist())
            self.row.extend(rows.flatten().tolist())
            self.B.extend(b.tolist())
        if self.shape == 'square':
            cc = np.einsum('ij,ik->ijk', c, c)  # np.dot(c,c.T)
            for k in range(len(elements)):
                for i in range(4):
                    if B != 0:
                        self.B[k] += B * c[k, i]
                    for j in range(4):
                        self.A.append(cc[k, i, j])
                        self.row.append(a[k, i])
                        self.col.append(a[k, j])
                        # don't need to do anything with b as its 0
    def add_elements_gradient_constraint(self, elements, normal, w):
        # normals=normals/np.einsum('ij,ij->i',normals,normals)[:,None]
        for element in elements:
            d_t = self.mesh.get_element_gradient(self.mesh.elements[element])  # calculate_d(t_points)
            norm = normal
            scalar_product = np.zeros(4)
            for j in range(3):  # loop through gradient compon
                leng = 0.
                for i, alpha in enumerate(self.mesh.elements[element]):
                    scalar_product[i] = np.dot(d_t[:, i], norm)  # 2d/3d compat
                    # leng+= scalar_product[i]*scalar_product[i]
                    leng += d_t[j, i] * d_t[j, i]  # np.dot(d_t[:,i],d_t[:,i])
                leng = np.sqrt(leng)
                for i, alpha in enumerate(self.mesh.elements[element]):

                    self.B[alpha] += (norm[j] / leng) * (d_t[j, i] / leng)
                    for k, alpha2 in enumerate(self.mesh.elements[element]):
                        self.A.append((d_t[j, i] / leng) * w[element] * (d_t[j, k] / leng) * w[element])
                        self.row.append(alpha2)
                        self.col.append(alpha)

