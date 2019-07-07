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

        self.B = []  # np.zeros((self.nx,1))
        if self.shape == 'square':
            self.B = np.zeros(self.nx)
        self.c_ = 0
        self.A = []  # sparse matrix storage coo format
        self.col = []
        self.row = []  # sparse matrix storage
        self.mesh.dinfo = {}  # [:] = False #intitialise dinfo to be 0

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
        A = np.array(A)
        # col = self.region_map[col]
        if self.shape == 'rectangular':
            A *= w
        if self.shape == 'square':
            A *= w * w
            # row = self.region_map[row]
        A = A.tolist()
        self.A.extend(A)
        if self.shape == 'rectangular':
            self.B.extend(B)
        self.row.extend(row)
        self.col.extend(col)
        self.c_ += c_
        return

    def add_ctr_pts(self, w=1.0):  # for now weight all value points the same
        """
        add value data to the interpolator
        :param w: weight per constraint
        :return:
        """
        for p in self.p_i:
            element, flag = self.mesh.get_element(p.pos)
            if flag == False:
                print('Could not find triangle for x:%f y:%f z:%f' % (p.pos[0], p.pos[1], p.pos[2]))
                continue
            if ~np.all(self.region[element]):
                print('Could not find element for x:%f y:%f z:%f inside region' % (p.pos[0], p.pos[1], p.pos[2]))

                continue
            points = self.mesh.nodes[element]
            # points -=points[0,:]

            vap = np.zeros(3)
            vbp = np.zeros(3)

            vcp = np.zeros(3)
            vdp = np.zeros(3)
            vab = np.zeros(3)
            vac = np.zeros(3)
            vad = np.zeros(3)
            vbc = np.zeros(3)
            vbd = np.zeros(3)
            temp = np.zeros(3)
            ##bp = 
            for i in range(3):
                vap[i] = p.pos[i] - points[0, i]
                vbp[i] = p.pos[i] - points[1, i]
                vcp[i] = p.pos[i] - points[2, i]
                vdp[i] = p.pos[i] - points[3, i]
                vab[i] = points[1, i] - points[0, i]
                vac[i] = points[2, i] - points[0, i]
                vad[i] = points[3, i] - points[0, i]
                vbc[i] = points[2, i] - points[1, i]
                vbd[i] = points[3, i] - points[1, i]

            va = np.dot(vbp, np.cross(vbd, vbc)) / 6.
            vb = np.dot(vap, np.cross(vac, vad)) / 6.
            vc = np.dot(vap, np.cross(vad, vab)) / 6.
            vd = np.dot(vap, np.cross(vab, vac)) / 6.
            v = np.dot(vab, np.cross(vac, vad)) / 6.
            c = np.zeros(4)
            c[0] = va / v
            c[1] = vb / v
            c[2] = vc / v
            c[3] = vd / v
            if self.shape == 'rectangular':
                for i in range(len(c)):
                    self.A.append(c[i] * w)
                    self.row.append(self.c_)
                    self.col.append(element[i])
                self.B.append(p.val * w)
                self.c_ += 1
            if self.shape == 'square':
                for i in range(4):
                    self.B[element[i]] += (p.val * c[i] * w)
                    for j in range(4):
                        self.A.append(c[i] * w * c[j] * w)
                        self.col.append(element[i])
                        self.row.append(element[j])

    def add_gradient_ctr_pts(self, w=1.0):  # for now weight all gradient points the same
        """
        add gradient norm constraints to the interpolator
        :param w: weight per constraint
        :return:
        """
        for p in self.p_g:

            element, flag = self.mesh.get_element(p.pos)
            if flag == False:
                print('Could not find triangle for %f %f %f' % (p.pos[0], p.pos[1], p.pos[2]))
                continue
            if ~np.all(self.region[element]):
                print('Could not find element for x:%f y:%f z:%f inside region' % (p.pos[0], p.pos[1], p.pos[2]))

                continue
            d_t = self.mesh.get_element_gradient(element)  # calculate_d(t_points)
            norm = p.dir_() / la.norm(p.dir_())
            scalar_product = np.zeros(len(element))
            if self.shape == 'rectangular':
                for j in range(len(p.pos)):  # loop through gradient compon
                    leng = 0.
                    for i, alpha in enumerate(element):
                        leng += d_t[j, i] * d_t[j, i]  # np.dot(d_t[:,i],d_t[:,i])
                    leng = np.sqrt(leng)
                    for i, alpha in enumerate(element):
                        self.A.append((d_t[j, i] / leng) * w)
                        self.row.append(self.c_)
                        self.col.append(alpha)
                    self.B.append(norm[j] / leng)
                    self.c_ = self.c_ + 1
            if self.shape == 'square':
                for j in range(len(p.pos)):  # loop through gradient compon
                    leng = 0.
                    for i, alpha in enumerate(element):
                        leng += d_t[j, i] * d_t[j, i]  # np.dot(d_t[:,i],d_t[:,i])
                    leng = np.sqrt(leng)
                    for i, alpha in enumerate(element):

                        self.B[element[i]] += (norm[j] / leng) * (d_t[j, i] / leng)
                        for k, alpha2 in enumerate(element):
                            self.A.append((d_t[j, i] / leng) * w * (d_t[j, k] / leng) * w)
                            self.row.append(alpha2)
                            self.col.append(alpha)
                    # self.B.append(norm[j]/leng)
                    # self.c_ = self.c_+1

    def add_tangent_ctr_pts(self, w=1.0):
        """
        Add tangent constraint to the scalar field
        :param w: weight per constraint
        :return:
        """
        for p in self.p_t:
            element = self.mesh.get_closest_triangle(p.pos)
            if element is None:
                print('Could not find triangle for %f %f %f' % (p.pos[0], p.pos[1], p.pos[2]))
                return False
            t_points = self.mesh.nodes[element]
            d_t = self.calculate_d(t_points)
            norm = p.dir_() / la.norm(g)
            scalar_product = np.zeros(3)
            n = 0.
            for i in range(len(t)):
                scalar_product[i] = np.dot(d_t[:, i], norm)
                n += scalar_product[i] * scalar_product[i]
            n = np.sqrt(n)
            scalar_product /= n
            for i, alpha in enumerate(t):
                self.A.append(scalar_product[i] * w)
                self.row.append(self.c_)
                self.col.append(alpha)
            self.B.append(0.0)
            self.c_ = self.c_ + 1

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

