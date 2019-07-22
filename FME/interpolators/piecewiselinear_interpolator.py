from FME.interpolators.discete_interpolator import DiscreteInterpolator
from FME.modelling.scalar_field import TetrahedralMeshScalarField
import numpy as np

class PiecewiseLinearInterpolator(DiscreteInterpolator):
    """
    Piecewise Linear Interpolator
    Approximates scalar field by finding coefficients to a piecewise linear
    equation on a tetrahedral mesh

    """

    def __init__(self, mesh, propertyname = 'InterpolatedProperty', **kwargs):
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
        self.propertyname = propertyname
        # whether to assemble a rectangular matrix or a square matrix
        self.shape = 'rectangular'
        if 'shape' in kwargs:
            self.shape = kwargs['shape']
        self.mesh = mesh

        self.interpolator_type = 'PLI'
        self.region = self.mesh.regions[region]
        self.region_map = np.zeros(mesh.n_nodes).astype(int)
        self.region_map[self.region] = np.array(range(0, len(self.region_map[self.region])))
        self.nx = len(self.mesh.nodes[self.region])

        DiscreteInterpolator.__init__(self)





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
        A, idc, B = self.mesh.get_constant_gradient(region=self.region, shape='rectangular')
        A = np.array(A)
        B = np.array(B)
        idc = np.array(idc)
        self.add_constraints_to_least_squares(A*w,B*w,idc)
        return

    def add_gradient_ctr_pts(self, w=1.0):  # for now weight all gradient points the same
        """
        add gradient norm constraints to the interpolator
        :param w: weight per constraint
        :return:
        """
        points = self.get_gradient_control()
        if points.shape[0] > 0:
            e, inside = self.mesh.elements_for_array(points[:,3:])
            d_t = self.mesh.get_elements_gradients(e)
            points[:,3:] /= np.linalg.norm(points[:,3:],axis=1)[:,None]
            #add in the element gradient matrix into the inte
            e=np.tile(e,(3,1)).T
            idc = self.mesh.elements[e]
            self.add_constraints_to_least_squares(d_t*w,points[:,3:]*w,idc)

    def add_tangent_ctr_pts(self, w=1.0):
        """
        Add tangent constraint to the scalar field
        :param w: weight per constraint
        :return:
        """
        return

    def add_ctr_pts(self, w=1.0):  # for now weight all value points the same
        """
        add value data to the interpolator
        :param w: weight per constraint
        :return:
        """
        #get elements for points
        points = self.get_control_points()
        if points.shape[0] > 1:
            e, inside = self.mesh.elements_for_array(points[:,:3])
            #get barycentric coordinates for points

            A = self.mesh.calc_bary_c(e,points[:,:3])
            idc = self.mesh.elements[e]
            self.add_constraints_to_least_squares(A.T*w,points[:,3]*w,idc)

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
        dot_p = np.einsum('ij,ij->i', normals, normals)[:, None]
        mask = np.abs(dot_p) > 0
        normals[mask[:,0] ,:] =  normals[mask[:,0],:] / dot_p[mask][:,None]
        normals = normals / np.einsum('ij,ij->i', normals, normals)[:, None]
        A = np.einsum('ij,ijk->ik', normals, d_t)
        idc = self.mesh.elements[elements]
        B = np.zeros(len(elements))
        self.add_constraints_to_least_squares(A,B,idc)

    def get_support(self):

        # TODO add check to see if interpolant is up to date
        # this requires adding solving parameters to interpolator object
        # if self.up_to_date == False:
        #     self.so

        return TetrahedralMeshScalarField.from_node_values(
            self.mesh,
            self.propertyname,
            self.node_values)