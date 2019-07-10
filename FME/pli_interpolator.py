from .discete_interpolator import DiscreteInterpolator
import numpy as np

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
        A, B, row, col, c_ = self.mesh.get_constant_gradient(region=self.region, shape='rectangular')
        A = np.array([A])
        B= np.array(B)
        col = np.array(col)
        self.add_constraints_to_least_squares(A*w,B*w,col)
        return

    def add_gradient_ctr_pts(self, w=1.0):  # for now weight all gradient points the same
        """
        add gradient norm constraints to the interpolator
        :param w: weight per constraint
        :return:
        """
        points = self.get_gradient_control()
        if points.shape[0] > 1:
            e,inside = self.mesh.elements_for_array(points[:,3:])
            d_t = self.mesh.get_element_gradients(e)
            points[3:] /= np.linalg.norm(points[:,3:],axis=1)
            #add in the element gradient matrix into the inte
            self.add_constraints_to_least_squares(d_t*w,points[:,3:]*w,e)

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
            self.add_constraints_to_least_squares(A*w,points[:,3]*w,idc)

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
        A = np.einsum('ij,ijk->ik', normals, d_t)
        idc = self.mesh.elements[elements]
        B = np.zeros(len(elements))
        self.add_constraints_to_least_squares(A,B,idc)
