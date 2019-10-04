from LoopStructural.interpolators.discete_interpolator import DiscreteInterpolator
import numpy as np


class PiecewiseLinearInterpolator(DiscreteInterpolator):
    """
    Piecewise Linear Interpolator
    Approximates scalar field by finding coefficients to a piecewise linear
    equation on a tetrahedral mesh

    """

    def __init__(self, mesh):
        """

        Parameters
        ----------
        mesh
        """

        # whether to assemble a rectangular matrix or a square matrix
        self.shape = 'rectangular'
        self.support = mesh
        self.interpolator_type = 'PLI'

        self.region = self.support.regions['everywhere']
        self.region_map = np.zeros(mesh.n_nodes).astype(int)
        self.region_map[self.region] = np.array(range(0,len(self.region_map[self.region])))
        self.nx = len(self.support.nodes[self.region])

        DiscreteInterpolator.__init__(self)
        # TODO need to fix this, constructor of DI is breaking support
        self.support = mesh

        self.interpolation_weights = {'cgw': 0.1, 'cpw' : 1., 'gpw':1., 'tpw':1.}

    def copy(self):
        return PiecewiseLinearInterpolator(self.support)

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
        # can't reset here, clears fold constraints
        #self.reset()
        for key in kwargs:
            self.up_to_date = False
            self.interpolation_weights[key] = kwargs[key]
        if self.interpolation_weights['cgw'] > 0.:
            self.up_to_date = False
            self.add_constant_gradient(self.interpolation_weights['cgw'])
        self.add_gradient_ctr_pts(self.interpolation_weights['gpw'])
        self.add_ctr_pts(self.interpolation_weights['cpw'])
        self.add_tangent_ctr_pts(self.interpolation_weights['tpw'])

    def add_constant_gradient(self, w=0.1):
        """
        Add the constant gradient regularisation to the system
        Parameters
        ----------
        w (double) - weighting of the cg parameter

        Returns
        -------

        """
        # iterate over all elements
        A, idc, B = self.support.get_constant_gradient(region=self.region, shape='rectangular')
        A = np.array(A)
        B = np.array(B)
        idc = np.array(idc)
        # w/=A.shape[0]
        self.add_constraints_to_least_squares(A*w,B*w,idc)
        return

    def add_gradient_ctr_pts(self, w=1.0):  # for now weight all gradient points the same
        """

        Parameters
        ----------
        w

        Returns
        -------

        """
        points = self.get_gradient_constraints()
        if points.shape[0] > 0:
            e, inside = self.support.elements_for_array(points[:, :3])
            nodes = self.support.nodes[self.support.elements[e]]
            vecs = nodes[:,1:,:] - nodes[:,0,None,:]
            vol = np.abs(np.linalg.det(vecs)) #/ 6
            d_t = self.support.get_elements_gradients(e)
            norm = np.linalg.norm(d_t,axis=2)
            d_t /= norm[:,:,None]
            # d_t *= vol[:,None,None]
            # w*=10^11

            points[:,3:] /= norm

            # add in the element gradient matrix into the inte
            e = np.tile(e,(3,1)).T
            idc = self.support.elements[e]
            w /= 3
            self.add_constraints_to_least_squares(d_t*w,points[:,3:]*w,idc)
    def add_norm_ctr_pts(self, w=1.0):
        """

        Parameters
        ----------
        w

        Returns
        -------

        """

        points = self.get_norm_constraints()
        if points.shape[0] > 0:
            e, inside = self.support.elements_for_array(points[:, :3])
            nodes = self.support.nodes[self.support.elements[e]]
            vecs = nodes[:,1:,:] - nodes[:,0,None,:]
            vol = np.abs(np.linalg.det(vecs)) #/ 6
            d_t = self.support.get_elements_gradients(e)

            d_t *= vol[:,None,None]
            # w*=10^11

            points[:,3:] /= np.linalg.norm(points[:,3:],axis=1)[:,None]

            # add in the element gradient matrix into the inte
            e=np.tile(e,(3,1)).T
            idc = self.support.elements[e]
            w /= 3
            self.add_constraints_to_least_squares(d_t*w,points[:,3:]*w*vol[:,None],idc)
    def add_tangent_ctr_pts(self, w=1.0):
        """

        Parameters
        ----------
        w

        Returns
        -------

        """
        return

    def add_ctr_pts(self, w=1.0):  # for now weight all value points the same
        """

        Parameters
        ----------
        w

        Returns
        -------

        """

        #get elements for points
        points = self.get_value_constraints()
        if points.shape[0] > 1:
            e, inside = self.support.elements_for_array(points[:, :3])
            # get barycentric coordinates for points
            nodes = self.support.nodes[self.support.elements[e]]
            vecs = nodes[:, 1:, :] - nodes[:, 0, None, :]
            vol = np.abs(np.linalg.det(vecs)) / 6
            A = self.support.calc_bary_c(e, points[:, :3])
            A *= vol[None,:]
            idc = self.support.elements[e]
            # w /= points.shape[0]
            self.add_constraints_to_least_squares(A.T*w, points[:, 3]*w*vol[None,:], idc)

    def add_gradient_orthogonal_constraint(self, elements, normals, w=1.0, B=0):
        """
        constraints scalar field to be orthogonal to a given vector
        Parameters
        ----------
        elements
        normals
        w
        B

        Returns
        -------

        """
        nodes = self.support.nodes[self.support.elements[elements]]
        vecs = nodes[:,1:,:] - nodes[:,0,None,:]
        vol = np.abs(np.linalg.det(vecs)) / 6
        d_t = self.support.get_elements_gradients(elements)
        dot_p = np.einsum('ij,ij->i', normals, normals)[:, None]
        mask = np.abs(dot_p) > 0
        normals[mask[:,0] ,:] =  normals[mask[:,0],:] / dot_p[mask][:,None]
        magnitude = np.einsum('ij,ij->i', normals, normals)
        normals[magnitude>0] = normals[magnitude>0] / magnitude[magnitude>0,None]
        A = np.einsum('ij,ijk->ik', normals, d_t)
        A *= vol[:,None]
        idc = self.support.elements[elements]
        B = np.zeros(len(elements))
        self.add_constraints_to_least_squares(A*w, B, idc)
