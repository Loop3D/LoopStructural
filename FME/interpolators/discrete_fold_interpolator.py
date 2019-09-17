from ..interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator
from ..cython.dsi_helper import fold_cg
import numpy as np

class DiscreteFoldInterpolator(PiecewiseLinearInterpolator):
    def __init__(self, mesh, fold ):
         PiecewiseLinearInterpolator.__init__(self, mesh)
         self.type = ['foldinterpolator']
         self.fold = fold

    def update_fold(self, fold):
        self.fold = fold

    def add_fold_constraints(self,**kwargs):
        """
        add the fold constraints to the interpolation matrix
        using the fold object
        :param kwargs:
        :return:
        """
        

        #get the gradient of all of the elements of the mesh
        eg = self.support.get_elements_gradients(np.arange(self.support.n_elements))
        # get array of all nodes for all elements N,4,3
        nodes = self.support.nodes[self.support.elements[np.arange(self.support.n_elements)]]
        #calculate the fold geometry for the elements barycentre
        deformed_orientation, fold_axis, dgz = \
            self.fold.get_deformed_orientation(self.support.barycentre)

        # calculate element volume for weighting
        vecs = nodes[:, 1:, :] - nodes[:, 0, None, :]
        vol = np.abs(np.linalg.det(vecs)) / 6
        if "fold_orientation" in kwargs:
            """
            dot product between vector in deformed ori plane = 0
            """
            A = np.einsum('ij,ijk->ik', deformed_orientation, eg)
            A*=vol[:,None]
            B = np.zeros(self.support.n_elements)
            idc = self.support.elements
            self.add_constraints_to_least_squares(A*kwargs['fold_orientation'], B, idc)
        if "fold_axis" in kwargs:
            """
            dot product between axis and gradient should be 0
            """
            A = np.einsum('ij,ijk->ik', fold_axis, eg)
            A*=vol[:,None]
            B = np.zeros(self.support.n_elements).tolist()
            self.add_constraints_to_least_squares(A * kwargs['fold_axis'], B, self.support.elements)

        if "fold_normalisation" in kwargs:
            """
            specify scalar norm in X direction
            """

            A  = np.einsum('ij,ijk->ik', dgz, eg)
            A*= kwargs['fold_normalisation']
            A*=vol[:,None]
            B = np.ones(self.support.n_elements)
            if "fold_norm" in kwargs:
                B[:] = kwargs['fold_norm']
            self.add_constraints_to_least_squares(A, B*vol*kwargs['fold_normalisation'], self.support.elements)

        if "fold_regularisation" in kwargs:
            """
            fold constant gradient  
            """
            idc, c, ncons = fold_cg(eg, dgz, self.support.neighbours, self.support.elements, self.support.nodes)
            A = np.array(c)
            A *= kwargs['fold_regularisation']
            B = np.zeros(len(c))
            idc = np.array(idc)
            self.add_constraints_to_least_squares(A, B, idc)

