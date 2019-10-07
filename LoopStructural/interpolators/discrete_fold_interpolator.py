from LoopStructural.interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator
from LoopStructural.cython.dsi_helper import fold_cg
import numpy as np

import logging
logger = logging.getLogger(__name__)


class DiscreteFoldInterpolator(PiecewiseLinearInterpolator):
    def __init__(self, mesh, fold ):
         PiecewiseLinearInterpolator.__init__(self, mesh)
         self.type = ['foldinterpolator']
         self.fold = fold

    @classmethod
    def from_piecewise_linear_and_fold(cls, pli, fold):

        # create a blank fold interpolator
        interpolator = cls(pli.support, fold)

        # copy the data and stuff from the existing interpolator
        interpolator.region = pli.region
        interpolator.shape = pli.shape
        interpolator.region_map = pli.region_map
        interpolator.p_i = pli.p_i
        interpolator.p_g = pli.p_g
        interpolator.p_t = pli.p_t
        interpolator.n_i = pli.n_i
        interpolator.n_g = pli.n_g
        interpolator.n_t = pli.n_t
        interpolator.propertyname = pli.propertyname
        return interpolator

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
            A*=kwargs['fold_orientation']
            B = np.zeros(self.support.n_elements)
            idc = self.support.elements
            self.add_constraints_to_least_squares(A, B, idc)

        if "fold_axis" in kwargs:
            """
            dot product between axis and gradient should be 0
            """
            A = np.einsum('ij,ijk->ik', fold_axis, eg)
            A*=vol[:,None]
            A*=kwargs['fold_axis']
            B = np.zeros(self.support.n_elements).tolist()
            self.add_constraints_to_least_squares(A , B, self.support.elements)

        if "fold_normalisation" in kwargs:
            """
            specify scalar norm in X direction
            """

            A  = np.einsum('ij,ijk->ik', dgz, eg)
            A*=vol[:,None]
            A*= kwargs['fold_normalisation']
            B = np.ones(self.support.n_elements)
            
            if "fold_norm" in kwargs:
                B[:] = kwargs['fold_norm']
            B*=kwargs['fold_normalisation']
            B*=vol
            self.add_constraints_to_least_squares(A, B, self.support.elements)

        if "fold_regularisation" in kwargs:
            """
            fold constant gradient  
            """
            idc, c, ncons = fold_cg(eg, dgz, self.support.neighbours, self.support.elements, self.support.nodes)
            A = np.array(c[:ncons,:])
            A *= kwargs['fold_regularisation']
            B = np.zeros(A.shape[0])
            idc = np.array(idc[:ncons,:])
            self.add_constraints_to_least_squares(A, B, idc)

