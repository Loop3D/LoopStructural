"""
Piecewise linear interpolator using folds
"""
import logging

import numpy as np
from LoopStructural.interpolators.cython.dsi_helper import fold_cg

from LoopStructural.interpolators.piecewiselinear_interpolator import \
    PiecewiseLinearInterpolator

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class DiscreteFoldInterpolator(PiecewiseLinearInterpolator):
    """

    """
    def __init__(self, support, fold):
        """
        A piecewise linear interpolator that can also use fold constraints defined in Laurent et al., 2016

        Parameters
        ----------
        support
            discrete support with nodes and elements etc
        fold FoldEvent
            a fold event with a valid geometry
        """

        PiecewiseLinearInterpolator.__init__(self, support)
        self.type = ['foldinterpolator']
        self.fold = fold


    @classmethod
    def from_piecewise_linear_and_fold(cls, pli, fold):
        """
        Constructor from an existing piecewise linear interpolation object and a fold object
        copies data from the PLI to the DFI

        Parameters
        ----------
        pli : PiecewiseLinearInterpolator
            existing interpolator
        fold : FoldEvent
            a fold event with a valid

        Returns
        -------
        DiscreteFoldInterpolator

        """
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
        """

        Parameters
        ----------
        fold : FoldEvent
            a fold that contrains the geometry we are trying to add

        Returns
        -------

        """
        logger.error('updating fold, this should be done by accessing the fold attribute')
        self.fold = fold

    def add_fold_constraints(self, fold_orientation=10., fold_axis_w=10., fold_regularisation=[.1,0.01,0.01],
                             fold_normalisation=1.,
                             fold_norm=1., step=2):
        """

        Parameters
        ----------
        fold_orientation : double
            weight for the fold direction/orientation in the least squares system
        fold_axis_w : double
            weight for the fold axis in the least squares system
        fold_regularisation : list
            weight for the fold regularisation in the least squares system
        fold_normalisation : double
            weight for the fold norm constraint in the least squares system
        fold_norm
            length of the interpolation norm in the least squares system
        step: int
            array step for adding constraints


        Returns
        -------

        Notes
        -----
        For more information about the fold weights see EPSL paper by Gautier Laurent 2016

        """
        # get the gradient of all of the elements of the mesh
        eg = self.support.get_element_gradients(np.arange(self.support.n_elements))
        # get array of all nodes for all elements N,4,3
        nodes = self.support.nodes[self.support.get_elements()[np.arange(self.support.n_elements)]]
        # calculate the fold geometry for the elements barycentre
        deformed_orientation, fold_axis, dgz = \
            self.fold.get_deformed_orientation(self.support.barycentre())
        element_idx = np.arange(self.support.n_elements)
        np.random.shuffle(element_idx)
        # calculate element volume for weighting
        vecs = nodes[:, 1:, :] - nodes[:, 0, None, :]
        vol = np.abs(np.linalg.det(vecs)) / 6
        if fold_orientation is not None:
            """
            dot product between vector in deformed ori plane = 0
            """
            np.random.shuffle(element_idx)

            logger.info("Adding fold orientation constraint to %s w = %f"%(self.propertyname, fold_orientation))
            A = np.einsum('ij,ijk->ik', deformed_orientation[element_idx[::step],:], eg[element_idx[::step],:,:])
            A *= vol[element_idx[::step], None]
            A *= fold_orientation
            B = np.zeros(A.shape[0])
            idc = self.support.get_elements()[element_idx[::step],:]
            self.add_constraints_to_least_squares(A, B, idc, name='fold orientation')

        if fold_axis_w is not None:
            """
            dot product between axis and gradient should be 0
            """
            np.random.shuffle(element_idx)

            logger.info("Adding fold axis constraint to %s w = %f"%(self.propertyname,fold_axis_w))
            A = np.einsum('ij,ijk->ik', fold_axis[element_idx[::step],:], eg[element_idx[::step],:,:])
            A *= vol[element_idx[::step], None]
            A *= fold_axis_w
            B = np.zeros(A.shape[0]).tolist()
            idc = self.support.get_elements()[element_idx[::step],:]

            self.add_constraints_to_least_squares(A, B, idc, name='fold axis')

        if fold_normalisation is not None:
            """
            specify scalar norm in X direction
            """
            np.random.shuffle(element_idx)

            logger.info("Adding fold normalisation constraint to %s w = %f"%(self.propertyname,fold_normalisation))
            A = np.einsum('ij,ijk->ik', dgz[element_idx[::step],:], eg[element_idx[::step],:,:])
            A *= vol[element_idx[::step], None]
            A *= fold_normalisation
            B = np.ones(A.shape[0])

            if fold_norm is not None:
                B[:] = fold_norm
            B *= fold_normalisation
            B *= vol[element_idx[::step]]
            idc = self.support.get_elements()[element_idx[::step],:]

            self.add_constraints_to_least_squares(A, B, idc, name='fold normalisation')

        if fold_regularisation is not None:
            """
            fold constant gradient  
            """
            logger.info("Adding fold regularisation constraint to {} w = {} {} {}".format(self.propertyname,
            fold_regularisation[0],fold_regularisation[1],fold_regularisation[1]))

            idc, c, ncons = fold_cg(eg, dgz, self.support.get_neighbours(), self.support.get_elements(), self.support.nodes)
            A = np.array(c[:ncons, :])
            A *= fold_regularisation[0]
            B = np.zeros(A.shape[0])
            idc = np.array(idc[:ncons, :])
            self.add_constraints_to_least_squares(A, B, idc, name='fold regularisation 1')

            idc, c, ncons = fold_cg(eg, deformed_orientation, self.support.get_neighbours(), self.support.get_elements(), self.support.nodes)
            A = np.array(c[:ncons, :])
            A *= fold_regularisation[1]
            B = np.zeros(A.shape[0])
            idc = np.array(idc[:ncons, :])
            self.add_constraints_to_least_squares(A, B, idc, name='fold regularisation 2')

            idc, c, ncons = fold_cg(eg, fold_axis, self.support.get_neighbours(), self.support.get_elements(), self.support.nodes)
            A = np.array(c[:ncons, :])
            A *= fold_regularisation[2]
            B = np.zeros(A.shape[0])
            idc = np.array(idc[:ncons, :])
            self.add_constraints_to_least_squares(A, B, idc, name='fold regularisation 3')
