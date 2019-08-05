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
        eg = self.mesh.get_elements_gradients(np.arange(self.mesh.n_elements))
        #calculate the fold geometry for the elements barycentre
        deformed_orientation, fold_axis, dgz = \
            self.fold.get_deformed_orientation(self.mesh.barycentre)
        if "fold_orientation" in kwargs:
            """
            dot product between vector in deformed ori plane = 0
            """
            A = np.einsum('ij,ijk->ik', deformed_orientation, eg)
            B = np.zeros(self.mesh.n_elements)
            idc = self.mesh.elements
            self.add_constraints_to_least_squares(A*kwargs['fold_orientation'], B, idc)
        if "fold_axis" in kwargs:
            """
            dot product between axis and gradient should be 0
            """
            A = np.einsum('ij,ijk->ik', fold_axis, eg)
            B = np.zeros(self.mesh.n_elements).tolist()
            self.add_constraints_to_least_squares(A*kwargs['fold_axis'], B, self.mesh.elements)

        if "fold_normalisation" in kwargs:
            """
            specify scalar norm in X direction
            """

            A  = np.einsum('ij,ijk->ik', dgz, eg)
            A*= kwargs['fold_normalisation']
            B = np.ones(self.mesh.n_elements)
            if "fold_norm" in kwargs:
                B[:] = kwargs['fold_norm']
            self.add_constraints_to_least_squares(A, B, self.mesh.elements)

        if "fold_regularisation" in kwargs:
            """
            fold constant gradient  
            """
            idc, c, ncons = fold_cg(eg, dgz, self.mesh.neighbours, self.mesh.elements, self.mesh.nodes)
            A = np.array(c)
            A*=kwargs['fold_regularisation']
            B = np.zeros(len(c))
            idc = np.array(idc)
            self.add_constraints_to_least_squares(A, B, idc)

