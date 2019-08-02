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
        
        A = []
        B = []
        col = []
        #get the gradient of all of the elements of the mesh
        eg = self.mesh.get_elements_gradients(np.arange(self.mesh.n_elements))
        #calculate the fold geometry for the elements barycentre
        deformed_orientation, fold_axis, dgz = \
            self.fold.get_deformed_orientation(self.mesh.barycentre)
        if "fold_orientation" in kwargs:
            """
            dot product between vector in deformed ori plane = 0
            """
            c = np.einsum('ij,ijk->ik', deformed_orientation, eg)
            A.extend(c*kwargs['fold_orientation'])
            B.extend(np.zeros(self.mesh.n_elements).tolist())
            col.extend(self.mesh.elements.tolist())
        if "fold_axis" in kwargs:
            """
            dot product between axis and gradient should be 0
            """
            c  = np.einsum('ij,ijk->ik', fold_axis, eg)
            A.extend(c*kwargs['fold_axis'])
            B.extend(np.zeros(self.mesh.n_elements).tolist())
            col.extend(self.mesh.elements.tolist())
        if "fold_normalisation" in kwargs:
            """
            specify scalar norm in X direction
            """

            c  = np.einsum('ij,ijk->ik', dgz, eg)
            A.extend(c*kwargs['fold_normalisation'])
            b = np.ones(self.mesh.n_elements)
            if "fold_norm" in kwargs:
                b[:] = kwargs['fold_norm']

            B.extend(b)
            col.extend(self.mesh.elements.tolist())
        if "fold_regularisation" in kwargs:
            """
            fold constant gradient  
            """
            idc, c, ncons = fold_cg(eg, dgz, self.mesh.neighbours, self.mesh.elements, self.mesh.nodes)
            c = np.array(c)
            c*=kwargs['fold_regularisation']
            A.extend(c.tolist())
            B.extend(np.zeros(ncons))
            col.extend(np.array(idc).tolist())
        self.add_constraints_to_least_squares(A, B, col)


