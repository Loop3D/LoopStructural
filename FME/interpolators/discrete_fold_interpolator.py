from ..interpolators.piecewiselinear_interpolator import PiecewiseLinearInterpolator


class DiscreteFoldInterpolator(PiecewiseLinearInterpolator):
    def __init__(self,**kwargs):
         PiecewiseLinearInterpolator.__init__(self,**kwargs)
         self.type = ['foldinterpolator']

         self.mesh = mesh
         self.fold = fold
         self.shape = shape


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
        eg = self.mesh.get_elements_gradients(np.arange(self.mesh.n_elements))
        deformed_orientation, fold_axis, dgz = self.fold.get_deformed_orientation(self.mesh.barycentre)
        if "fold_orientation" in kwargs:
            c = np.einsum('ij,ijk->ik', deformed_orientation, eg)
            A.extend(c*kwargs['fold_orienation'])
            B.extend(np.zeros(mesh.n_elements).tolist())
            col.extend(mesh.elements.tolist())
        if "fold_axis" in kwargs:
            c  = np.einsum('ij,ijk->ik', fold_axis, eg)
            A.extend(c*kwargs['fold_axis'])
            B.extend(np.zeros(mesh.n_elements).tolist())
            col.extend(mesh.elements.tolist())
        if "fold_normalisation" in kwargs:
            c  = np.einsum('ij,ijk->ik', dgz, eg)
            A.extend(c*kwargs['fold_normalisation'])
            b = np.ones(mesh.n_elements)
            if "fold_norm" in kwargs:
                b[:] = kwargs['fold_norm']

            B.extend(b)
            col.extend(mesh.elements.tolist())
        if "fold_regularisation" in kwargs:
            idc, c, ncons = fold_cg(eg, dgz, self.mesh.neighbours, self.mesh.elements, self.mesh.nodes)
            c = np.array(c)
            c*=kwargs['fold_regularisation']
            A.extend(c.tolist)
            B.extend(np.zeros(ncons))
            col.extend(np.array(idc).tolist())
        self.add_constraints_to_least_squares(A, B, col)


