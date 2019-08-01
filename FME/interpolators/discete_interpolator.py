from .geological_interpolator import GeologicalInterpolator
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse import linalg as sla
import timeit

class DiscreteInterpolator(GeologicalInterpolator):
    """
    Base class for a discrete interpolator e.g. piecewise linear or finite difference which is
    any interpolator that solves the system using least squares approximation
    """
    def __init__(self):
        GeologicalInterpolator.__init__(self)
        self.B = []
        if self.shape == 'square':
            self.B = np.zeros(self.nx)
        self.c_ = 0
        self.A = []  # sparse matrix storage coo format
        self.col = []
        self.row = []  # sparse matrix storage
        self.support = None
        self.solver = None
    def add_constraints_to_least_squares(self, A, B, idc):
        """
        Adds constraints to the least squares system
        :param A: A N x C block to add
        :param B: B N vector
        :param idc: N x C node indices
        :return:
        """

        nr = A.shape[0]
        if len(A.shape) >2:
            nr = A.shape[0]*A.shape[1]
        rows = np.arange(0,nr)
        rows = np.tile(rows, (A.shape[-1],1)).T
        rows+= self.c_
        self.c_+=nr

        if self.shape == 'rectangular':
            self.A.extend(A.flatten().tolist())
            self.row.extend(rows.flatten().tolist())
            self.col.extend(idc.flatten().tolist())
            self.B.extend(B.flatten().tolist())
        if self.shape == 'square':
            print('square')
            for i in range(4):
                self.B[element[i]] += (p.val * c[i] * w)
                for j in range(4):
                    self.A.append(c[i] * w * c[j] * w)
                    self.col.append(element[i])
                    self.row.append(element[j])

    def _solve(self, solver='spqr', clear=True):
        """
        Solve the least squares problem with specified solver
        :param solver: string for solver
        :param clear:
        :return:
        """
        # save the solver so we can rerun the interpolation at a later stage
        self.solver = solver
        # map node indicies from global to region
        if self.shape == 'rectangular':
            # print("Building rectangular sparse matrix")
            start = timeit.default_timer()
            cols = self.region_map[np.array(self.col)]

            self.AA = coo_matrix((np.array(self.A), (np.array(self.row), \
                                                     cols)), shape=(self.c_, self.nx), dtype=float)  # .tocsr()
            # print("Sparse matrix built in %f seconds"%(timeit.default_timer()-start))
        if self.shape == 'square':
            cols = np.array(self.col)
            rows = np.array(self.row)
            self.AA = coo_matrix((np.array(self.A), (np.array(rows).astype(np.int64), \
                                                     np.array(cols).astype(np.int64))), dtype=float)
            d = np.zeros(self.nx)
            d += np.finfo('float').eps
            self.AA += spdiags(d, 0, self.nx, self.nx)
        B = np.array(self.B)
        self.cc_ = [0, 0, 0, 0]
        self.c = np.zeros(self.mesh.n_nodes)
        self.c[:] = np.nan
        if solver == 'lsqr':
            # print("Solving using scipy lsqr")
            start = timeit.default_timer()
            self.cc_ = sla.lsqr(self.AA, B)
            self.up_to_date = True
            # print("Solving took %f seconds"%(timeit.default_timer()-start))
        elif solver == 'lsmr' and self.shape == 'rectangular':
            # print("Solving using scipy lsmr")
            start = timeit.default_timer()
            self.cc_ = sla.lsmr(self.AA, B)
            self.up_to_date = True

            # print("Solving took %f seconds" %(timeit.default_timer()-start))
        elif solver == 'eigenlsqr' and self.shape == 'rectangular':
            try:
                import eigensparse
            except ImportError:
                print("eigen sparse not installed")
            self.c[self.region] = eigensparse.lsqrcg(self.AA.tocsr(), self.B)
            self.up_to_date = True

            return
        elif solver == 'spqr' and self.shape == 'rectangular':
            sys.path.insert(0, '/home/lgrose/dev/cpp/PyEigen/build')
            import eigensparse
            self.c[self.region] = eigensparse.lsqspqr(self.AA.tocsr(), self.B)
            self.up_to_date = True

            return
        if solver == 'lu':
            # print("Solving using scipy LU decomposition")
            start = timeit.default_timer()
            if self.shape == 'rectangular':
                # print("Performing A.T @ A")
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            lu = sla.splu(A.tocsc())
            b = B  # np.array([1, 2, 3, 4])
            self.c[self.region] = lu.solve(b)
            self.up_to_date = True

            # print("Solving took %f seconds"%(timeit.default_timer()-start))

            return
        if solver == 'chol':
            try:
                from sksparse.cholmod import cholesky
            except ImportError:
                print("Scikit Sparse not installed try another solver e.g. lu")
                return
            # print("Solving using cholmod's cholesky decomposition")
            start = timeit.default_timer()
            if self.shape == 'rectangular':
                # print("Performing A.T @ A")
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            factor = cholesky(A.tocsc())
            self.c = factor(B)
            self.up_to_date = True

            # print("Solving took %f seconds"%(timeit.default_timer()-start))
            return
        if solver == 'cgp':
            num_iter=0
            def call(xk):
                nonlocal num_iter
                num_iter+=1
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B

            # precon = sla.spilu(A)
            # M2 = sla.LinearOperator(A.shape, precon.solve)
            # #print(precon)
            self.cc_ = sla.cg(A, B, M=diags(1/A.diagonal()),callback=call)
            print("num",num_iter)
            self.up_to_date = True
        if solver == 'eigencg':
            import sys
            sys.path.insert(0, '/home/lgrose/dev/cpp/PyEigen/build')
            import eigensparse
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.c = eigensparse.cg(A, B)
            self.up_to_date = True
            return
        if solver == 'cg':

            num_iter = 0

            def call(xk):
                nonlocal num_iter
                num_iter += 1
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            # precon = sla.spilu(A)
            # M2 = sla.LinearOperator(A.shape, precon.solve)
            #print(precon)
            self.cc_ = sla.cg(A,B,callback=call)#,M=M2)
            print("num",num_iter)
            if self.cc_[1] == 0:
                print("Conjugate gradient converged")
            if self.cc_[1] > 0:
                print("CG used %i iterations and didn't converge"%i)
            self.up_to_date = True
        if solver == 'cgs':
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.cc_ = sla.cgs(A,B)
            self.up_to_date = True

        if solver == 'bicg':
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.cc_ = sla.bicg(A,B)
            self.up_to_date = True

        if solver == 'qmr':
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.cc_ = sla.qmr(A,B)
            self.up_to_date = True

        if solver == 'gmres':
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.cc_ = sla.gmres(A,B)
            self.up_to_date = True

        if solver == 'lgmres':
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.cc_ = sla.lgmres(A,B)
            self.up_to_date = True

        if solver == 'minres':
            if self.shape == 'rectangular':
                A = self.AA.T.dot(self.AA)
                B = self.AA.T.dot(self.B)
            if self.shape == 'square':
                A = self.AA
                B = self.B
            self.cc_ = sla.minres(A,B)
            self.up_to_date = True

        self.c[self.region] = self.cc_[0]
        self.node_values = self.c

    def update(self):
        if self.solver is None:
            print("Cannot rerun interpolator")
            return
        if not self.up_to_date:
            self._solve(self.solver)