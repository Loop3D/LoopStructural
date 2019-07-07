class DiscreteInterpolator(GeologicalInterpolator):
    """
    Base class for a discrete interpolator e.g. piecewise linear or finite difference
    """
    def __init__(self):


    def add_constraints_to_least_squares(self, A, B, idc):
        """
        Adds constraints to the least squares system
        :param A: A N x C block to add
        :param B: B N vector
        :param idc: N x C node indices
        :return:
        """
        if self.shape == 'rectangular':
            self.A.append(c[i] * w)
            self.row.append(self.c_)
            self.col.append(element[i])
        if self.shape == 'square':
            for i in range(4):
                self.B[element[i]] += (p.val * c[i] * w)
                for j in range(4):
                    self.A.append(c[i] * w * c[j] * w)
                    self.col.append(element[i])
                    self.row.append(element[j])

    def _solve(self, solver='spqr', clear=True):
        """
        Solve the least squares problem definied by PLI
        :param solver: string for solver
        :param clear:
        :return:
        """
        # map node indicies from global to region
        if self.shape == 'rectangular':
            cols = self.region_map[np.array(self.col)]

            self.AA = coo_matrix((np.array(self.A), (np.array(self.row), \
                                                     cols)), shape=(self.c_, self.nx), dtype=float)  # .tocsr()
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
            self.cc_ = sla.lsqr(self.AA, B)
        elif solver == 'lsmr' and self.shape == 'rectangular':
            self.cc_ = sla.lsmr(self.AA, B)
        elif solver == 'eigenlsqr' and self.shape == 'rectangular':
            try:
                import eigensparse
            except ImportError:
                print("eigen sparse not installed")
            self.c[self.region] = eigensparse.lsqrcg(self.AA.tocsr(), self.B)
            return
        elif solver == 'spqr' and self.shape == 'rectangular':
            sys.path.insert(0, '/home/lgrose/dev/cpp/PyEigen/build')
            import eigensparse
            self.c[self.region] = eigensparse.lsqspqr(self.AA.tocsr(), self.B)
            return
        if solver == 'lu' and self.shape == 'square':
            lu = sla.splu(self.AA.tocsr())
            b = self.B  # np.array([1, 2, 3, 4])
            self.c[self.region] = lu.solve(b)
            if clear:
                self.AA = None
                self.A = []
                self.col = []
                self.row = []
                lu = None
            return
        if solver == 'chol' and self.shape == 'square':
            try:
                from sksparse.cholmod import cholesky
            except ImportError:
                print("Scikit Sparse not installed try another solver e.g. lu")
                return
            factor = cholesky(self.AA.tocsc())
            self.c = factor(self.B)
            if clear:
                self.AA = None
                self.A = []
                self.col = []
                self.row = []
                factor = None
            return
        # M = self.AA.diagonal()
        # z = sla.lsqr(self.AA*sla.inv(M),B)
        self.c[self.region] = self.cc_[0]