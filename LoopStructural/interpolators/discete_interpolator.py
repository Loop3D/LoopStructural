from LoopStructural.interpolators.geological_interpolator import GeologicalInterpolator
import numpy as np
from scipy.sparse import coo_matrix, diags, bmat, eye
from scipy.sparse import linalg as sla

import logging
logger = logging.getLogger(__name__)

class DiscreteInterpolator(GeologicalInterpolator):
    """
    Base class for a discrete interpolator e.g. piecewise linear or finite difference which is
    any interpolator that solves the system using least squares approximation
    """
    def __init__(self, support):
        GeologicalInterpolator.__init__(self)
        self.B = []
        self.support = support

        region = 'everywhere'
        self.region = self.support.regions[region]
        self.region_map = np.zeros(support.n_nodes).astype(int)
        self.region_map[self.region] = np.array(range(0,len(self.region_map[self.region])))
        self.nx = len(self.support.nodes[self.region])
        if self.shape == 'square':
            self.B = np.zeros(self.nx)
        self.c_ = 0
        self.A = []  # sparse matrix storage coo format
        self.col = []
        self.row = []  # sparse matrix storage
        self.solver = None
        self.eq_const_C = []
        self.eq_const_row = []
        self.eq_const_col = []
        self.eq_const_d = []
        self.eq_const_c_ = 0


    def set_property_name(self, propertyname):
        """
        Set the property name attribute, this is usually used to
        save the property on the support
        Parameters
        ----------
        propertyname

        Returns
        -------

        """
        self.propertyname = propertyname

    def set_region(self, regionname=None, region=None):
        """
        Set the region of the support the interpolator is working on
        Parameters
        ----------
        regionname - string with name of region
        region - numpy mask for region

        Returns
        -------

        """
        if region is not None:
            self.region = region
        if regionname is not None:
            self.region = self.support.regions[regionname]

        self.region_map = np.zeros(self.support.n_nodes).astype(int)
        self.region_map[self.region] = np.array(range(0,len(self.region_map[self.region])))
        self.nx = len(self.support.nodes[self.region])

    def set_interpolation_weights(self, weights):
        """
        Set the interpolation weights dictionary
        Parameters
        ----------
        weights - dictionary with the interpolation weights

        Returns
        -------

        """
        for key in weights:
            self.up_to_date = False
            self.interpolation_weights[key] = weights[key]

    def reset(self):
        """
        Reset the interpolation constraints
        """
        logger.debug("Resetting interpolation constraints")
        self.c_ = 0
        self.A = []  # sparse matrix storage coo format
        self.col = []
        self.row = []  # sparse matrix storage
        self.eq_const_C = []
        self.eq_const_row = []
        self.eq_const_col = []
        self.eq_const_d = []
        self.eq_const_c_ = 0
        self.B = []

    def add_constraints_to_least_squares(self, A, B, idc):
        """
        Adds constraints to the least squares system. Automatically works out the row
        index given the shape of the input arrays
        Parameters
        ----------
        A - RxC numpy array of constraints where C is number of columns,R rows
        B - B values array length R
        idc - RxC column index

        Returns
        -------

        """
        A = np.array(A)
        B = np.array(B)
        idc = np.array(idc)
        if np.any(np.isnan(idc)) or np.any(np.isnan(A)) or np.any(np.isnan(B)):
            logger.warning("Constraints contain nan not adding constraints")
            return
        nr = A.shape[0]
        if len(A.shape) >2:
            nr = A.shape[0]*A.shape[1]
        rows = np.arange(0,nr).astype(int)
        rows = np.tile(rows, (A.shape[-1],1)).T
        rows+= self.c_
        self.c_+=nr
        if self.shape == 'rectangular':
            # don't add operator where it is = 0 to the sparse matrix!
            A = A.flatten()
            rows = rows.flatten()
            idc = idc.flatten()
            B = B.flatten()
            mask = A == 0
            self.A.extend(A[~mask].tolist())
            self.row.extend(rows[~mask].tolist())
            self.col.extend(idc[~mask].tolist())
            self.B.extend(B.tolist())

    def add_equality_constraints(self, node_idx, values):
        """
        Adds hard constraints to the least squares system. For now this just sets
        the node values to be fixed using a lagrangian.
        Parameters
        ----------
        node_idx - int array of node indexes
        values - array of node values

        Returns
        -------

        """
        # map from mesh node index to region node index
        gi = np.zeros(self.support.n_nodes)
        gi[:] = -1
        gi[self.region] = np.arange(0,self.nx)
        idc = gi[node_idx]
        # check
        outside = ~(idc == -1)
        self.eq_const_C.extend(np.ones(idc[outside].shape[0]).tolist())
        self.eq_const_col.extend(idc[outside].tolist())
        self.eq_const_row.extend((np.arange(0,idc[outside].shape[0])))
        self.eq_const_d.extend(values[outside].tolist())
        self.eq_const_c_ += node_idx[outside].shape[0]

    def build_matrix(self, damp=True):
        """
        Assemble constraints into interpolation matrix. Adds equaltiy constraints
        using lagrange modifiers if necessary
        Parameters
        ----------
        damp: bool
            Flag whether damping should be added to the diagonal of the matrix
        Returns
        -------
        Interpolation matrix and B
        """
        cols = np.array(self.col)
        A = coo_matrix((np.array(self.A), (np.array(self.row), \
                                                 cols)), shape=(self.c_, self.nx), dtype=float)  # .tocsr()
        B = np.array(self.B)
        AAT = A.T.dot(A)
        BT = A.T.dot(B)
        # add a small number to the matrix diagonal to smooth the results
        # can help speed up solving, but might also introduce some errors
        if damp:
            AAT += eye(AAT.shape[0])*np.finfo('float').eps
        if self.eq_const_c_ > 0:
            # solving constrained least squares using
            # | ATA CT | |c| = b
            # | C   0  | |y|   d
            # where A is the interpoaltion matrix
            # C is the equality constraint matrix
            # b is the interpolation constraints to be honoured
            # in a least squares sense
            # and d are the equality constraints
            # c are the node values and y are the
            # lagrange multipliers#
            C = coo_matrix((np.array(self.eq_const_C), (np.array(self.eq_const_row),
                                                             np.array(self.eq_const_col))),
                                shape=(self.eq_const_c_, self.nx))
            d = np.array(self.eq_const_d)
            AAT = bmat([[AAT, C.T], [C, None]])
            BT = np.hstack([BT,d])
        return AAT, BT

    def _solve_lu(self, A, B):
        """
        Call scipy LU decomoposition
        Parameters
        ----------
        A - square sparse matrix
        B

        Returns
        -------

        """
        lu = sla.splu(A.tocsc())
        # b = np.hstack([B, d])  # np.array([1, 2, 3, 4])
        sol = lu.solve(B)
        return sol[:self.nx]

    def _solve_chol(self, A ,B):
        """
        Call suitesparse cholmod through scikitsparse
        LINUX ONLY!
        Parameters
        ----------
        A - square sparse matrix
        B - numpy vector

        Returns
        -------

        """
        try:
            from sksparse.cholmod import cholesky
            factor = cholesky(A.tocsc())
            return factor(B)[:self.nx]
        except ImportError:
            logger.warning("Scikit Sparse not installed try using cg instead")
            return False

    def _solve_cg(self, A, B, precon=None, **kwargs):
        """
        Call scipy conjugate gradient
        Parameters
        ----------
        A - square sparse matrix
        B - numpy vector
        precon
        kwargs

        Returns
        -------

        """
        cgargs={}
        if 'maxiter' in kwargs:
            cgargs['maxiter'] = kwargs['maxiter']
        if 'x0' in kwargs:
            cgargs['x0'] = kwargs['x0']
        if 'tol' in kwargs:
            cgargs['tol'] = kwargs['tol']
        if 'atol' in kwargs:
            cgargs['atol'] = kwargs['atol']
        if 'callback' in kwargs:
            cgargs['callback'] = kwargs['callback']
        if precon is not None:
            cgargs['M'] = precon(A)
        return sla.cg(A,B,**cgargs)[0][:self.nx]

    def _solve(self, solver, **kwargs):
        """
        Main entry point to run the solver and update the node value attribute for the
        discreteinterpolator class
        Parameters
        ----------
        solver - string of solver e.g. cg, lu, chol, custom
        kwargs - kwargs for solver e.g. maxiter, preconditioner etc, damping for

        Returns
        -------
        True if the interpolation is run

        """
        self.c = np.zeros(self.support.n_nodes)
        self.c[:] = np.nan
        damp = False
        if 'damp' in kwargs:
            damp = kwargs['damp']
        A, B = self.build_matrix(damp=damp)

        # run the chosen solver
        if solver == 'cg':
            logger.info("Solving using conjugate gradient")
            self.c[self.region] = self._solve_cg(A, B, **kwargs)
        if solver == 'chol':
            self.c[self.region] = self._solve_chol(A, B)
        if solver == 'lu':
            logger.info("Solving using scipy LU")
            self.c[self.region] = self._solve_lu(A, B)

        if solver == 'external':
            logger.warning("Using external solver")
            self.c[self.region] = kwargs['external'](A, B)[:self.nx]
        # check solution is not nan
        if np.all(self.c == np.nan):
            logger.warning("Solver not run, no scalar field")
        # if solution is all 0, probably didn't work
        if np.all(self.c[self.region] == 0):
            logger.warning("No solution, scalar field 0. Add more data.")

    def update(self):
        """
        Check if the solver is up to date, if not rerun interpolation using
        the previously used solver. If the interpolation has not been run before it will
        return False
        Returns
        -------

        """
        if self.solver is None:
            logging.debug("Cannot rerun interpolator")
            return False
        if not self.up_to_date:
            self.setup_interpolator()
            return self._solve(self.solver)
