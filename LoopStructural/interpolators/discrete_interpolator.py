"""
Discrete interpolator base for least squares
"""
import logging

import numpy as np
from scipy.sparse import coo_matrix, bmat, eye
from scipy.sparse import linalg as sla

from LoopStructural.interpolators.geological_interpolator import \
    GeologicalInterpolator

from LoopStructural.utils import getLogger
logger = getLogger(__name__)


class DiscreteInterpolator(GeologicalInterpolator):
    """

    """
    def __init__(self, support):
        """
        Base class for a discrete interpolator e.g. piecewise linear or finite difference which is
        any interpolator that solves the system using least squares approximation

        Parameters
        ----------
        support 
            A discrete mesh with, nodes, elements, etc
        """
        GeologicalInterpolator.__init__(self)
        self.B = []
        self.support = support
        self.region_function = lambda xyz : np.ones(xyz.shape[0],dtype=int)
        # self.region_map[self.region] = np.array(range(0,
        # len(self.region_map[self.region])))
        self.shape = 'rectangular'
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
        self.constraints = {}
        self.interpolation_weights= {}
        logger.info("Creating discrete interpolator with {} degrees of freedom".format(self.nx))
        self.type = 'discrete'
    @property
    def nx(self):
        return len(self.support.nodes[self.region])

    @property
    def region(self):
        return self.region_function(self.support.nodes)

    @property
    def region_map(self):
        region_map = np.zeros(self.support.n_nodes).astype(int)
        region_map[self.region] = np.array(
            range(0, len(region_map[self.region])))
        return region_map
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

    def set_region(self, region=None):
        """
        Set the region of the support the interpolator is working on

        Parameters
        ----------
        region - function(position)
            return true when in region, false when out

        Returns
        -------

        """
        # evaluate the region function on the support to determine
        # which nodes are inside update region map and degrees of freedom
        self.region_function = region
        logger.info("Interpolation now uses region and has {} degrees of freedom".format(self.nx))

    def set_interpolation_weights(self, weights):
        """
        Set the interpolation weights dictionary

        Parameters
        ----------
        weights - dictionary
            Entry of new weights to assign to self.interpolation_weights

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
        self.n_constraints = 0

    def add_constraints_to_least_squares(self, A, B, idc, name='undefined'):
        """
        Adds constraints to the least squares system. Automatically works
        out the row
        index given the shape of the input arrays

        Parameters
        ----------
        A : numpy array / list
            RxC numpy array of constraints where C is number of columns,R rows
        B : numpy array /list
            B values array length R
        idc : numpy array/list
            RxC column index

        Returns
        -------
        list of constraint ids

        """
        A = np.array(A)
        B = np.array(B)
        idc = np.array(idc)
        nr = A.shape[0]
        #logger.debug('Adding constraints to interpolator: {} {} {}'.format(A.shape[0]))
        # print(A.shape,B.shape,idc.shape)
        if A.shape != idc.shape:
            return
        
        if len(A.shape) > 2:
            nr = A.shape[0] * A.shape[1]
            A = A.reshape((A.shape[0]*A.shape[1],A.shape[2]))
            idc = idc.reshape((idc.shape[0]*idc.shape[1],idc.shape[2]))
        # going to assume if any are nan they are all nan
        mask = np.any(np.isnan(A),axis=1)
        A[mask,:] = 0
        if np.any(np.isnan(idc)) or np.any(np.isnan(A)) or np.any(np.isnan(B)):
            logger.warning("Constraints contain nan not adding constraints: {}".format(name))
            # return
        
        rows = np.arange(0, nr).astype(int)
        rows += self.c_
        constraint_ids = rows.copy()

        if name in self.constraints:            
            
            self.constraints[name]['A'] =  np.vstack([self.constraints[name]['A'],A])
            self.constraints[name]['B'] =  np.hstack([self.constraints[name]['B'], B])
            self.constraints[name]['idc'] = np.vstack([self.constraints[name]['idc'],
                                                idc])
                                   
        if name not in self.constraints:
            self.constraints[name] = {'node_indexes':constraint_ids,'A':A,'B':B.flatten(),'idc':idc}
        rows = np.tile(rows, (A.shape[-1], 1)).T

        self.c_ += nr
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
    
    def calculate_residual_for_constraints(self):
        residuals = {}
        for constraint_name, constraint in self.constraints:
            residuals[constraint_name] = np.einsum('ij,ij->i',constraint['A'],self.c[constraint['idc'].astype(int)]) - constraint['B'].flatten()
        return residuals
    def remove_constraints_from_least_squares(self, name='undefined',
                                              constraint_ids=None):
        """
        Remove constraints from the least squares system using the constraint ids
        which corresponds to the rows in the interpolation matrix.

        Parameters
        ----------
        constraint_ids : np.array(dtype=int)
            id of constraints to remove

        Returns
        -------

        """

        if constraint_ids is None:
            constraint_ids = self.constraints[name]
        print("Removing {} {} constraints from least squares".format(len(constraint_ids), name))
        A = np.array(self.A)
        B = np.array(self.B)
        col = np.array(self.col)
        row = np.array(self.row)
        mask = np.any((row[:,None] == constraint_ids[None,:]) == True,
                      axis=1)
        # np.any((numbers[:, None] == np.array([0, 10, 30])[None, :]) == True,
        #        axis=1)
        bmask = np.ones(B.shape,dtype=bool)
        bmask[constraint_ids] = 0
        self.A = A[~mask].tolist()
        self.B = B[bmask]
        self.col = col[~mask].tolist()
        rowmax = np.max(row[mask])
        rowrange = rowmax-np.min(row[mask])
        # row[np.logical_and(~mask,row>rowmax)] -= rowrange
        return row[~mask]

    def add_equality_constraints(self, node_idx, values):
        """
        Adds hard constraints to the least squares system. For now this just
        sets
        the node values to be fixed using a lagrangian.

        Parameters
        ----------
        node_idx : numpy array/list
            int array of node indexes
        values : numpy array/list
            array of node values

        Returns
        -------

        """
        # map from mesh node index to region node index
        gi = np.zeros(self.support.n_nodes)
        gi[:] = -1
        gi[self.region] = np.arange(0, self.nx)
        idc = gi[node_idx]
        outside = ~(idc == -1)

        self.eq_const_C.extend(np.ones(idc[outside].shape[0]).tolist())
        self.eq_const_col.extend(idc[outside].tolist())
        self.eq_const_row.extend((np.arange(0, idc[outside].shape[0])))
        self.eq_const_d.extend(values[outside].tolist())
        self.eq_const_c_ += idc[outside].shape[0]

    def add_tangent_ctr_pts(self, w=1.0):
        """

        Parameters
        ----------
        w : double


        Returns
        -------

        """
        points = self.get_tangent_constraints()
        if points.shape[0] > 1:
            self.add_gradient_orthogonal_constraint(points[:,:3],points[:,3:6],w)

    def build_matrix(self, square=True, damp=True):
        """
        Assemble constraints into interpolation matrix. Adds equaltiy
        constraints
        using lagrange modifiers if necessary

        Parameters
        ----------
        damp: bool
            Flag whether damping should be added to the diagonal of the matrix
        Returns
        -------
        Interpolation matrix and B
        """

        logger.info("Interpolation matrix is %i x %i"%(self.c_,self.nx))
        cols = np.array(self.col)
        A = coo_matrix((np.array(self.A), (np.array(self.row), \
                                           cols)), shape=(self.c_, self.nx),
                       dtype=float)  # .tocsr()
        B = np.array(self.B)
        if not square:
            logger.info("Using rectangular matrix, equality constraints are not used")
            return A, B
        AAT = A.T.dot(A)
        BT = A.T.dot(B)
        # add a small number to the matrix diagonal to smooth the results
        # can help speed up solving, but might also introduce some errors

        if self.eq_const_c_ > 0:
            logger.info("Equality block is %i x %i"%(self.eq_const_c_,self.nx))
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
            C = coo_matrix(
                (np.array(self.eq_const_C), (np.array(self.eq_const_row),
                                             np.array(self.eq_const_col))),
                shape=(self.eq_const_c_, self.nx))
            d = np.array(self.eq_const_d)
            AAT = bmat([[AAT, C.T], [C, None]])
            BT = np.hstack([BT, d])
        if damp:
            logger.info("Adding eps to matrix diagonal")
            AAT += eye(AAT.shape[0]) * np.finfo('float').eps
        return AAT, BT

    def _solve_lu(self, A, B):
        """
        Call scipy LU decomoposition

        Parameters
        ----------
        A : scipy square sparse matrix
        B : numpy vector

        Returns
        -------

        """
        lu = sla.splu(A.tocsc())
        sol = lu.solve(B)
        return sol[:self.nx]

    def _solve_lsqr(self, A, B, **kwargs):
        """
        Call scipy lsqr

        Parameters
        ----------
        A : rectangular sparse matrix
        B : vector

        Returns
        -------

        """

        lsqrargs = {}
        lsqrargs['btol'] = 1e-12
        lsqrargs['atol'] = 0
        if 'iter_lim' in kwargs:
            logger.info("Using %i maximum iterations" % kwargs['iter_lim'])
            lsqrargs['iter_lim'] = kwargs['iter_lim']
        if 'damp' in kwargs:
            logger.info("Using damping coefficient")
            lsqrargs['damp'] = kwargs['damp']
        if 'atol' in kwargs:
            logger.info('Using a tolerance of %f' % kwargs['atol'])
            lsqrargs['atol'] = kwargs['atol']
        if 'btol' in kwargs:
            logger.info('Using btol of %f' % kwargs['btol'])
            lsqrargs['btol'] = kwargs['btol']
        if 'show' in kwargs:
            lsqrargs['show'] = kwargs['show']
        if 'conlim' in kwargs:
            lsqrargs['conlim'] = kwargs['conlim']
        return sla.lsqr(A,B, **lsqrargs)[0]

    def _solve_chol(self, A, B):
        """
        Call suitesparse cholmod through scikitsparse
        LINUX ONLY!

        Parameters
        ----------
        A : scipy.sparse.matrix
            square sparse matrix
        B : numpy array
            RHS of equation

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
        A : scipy.sparse.matrix
            square sparse matrix
        B : numpy vector
        precon : scipy.sparse.matrix
            a preconditioner for the conjugate gradient system
        kwargs
            kwargs to pass to scipy solve e.g. atol, btol, callback etc

        Returns
        -------
        numpy array
        """
        cgargs = {}
        cgargs['tol'] = 1e-12
        cgargs['atol'] = 0
        if 'maxiter' in kwargs:
            logger.info("Using %i maximum iterations"%kwargs['maxiter'])
            cgargs['maxiter'] = kwargs['maxiter']
        if 'x0' in kwargs:
            logger.info("Using starting guess")
            cgargs['x0'] = kwargs['x0']
        if 'tol' in kwargs:
            logger.info('Using tolerance of %f'%kwargs['tol'])
            cgargs['tol'] = kwargs['tol']
        if 'atol' in kwargs:
            logger.info('Using atol of %f'%kwargs['atol'])
            cgargs['atol'] = kwargs['atol']
        if 'callback' in kwargs:
            cgargs['callback'] = kwargs['callback']
        if precon is not None:
            cgargs['M'] = precon(A)
        return sla.cg(A, B, **cgargs)[0][:self.nx]

    def _solve_pyamg(self, A, B, tol=1e-12,x0=None,**kwargs):
        """
        Solve least squares system using pyamg algorithmic multigrid solver

        Parameters
        ----------
        A :  scipy.sparse.matrix
        B : numpy array

        Returns
        -------

        """
        import pyamg
        logger.info("Solving using pyamg: tol {}".format(tol))
        return pyamg.solve(A, B, tol=tol, x0=x0, verb=False)[:self.nx]

    def _solve(self, solver='cg', **kwargs):
        """
        Main entry point to run the solver and update the node value
        attribute for the
        discreteinterpolator class

        Parameters
        ----------
        solver : string
            solver e.g. cg, lu, chol, custom
        kwargs
            kwargs for solver e.g. maxiter, preconditioner etc, damping for
        
        Returns
        -------
        bool
            True if the interpolation is run

        """
        logger.info("Solving interpolation for {}".format(self.propertyname))
        self.c = np.zeros(self.support.n_nodes)
        self.c[:] = np.nan
        damp = True
        if 'damp' in kwargs:
            damp = kwargs['damp']
        if solver == 'lu':
            logger.info("Forcing matrix damping for LU")
            damp = True
        if solver == 'lsqr':
            A, B =  self.build_matrix(False)
        else:
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
        if solver == 'pyamg':
            try:
                logger.info("Solving with pyamg solve")
                self.c[self.region] = self._solve_pyamg(A, B,**kwargs)
            except ImportError:
                logger.warn("Pyamg not installed using cg instead")
                self.c[self.region] = self._solve_cg(A, B)
        if solver == 'lsqr':
            self.c[self.region] = self._solve_lsqr(A, B, **kwargs)
        if solver == 'external':
            logger.warning("Using external solver")
            self.c[self.region] = kwargs['external'](A, B)[:self.nx]
        # check solution is not nan
        # self.support.properties[self.propertyname] = self.c
        if np.all(self.c == np.nan):
            logger.warning("Solver not run, no scalar field")
        # if solution is all 0, probably didn't work
        if np.all(self.c[self.region] == 0):
            logger.warning("No solution, {} scalar field 0. Add more data.".format(self.propertyname))

    def update(self):
        """
        Check if the solver is up to date, if not rerun interpolation using
        the previously used solver. If the interpolation has not been run
        before it will
        return False

        Returns
        -------
        bool

        """
        if self.solver is None:
            logging.debug("Cannot rerun interpolator")
            return False
        if not self.up_to_date:
            self.setup_interpolator()
            return self._solve(self.solver)

    def evaluate_value(self, evaluation_points):
        evaluation_points = np.array(evaluation_points)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(evaluation_points == np.nan, axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.support.evaluate_value(
                evaluation_points[~mask], self.c)
        return evaluated

    def evaluate_gradient(self, evaluation_points):
        """
        Evaluate the gradient of the scalar field at the evaluation points
        Parameters
        ----------
        evaluation_points : np.array
            xyz locations to evaluate the gradient

        Returns
        -------

        """
        if evaluation_points.shape[0] > 0:
            return self.support.evaluate_gradient(evaluation_points,
                                                  self.c)
        return np.zeros((0, 3))