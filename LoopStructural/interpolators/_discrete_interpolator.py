"""
Discrete interpolator base for least squares
"""
import logging

import numpy as np
from scipy.sparse import coo_matrix, bmat, eye
from scipy.sparse import linalg as sla
from scipy.sparse.linalg import norm
from sklearn.preprocessing import normalize
from LoopStructural.interpolators import InterpolatorType

from LoopStructural.interpolators import GeologicalInterpolator
from LoopStructural.utils import getLogger
from LoopStructural.utils.exceptions import LoopImportError

logger = getLogger(__name__)

from ._geological_interpolator import GeologicalInterpolator


class DiscreteInterpolator(GeologicalInterpolator):
    """ """

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
        self.region_function = lambda xyz: np.ones(xyz.shape[0], dtype=bool)
        # self.region_map[self.region] = np.array(range(0,
        # len(self.region_map[self.region])))
        self.shape = "rectangular"
        if self.shape == "square":
            self.B = np.zeros(self.nx)
        self.c_ = 0
        # self.A = []  # sparse matrix storage coo format
        # self.col = []
        # self.row = []  # sparse matrix storage
        # self.w = []
        self.solver = None

        self.eq_const_C = []
        self.eq_const_row = []
        self.eq_const_col = []
        self.eq_const_d = []

        self.equal_constraints = {}
        self.eq_const_c = 0
        self.ineq_constraints = {}
        self.ineq_const_c = 0

        self.non_linear_constraints = []
        self.constraints = {}
        self.interpolation_weights = {}
        logger.info(
            "Creating discrete interpolator with {} degrees of freedom".format(self.nx)
        )
        self.type = InterpolatorType.BASE_DISCRETE

    @property
    def nx(self):
        return len(self.support.nodes[self.region])

    @property
    def region(self):
        return self.region_function(self.support.nodes).astype(bool)

    @property
    def region_map(self):
        region_map = np.zeros(self.support.n_nodes).astype(int)
        region_map[self.region] = np.array(range(0, len(region_map[self.region])))
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
        logger.info(
            "Interpolation now uses region and has {} degrees of freedom".format(
                self.nx
            )
        )

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

    def add_constraints_to_least_squares(self, A, B, idc, w=1.0, name="undefined"):
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
        # logger.debug('Adding constraints to interpolator: {} {} {}'.format(A.shape[0]))
        # print(A.shape,B.shape,idc.shape)
        if A.shape != idc.shape:
            logger.error(
                f"Cannot add constraints: A and indexes have different shape : {name}"
            )
            return

        if len(A.shape) > 2:
            nr = A.shape[0] * A.shape[1]
            if isinstance(w, np.ndarray):
                w = np.tile(w, (A.shape[1]))
            A = A.reshape((A.shape[0] * A.shape[1], A.shape[2]))
            idc = idc.reshape((idc.shape[0] * idc.shape[1], idc.shape[2]))
            B = B.reshape((A.shape[0]))
            # w = w.reshape((A.shape[0]))
        # normalise by rows of A
        length = np.linalg.norm(A, axis=1)  # .getcol(0).norm()
        B[length > 0] /= length[length > 0]
        # going to assume if any are nan they are all nan
        mask = np.any(np.isnan(A), axis=1)
        A[mask, :] = 0
        A[length > 0, :] /= length[length > 0, None]
        if isinstance(w, (float, int)):
            w = np.ones(A.shape[0]) * w
        if isinstance(w, np.ndarray) == False:
            raise BaseException("w must be a numpy array")

        if w.shape[0] != A.shape[0]:
            #     # make w the same size as A
            #     w = np.tile(w,(A.shape[1],1)).T
            # else:
            raise BaseException("Weight array does not match number of constraints")
        if np.any(np.isnan(idc)) or np.any(np.isnan(A)) or np.any(np.isnan(B)):
            logger.warning(
                "Constraints contain nan not adding constraints: {}".format(name)
            )
            # return
        rows = np.arange(0, nr).astype(int)
        rows += self.c_
        constraint_ids = rows.copy()
        base_name = name
        while name in self.constraints:
            count = 0
            if "_" in name:
                count = int(name.split("_")[1]) + 1
            name = base_name + "_{}".format(count)

            # self.constraints[name]['A'] =  A#np.vstack([self.constraints[name]['A'],A])
            # self.constraints[name]['B'] =  B#np.hstack([self.constraints[name]['B'], B])
            # self.constraints[name]['idc'] = idc#np.vstack([self.constraints[name]['idc'],
            #                                     idc])
        rows = np.tile(rows, (A.shape[-1], 1)).T
        self.constraints[name] = {
            "node_indexes": constraint_ids,
            "A": A,
            "B": B.flatten(),
            "col": idc,
            "w": w,
            "row": rows,
        }

        self.c_ += nr

    def calculate_residual_for_constraints(self):
        residuals = {}
        for constraint_name, constraint in self.constraints:
            residuals[constraint_name] = (
                np.einsum(
                    "ij,ij->i", constraint["A"], self.c[constraint["idc"].astype(int)]
                )
                - constraint["B"].flatten()
            )
        return residuals

    def remove_constraints_from_least_squares(
        self, name="undefined", constraint_ids=None
    ):
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

        pass

    def add_equality_constraints(self, node_idx, values, name="undefined"):
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

        self.equal_constraints[name] = {
            "A": np.ones(idc[outside].shape[0]),
            "B": values[outside],
            "col": idc[outside],
            # "w": w,
            "row": np.arange(self.eq_const_c, self.eq_const_c + idc[outside].shape[0]),
        }
        self.eq_const_c += idc[outside].shape[0]

    def add_non_linear_constraints(self, nonlinear_constraint):
        self.non_linear_constraints.append(nonlinear_constraint)

    def add_inequality_constraints_to_matrix(self, A, l, u, idc, name="undefined"):
        """Adds constraints for a matrix where the linear function
        l < Ax > u constrains the objective function


        Parameters
        ----------
        A : numpy array
            matrix of coefficients
        l : numpy array
            lower bounds
        u : numpy array
            upper bounds
        idc : numpy array
            index of constraints
        Returns
        -------

        """
        # map from mesh node index to region node index
        gi = np.zeros(self.support.n_nodes, dtype=int)
        gi[:] = -1
        gi[self.region] = np.arange(0, self.nx, dtype=int)
        idc = gi[idc]
        rows = np.arange(self.ineq_const_c, self.ineq_const_c + idc.shape[0])
        rows = np.tile(rows, (A.shape[-1], 1)).T
        self.ineq_constraints[name] = {"A": A, "l": l, "col": idc, "u": u, "row": rows}
        self.ineq_const_c += idc.shape[0]

    def add_inequality_feature(self, feature, lower=True, mask=None):

        # add inequality value for the nodes of the mesh
        # flag lower determines whether the feature is a lower bound or upper bound
        # mask is just a boolean array determining which nodes to apply it to

        value = feature(self.support.nodes)
        if mask is None:
            mask = np.ones(value.shape[0], dtype=bool)
        l = np.zeros(value.shape[0]) - np.inf
        u = np.zeros(value.shape[0]) + np.inf
        mask = np.logical_and(mask, ~np.isnan(value))
        if lower:
            l[mask] = value[mask]
        if lower == False:
            u[mask] = value[mask]

        self.add_inequality_constraints_to_matrix(
            np.ones((value.shape[0], 1)),
            l,
            u,
            np.arange(0, self.nx, dtype=int),
        )

    def add_tangent_constraints(self, w=1.0):
        """

        Parameters
        ----------
        w : double


        Returns
        -------

        """
        points = self.get_tangent_constraints()
        if points.shape[0] > 1:
            self.add_gradient_orthogonal_constraints(points[:, :3], points[:, 3:6], w)

    def build_matrix(self, square=True, damp=0.0, ie=False):
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

        logger.info("Interpolation matrix is %i x %i" % (self.c_, self.nx))
        # To keep the solvers consistent for different model scales the range of the constraints should be similar.
        # We normalise the row vectors for the interpolation matrix
        # Each constraint can then be weighted separately for the least squares problem
        # The weights are normalised so that the max weight is 1.0
        # This means that the tolerance and other parameters for the solver
        # are kept the same between iterations.
        # #TODO currently the element size is not incorporated into the weighting.
        # For cartesian grids this is probably ok but for tetrahedron could be more problematic if
        # the tetras have different volumes. Would expect for the size of the element to influence
        # how much it contributes to the system.
        # It could be implemented by multiplying the weight array by the element size.
        # I am not sure how to integrate regularisation into this framework as my gut feeling is the regularisation
        # should be weighted by the area of the element face and not element volume, but this means the weight decreases with model scale
        # which is not ideal.
        max_weight = 0
        for c in self.constraints.values():
            if len(c["w"]) == 0:
                continue
            if c["w"].max() > max_weight:
                max_weight = c["w"].max()
        a = []
        b = []
        rows = []
        cols = []
        for c in self.constraints.values():
            if len(c["w"]) == 0:
                continue
            aa = (c["A"] * c["w"][:, None] / max_weight).flatten()
            b.extend((c["B"] * c["w"] / max_weight).tolist())
            mask = aa == 0
            a.extend(aa[~mask].tolist())
            rows.extend(c["row"].flatten()[~mask].tolist())
            cols.extend(c["col"].flatten()[~mask].tolist())

        A = coo_matrix(
            (np.array(a), (np.array(rows), cols)), shape=(self.c_, self.nx), dtype=float
        ).tocsc()  # .tocsr()

        B = np.array(b)
        if not square:
            logger.info("Using rectangular matrix, equality constraints are not used")
            return A, B
        ATA = A.T.dot(A)
        ATB = A.T.dot(B)
        # add a small number to the matrix diagonal to smooth the results
        # can help speed up solving, but might also introduce some errors

        if len(self.equal_constraints) > 0:
            logger.info(f"Equality block is {self.eq_const_c} x {self.nx}")
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
            nc = 0
            a = []
            rows = []
            cols = []
            b = []
            for c in self.equal_constraints.values():
                b.extend((c["B"]).tolist())
                aa = c["A"].flatten()
                mask = aa == 0
                a.extend(aa[~mask].tolist())
                rows.extend(c["row"].flatten()[~mask].tolist())
                cols.extend(c["col"].flatten()[~mask].tolist())

            C = coo_matrix(
                (np.array(a), (np.array(rows), cols)),
                shape=(self.eq_const_c, self.nx),
                dtype=float,
            ).tocsr()

            d = np.array(b)
            ATA = bmat([[ATA, C.T], [C, None]])
            ATB = np.hstack([ATB, d])

        if isinstance(damp, bool):
            if damp == True:
                damp = np.finfo("float").eps
            if damp == False:
                damp = 0.0
        if isinstance(damp, float):
            logger.info("Adding eps to matrix diagonal")
            ATA += eye(ATA.shape[0]) * damp
        if len(self.ineq_constraints) > 0 and ie:
            print("using inequality constraints")
            a = []
            l = []
            u = []
            rows = []
            cols = []
            for c in self.ineq_constraints.values():
                aa = (c["A"]).flatten()
                l.extend((c["l"]).tolist())
                u.extend((c["u"]).tolist())

                mask = aa == 0
                a.extend(aa[~mask].tolist())
                rows.extend(c["row"].flatten()[~mask].tolist())
                cols.extend(c["col"].flatten()[~mask].tolist())
            Aie = coo_matrix(
                (np.array(a), (np.array(rows), cols)),
                shape=(self.ineq_const_c, self.nx),
                dtype=float,
            ).tocsc()  # .tocsr()

            uie = np.array(u)
            lie = np.array(l)

            return ATA, ATB, Aie.T.dot(Aie), Aie.T.dot(uie), Aie.T.dot(lie)
        return ATA, ATB

    def _solve_osqp(self, P, A, q, l, u, mkl=False):

        try:
            import osqp
        except ImportError:
            raise LoopImportError("Missing osqp pip install osqp")
        prob = osqp.OSQP()

        # Setup workspace
        # osqp likes csc matrices
        linsys_solver = "qdldl"
        if mkl:
            linsys_solver = "mkl pardiso"

        try:
            prob.setup(
                P.tocsc(),
                np.array(q),
                A.tocsc(),
                np.array(u),
                np.array(l),
                linsys_solver=linsys_solver,
            )
        except ValueError:
            if mkl:
                logger.error(
                    "MKL solver library path not correct. Please add to LD_LIBRARY_PATH"
                )
                raise LoopImportError("Cannot import MKL pardiso")
        res = prob.solve()
        return res.x

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
        return sol[: self.nx]

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
        lsqrargs["btol"] = 1e-12
        lsqrargs["atol"] = 0
        if "iter_lim" in kwargs:
            logger.info("Using %i maximum iterations" % kwargs["iter_lim"])
            lsqrargs["iter_lim"] = kwargs["iter_lim"]
        if "damp" in kwargs:
            logger.info("Using damping coefficient")
            lsqrargs["damp"] = kwargs["damp"]
        if "atol" in kwargs:
            logger.info("Using a tolerance of %f" % kwargs["atol"])
            lsqrargs["atol"] = kwargs["atol"]
        if "btol" in kwargs:
            logger.info("Using btol of %f" % kwargs["btol"])
            lsqrargs["btol"] = kwargs["btol"]
        if "show" in kwargs:
            lsqrargs["show"] = kwargs["show"]
        if "conlim" in kwargs:
            lsqrargs["conlim"] = kwargs["conlim"]
        return sla.lsqr(A, B, **lsqrargs)[0]

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
            return factor(B)[: self.nx]
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
        cgargs["tol"] = 1e-12
        cgargs["atol"] = 1e-10
        if "maxiter" in kwargs:
            logger.info("Using %i maximum iterations" % kwargs["maxiter"])
            cgargs["maxiter"] = kwargs["maxiter"]
        if "x0" in kwargs:
            logger.info("Using starting guess")
            cgargs["x0"] = kwargs["x0"]
        if "tol" in kwargs:
            logger.info("Using tolerance of %f" % kwargs["tol"])
            cgargs["tol"] = kwargs["tol"]
        if "atol" in kwargs:
            logger.info("Using atol of %f" % kwargs["atol"])
            cgargs["atol"] = kwargs["atol"]
        if "callback" in kwargs:
            cgargs["callback"] = kwargs["callback"]
        if precon is not None:
            cgargs["M"] = precon(A)
        return sla.cg(A, B, **cgargs)[0][: self.nx]

    def _solve_pyamg(self, A, B, tol=1e-12, x0=None, verb=False, **kwargs):
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
        return pyamg.solve(A, B, tol=tol, x0=x0, verb=verb)[: self.nx]

    def _solve(self, solver="cg", **kwargs):
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
        if "damp" in kwargs:
            damp = kwargs["damp"]
        if solver == "lu":
            logger.info("Forcing matrix damping for LU")
            damp = True
        if solver == "lsqr":
            A, B = self.build_matrix(False)
        elif solver == "osqp":
            P, q, A, l, u = self.build_matrix(True, ie=True)
        else:
            A, B = self.build_matrix(damp=damp)
        # run the chosen solver
        if solver == "cg":
            logger.info("Solving using conjugate gradient")
            self.c[self.region] = self._solve_cg(A, B, **kwargs)
        if solver == "chol":
            self.c[self.region] = self._solve_chol(A, B)
        if solver == "lu":
            logger.info("Solving using scipy LU")
            self.c[self.region] = self._solve_lu(A, B)
        if solver == "pyamg":
            try:
                logger.info("Solving with pyamg solve")
                self.c[self.region] = self._solve_pyamg(A, B, **kwargs)
            except ImportError:
                logger.warn("Pyamg not installed using cg instead")
                self.c[self.region] = self._solve_cg(A, B)
        if solver == "lsqr":
            self.c[self.region] = self._solve_lsqr(A, B, **kwargs)
        if solver == "external":
            logger.warning("Using external solver")
            self.c[self.region] = kwargs["external"](A, B)[: self.nx]
        if solver == "osqp":
            self.c[self.region] = self._solve_osqp(
                P, A, q, l, u, mkl=kwargs.get("mkl", False)
            )  # , **kwargs)
        # check solution is not nan
        # self.support.properties[self.propertyname] = self.c
        if np.all(self.c == np.nan):
            self.valid = False
            logger.warning("Solver not run, no scalar field")
            return
        # if solution is all 0, probably didn't work
        if np.all(self.c[self.region] == 0):
            self.valid = False
            logger.warning(
                "No solution, {} scalar field 0. Add more data.".format(
                    self.propertyname
                )
            )
            return
        self.valid = True

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
                evaluation_points[~mask], self.c
            )
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
            return self.support.evaluate_gradient(evaluation_points, self.c)
        return np.zeros((0, 3))
