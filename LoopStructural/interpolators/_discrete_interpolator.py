"""
Discrete interpolator base for least squares
"""

from abc import abstractmethod
import logging

from time import time
import numpy as np
from scipy import sparse  # import sparse.coo_matrix, sparse.bmat, sparse.eye
from scipy.sparse import linalg as sla
from ..interpolators import InterpolatorType

from ..interpolators import GeologicalInterpolator
from ..utils import getLogger
from ..utils.exceptions import LoopImportError

logger = getLogger(__name__)


class DiscreteInterpolator(GeologicalInterpolator):
    """ """

    def __init__(self, support, data={}, c=None, up_to_date=False):
        """
        Base class for a discrete interpolator e.g. piecewise linear or finite difference which is
        any interpolator that solves the system using least squares approximation

        Parameters
        ----------
        support
            A discrete mesh with, nodes, elements, etc
        """
        GeologicalInterpolator.__init__(self, data=data, up_to_date=up_to_date)
        self.B = []
        self.support = support
        self.c = (
            np.array(c)
            if c is not None and np.array(c).shape[0] == self.support.n_nodes
            else np.zeros(self.support.n_nodes)
        )
        self.region_function = lambda xyz: np.ones(xyz.shape[0], dtype=bool)

        self.shape = "rectangular"
        if self.shape == "square":
            self.B = np.zeros(self.nx)
        self.c_ = 0

        self.solver = "cg"

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
        logger.info("Creating discrete interpolator with {} degrees of freedom".format(self.nx))
        self.type = InterpolatorType.BASE_DISCRETE
        self.c = np.zeros(self.support.n_nodes)

    @property
    def nx(self) -> int:
        """Number of degrees of freedom for the interpolator

        Returns
        -------
        int
            number of degrees of freedom, positve
        """
        return len(self.support.nodes[self.region])

    @property
    def region(self) -> np.ndarray:
        """The active region of the interpolator. A boolean
        mask for all elements that are interpolated

        Returns
        -------
        np.ndarray

        """

        return self.region_function(self.support.nodes).astype(bool)

    @property
    def region_map(self):
        region_map = np.zeros(self.support.n_nodes).astype(int)
        region_map[self.region] = np.array(range(0, len(region_map[self.region])))
        return region_map

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
        # self.region_function = region
        logger.info(
            "Cannot use region at the moment. Interpolation now uses region and has {} degrees of freedom".format(
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
        self.constraints = {}
        self.c_ = 0
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
        n_rows = A.shape[0]
        # logger.debug('Adding constraints to interpolator: {} {} {}'.format(A.shape[0]))
        # print(A.shape,B.shape,idc.shape)
        if A.shape != idc.shape:
            logger.error(f"Cannot add constraints: A and indexes have different shape : {name}")
            return

        if len(A.shape) > 2:
            n_rows = A.shape[0] * A.shape[1]
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
        if not isinstance(w, np.ndarray):
            raise BaseException("w must be a numpy array")

        if w.shape[0] != A.shape[0]:
            #     # make w the same size as A
            #     w = np.tile(w,(A.shape[1],1)).T
            # else:
            raise BaseException("Weight array does not match number of constraints")
        if np.any(np.isnan(idc)) or np.any(np.isnan(A)) or np.any(np.isnan(B)):
            logger.warning("Constraints contain nan not adding constraints: {}".format(name))
            # return
        rows = np.arange(0, n_rows).astype(int)
        base_name = name
        while name in self.constraints:
            count = 0
            if "_" in name:
                count = int(name.split("_")[1]) + 1
            name = base_name + "_{}".format(count)

        rows = np.tile(rows, (A.shape[-1], 1)).T
        self.constraints[name] = {
            'matrix': sparse.coo_matrix(
                (A.flatten(), (rows.flatten(), idc.flatten())), shape=(n_rows, self.nx)
            ).tocsc(),
            'b': B.flatten(),
            'w': w,
        }

    @abstractmethod
    def add_gradient_orthogonal_constraints(
        self, points: np.ndarray, vectors: np.ndarray, w: float = 1.0
    ):
        pass

    def calculate_residual_for_constraints(self):
        """Calculates Ax-B for all constraints added to the interpolator
        This could be a proxy to identify which constraints are controlling the model

        Returns
        -------
        np.ndarray
            vector of Ax-B
        """
        residuals = {}
        for constraint_name, constraint in self.constraints:
            residuals[constraint_name] = (
                np.einsum("ij,ij->i", constraint["A"], self.c[constraint["idc"].astype(int)])
                - constraint["B"].flatten()
            )
        return residuals

    def remove_constraints_from_least_squares(self, name="undefined", constraint_ids=None):
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

    def add_tangent_constraints(self, w=1.0):
        """Adds the constraints :math:`f(X)\cdotT=0`

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
        # max_weight = 0
        # for c in self.constraints.values():
        #     if len(c["w"]) == 0:
        #         continue
        #     if c["w"].max() > max_weight:
        #         max_weight = c["w"].max()
        # a = []
        # b = []
        # rows = []
        # cols = []
        # for c in self.constraints.values():
        #     if len(c["w"]) == 0:
        #         continue
        #     aa = (c["A"] * c["w"][:, None] / max_weight).flatten()
        #     b.extend((c["B"] * c["w"] / max_weight).tolist())
        #     mask = aa == 0
        #     a.extend(aa[~mask].tolist())
        #     rows.extend(c["row"].flatten()[~mask].tolist())
        #     cols.extend(c["col"].flatten()[~mask].tolist())
        mats = []
        bs = []
        for c in self.constraints.values():
            if len(c["w"]) == 0:
                continue
            mats.append(c['matrix'].multiply(c['w'][:, None]))
            bs.append(c['b'] * c['w'])
        A = sparse.vstack(mats)
        B = np.hstack(bs)

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

            C = sparse.coo_matrix(
                (np.array(a), (np.array(rows), cols)),
                shape=(self.eq_const_c, self.nx),
                dtype=float,
            ).tocsr()

            d = np.array(b)
            ATA = sparse.bmat([[ATA, C.T], [C, None]])
            ATB = np.hstack([ATB, d])

        if isinstance(damp, bool):
            if damp:
                damp = np.finfo("float").eps
            if not damp:
                damp = 0.0
        if isinstance(damp, float):
            logger.info("Adding eps to matrix diagonal")
            ATA += sparse.eye(ATA.shape[0]) * damp
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
            Aie = sparse.coo_matrix(
                (np.array(a), (np.array(rows), cols)),
                shape=(self.ineq_const_c, self.nx),
                dtype=float,
            ).tocsc()  # .tocsr()

            uie = np.array(u)
            lie = np.array(l)

            return ATA, ATB, Aie.T.dot(Aie), Aie.T.dot(uie), Aie.T.dot(lie)
        return ATA, ATB

    def _solve_osqp(self, P, A, q, l, u, mkl=False):
        """Wrapper to use osqp solver

        Parameters
        ----------
        P : _type_
            _description_
        A : _type_
            _description_
        q : _type_
            _description_
        l : _type_
            _description_
        u : _type_
            _description_
        mkl : bool, optional
            _description_, by default False

        Returns
        -------
        _type_
            _description_

        Raises
        ------
        LoopImportError
            _description_
        LoopImportError
            _description_
        """
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
                logger.error("MKL solver library path not correct. Please add to LD_LIBRARY_PATH")
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

    def solve_system(self, solver="cg", **kwargs):
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
        starttime = time()
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
        if np.all(self.c == np.nan):
            self.valid = False
            logger.warning("Solver not run, no scalar field")
            return
        # if solution is all 0, probably didn't work
        if np.all(self.c[self.region] == 0):
            self.valid = False
            logger.warning("No solution, scalar field 0. Add more data.")
            self.up_to_date = True
            return
        self.valid = True
        logging.info(f"Solving interpolation took: {time()-starttime}")
        self.up_to_date = True

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
            return self.solve_system(self.solver)

    def evaluate_value(self, locations: np.ndarray) -> np.ndarray:
        """Evaluate the value of the interpolator at location

        Parameters
        ----------
        evaluation_points : np.ndarray
            location to evaluate the interpolator

        Returns
        -------
        np.ndarray
            value of the interpolator
        """
        self.update()
        evaluation_points = np.array(locations)
        evaluated = np.zeros(evaluation_points.shape[0])
        mask = np.any(evaluation_points == np.nan, axis=1)

        if evaluation_points[~mask, :].shape[0] > 0:
            evaluated[~mask] = self.support.evaluate_value(evaluation_points[~mask], self.c)
        return evaluated

    def evaluate_gradient(self, locations: np.ndarray) -> np.ndarray:
        """
        Evaluate the gradient of the scalar field at the evaluation points
        Parameters
        ----------
        evaluation_points : np.array
            xyz locations to evaluate the gradient

        Returns
        -------

        """
        self.update()
        if locations.shape[0] > 0:
            return self.support.evaluate_gradient(locations, self.c)
        return np.zeros((0, 3))

    def to_dict(self):
        return {
            "type": self.type.name,
            "support": self.support.to_dict(),
            "c": self.c,
            **super().to_dict(),
            # 'region_function':self.region_function,
        }
