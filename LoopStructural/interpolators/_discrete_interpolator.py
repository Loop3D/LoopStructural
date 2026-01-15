"""
Discrete interpolator base for least squares
"""

from abc import abstractmethod
from typing import Callable, Optional, Union
import logging

import numpy as np
from scipy import sparse  # import sparse.coo_matrix, sparse.bmat, sparse.eye
from ..interpolators import InterpolatorType

from ..interpolators import GeologicalInterpolator
from ..utils import getLogger

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
        self.dimensions = support.dimension
        self.c = (
            np.array(c)
            if c is not None and np.array(c).shape[0] == self.support.n_nodes
            else np.zeros(self.support.n_nodes)
        )
        self.region_function = lambda xyz: np.ones(xyz.shape[0], dtype=bool)

        self.shape = "rectangular"
        if self.shape == "square":
            self.B = np.zeros(self.dof)
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
        logger.info("Creating discrete interpolator with {} degrees of freedom".format(self.dof))
        self.type = InterpolatorType.BASE_DISCRETE
        self.apply_scaling_matrix = True
        self.add_ridge_regulatisation = True
        self.ridge_factor = 1e-8
    def set_nelements(self, nelements: int) -> int:
        return self.support.set_nelements(nelements)

    @property
    def n_elements(self) -> int:
        """Number of elements in the interpolator

        Returns
        -------
        int
            number of elements, positive
        """
        return self.support.n_elements

    @property
    def dof(self) -> int:
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
                self.dof
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

    def _pre_solve(self):
        """
        Pre solve function to be run before solving the interpolation
        """
        self.c = np.zeros(self.support.n_nodes)
        self.c[:] = np.nan
        return True

    def _post_solve(self):
        """Post solve function(s) to be run after the solver has been called"""
        self.clear_constraints()
        return True

    def clear_constraints(self):
        """
        Clear the constraints from the interpolator, this makes sure we are not storing
        the constraints after the solver has been run
        """
        self.constraints = {}
        self.ineq_constraints = {}
        self.equal_constraints = {}

    def reset(self):
        """
        Reset the interpolation constraints

        """
        self.constraints = {}
        self.c_ = 0
        self.regularisation_scale = np.ones(self.dof)
        logger.info("Resetting interpolation constraints")

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
        # Should this be done? It should make the solution more stable
        length = np.linalg.norm(A, axis=1)
        # length[length>0] = 1.
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
                (A.flatten(), (rows.flatten(), idc.flatten())), shape=(n_rows, self.dof)
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

    def add_inequality_constraints_to_matrix(
        self, A: np.ndarray, bounds: np.ndarray, idc: np.ndarray, name: str = "undefined"
    ):
        """Adds constraints for a matrix where the linear function
        l < Ax > u constrains the objective function


        Parameters
        ----------
        A : numpy array
            matrix of coefficients
        bounds : numpy array
            n*3 lower, upper, 1
        idc : numpy array
            index of constraints in the matrix
        Returns
        -------

        """
        # map from mesh node index to region node index
        gi = np.zeros(self.support.n_nodes, dtype=int)
        gi[:] = -1
        gi[self.region] = np.arange(0, self.dof, dtype=int)
        idc = gi[idc]
        rows = np.arange(0, idc.shape[0])
        rows = np.tile(rows, (A.shape[-1], 1)).T

        self.ineq_constraints[name] = {
            'matrix': sparse.coo_matrix(
                (A.flatten(), (rows.flatten(), idc.flatten())), shape=(rows.shape[0], self.dof)
            ).tocsc(),
            "bounds": bounds,
        }

    def add_value_inequality_constraints(self, w: float = 1.0):
        points = self.get_inequality_value_constraints()
        # check that we have added some points
        if points.shape[0] > 0:
            vertices, a, element, inside = self.support.get_element_for_location(points)
            rows = np.arange(0, points[inside, :].shape[0], dtype=int)
            rows = np.tile(rows, (a.shape[-1], 1)).T
            a = a[inside]
            cols = self.support.elements[element[inside]]
            self.add_inequality_constraints_to_matrix(a, points[:, 3:5], cols, 'inequality_value')

    def add_inequality_pairs_constraints(
        self,
        w: float = 1.0,
        upper_bound=np.finfo(float).eps,
        lower_bound=-np.inf,
        pairs: Optional[list] = None,
    ):

        points = self.get_inequality_pairs_constraints()
        if points.shape[0] > 0:
            # assemble a list of pairs in the model
            # this will make pairs even across stratigraphic boundaries
            # TODO add option to only add stratigraphic pairs
            if not pairs:
                pairs = {}
                k = 0
                for i in np.unique(points[:, self.support.dimension]):
                    for j in np.unique(points[:, self.support.dimension]):
                        if i == j:
                            continue
                        if tuple(sorted([i, j])) not in pairs:
                            pairs[tuple(sorted([i, j]))] = k
                            k += 1
                pairs = list(pairs.keys())
            for pair in pairs:
                upper_points = points[points[:, self.support.dimension] == pair[0]]
                lower_points = points[points[:, self.support.dimension] == pair[1]]

                upper_interpolation = self.support.get_element_for_location(upper_points)
                lower_interpolation = self.support.get_element_for_location(lower_points)
                if (~upper_interpolation[3]).sum() > 0:
                    logger.warning(
                        f'Upper points not in mesh {upper_points[~upper_interpolation[3]]}'
                    )
                if (~lower_interpolation[3]).sum() > 0:
                    logger.warning(
                        f'Lower points not in mesh {lower_points[~lower_interpolation[3]]}'
                    )
                ij = np.array(
                    [
                        *np.meshgrid(
                            np.arange(0, int(upper_interpolation[3].sum()), dtype=int),
                            np.arange(0, int(lower_interpolation[3].sum()), dtype=int),
                        )
                    ],
                    dtype=int,
                )

                ij = ij.reshape(2, -1).T
                rows = np.arange(0, ij.shape[0], dtype=int)
                rows = np.tile(rows, (upper_interpolation[1].shape[-1], 1)).T
                rows = np.hstack([rows, rows])
                a = upper_interpolation[1][upper_interpolation[3]][ij[:, 0]]
                a = np.hstack([a, -lower_interpolation[1][lower_interpolation[3]][ij[:, 1]]])
                cols = np.hstack(
                    [
                        self.support.elements[
                            upper_interpolation[2][upper_interpolation[3]][ij[:, 0]]
                        ],
                        self.support.elements[
                            lower_interpolation[2][lower_interpolation[3]][ij[:, 1]]
                        ],
                    ]
                )

                bounds = np.zeros((ij.shape[0], 2))
                bounds[:, 0] = lower_bound
                bounds[:, 1] = upper_bound

                self.add_inequality_constraints_to_matrix(
                    a, bounds, cols, f'inequality_pairs_{pair[0]}_{pair[1]}'
                )

    def add_inequality_feature(
        self,
        feature: Callable[[np.ndarray], np.ndarray],
        lower: bool = True,
        mask: Optional[np.ndarray] = None,
    ):
        """Add an inequality constraint to the interpolator using an existing feature.
        This will make the interpolator greater than or less than the exising feature.
        Evaluate the feature at the interpolation nodes.
        Can provide a boolean mask to restrict to only some parts

        Parameters
        ----------
        feature : BaseFeature
            the feature that will be used to constraint the interpolator
        lower : bool, optional
            lower or upper constraint, by default True
        mask : np.ndarray, optional
            restrict the nodes to evaluate on, by default None
        """
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
        if not lower:
            u[mask] = value[mask]

        self.add_inequality_constraints_to_matrix(
            np.ones((value.shape[0], 1)),
            l,
            u,
            np.arange(0, self.dof, dtype=int),
        )

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
        gi[self.region] = np.arange(0, self.dof)
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

    def build_matrix(self):
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

        mats = []
        bs = []
        for c in self.constraints.values():
            if len(c["w"]) == 0:
                continue
            mats.append(c['matrix'].multiply(c['w'][:, None]))
            bs.append(c['b'] * c['w'])
        A = sparse.vstack(mats)
        logger.info(f"Interpolation matrix is {A.shape[0]} x {A.shape[1]}")

        B = np.hstack(bs)
        return A, B
    def compute_column_scaling_matrix(self, A: sparse.csr_matrix) -> sparse.dia_matrix:
        """Compute column scaling matrix S for matrix A so that A @ S has columns with unit norm.

        Parameters
        ----------
        A : sparse.csr_matrix
            interpolation matrix

        Returns
        -------
        scipy.sparse.dia_matrix
            diagonal scaling matrix S
        """
        col_norms = sparse.linalg.norm(A, axis=0)
        scaling_factors = np.ones(A.shape[1])
        mask = col_norms > 0
        scaling_factors[mask] = 1.0 / col_norms[mask]
        S = sparse.diags(scaling_factors)
        return S

    def add_equality_block(self, A, B):
        if len(self.equal_constraints) > 0:
            ATA = A.T.dot(A)
            ATB = A.T.dot(B)
            logger.info(f"Equality block is {self.eq_const_c} x {self.dof}")
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
                shape=(self.eq_const_c, self.dof),
                dtype=float,
            ).tocsr()

            d = np.array(b)
            ATA = sparse.bmat([[ATA, C.T], [C, None]])
            ATB = np.hstack([ATB, d])

            return ATA, ATB

    def build_inequality_matrix(self):
        mats = []
        bounds = []
        for c in self.ineq_constraints.values():
            mats.append(c['matrix'])
            bounds.append(c['bounds'])
        if len(mats) == 0:
            return sparse.csr_matrix((0, self.dof), dtype=float), np.zeros((0, 3))
        Q = sparse.vstack(mats)
        bounds = np.vstack(bounds)
        return Q, bounds

    def solve_system(
        self,
        solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]] = None,
        tol: Optional[float] = None,
        solver_kwargs: dict = {},
    ) -> bool:
        """
        Main entry point to run the solver and update the node value
        attribute for the
        discreteinterpolator class

        Parameters
        ----------
        solver : string/callable
            solver 'cg' conjugate gradient, 'lsmr' or callable function
        solver_kwargs
            kwargs for solver check scipy documentation for more information

        Returns
        -------
        bool
            True if the interpolation is run

        """
        if not self._pre_solve():
            raise ValueError("Pre solve failed")

        A, b = self.build_matrix()
        if self.add_ridge_regulatisation:
            ridge = sparse.eye(A.shape[1]) * self.ridge_factor
            A = sparse.vstack([A, ridge])
            b = np.hstack([b, np.zeros(A.shape[1])])
            logger.info("Adding ridge regularisation to interpolation matrix")
        if self.apply_scaling_matrix:
            S = self.compute_column_scaling_matrix(A)
            A = A @ S

        Q, bounds = self.build_inequality_matrix()
        if callable(solver):
            logger.warning('Using custom solver')
            self.c = solver(A.tocsr(), b)
            self.up_to_date = True
        elif isinstance(solver, str) or solver is None:
            if solver not in ['cg', 'lsmr', 'admm']:
                logger.warning(
                    f'Unknown solver {solver} using cg. \n Available solvers are cg and lsmr or a custom solver as a callable function'
                )
                solver = 'cg'
        if solver == 'cg':
            logger.info("Solving using cg")
            if 'atol' not in solver_kwargs or 'rtol' not in solver_kwargs:
                if tol is not None:
                    solver_kwargs['atol'] = tol

            logger.info(f"Solver kwargs: {solver_kwargs}")

            res = sparse.linalg.cg(A.T @ A, A.T @ b, **solver_kwargs)
            if res[1] > 0:
                logger.warning(
                    f'CG reached iteration limit ({res[1]})and did not converge, check input data. Setting solution to last iteration'
                )
            self.c = res[0]
            self.up_to_date = True

        elif solver == 'lsmr':
            logger.info("Solving using lsmr")
            # if 'atol' not in solver_kwargs:
            #     if tol is not None:
            #         solver_kwargs['atol'] = tol
            if 'btol' not in solver_kwargs:
                if tol is not None:
                    solver_kwargs['btol'] = tol
                    solver_kwargs['atol'] = 0.
                    logger.info(f"Setting lsmr btol to {tol}")
            logger.info(f"Solver kwargs: {solver_kwargs}")
            res = sparse.linalg.lsmr(A, b, **solver_kwargs)
            if res[1] == 1 or res[1] == 4 or res[1] == 2 or res[1] == 5:
                self.c = res[0]
            elif res[1] == 0:
                logger.warning("Solution to least squares problem is all zeros, check input data")
            elif res[1] == 3 or res[1] == 6:
                logger.warning("COND(A) seems to be greater than CONLIM, check input data")
                # self.c = res[0]
            elif res[1] == 7:
                logger.warning(
                    "LSMR reached iteration limit and did not converge, check input data. Setting solution to last iteration"
                )
                self.c = res[0]
            self.up_to_date = True

        elif solver == 'admm':
            logger.info("Solving using admm")

            if 'x0' in solver_kwargs:
                x0 = solver_kwargs['x0'](self.support)
            else:
                x0 = np.zeros(A.shape[1])
            solver_kwargs.pop('x0', None)
            if Q is None:
                logger.warning("No inequality constraints, using lsmr")
                return self.solve_system('lsmr', solver_kwargs)

            try:
                from loopsolver import admm_solve

                try:
                    linsys_solver = solver_kwargs.pop('linsys_solver', 'lsmr')
                    res = admm_solve(
                        A,
                        b,
                        Q,
                        bounds,
                        x0=x0,
                        admm_weight=solver_kwargs.pop('admm_weight', 0.01),
                        nmajor=solver_kwargs.pop('nmajor', 200),
                        linsys_solver_kwargs=solver_kwargs,
                        linsys_solver=linsys_solver,
                    )
                    self.c = res
                    self.up_to_date = True
                except ValueError as e:
                    logger.error(f"ADMM solver failed: {e}")
                    self.up_to_date = False
            except ImportError:
                logger.warning(
                    "Cannot import admm solver. Please install loopsolver or use lsmr or cg"
                )
                self.up_to_date = False
        else:
            logger.error(f"Unknown solver {solver}")
            self.up_to_date = False
        # self._post_solve()
        # apply scaling matrix to solution
        if self.apply_scaling_matrix:
            self.c = S @ self.c
        return self.up_to_date

    def update(self) -> bool:
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
            self.up_to_date = self.solve_system(self.solver)
            return self.up_to_date

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

    def vtk(self):
        if self.up_to_date is False:
            self.update()
        return self.support.vtk({'c': self.c})
