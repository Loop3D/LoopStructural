import numpy as np

from LoopStructural.interpolators._discrete_interpolator import DiscreteInterpolator
from LoopStructural.interpolators._finite_difference_interpolator import FiniteDifferenceInterpolator
from ._p1interpolator import P1Interpolator
from typing import Optional, Union, Callable
from scipy import sparse
from LoopStructural.utils import rng

class ConstantNormInterpolator:
    """Adds a non linear constraint to an interpolator to constrain
    the norm of the gradient to be a set value.

    Returns
    -------
    _type_
        _description_
    """
    def __init__(self, interpolator: DiscreteInterpolator,basetype):
        """Initialise the constant norm inteprolator
        with a discrete interpolator.

        Parameters
        ----------
        interpolator : DiscreteInterpolator
            The discrete interpolator to add constant norm to.
        """
        self.basetype = basetype
        self.interpolator = interpolator
        self.support = interpolator.support
        self.random_subset = False
        self.norm_length = 1.0
        self.n_iterations = 20
        self.store_solution_history = False
        self.solution_history = []#np.zeros((self.n_iterations, self.support.n_nodes))
        self.gradient_constraint_store = []
    def add_constant_norm(self, w:float):
        """Add a constraint to the interpolator to constrain the norm of the gradient
        to be a set value

        Parameters
        ----------
        w : float
            weighting of the constraint
        """
        if "constant norm" in self.interpolator.constraints:
            _ = self.interpolator.constraints.pop("constant norm")
        
        element_indices = np.arange(self.support.elements.shape[0])
        if self.random_subset:
            rng.shuffle(element_indices)
            element_indices = element_indices[: int(0.1 * self.support.elements.shape[0])]
        vertices, gradient, elements, inside = self.support.get_element_gradient_for_location(
            self.support.barycentre[element_indices]
        )

        t_g = gradient[:, :, :]
        # t_n = gradient[self.support.shared_element_relationships[:, 1], :, :]
        v_t = np.einsum(
            "ijk,ik->ij",
            t_g,
            self.interpolator.c[self.support.elements[elements]],
        )

        v_t = v_t / np.linalg.norm(v_t, axis=1)[:, np.newaxis]
        self.gradient_constraint_store.append(np.hstack([self.support.barycentre[element_indices],v_t]))
        A1 = np.einsum("ij,ijk->ik", v_t, t_g)
        volume = self.support.element_size[element_indices]
        A1 = A1 / volume[:, np.newaxis]  # normalise by element size
        
        b = np.zeros(A1.shape[0]) + self.norm_length
        b = b / volume  # normalise by element size
        idc = np.hstack(
            [
                self.support.elements[elements],
            ]
        )
        self.interpolator.add_constraints_to_least_squares(A1, b, idc, w=w, name="constant norm")

    def solve_system(
        self,
        solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]] = None,
        tol: Optional[float] = None,
        solver_kwargs: dict = {},
    ) -> bool:
        """Solve the system of equations iteratively for the constant norm interpolator.

        Parameters
        ----------
        solver : Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]], optional
            Solver function or name, by default None
        tol : Optional[float], optional
            Tolerance for the solver, by default None
        solver_kwargs : dict, optional
            Additional arguments for the solver, by default {}

        Returns
        -------
        bool
            Success status of the solver
        """
        success = True
        for i in range(self.n_iterations):
            if i > 0:
                self.add_constant_norm(w=(0.1 * i) ** 2 + 0.01)
            # Ensure the interpolator is cast to P1Interpolator before calling solve_system
            if isinstance(self.interpolator, self.basetype):
                success = self.basetype.solve_system(self.interpolator, solver=solver, tol=tol, solver_kwargs=solver_kwargs)
                if self.store_solution_history:

                    self.solution_history.append(self.interpolator.c)
            else:
                raise TypeError("self.interpolator is not an instance of P1Interpolator")
            if not success:
                break
        return success

class ConstantNormP1Interpolator(P1Interpolator, ConstantNormInterpolator):
    """Constant norm interpolator using P1 base interpolator

    Parameters
    ----------
    P1Interpolator : class
        The P1Interpolator class.
    ConstantNormInterpolator : class
        The ConstantNormInterpolator class.
    """
    def __init__(self, support):
        """Initialise the constant norm P1 interpolator.

        Parameters
        ----------
        support : _type_
            _description_
        """
        P1Interpolator.__init__(self, support)
        ConstantNormInterpolator.__init__(self, self, P1Interpolator)

    def solve_system(
        self,
        solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]] = None,
        tol: Optional[float] = None,
        solver_kwargs: dict = {},
    ) -> bool:
        """Solve the system of equations for the constant norm P1 interpolator.

        Parameters
        ----------
        solver : Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]], optional
            Solver function or name, by default None
        tol : Optional[float], optional
            Tolerance for the solver, by default None
        solver_kwargs : dict, optional
            Additional arguments for the solver, by default {}

        Returns
        -------
        bool
            Success status of the solver
        """
        return ConstantNormInterpolator.solve_system(self, solver=solver, tol=tol, solver_kwargs=solver_kwargs)

class ConstantNormFDIInterpolator(FiniteDifferenceInterpolator, ConstantNormInterpolator):
    """Constant norm interpolator using finite difference base interpolator

    Parameters
    ----------
    FiniteDifferenceInterpolator : class
        The FiniteDifferenceInterpolator class.
    ConstantNormInterpolator : class
        The ConstantNormInterpolator class.
    """
    def __init__(self, support):
        """Initialise the constant norm finite difference interpolator.

        Parameters
        ----------
        support : _type_
            _description_
        """
        FiniteDifferenceInterpolator.__init__(self, support)
        ConstantNormInterpolator.__init__(self, self, FiniteDifferenceInterpolator)
    def solve_system(
        self,
        solver: Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]] = None,
        tol: Optional[float] = None,
        solver_kwargs: dict = {},
    ) -> bool:
        """Solve the system of equations for the constant norm finite difference interpolator.

        Parameters
        ----------
        solver : Optional[Union[Callable[[sparse.csr_matrix, np.ndarray], np.ndarray], str]], optional
            Solver function or name, by default None
        tol : Optional[float], optional
            Tolerance for the solver, by default None
        solver_kwargs : dict, optional
            Additional arguments for the solver, by default {}

        Returns
        -------
        bool
            Success status of the solver
        """
        return ConstantNormInterpolator.solve_system(self, solver=solver, tol=tol, solver_kwargs=solver_kwargs)