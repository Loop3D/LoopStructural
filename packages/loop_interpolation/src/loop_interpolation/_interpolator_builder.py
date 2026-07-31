"""Fluent builder for configuring interpolators.

This class intentionally stays thin: it delegates object creation to
InterpolatorFactory and provides a chainable API for adding constraints,
configuring setup options, and solving.
"""
from __future__ import annotations

import numpy as np
from loop_common.geometry import BoundingBox

from ._interpolatortype import InterpolatorType
from ._interpolator_factory import InterpolatorFactory


class InterpolatorBuilder:
    def __init__(
        self,
        interpolatortype: str | InterpolatorType = InterpolatorType.FINITE_DIFFERENCE,
        bounding_box: BoundingBox | None = None,
        nelements: int | None = None,
        buffer: float | None = None,
        **kwargs,
    ):
        """This class helps initialise and setup a geological interpolator.

        Parameters
        ----------
        interpolatortype : Union[str, InterpolatorType]
            type of interpolator
        bounding_box : BoundingBox
            bounding box of the area to interpolate
        nelements : int, optional
            degrees of freedom of the interpolator, by default 1000
        buffer : float, optional
            how much of a buffer around the bounding box should be used, by default 0.2
        """
        self.interpolatortype = interpolatortype
        if bounding_box is None:
            bounding_box = BoundingBox()
        self.bounding_box = bounding_box
        self.nelements = nelements
        self.buffer = buffer
        self.solver = kwargs.pop("solver", None)
        self.kwargs = kwargs
        self.setup_kwargs = {}
        self.solver_kwargs = {}
        self.interpolator = InterpolatorFactory.create_interpolator(
            interpolatortype=self.interpolatortype,
            boundingbox=self.bounding_box,
            nelements=self.nelements,
            buffer=self.buffer,
            solver=self.solver,
            **self.kwargs,
        )

    def use_solver(self, solver: str, **solver_kwargs) -> InterpolatorBuilder:
        """Configure the solver used when calling solve().

        Parameters
        ----------
        solver : str
            Solver name, e.g. ``"cg"``, ``"lsmr"``, or ``"admm"``.
        **solver_kwargs
            Keyword arguments forwarded to ``solve_system``.

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        self.solver = solver
        self.solver_kwargs = dict(solver_kwargs)
        if self.interpolator is not None:
            self.interpolator.solver = solver
        return self

    def solve(
        self,
        solver: str | None = None,
        tol: float | None = None,
        **solver_kwargs,
    ) -> InterpolatorBuilder:
        """Solve the configured interpolator system.

        Parameters
        ----------
        solver : Optional[str], optional
            Override solver name for this call.
        tol : Optional[float], optional
            Optional solver tolerance.
        **solver_kwargs
            Additional arguments passed to ``solve_system``.

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        if self.interpolator:
            selected_solver = solver if solver is not None else self.solver
            merged_solver_kwargs = dict(self.solver_kwargs)
            merged_solver_kwargs.update(solver_kwargs)
            self.interpolator.solve_system(
                solver=selected_solver,
                tol=tol,
                solver_kwargs=merged_solver_kwargs,
            )
        return self

    def use_regularisation_weight_scale(self, enabled: bool = True) -> InterpolatorBuilder:
        """Configure whether regularisation terms use spatial weighting.

        Parameters
        ----------
        enabled : bool, optional
            If True, enable regularisation-weight scaling in interpolators
            that support it.

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        self.setup_kwargs["use_regularisation_weight_scale"] = bool(enabled)
        return self

    def regularisation_weight_sigma(self, sigma: float) -> InterpolatorBuilder:
        """Configure spatial decay sigma for regularisation weight scaling.

        Parameters
        ----------
        sigma : float
            Gaussian decay width used by interpolators that support
            regularisation-weight scaling.

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        self.setup_kwargs["regularisation_weight_sigma"] = float(sigma)
        return self

    def _set_constraint(self, setter_name: str, values: np.ndarray) -> InterpolatorBuilder:
        """Forward constraint arrays to the underlying interpolator."""
        if self.interpolator:
            getattr(self.interpolator, setter_name)(values)
        return self

    def add_value_constraints(self, value_constraints: np.ndarray) -> InterpolatorBuilder:
        """Add value constraints to the interpolator

        Parameters
        ----------
        value_constraints : np.ndarray
            x,y,z,value of the constraints

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        return self._set_constraint("set_value_constraints", value_constraints)

    def add_gradient_constraints(self, gradient_constraints: np.ndarray) -> InterpolatorBuilder:
        """Add gradient constraints to the interpolator.

        Where g1 and g2 are two vectors that are orthogonal to the gradient:
        f'(X) · g1 = 0 and f'(X) · g2 = 0

        Parameters
        ----------
        gradient_constraints : np.ndarray
            Array with columns [x, y, z, gradient_x, gradient_y, gradient_z] of the constraints

        Returns
        -------
        bool
            True if constraints were added successfully
        """

        return self._set_constraint("set_gradient_constraints", gradient_constraints)

    def add_normal_constraints(self, normal_constraints: np.ndarray) -> InterpolatorBuilder:
        """Add normal constraints to the interpolator
        Where n is the normal vector to the surface
        $f'(X).dx = nx$
        $f'(X).dy = ny$
        $f'(X).dz = nz$
        Parameters
        ----------
        normal_constraints : np.ndarray
            x,y,z,nx,ny,nz of the constraints

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        return self._set_constraint("set_normal_constraints", normal_constraints)

    def add_tangent_constraints(self, tangent_constraints: np.ndarray) -> InterpolatorBuilder:
        """Add tangent constraints to the interpolator.

        Parameters
        ----------
        tangent_constraints : np.ndarray
            Array with columns [x, y, z, tx, ty, tz] where (tx, ty, tz) is a
            vector that lies *in* the surface (orthogonal to the gradient).

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        return self._set_constraint("set_tangent_constraints", tangent_constraints)

    def add_inequality_constraints(
        self, inequality_constraints: np.ndarray
    ) -> InterpolatorBuilder:
        return self._set_constraint("set_value_inequality_constraints", inequality_constraints)

    def add_inequality_pair_constraints(
        self, inequality_pair_constraints: np.ndarray
    ) -> InterpolatorBuilder:
        return self._set_constraint("set_inequality_pairs_constraints", inequality_pair_constraints)

    def setup_interpolator(self, **kwargs) -> InterpolatorBuilder:
        """This adds all of the constraints to the interpolator and
        sets the regularisation constraints

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        if self.interpolator:
            setup_kwargs = {**self.setup_kwargs, **kwargs}
            self.interpolator.setup(**setup_kwargs)
        return self

    def build(self) -> GeologicalInterpolator:
        """Builds the interpolator and returns it

        Returns
        -------
        GeologicalInterpolator
            The interpolator fitting all of the constraints provided
        """
        return self.interpolator
