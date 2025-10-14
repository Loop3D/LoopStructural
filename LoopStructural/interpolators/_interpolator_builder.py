from LoopStructural.interpolators import (
    InterpolatorFactory,
    InterpolatorType,
)
from LoopStructural.datatypes import BoundingBox
from typing import  Union, Optional
import numpy as np

from LoopStructural.interpolators._geological_interpolator import GeologicalInterpolator


class InterpolatorBuilder:
    def __init__(
        self,
        interpolatortype: Union[str, InterpolatorType],
        bounding_box: BoundingBox,
        nelements: Optional[int] = None,
        buffer: Optional[float] = None,
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
        self.bounding_box = bounding_box
        self.nelements = nelements
        self.buffer = buffer
        self.kwargs = kwargs
        self.interpolator = InterpolatorFactory.create_interpolator(
                    interpolatortype=self.interpolatortype,
                    boundingbox=self.bounding_box,
                    nelements=self.nelements,
                    buffer=self.buffer,
                    **self.kwargs,
                )

    def add_value_constraints(self, value_constraints: np.ndarray) -> 'InterpolatorBuilder':
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
        if self.interpolator:
            self.interpolator.set_value_constraints(value_constraints)
        return self

    def add_gradient_constraints(self, gradient_constraints: np.ndarray) -> 'InterpolatorBuilder':
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

        if self.interpolator:
            self.interpolator.set_gradient_constraints(gradient_constraints)
        return self

    def add_normal_constraints(self, normal_constraints: np.ndarray) -> 'InterpolatorBuilder':
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
        if self.interpolator:
            self.interpolator.set_normal_constraints(normal_constraints)
        return self
    def add_inequality_constraints(self, inequality_constraints: np.ndarray) -> 'InterpolatorBuilder':
        if self.interpolator:
            self.interpolator.set_value_inequality_constraints(inequality_constraints)
        return self
    def add_inequality_pair_constraints(self, inequality_pair_constraints: np.ndarray) -> 'InterpolatorBuilder':
        if self.interpolator:
            self.interpolator.set_inequality_pairs_constraints(inequality_pair_constraints)
        return self

    def setup_interpolator(self, **kwargs) -> 'InterpolatorBuilder':
        """This adds all of the constraints to the interpolator and
        sets the regularisation constraints

        Returns
        -------
        InterpolatorBuilder
            reference to the builder
        """
        if self.interpolator:
            self.interpolator.setup(**kwargs)
        return self

    def build(self)->GeologicalInterpolator:
        """Builds the interpolator and returns it

        Returns
        -------
        GeologicalInterpolator
            The interpolator fitting all of the constraints provided
        """
        return self.interpolator
