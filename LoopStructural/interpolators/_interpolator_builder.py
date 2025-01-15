from LoopStructural.interpolators import (
    GeologicalInterpolator,
    InterpolatorFactory,
    InterpolatorType,
)
from LoopStructural.datatypes import BoundingBox
from typing import Optional, Union
import numpy as np


class InterpolatorBuilder:
    def __init__(
        self,
        interpolatortype: Union[str, InterpolatorType],
        bounding_box: BoundingBox,
        nelements: int = 1000,
        buffer: float = 0.2,
        **kwargs,
    ):
        self.interpolatortype = interpolatortype
        self.bounding_box = bounding_box
        self.nelements = nelements
        self.buffer = buffer
        self.kwargs = kwargs
        self.interpolator : Optional[GeologicalInterpolator]= None

    def create_interpolator(self) -> 'InterpolatorBuilder':
        self.interpolator = InterpolatorFactory.create_interpolator(
            interpolatortype=self.interpolatortype,
            boundingbox=self.bounding_box,
            nelements=self.nelements,
            buffer=self.buffer,
            **self.kwargs,
        )
        return self

    def set_value_constraints(self, value_constraints: np.ndarray) -> 'InterpolatorBuilder':
        if self.interpolator:
            self.interpolator.set_value_constraints(value_constraints)
        return self

    def set_gradient_constraints(self, gradient_constraints: np.ndarray) -> 'InterpolatorBuilder':
        if self.interpolator:
            self.interpolator.set_gradient_constraints(gradient_constraints)
        return self

    def set_normal_constraints(self, normal_constraints: np.ndarray) -> 'InterpolatorBuilder':
        if self.interpolator:
            self.interpolator.set_normal_constraints(normal_constraints)
        return self

    def setup_interpolator(self, **kwargs) -> Optional[GeologicalInterpolator]:
        if self.interpolator:
            self.interpolator.setup(**kwargs)
            return self.interpolator
