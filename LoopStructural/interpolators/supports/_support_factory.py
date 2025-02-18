from LoopStructural.interpolators.supports import support_map, SupportType
import numpy as np
from typing import Optional


class SupportFactory:
    @staticmethod
    def create_support(support_type, **kwargs):
        if support_type is None:
            raise ValueError("No support type specified")
        if isinstance(support_type, str):
            support_type = SupportType._member_map_[support_type].numerator
        return support_map[support_type](**kwargs)

    @staticmethod
    def from_dict(d):
        d = d.copy()
        support_type = d.pop("type", None)
        if support_type is None:
            raise ValueError("No support type specified")
        return SupportFactory.create_support(support_type, **d)

    @staticmethod
    def create_support_from_bbox(
        support_type, bounding_box, nelements, element_volume=None, buffer: Optional[float] = None
    ):
        if isinstance(support_type, str):
            support_type = SupportType._member_map_[support_type].numerator
        if buffer is not None:
            bounding_box = bounding_box.with_buffer(buffer=buffer)
        if element_volume is not None:
            nelements = int(np.prod(bounding_box.length) / element_volume)
        if nelements is not None:
            bounding_box.nelements = nelements

        return support_map[support_type](
            origin=bounding_box.origin,
            step_vector=bounding_box.step_vector,
            nsteps=bounding_box.nsteps,
        )
