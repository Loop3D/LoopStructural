from loop_common.supports import support_map, SupportType
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

    # Support types whose constructor takes nsteps as a *cell* count
    # (translated internally to a node count via BaseStructuredSupport).
    _CELL_COUNT_SUPPORT_TYPES = {
        SupportType.StructuredGrid,
        SupportType.TetMesh,
        SupportType.P2UnstructuredTetMesh,
    }

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

        nsteps_kwarg = (
            "nsteps_cells" if support_type in SupportFactory._CELL_COUNT_SUPPORT_TYPES else "nsteps"
        )
        return support_map[support_type](
            origin=bounding_box.origin,
            step_vector=bounding_box.step_vector,
            **{nsteps_kwarg: bounding_box.nsteps},
        )
