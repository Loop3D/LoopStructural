from __future__ import annotations

from typing import ClassVar

import numpy as np

from loop_common.supports import SupportType, support_map


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
    _CELL_COUNT_SUPPORT_TYPES: ClassVar[frozenset[SupportType]] = frozenset(
        {
            SupportType.StructuredGrid,
            SupportType.TetMesh,
            SupportType.P2UnstructuredTetMesh,
        }
    )

    @staticmethod
    def create_support_from_bbox(
        support_type,
        bounding_box,
        nelements,
        element_volume=None,
        buffer: float | None = None,
        local_coordinates: bool = True,
    ):
        if isinstance(support_type, str):
            support_type = SupportType._member_map_[support_type].numerator
        if buffer is not None:
            bounding_box = bounding_box.with_buffer(buffer=buffer)
        if element_volume is not None:
            nelements = int(np.prod(bounding_box.length) / element_volume)
        if nelements is not None:
            bounding_box.nelements = nelements

        if local_coordinates:
            # Build the mesh in the bounding box's local (interpolation) frame
            # rather than raw world coordinates -- keeps node coordinates
            # numerically well-conditioned. The local box encloses every
            # projected corner, so a rotated frame gives a full support.
            origin, local_maximum = bounding_box.local_enclosure()
            # A rotated enclosure is longer than the world box: add cells so each
            # axis keeps the world cell size (an unrotated box keeps its nsteps).
            ratio = (local_maximum - origin) / (bounding_box.maximum - bounding_box.origin)
            nsteps = np.maximum(np.ceil(bounding_box.nsteps * ratio - 1e-9).astype(int), 1)
            step_vector = (local_maximum - origin) / nsteps
        else:
            origin = bounding_box.origin
            step_vector = bounding_box.step_vector
            nsteps = bounding_box.nsteps

        nsteps_kwarg = (
            "nsteps_cells" if support_type in SupportFactory._CELL_COUNT_SUPPORT_TYPES else "nsteps"
        )
        return support_map[support_type](
            origin=origin,
            step_vector=step_vector,
            **{nsteps_kwarg: nsteps},
        )
