from __future__ import annotations

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
    _CELL_COUNT_SUPPORT_TYPES = {
        SupportType.StructuredGrid,
        SupportType.TetMesh,
        SupportType.P2UnstructuredTetMesh,
    }

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
            # numerically well-conditioned. Project origin/maximum directly
            # (rather than all corners, which BoundingBox.corners only
            # supports in 3D) -- exact for translation-only transforms, which
            # is the only kind in use today (no code sets a non-identity
            # rotation on a bounding box).
            local_points = bounding_box.project(np.array([bounding_box.origin, bounding_box.maximum]))
            origin = np.min(local_points, axis=0)
            local_maximum = np.max(local_points, axis=0)
            step_vector = (local_maximum - origin) / bounding_box.nsteps
        else:
            origin = bounding_box.origin
            step_vector = bounding_box.step_vector

        nsteps_kwarg = (
            "nsteps_cells" if support_type in SupportFactory._CELL_COUNT_SUPPORT_TYPES else "nsteps"
        )
        return support_map[support_type](
            origin=origin,
            step_vector=step_vector,
            **{nsteps_kwarg: bounding_box.nsteps},
        )
