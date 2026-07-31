
from loop_common.base import LoopEntity, NumpyArray


class PointSet(LoopEntity):
    """A set of XYZ points representing a contact or fault trace."""

    coords: NumpyArray  # Shape (3,) or (N, 3)
