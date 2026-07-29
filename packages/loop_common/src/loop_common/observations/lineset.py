from .orientation import Orientation, OrientationType
from .pointset import PointSet
import numpy as np
from typing import List
from loop_common.base import LoopEntity, NumpyArray


class LineSet(LoopEntity):
    """A set of lines representing geological features like faults or horizons."""

    vertices: NumpyArray  # Shape (N, 3) for N points along the line
    # Indices that mark the START of each new line segment
    offsets: NumpyArray  # Shape (M,) - e.g., [0, 5, 12]

    def to_tangent_vectors(self) -> List[Orientation]:
        """Compute tangent vectors for each line segment."""
        tangents = []
        for start, end in zip(self.offsets[:-1], self.offsets[1:]):
            # Compute tangent as the difference between consecutive points
            segment_tangents = self.vertices[start + 1 : end] - self.vertices[start : end - 1]
            segment_centres = (self.vertices[start : end - 1] + self.vertices[start + 1 : end]) / 2
            segment_tangents = segment_tangents / np.linalg.norm(
                segment_tangents, axis=1, keepdims=True
            )
            tangents.append(
                Orientation(
                    coords=segment_centres,
                    vector=segment_tangents,
                    magnitude=np.ones(segment_tangents.shape[0]),
                    polarity=np.ones(segment_tangents.shape[0]),
                    type=OrientationType.TANGENT,
                )
            )

        return tangents

    def to_point_set(self) -> PointSet:
        """Convert LineSet to PointSet by taking the vertices."""
        return PointSet(coords=self.vertices)
