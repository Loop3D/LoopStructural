from typing import Optional, Union
import numpy as np
from LoopStructural.utils import getLogger

logger = getLogger(__name__)
try:
    from skimage.measure import marching_cubes
except ImportError:
    logger.warning("Using deprecated version of scikit-image")
    from skimage.measure import marching_cubes_lewiner as marching_cubes

from LoopStructural.interpolators import GeologicalInterpolator
from LoopStructural.utils import BoundingBox
from LoopStructural.datatypes import Surface

surface_list = dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]


class LoopIsosurfacer:
    def __init__(self, bounding_box: BoundingBox, interpolator: GeologicalInterpolator):
        self.bounding_box = bounding_box
        self.interpolator = interpolator

    def fit(self, values: Union[list, int, float]) -> surface_list:
        surfaces = {}
        all_values = self.interpolator.evaluate_value(self.bounding_box.regular_grid())
        if isinstance(values, list):
            isovalues = values
        elif isinstance(values, float):
            isovalues = [values]
        elif isinstance(values, int):
            isovalues = np.linspace(np.min(all_values), np.max(all_values), values)
        for isovalue in isovalues:
            verts, faces, normals, values = marching_cubes(
                all_values.reshape(self.bounding_box.nsteps, order="C"),
                isovalue,
                spacing=self.bounding_box.step_vector,
            )
            values = np.zeros(verts.shape[0]) + isovalue
            surfaces[f"surface_{isovalue}"] = Surface(
                vertices=verts,
                triangles=faces,
                normals=normals,
                name=f"surface_{isovalue}",
                values=values,
            )
        return surfaces
