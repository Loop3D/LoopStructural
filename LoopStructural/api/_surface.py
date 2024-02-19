from typing import Optional, Union, Callable
import numpy as np
import numpy.typing as npt
from LoopStructural.utils import getLogger

logger = getLogger(__name__)
try:
    from skimage.measure import marching_cubes
except ImportError:
    logger.warning("Using deprecated version of scikit-image")
    from skimage.measure import marching_cubes_lewiner as marching_cubes

from LoopStructural.interpolators import GeologicalInterpolator
from LoopStructural.datatypes import BoundingBox
from LoopStructural.datatypes import Surface

surface_list = dict[str, tuple[npt.ArrayLike, npt.ArrayLike, npt.ArrayLike, npt.ArrayLike]]


class LoopIsosurfacer:
    def __init__(
        self,
        bounding_box: BoundingBox,
        interpolator: Optional[GeologicalInterpolator] = None,
        callable: Optional[Callable[[npt.ArrayLike], npt.ArrayLike]] = None,
    ):
        """Extract isosurfaces from a geological interpolator or a callable function.


        Parameters
        ----------
        bounding_box : BoundingBox
            _description_
        interpolator : Optional[GeologicalInterpolator], optional
            interpolator object, by default None
        callable : Optional[Callable[[npt.ArrayLike], npt.ArrayLike]], optional
            callable object, by default None

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        ValueError
            _description_
        """
        self.bounding_box = bounding_box
        self.callable = callable
        if interpolator is None and callable is None:
            raise ValueError("Must specify either interpolator or callable")
        if interpolator is not None and self.callable is not None:
            raise ValueError("Must specify either interpolator or callable")

        if interpolator is not None:
            self.callable = interpolator.evaluate_value
        if self.callable is None:
            raise ValueError("Must specify either interpolator or callable")

    def fit(self, values: Union[list, int, float]) -> surface_list:
        """Extract isosurfaces from the interpolator

        Parameters
        ----------
        values : Union[list, int, float]
            Either a list of values to extract isosurfaces for, or a single value
            to extract a single isosurface for, or an integer to extract that many
            isosurfaces evenly spaced between the minimum and maximum values of the
            interpolator.

        Returns
        -------
        surface_list
            a dictionary containing the extracted isosurfaces
        """
        if not callable(self.callable):
            raise ValueError("No interpolator of callable function set")
        surfaces = {}
        all_values = self.callable(self.bounding_box.regular_grid())
        if isinstance(values, list):
            isovalues = values
        elif isinstance(values, float):
            isovalues = [values]
        elif isinstance(values, int):
            isovalues = np.linspace(
                np.min(all_values) + np.finfo(float).eps,
                np.max(all_values) + np.finfo(float).eps,
                values,
            )
        for isovalue in isovalues:
            verts, faces, normals, values = marching_cubes(
                all_values.reshape(self.bounding_box.nsteps, order="C"),
                isovalue,
                spacing=self.bounding_box.step_vector,
            )
            values = np.zeros(verts.shape[0]) + isovalue
            surfaces[f"surface_{isovalue}"] = Surface(
                vertices=verts + self.bounding_box.origin,
                triangles=faces,
                normals=normals,
                name=f"surface_{isovalue}",
                values=values,
            )
        return surfaces
