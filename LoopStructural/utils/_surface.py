from __future__ import annotations

from typing import Optional, Union, Callable, List
from collections.abc import Iterable
import numpy as np
import numpy.typing as npt
from LoopStructural.utils.logging import getLogger

logger = getLogger(__name__)
try:
    from skimage.measure import marching_cubes
except ImportError:
    logger.warning("Using deprecated version of scikit-image")
    from skimage.measure import marching_cubes_lewiner as marching_cubes

# from LoopStructural.interpolators._geological_interpolator import GeologicalInterpolator
from LoopStructural.datatypes import Surface, BoundingBox

surface_list = List[Surface]


class LoopIsosurfacer:
    def __init__(
        self,
        bounding_box: BoundingBox,
        interpolator=None,
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

    def fit(
        self,
        values: Optional[Union[list, int, float]],
        name: Optional[Union[List[str], str]] = None,
        local=False,
        colours: Optional[List] = None,
    ) -> surface_list:
        """Extract isosurfaces from the interpolator

        Parameters
        ----------
        values : Union[list, int, float]
            Either a list of values to extract isosurfaces for, or a single value
            to extract a single isosurface for, or an integer to extract that many
            isosurfaces evenly spaced between the minimum and maximum values of the
            interpolator.
        name : Optional[Union[List[str], str]], optional
            name of the isosurface, by default None
        local : bool, optional
            whether to use the local regular grid for the bounding box or global

        Returns
        -------
        surface_list
            a dictionary containing the extracted isosurfaces
        """

        if not callable(self.callable):
            raise ValueError("No interpolator of callable function set")

        surfaces = []
        all_values = self.callable(self.bounding_box.regular_grid(local=local, order='C'))
        ## set value to mean value if its not specified
        if values is None:
            values = [((np.nanmax(all_values) - np.nanmin(all_values)) / 2) + np.nanmin(all_values)]
        if isinstance(values, Iterable):
            isovalues = values
        elif isinstance(values, float):
            isovalues = [values]
        if isinstance(values, int) and values == 0:
            values = 0.0  # assume 0 isosurface is meant to be a float

        elif isinstance(values, int) and values < 1:
            raise ValueError(
                "Number of isosurfaces must be greater than 1. Either use a positive integer or provide a list or float for a specific isovalue."
            )
        elif isinstance(values, int):
            var = np.nanmax(all_values) - np.nanmin(all_values)
            # buffer slices by 5% to make sure that we don't get isosurface does't exist issues
            buffer = var * 0.05
            isovalues = np.linspace(
                np.nanmin(all_values) + buffer,
                np.nanmax(all_values) - buffer,
                values,
            )
        logger.info(f'Isosurfacing at values: {isovalues}')
        individual_names = False
        if name is None:
            names = ["surface"] * len(isovalues)
        if isinstance(name, str):
            names = [name] * len(isovalues)
            if len(isovalues) == 1:
                individual_names = True
        if isinstance(name, list):
            names = name
            if len(names) == len(isovalues):
                individual_names = True
        if colours is None:
            colours = [None] * len(isovalues)
        for name, isovalue, colour in zip(names, isovalues, colours):
            try:
                step_vector = (self.bounding_box.maximum - self.bounding_box.origin) / (
                    np.array(self.bounding_box.nsteps) - 1
                )
                verts, faces, normals, values = marching_cubes(
                    # np.rot90(
                    all_values.reshape(self.bounding_box.nsteps, order="C"),  # k=2, axes=(0, 1)
                    # ),
                    isovalue,
                    spacing=step_vector,
                    mask=~np.isnan(all_values.reshape(self.bounding_box.nsteps, order="C")),
                )
            except RuntimeError:
                logger.warning(f"Failed to extract isosurface for {isovalue}")
                continue
            except ValueError:
                logger.warning(f"Failed to extract isosurface for {isovalue}")
                continue
            values = np.zeros(verts.shape[0]) + isovalue
            # need to add both global and local origin. If the bb is a buffer the local
            # origin may not be 0
            verts += self.bounding_box.global_origin+self.bounding_box.origin
            surfaces.append(
                Surface(
                    vertices=verts,
                    triangles=faces,
                    normals=normals,
                    name=name if individual_names else f"{name}_{isovalue}",
                    values=values,
                    colour=colour,
                )
            )
        return surfaces
