from abc import abstractmethod

import numpy as np

from loop_common.geometry import StructuredGrid3DGeometry
from loop_common.logging import get_logger as getLogger

from . import SupportType

logger = getLogger(__name__)

from ._base_support import BaseSupport


class LoopException(Exception):
    """Custom exception for LoopStructural errors."""


class BaseStructuredSupport(StructuredGrid3DGeometry, BaseSupport):
    """A structured grid used as an interpolation support.

    The grid geometry (origin, step vector, node counts, rotation, node and cell indexing)
    comes from :class:`loop_common.geometry.StructuredGrid3DGeometry`; this class adds the
    interpolation machinery. Note the constructor takes the number of **cells**
    (``nsteps``) whereas the geometry class takes the number of nodes.
    """

    dimension = 3

    def __init__(
        self,
        origin=None,
        nsteps=None,
        step_vector=None,
        rotation_xy=None,
    ):
        """

        Parameters
        ----------
        origin - 3d list or numpy array
        nsteps - 3d list or numpy array of ints, the number of cells in each direction
        step_vector - 3d list or numpy array of int
        """
        if nsteps is None:
            nsteps = np.array([10, 10, 10])
        nsteps = np.array(nsteps)
        if np.any(nsteps == 0):
            raise LoopException("nsteps cannot be zero")
        if np.any(nsteps < 0):
            raise LoopException("nsteps cannot be negative")
        StructuredGrid3DGeometry.__init__(
            self,
            origin=origin,
            nsteps=nsteps.astype(int) + 1,
            step_vector=step_vector,
            rotation_xy=rotation_xy,
        )
        self.type = SupportType.BaseStructured
        self.supporttype = "Base"
        self.interpolator = None

    def _on_geometry_change(self):
        self.onGeometryChange()

    @abstractmethod
    def onGeometryChange(self):
        """Function to be called when the geometry of the support changes"""

    def associateInterpolator(self, interpolator):
        self.interpolator = interpolator

    def __str__(self):
        return (
            f"LoopStructural interpolation support:  {self.supporttype} \n"
            f"Origin: {self.origin[0]} {self.origin[1]} {self.origin[2]} \n"
            f"Maximum: {self.maximum[0]} {self.maximum[1]} {self.maximum[2]} \n"
            f"Step Vector: {self.step_vector[0]} {self.step_vector[1]} {self.step_vector[2]} \n"
            f"Number of Steps: {self.nsteps[0]} {self.nsteps[1]} {self.nsteps[2]} \n"
            f"Degrees of freedon {self.n_nodes}"
        )

    def position_to_cell_global_index(self, pos):
        _ix, _iy, _iz = self.position_to_cell_index(pos)

    def check_position(self, pos: np.ndarray) -> np.ndarray:
        """[summary]

        [extended_summary]

        Parameters
        ----------
        pos : [type]
            [description]

        Returns
        -------
        [type]
            [description]
        """
        if not isinstance(pos, np.ndarray):
            try:
                pos = np.array(pos, dtype=float)
            except Exception as e:
                logger.error(
                    f"Position array should be a numpy array or list of points, not {type(pos)}"
                )
                raise ValueError(
                    f"Position array should be a numpy array or list of points, not {type(pos)}"
                ) from e

        if len(pos.shape) == 1:
            pos = np.array([pos])
        if len(pos.shape) != 2:
            logger.error("Position array needs to be a list of points or a point")
            raise ValueError("Position array needs to be a list of points or a point")
        return pos

    def vtk(self, node_properties=None, cell_properties=None):
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for vtk support")

        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}
        from pyvista import CellType

        celltype = np.full(self.n_elements, CellType.VOXEL, dtype=np.uint8)
        elements = np.hstack(
            [np.zeros(self.elements.shape[0], dtype=int)[:, None] + 8, self.elements]
        )
        elements = elements.flatten()
        grid = pv.UnstructuredGrid(elements, celltype, self.nodes)
        for key, value in node_properties.items():
            grid[key] = value
        for key, value in cell_properties.items():
            grid.cell_arrays[key] = value
        return grid
