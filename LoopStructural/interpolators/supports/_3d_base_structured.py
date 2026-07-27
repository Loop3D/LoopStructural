from LoopStructural.utils.exceptions import LoopException
from abc import abstractmethod
import numpy as np
from LoopStructural.utils import getLogger
from LoopStructural.geometry import StructuredGrid3DGeometry
from . import SupportType
from typing import Tuple

logger = getLogger(__name__)

from ._base_support import BaseSupport


class BaseStructuredSupport(BaseSupport):
    """ """

    dimension = 3

    def __init__(
        self,
        origin=None,
        nsteps_cells=None,
        step_vector=None,
        rotation_xy=None,
    ):
        """

        Parameters
        ----------
        origin - 3d list or numpy array
        nsteps_cells - 3d list or numpy array of ints, number of cells in each direction
        step_vector - 3d list or numpy array of int
        """
        if origin is None:
            origin = np.zeros(3)
        if nsteps_cells is None:
            nsteps_cells = np.array([10, 10, 10])
        if step_vector is None:
            step_vector = np.ones(3)
        # cast to numpy array, to allow list like input
        nsteps_cells = np.array(nsteps_cells)

        self.type = SupportType.BaseStructured
        if np.any(nsteps_cells == 0):
            raise LoopException("nsteps cannot be zero")
        if np.any(nsteps_cells < 0):
            raise LoopException("nsteps cannot be negative")
        # BaseStructuredSupport's constructor takes nsteps_cells as a *cell* count,
        # while StructuredGrid3DGeometry (like geometry.StructuredGrid/BoundingBox)
        # takes nsteps as a *node* count -- translate here, at the support boundary.
        nsteps_nodes = np.array(nsteps_cells, dtype=int) + 1
        self._geom = StructuredGrid3DGeometry(
            origin=origin, nsteps=nsteps_nodes, step_vector=step_vector, rotation_xy=rotation_xy
        )
        self.supporttype = "Base"
        self.interpolator = None

    @property
    def volume(self):
        return self._geom.volume

    def set_nelements(self, nelements) -> int:
        result = self._geom.set_nelements(nelements)
        self.onGeometryChange()
        return result

    def to_dict(self):
        return {
            "origin": self.origin,
            "nsteps": self.nsteps,
            "step_vector": self.step_vector,
            "rotation_xy": self.rotation_xy,
        }

    @abstractmethod
    def onGeometryChange(self):
        """Function to be called when the geometry of the support changes"""
        pass

    def associateInterpolator(self, interpolator):
        self.interpolator = interpolator

    @property
    def nsteps(self):
        return self._geom.nsteps

    @nsteps.setter
    def nsteps(self, nsteps):
        # if nsteps changes we need to change the step vector
        self._geom.nsteps = nsteps
        self.onGeometryChange()

    @property
    def nsteps_cells(self):
        return self._geom.nsteps_cells

    @property
    def rotation_xy(self):
        return self._geom.rotation_xy

    @rotation_xy.setter
    def rotation_xy(self, rotation_xy):
        self._geom.rotation_xy = rotation_xy

    @property
    def step_vector(self):
        return self._geom.step_vector

    @step_vector.setter
    def step_vector(self, step_vector):
        self._geom.step_vector = step_vector
        self.onGeometryChange()

    @property
    def origin(self):
        return self._geom.origin

    @origin.setter
    def origin(self, origin):
        self._geom.origin = origin
        self.onGeometryChange()

    @property
    def maximum(self):
        return self._geom.maximum

    @maximum.setter
    def maximum(self, maximum):
        """
        update the number of steps to fit new boundary
        """
        self._geom.maximum = maximum
        self.onGeometryChange()

    @property
    def n_nodes(self):
        return self._geom.n_nodes

    @property
    def n_elements(self):
        return self._geom.n_elements

    @property
    def elements(self):
        return self._geom.elements

    def __str__(self):
        return (
            "LoopStructural interpolation support:  {} \n"
            "Origin: {} {} {} \n"
            "Maximum: {} {} {} \n"
            "Step Vector: {} {} {} \n"
            "Number of Steps: {} {} {} \n"
            "Degrees of freedon {}".format(
                self.supporttype,
                self.origin[0],
                self.origin[1],
                self.origin[2],
                self.maximum[0],
                self.maximum[1],
                self.maximum[2],
                self.step_vector[0],
                self.step_vector[1],
                self.step_vector[2],
                self.nsteps[0],
                self.nsteps[1],
                self.nsteps[2],
                self.n_nodes,
            )
        )

    @property
    def nodes(self):
        return self._geom.nodes

    def rotate(self, pos):
        """ """
        return self._geom.rotate(pos)

    def position_to_cell_index(self, pos: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Get the indexes (i,j,k) of a cell
        that a point is inside


        Parameters
        ----------
        pos : np.array
            Nx3 array of xyz locations

        Returns
        -------
        np.ndarray
            N,3 i,j,k indexes of the cell that the point is in
        """
        return self._geom.position_to_cell_index(pos)

    def position_to_cell_global_index(self, pos):
        return self._geom.position_to_cell_global_index(pos)

    def inside(self, pos):
        return self._geom.inside(pos)

    def check_position(self, pos: np.ndarray) -> np.ndarray:
        return self._geom.check_position(pos)

    def cell_corner_indexes(self, cell_indexes: np.ndarray) -> np.ndarray:
        """
        Returns the indexes of the corners of a cell given its location xi,
        yi, zi

        Parameters
        ----------
        x_cell_index
        y_cell_index
        z_cell_index

        Returns
        -------

        """
        return self._geom.cell_corner_indexes(cell_indexes)

    def position_to_cell_corners(self, pos):
        """Get the global indices of the vertices (corners) of the cell containing each point.

        Parameters
        ----------
        pos : np.array
            (N, 3) array of xyz coordinates representing the positions of N points.

        Returns
        -------
        globalidx : np.array
            (N, 8) array of global indices corresponding to the 8 corner nodes of the cell
            each point lies in. If a point lies outside the support, its corresponding entry
            will be set to -1.
        inside : np.array
            (N,) boolean array indicating whether each point is inside the support domain.
        """
        return self._geom.position_to_cell_corners(pos)

    def position_to_cell_vertices(self, pos):
        """Get the vertices of the cell a point is in

        Parameters
        ----------
        pos : np.array
            Nx3 array of xyz locations

        Returns
        -------
        np.array((N,3),dtype=float), np.array(N,dtype=int)
            vertices, inside
        """
        return self._geom.position_to_cell_vertices(pos)

    def node_indexes_to_position(self, node_indexes: np.ndarray) -> np.ndarray:
        return self._geom.node_indexes_to_position(node_indexes)

    def global_index_to_cell_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        return self._geom.global_index_to_cell_index(global_index)

    def global_index_to_node_index(self, global_index):
        """
        Convert from global indexes to xi,yi,zi

        Parameters
        ----------
        global_index

        Returns
        -------

        """
        return self._geom.global_index_to_node_index(global_index)

    def global_node_indices(self, indexes) -> np.ndarray:
        """
        Convert from node indexes to global node index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        return self._geom.global_node_indices(indexes)

    def global_cell_indices(self, indexes) -> np.ndarray:
        """
        Convert from cell indexes to global cell index

        Parameters
        ----------
        indexes

        Returns
        -------

        """
        return self._geom.global_cell_indices(indexes)

    @property
    def element_size(self):
        return self._geom.element_size

    @property
    def element_scale(self):
        # all elements are the same size
        return self._geom.element_scale

    def vtk(self, node_properties=None, cell_properties=None):
        if node_properties is None:
            node_properties = {}
        if cell_properties is None:
            cell_properties = {}
        return self._geom.vtk(node_properties=node_properties, cell_properties=cell_properties)
