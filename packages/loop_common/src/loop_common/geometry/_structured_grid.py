from dataclasses import dataclass, field

import numpy as np

from loop_common.logging import get_logger as getLogger

from ._structured_grid_3d import StructuredGrid3DGeometry

logger = getLogger(__name__)


@dataclass
class PropertyGrid:
    """A regular 3D grid for storing property values (the result container).

    This is a data container, not an interpolation support: it holds
    node-based ``properties`` and ``cell_properties`` on a regular grid and can
    export them (``vtk``, ``save``). For the class used to interpolate on a
    grid see :class:`loop_common.supports.StructuredGrid`.

    The grid geometry (node and cell positions, point-in-cell lookup) is defined once, by
    :class:`StructuredGrid3DGeometry`, which is built from ``origin``, ``step_vector`` and
    ``nsteps`` (the number of **nodes**) whenever it is needed; see :attr:`geometry`. Points
    passed to the lookup methods are in the grid's own frame (the grid is not rotated).

    It was previously named ``StructuredGrid``; that name remains as an alias.

    This class represents a regular 3D grid with properties and cell properties
    that can be used for geological modelling and visualisation.

    Parameters
    ----------
    origin : np.ndarray, optional
        Origin point of the grid, by default [0, 0, 0]
    step_vector : np.ndarray, optional
        Step size in each direction, by default [1, 1, 1]
    nsteps : np.ndarray, optional
        Number of steps in each direction, by default [10, 10, 10]
    cell_properties : Dict[str, np.ndarray], optional
        Properties defined at cell centres, by default empty dict
    properties : Dict[str, np.ndarray], optional
        Properties defined at grid nodes, by default empty dict
    name : str, optional
        Name of the grid, by default "default_grid"
    """

    origin: np.ndarray = field(default_factory=lambda: np.array([0, 0, 0]))
    step_vector: np.ndarray = field(default_factory=lambda: np.array([1, 1, 1]))
    nsteps: np.ndarray = field(default_factory=lambda: np.array([10, 10, 10]))
    cell_properties: dict[str, np.ndarray] = field(default_factory=dict)
    properties: dict[str, np.ndarray] = field(default_factory=dict)
    name: str = "default_grid"

    def to_dict(self):
        """Convert the structured grid to a dictionary representation.

        Returns
        -------
        dict
            Dictionary containing all grid properties and metadata
        """
        return {
            "origin": self.origin,
            "maximum": self.maximum,
            "step_vector": self.step_vector,
            "nsteps": self.nsteps,
            "cell_properties": self.cell_properties,
            "properties": self.properties,
            "name": self.name,
        }

    @property
    def geometry(self) -> StructuredGrid3DGeometry:
        """The grid geometry, built from the current ``origin``, ``step_vector`` and ``nsteps``.

        It is a fresh object each time, so it always reflects the current fields.
        """
        return StructuredGrid3DGeometry(
            origin=self.origin, nsteps=self.nsteps, step_vector=self.step_vector
        )

    @property
    def maximum(self):
        """Calculate the maximum coordinates of the grid.

        Returns
        -------
        np.ndarray
            Maximum coordinates (origin + (nsteps - 1) * step_vector)
        """
        return self.geometry.maximum

    @property
    def nsteps_cells(self):
        """Number of cells in each direction (``nsteps - 1``)."""
        return self.geometry.nsteps_cells

    @property
    def n_nodes(self):
        """Total number of nodes."""
        return self.geometry.n_nodes

    @property
    def n_cells(self):
        """Total number of cells."""
        return self.geometry.n_elements

    def inside(self, points):
        """Boolean mask of the points inside the grid (boundary included)."""
        return self.geometry.inside(np.asarray(points, dtype=float))

    def cell_index_at(self, points):
        """Flat index of the cell each point falls in.

        The flat index addresses ``cell_properties`` arrays and the rows of ``cell_centres``
        (x varies fastest). Points on the maximum boundary belong to the last cell.

        Parameters
        ----------
        points : array-like
            N x 3 locations (or one location) in the grid's own frame

        Returns
        -------
        np.ndarray
            N integer indexes, ``-1`` for points outside the grid
        """
        geometry = self.geometry
        cell_ijk, inside = geometry.position_to_cell_index(np.asarray(points, dtype=float))
        index = np.full(len(inside), -1, dtype=np.int64)
        index[inside] = geometry.global_cell_indices(cell_ijk[inside])
        return index

    def evaluate_value(self, name, points):
        """Trilinear interpolation of the node property ``name`` at the points.

        This is what is needed to transfer a property onto another support or mesh: pass the
        other support's node coordinates as ``points``.

        Parameters
        ----------
        name : str
            key of ``properties``; values are per node, in node order (x varies fastest), either
            ``(n_nodes,)`` or ``(n_nodes, k)``
        points : array-like
            N x 3 locations (or one location) in the grid's own frame

        Returns
        -------
        np.ndarray
            ``(N,)`` (or ``(N, k)``) values; NaN for points outside the grid
        """
        return self.geometry.evaluate_value(np.asarray(points, dtype=float), self.properties[name])

    def cell_value_at(self, name, points):
        """Value of the cell property ``name`` in the cell each point falls in.

        Parameters
        ----------
        name : str
            key of ``cell_properties``
        points : array-like
            N x 3 locations (or one location) in the grid's own frame

        Returns
        -------
        np.ndarray
            one entry per point; NaN for points outside the grid
        """
        values = np.asarray(self.cell_properties[name])
        index = self.cell_index_at(points)
        found = index >= 0
        result = np.full((len(index),) + values.shape[1:], np.nan)
        result[found] = values[index[found]]
        return result

    def vtk(self):
        """Convert the structured grid to a PyVista RectilinearGrid.

        Returns
        -------
        pv.RectilinearGrid
            PyVista grid object with all properties attached

        Raises
        ------
        ImportError
            If PyVista is not installed
        """
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for vtk support")
        x = np.linspace(self.origin[0], self.maximum[0], self.nsteps[0])
        y = np.linspace(self.origin[1], self.maximum[1], self.nsteps[1])
        z = np.linspace(self.origin[2], self.maximum[2], self.nsteps[2])
        grid = pv.RectilinearGrid(
            x,
            y,
            z,
        )
        for name, data in self.properties.items():
            grid[name] = data.reshape((grid.n_points, -1), order="F")
        for name, data in self.cell_properties.items():
            grid.cell_data[name] = data.reshape((grid.n_cells, -1), order="F")
        return grid

    def plot(self, pyvista_kwargs=None):
        """Calls pyvista plot on the vtk object

        Parameters
        ----------
        pyvista_kwargs : dict, optional
            kwargs passed to pyvista.DataSet.plot(), by default {}
        """
        if pyvista_kwargs is None:
            pyvista_kwargs = {}
        try:
            self.vtk().plot(**pyvista_kwargs)
            return
        except ImportError:
            logger.error("pyvista is required for vtk")

    @property
    def cell_centres(self):
        """Coordinates of the centres of the ``nsteps - 1`` cells in each direction.

        Returns
        -------
        np.ndarray
            (n_cells, 3) array in the same order as ``cell_properties`` (x varies fastest)
        """
        geometry = self.geometry
        cell_ijk = geometry.global_index_to_cell_index(np.arange(geometry.n_elements))
        return geometry.origin + (cell_ijk + 0.5) * geometry.step_vector

    @property
    def nodes(self):
        """Coordinates of the nodes, (n_nodes, 3), in the same order as ``properties``."""
        return self.geometry.nodes

    def merge(self, other):
        if not np.all(np.isclose(self.origin, other.origin)):
            raise ValueError("Origin of grids must be the same")
        if not np.all(np.isclose(self.step_vector, other.step_vector)):
            raise ValueError("Step vector of grids must be the same")
        if not np.all(np.isclose(self.nsteps, other.nsteps)):
            raise ValueError("Number of steps of grids must be the same")

        for name, data in other.cell_properties.items():
            self.cell_properties[name] = data
        for name, data in other.properties.items():
            self.properties[name] = data

    def save(self, filename, *, group="Loop"):
        filename = str(filename)
        ext = filename.split(".")[-1].lower()
        if ext == "json":
            import json

            with open(filename, "w") as f:
                json.dump(self.to_dict(), f)
        elif ext == "vtk":
            self.vtk().save(filename)

        elif ext == "geoh5":
            from LoopStructural.export.geoh5 import add_structured_grid_to_geoh5

            add_structured_grid_to_geoh5(filename, self, groupname=group)
        elif ext == "pkl":
            import pickle

            with open(filename, "wb") as f:
                pickle.dump(self, f)
        elif ext == "omf":
            from LoopStructural.export.omf_wrapper import add_structured_grid_to_omf

            add_structured_grid_to_omf(self, filename)
        elif ext == "vo" or ext == "gocad":
            from LoopStructural.export.gocad import _write_structured_grid_gocad

            _write_structured_grid_gocad(self, filename)
        else:
            raise ValueError(f"Unknown file extension {ext}")


# Backwards-compatible name (1.x); also keeps old pickles resolvable.
StructuredGrid = PropertyGrid
