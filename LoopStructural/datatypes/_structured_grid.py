from typing import Dict
import numpy as np
from dataclasses import dataclass, field
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


@dataclass
class StructuredGrid:
    """A structured grid for storing 3D geological data.

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
    cell_properties: Dict[str, np.ndarray] = field(default_factory=dict)
    properties: Dict[str, np.ndarray] = field(default_factory=dict)
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
    def maximum(self):
        """Calculate the maximum coordinates of the grid.

        Returns
        -------
        np.ndarray
            Maximum coordinates (origin + nsteps * step_vector)
        """
        return self.origin + self.nsteps * self.step_vector

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

    def plot(self, pyvista_kwargs={}):
        """Calls pyvista plot on the vtk object

        Parameters
        ----------
        pyvista_kwargs : dict, optional
            kwargs passed to pyvista.DataSet.plot(), by default {}
        """
        try:
            self.vtk().plot(**pyvista_kwargs)
            return
        except ImportError:
            logger.error("pyvista is required for vtk")

    @property
    def cell_centres(self):
        """Calculate the coordinates of cell centres.

        Returns
        -------
        tuple of np.ndarray
            X, Y, Z coordinates of all cell centres
        """
        x = np.linspace(
            self.origin[0] + self.step_vector[0] * 0.5,
            self.maximum[0] + self.step_vector[0] * 0.5,
            self.nsteps[0] - 1,
        )
        y = np.linspace(
            self.origin[1] + self.step_vector[1] * 0.5,
            self.maximum[1] - self.step_vector[1] * 0.5,
            self.nsteps[1] - 1,
        )
        z = np.linspace(
            self.origin[2] + self.step_vector[2] * 0.5,
            self.maximum[2] - self.step_vector[2] * 0.5,
            self.nsteps[2] - 1,
        )
        x, y, z = np.meshgrid(x, y, z, indexing="ij")
        return np.vstack([x.flatten(order='f'), y.flatten(order='f'), z.flatten(order='f')]).T

    @property
    def nodes(self):
        x = np.linspace(self.origin[0], self.maximum[0], self.nsteps[0])
        y = np.linspace(self.origin[1], self.maximum[1], self.nsteps[1])
        z = np.linspace(self.origin[2], self.maximum[2], self.nsteps[2])
        x, y, z = np.meshgrid(x, y, z, indexing="ij")
        return np.vstack([x.flatten(order='f'), y.flatten(order='f'), z.flatten(order='f')]).T

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

    def save(self, filename, *,group='Loop'):
        filename = str(filename)
        ext = filename.split('.')[-1]
        if ext == 'json':
            import json

            with open(filename, 'w') as f:
                json.dump(self.to_dict(), f)
        elif ext == 'vtk':
            self.vtk().save(filename)

        elif ext == 'geoh5':
            from LoopStructural.export.geoh5 import add_structured_grid_to_geoh5

            add_structured_grid_to_geoh5(filename, self, groupname=group)
        elif ext == 'pkl':
            import pickle

            with open(filename, 'wb') as f:
                pickle.dump(self, f)
        elif ext == 'omf':
            from LoopStructural.export.omf_wrapper import add_structured_grid_to_omf

            add_structured_grid_to_omf(self, filename)
        elif ext == 'vs':
            raise NotImplementedError(
                "Saving structured grids in gocad format is not yet implemented"
            )
            # from LoopStructural.export.gocad import _write_structued_grid

            # _write_pointset(self, filename)
        else:
            raise ValueError(f'Unknown file extension {ext}')
