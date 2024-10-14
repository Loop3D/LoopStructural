from typing import Dict
import numpy as np
from dataclasses import dataclass
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


@dataclass
class StructuredGrid:
    origin: np.ndarray
    step_vector: np.ndarray
    nsteps: np.ndarray
    cell_properties: Dict[str, np.ndarray]
    properties: Dict[str, np.ndarray]
    name: str

    def to_dict(self):
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
        return self.origin + self.nsteps * self.step_vector

    def vtk(self):
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
            grid[name] = data.flatten(order="F")
        for name, data in self.cell_properties.items():
            grid.cell_data[name] = data.flatten(order="F")
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

    def save(self, filename):
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

            add_structured_grid_to_geoh5(filename, self)
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
