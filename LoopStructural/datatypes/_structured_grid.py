from typing import Dict
import numpy as np
from dataclasses import dataclass


@dataclass
class StructuredGrid:
    origin: np.ndarray
    step_vector: np.ndarray
    nsteps: np.ndarray
    properties_cell: Dict[str, np.ndarray]
    properties_vertex: Dict[str, np.ndarray]
    name: str

    def to_dict(self):
        return {
            "origin": self.origin,
            "maximum": self.maximum,
            "step_vector": self.step_vector,
            "nsteps": self.nsteps,
            "properties_cell": self.properties_cell,
            "properties_vertex": self.properties_vertex,
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
        for name, data in self.properties_vertex.items():
            grid[name] = data.flatten(order="F")
        for name, data in self.properties_cell.items():
            grid.cell_data[name] = data.flatten(order="F")
        return grid

    def merge(self, other):
        if not np.all(np.isclose(self.origin, other.origin)):
            raise ValueError("Origin of grids must be the same")
        if not np.all(np.isclose(self.step_vector, other.step_vector)):
            raise ValueError("Step vector of grids must be the same")
        if not np.all(np.isclose(self.nsteps, other.nsteps)):
            raise ValueError("Number of steps of grids must be the same")

        for name, data in other.properties_cell.items():
            self.properties_cell[name] = data
        for name, data in other.properties_vertex.items():
            self.properties_vertex[name] = data

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
        elif ext == 'vs':
            raise NotImplementedError(
                "Saving structured grids in gocad format is not yet implemented"
            )
            # from LoopStructural.export.gocad import _write_structued_grid

            # _write_pointset(self, filename)
        else:
            raise ValueError(f'Unknown file extension {ext}')
