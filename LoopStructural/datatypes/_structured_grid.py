import numpy as np
from dataclasses import dataclass


@dataclass
class StructuredGrid:
    origin: np.ndarray
    step_vector: np.ndarray
    nsteps: np.ndarray
    data: np.ndarray
    name: str

    def to_dict(self):
        return {
            "origin": self.origin,
            "maximum": self.maximum,
            "step_vector": self.step_vector,
            "nsteps": self.nsteps,
            "data": self.data,
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
        grid[self.name] = self.data.flatten(order="F")
        return grid

    def save(self, filename):
        raise NotImplementedError("Saving structured grids is not yet implemented")
