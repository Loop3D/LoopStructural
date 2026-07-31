from dataclasses import dataclass, field

import numpy as np

from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)


@dataclass
class StructuredGrid:
    """A structured grid for storing 3D geological data."""

    origin: np.ndarray = field(default_factory=lambda: np.array([0, 0, 0]))
    step_vector: np.ndarray = field(default_factory=lambda: np.array([1, 1, 1]))
    nsteps: np.ndarray = field(default_factory=lambda: np.array([10, 10, 10]))
    cell_properties: dict[str, np.ndarray] = field(default_factory=dict)
    properties: dict[str, np.ndarray] = field(default_factory=dict)
    name: str = "default_grid"

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
        return self.origin + (self.nsteps - 1) * self.step_vector

    def vtk(self):
        try:
            import pyvista as pv
        except ImportError as exc:
            raise ImportError("pyvista is required for vtk support") from exc
        x = np.linspace(self.origin[0], self.maximum[0], self.nsteps[0])
        y = np.linspace(self.origin[1], self.maximum[1], self.nsteps[1])
        z = np.linspace(self.origin[2], self.maximum[2], self.nsteps[2])
        grid = pv.RectilinearGrid(x, y, z)
        for name, data in self.properties.items():
            grid[name] = data.reshape((grid.n_points, -1), order="F")
        for name, data in self.cell_properties.items():
            grid.cell_data[name] = data.reshape((grid.n_cells, -1), order="F")
        return grid

    def plot(self, pyvista_kwargs=None):
        if pyvista_kwargs is None:
            pyvista_kwargs = {}
        try:
            self.vtk().plot(**pyvista_kwargs)
            return
        except ImportError:
            logger.error("pyvista is required for vtk")

    @property
    def cell_centres(self):
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
        return np.vstack([x.flatten(order="f"), y.flatten(order="f"), z.flatten(order="f")]).T

    @property
    def nodes(self):
        x = np.linspace(self.origin[0], self.maximum[0], self.nsteps[0])
        y = np.linspace(self.origin[1], self.maximum[1], self.nsteps[1])
        z = np.linspace(self.origin[2], self.maximum[2], self.nsteps[2])
        x, y, z = np.meshgrid(x, y, z, indexing="ij")
        return np.vstack([x.flatten(order="f"), y.flatten(order="f"), z.flatten(order="f")]).T

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
        elif ext == "pkl":
            import pickle
            with open(filename, "wb") as f:
                pickle.dump(self, f)
        else:
            raise ValueError(f"Unknown file extension {ext}")
