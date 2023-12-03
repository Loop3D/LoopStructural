from dataclasses import dataclass
from typing import Optional, Union
import numpy as np


@dataclass
class Surface:
    vertices: np.ndarray
    triangles: np.ndarray
    normals: np.ndarray
    name: str
    values: Optional[np.ndarray] = None

    @property
    def pyvista(self):
        import pyvista as pv

        surface = pv.PolyData.from_regular_faces(self.vertices, self.triangles)
        if self.values is not None:
            surface["values"] = self.values
        return surface

    def to_dict(self):
        return {
            "vertices": self.vertices,
            "triangles": self.triangles,
            "normals": self.normals,
            "name": self.name,
            "values": self.values,
        }
