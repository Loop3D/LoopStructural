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

        return pv.pv.PolyData.from_regular_faces(self.vertices, self.triangles)
