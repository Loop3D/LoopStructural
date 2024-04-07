from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class Surface:
    vertices: np.ndarray
    triangles: np.ndarray
    normals: np.ndarray
    name: str
    values: Optional[np.ndarray] = None

    @property
    def vtk(self):
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

    def save(self, filename):
        filename = str(filename)
        ext = filename.split('.')[-1]
        if ext == 'json':
            import json

            with open(filename, 'w') as f:
                json.dump(self.to_dict(), f)
        elif ext == 'vtk':
            self.vtk.save(filename)
        elif ext == 'obj':
            import meshio

            meshio.write_points_cells(
                filename,
                self.vertices,
                [("triangle", self.triangles)],
                point_data={"normals": self.normals},
            )
        elif ext == 'ts':
            from LoopStructural.export.exporters import _write_feat_surfs_gocad

            _write_feat_surfs_gocad(self, filename)
        elif ext == 'pkl':
            import pickle

            with open(filename, 'wb') as f:
                pickle.dump(self, f)
        else:
            raise ValueError(f"Extension {ext} not supported")
