from dataclasses import dataclass
import numpy as np


@dataclass
class ValuePoints:
    locations: np.ndarray
    values: np.ndarray
    name: str

    def to_dict(self):
        return {
            "locations": self.locations,
            "values": self.values,
            "name": self.name,
        }

    def vtk(self):
        import pyvista as pv

        points = pv.PolyData(self.locations)
        points["values"] = self.values
        return points


@dataclass
class VectorPoints:
    locations: np.ndarray
    vectors: np.ndarray
    name: str

    def to_dict(self):
        return {
            "locations": self.locations,
            "vectors": self.vectors,
            "name": self.name,
        }

    def vtk(self):
        import pyvista as pv

        points = pv.PolyData(self.locations)
        points.point_data.set_vectors(self.vectors, 'vectors')
        geom = pv.Arrow()

        # Perform the glyph
        return points.glyph(orient="vectors", geom=geom)
