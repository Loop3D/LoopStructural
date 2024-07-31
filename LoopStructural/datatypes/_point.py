from dataclasses import dataclass
import numpy as np

from typing import Optional, Union
import io


@dataclass
class ValuePoints:
    locations: np.ndarray
    values: np.ndarray
    name: str
    properties: Optional[dict] = None

    def to_dict(self):
        return {
            "locations": self.locations,
            "values": self.values,
            "name": self.name,
            "properties": (
                {k: p.tolist() for k, p in self.properties.items()} if self.properties else None
            ),
        }

    def vtk(self):
        import pyvista as pv

        points = pv.PolyData(self.locations)
        points["values"] = self.values
        return points

    def save(self, filename: Union[str, io.StringIO], ext=None):
        if isinstance(filename, io.StringIO):
            if ext is None:
                raise ValueError('Please provide an extension for StringIO')
            ext = ext.lower()
        else:
            ext = filename.split('.')[-1].lower()
            filename = str(filename)
        if ext == 'json':
            import json

            with open(filename, 'w') as f:
                json.dump(self.to_dict(), f)
        elif ext == 'vtk':
            self.vtk().save(filename)

        elif ext == 'geoh5':
            from LoopStructural.export.geoh5 import add_points_to_geoh5

            add_points_to_geoh5(filename, self)
        elif ext == 'pkl':
            import pickle

            with open(filename, 'wb') as f:
                pickle.dump(self, f)
        elif ext == 'vs':
            from LoopStructural.export.gocad import _write_pointset

            _write_pointset(self, filename)
        elif ext == 'csv':
            import pandas as pd

            df = pd.DataFrame(self.locations, columns=['x', 'y', 'z'])
            df['value'] = self.values
            if self.properties is not None:
                for k, v in self.properties.items():
                    df[k] = v
            df.to_csv(filename, index=False)
        elif ext == 'omf':
            from LoopStructural.export.omf_wrapper import add_pointset_to_omf

            add_pointset_to_omf(self, filename)
        else:
            raise ValueError(f'Unknown file extension {ext}')

    @classmethod
    def from_dict(cls, d, flatten=False):
        if 'locations' not in d:
            raise ValueError('locations not in dictionary')
        locations = np.array(d['locations'])
        if flatten:
            locations = locations.reshape((-1, 3))
        return ValuePoints(
            locations, d.get('values', None), d.get('name', 'unnamed'), d.get('properties', None)
        )


@dataclass
class VectorPoints:
    locations: np.ndarray
    vectors: np.ndarray
    name: str
    properties: Optional[dict] = None

    def to_dict(self):
        return {
            "locations": self.locations,
            "vectors": self.vectors,
            "name": self.name,
            "properties": (
                {k: p.tolist() for k, p in self.properties.items()} if self.properties else None
            ),
        }

    def from_dict(self, d):
        return VectorPoints(d['locations'], d['vectors'], d['name'], d.get('properties', None))

    def vtk(self, geom='arrow', scale=1.0, scale_function=None, tolerance=0.05):
        import pyvista as pv

        vectors = np.copy(self.vectors)
        if scale_function is not None:
            vectors /= np.linalg.norm(vectors, axis=1)[:, None]
            vectors *= scale_function(self.locations)[:, None]
        points = pv.PolyData(self.locations)
        points.point_data.set_vectors(vectors, 'vectors')
        if geom == 'arrow':
            geom = pv.Arrow(scale=scale)
        elif geom == 'disc':
            geom = pv.Disc(inner=0, outer=scale)

        # Perform the glyph
        return points.glyph(orient="vectors", geom=geom, tolerance=tolerance)

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
            from LoopStructural.export.geoh5 import add_points_to_geoh5

            add_points_to_geoh5(filename, self)
        elif ext == 'pkl':
            import pickle

            with open(filename, 'wb') as f:
                pickle.dump(self, f)
        elif ext == 'vs':
            from LoopStructural.export.gocad import _write_pointset

            _write_pointset(self, filename)
        elif ext == 'csv':
            import pandas as pd

            df = pd.DataFrame(self.locations, columns=['x', 'y', 'z'])
            df['vx'] = self.vectors[:, 0]
            df['vy'] = self.vectors[:, 1]
            df['vz'] = self.vectors[:, 2]
            if self.properties is not None:
                for k, v in self.properties.items():
                    df[k] = v
            df.to_csv(filename)
        elif ext == 'omf':
            from LoopStructural.export.omf_wrapper import add_pointset_to_omf

            add_pointset_to_omf(self, filename)
        else:
            raise ValueError(f'Unknown file extension {ext}')
