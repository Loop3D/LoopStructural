from dataclasses import dataclass, field
import numpy as np

from typing import Optional, Union
import io
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


@dataclass
class ValuePoints:
    locations: np.ndarray = field(default_factory=lambda: np.array([[0, 0, 0]]))
    values: np.ndarray = field(default_factory=lambda: np.array([0]))
    name: str = "unnamed"
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

    def vtk(self, scalars=None):
        import pyvista as pv

        points = pv.PolyData(self.locations)
        if scalars is not None and len(scalars) == len(self.locations):
            points.point_data['scalars'] = scalars
        else:
            points["values"] = self.values
        return points

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

    def save(self, filename: Union[str, io.StringIO], *, group='Loop',ext=None):
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

            add_points_to_geoh5(filename, self, groupname=group)
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
    locations: np.ndarray = field(default_factory=lambda: np.array([[0, 0, 0]]))
    vectors: np.ndarray = field(default_factory=lambda: np.array([[0, 0, 0]]))
    name: str = "unnamed"
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

    def vtk(
        self,
        geom='arrow',
        scale=1.0,
        scale_function=None,
        normalise=False,
        tolerance=0.05,
        bb=None,
        scalars=None,
    ):
        import pyvista as pv

        _projected = False
        vectors = np.copy(self.vectors)

        if normalise:
            norm = np.linalg.norm(vectors, axis=1)
            vectors[norm > 0, :] /= norm[norm > 0][:, None]
        else:
            norm = np.linalg.norm(vectors, axis=1)
            vectors[norm > 0, :] /= norm[norm > 0][:, None]
            norm[norm > 0] = norm[norm > 0] / norm[norm > 0].max()
            vectors[norm > 0, :] *= norm[norm > 0, None]
        if scale_function is not None:
            # vectors /= np.linalg.norm(vectors, axis=1)[:, None]
            vectors *= scale_function(self.locations)[:, None]
        locations = self.locations
        if bb is not None:
            try:
                locations = bb.project(locations)
                _projected = True
            except Exception as e:
                logger.error(f'Failed to project points to bounding box: {e}')
                logger.error('Using unprojected points, this may cause issues with the glyphing')
        points = pv.PolyData(locations)
        if scalars is not None and len(scalars) == len(self.locations):
            points['scalars'] = scalars
        points.point_data.set_vectors(vectors, 'vectors')
        if geom == 'arrow':
            geom = pv.Arrow(scale=scale)
        elif geom == 'disc':
            geom = pv.Disc(inner=0, outer=scale * 0.5, c_res=50).rotate_y(90)

        # Perform the glyph
        glyphed = points.glyph(orient="vectors", geom=geom, tolerance=tolerance)
        if _projected:
            glyphed.points = bb.reproject(glyphed.points)
        return glyphed

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

    def save(self, filename,*, group='Loop'):
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

            add_points_to_geoh5(filename, self, groupname=group)
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
