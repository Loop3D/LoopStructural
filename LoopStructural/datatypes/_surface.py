from dataclasses import dataclass, field
from typing import Optional, Union
import numpy as np
import io
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


@dataclass
class Surface:
    vertices: np.ndarray = field(default_factory=lambda: np.array([[0, 0, 0]]))
    triangles: np.ndarray = field(default_factory=lambda: np.array([[0, 0, 0]]))
    colour: Optional[Union[str, np.ndarray]] = field(default_factory=lambda: None)
    normals: Optional[np.ndarray] = None
    name: str = 'surface'
    values: Optional[np.ndarray] = None
    properties: Optional[dict] = None
    cell_properties: Optional[dict] = None
    def __post_init__(self):
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3:
            raise ValueError("vertices must be a Nx3 numpy array")
        if self.triangles.ndim != 2 or self.triangles.shape[1] != 3:
            raise ValueError("triangles must be a Mx3 numpy array")
        if self.normals is not None:
            if (self.normals.shape[1] != 3 or 
                (self.normals.shape[0] != self.vertices.shape[0] and self.normals.shape[0] != self.triangles.shape[0])):
                raise ValueError("normals must be a Nx3 numpy array where N is the number of vertices or triangles")
        if self.values is not None:
            if self.values.shape[0] != self.vertices.shape[0]:
                raise ValueError("values must be a N numpy array where N is the number of vertices")
        if self.properties is not None:
            for k, v in self.properties.items():
                if len(v) != self.vertices.shape[0]:
                    raise ValueError(f"property {k} must be a list or array of length {self.vertices.shape[0]}")
        if self.cell_properties is not None:
            for k, v in self.cell_properties.items():
                if len(v) != self.triangles.shape[0]:
                    raise ValueError(f"cell property {k} must be a list or array of length {self.triangles.shape[0]}")
        if np.isnan(self.vertices).any():
            self.remove_nan_vertices()
    def remove_nan_vertices(self):
        """Remove vertices with NaN values from the surface. Also removes any triangles that reference these vertices.
        This modifies the vertices and triangles in place. Any associated properties are also updated.
        """
        vertex_index = np.arange(0, self.vertices.shape[0])
        not_nan_indexes = np.where(~np.isnan(self.vertices).any(axis=1))[0]
        new_vertex_map = np.zeros(self.vertices.shape[0], dtype=int) - 1
        new_vertex_map[not_nan_indexes] = np.arange(0, not_nan_indexes.shape[0])
        nan_indexes = np.setdiff1d(vertex_index, not_nan_indexes)
        triangles_with_nan = np.any(np.isin(self.triangles, nan_indexes), axis=1)
        nan_triangle_indexes = np.where(triangles_with_nan)[0]
        new_triangles = np.delete(self.triangles, nan_triangle_indexes, axis=0)
        vertices = self.vertices[not_nan_indexes]
        self.triangles = new_vertex_map[new_triangles]
        self.vertices = vertices
        if self.normals is not None:
            self.normals = self.normals[not_nan_indexes]
        if self.values is not None:
            self.values = self.values[not_nan_indexes]
        if self.properties is not None:
            for k, v in self.properties.items():
                self.properties[k] = np.array(v)[not_nan_indexes]
        if self.cell_properties is not None:
            for k, v in self.cell_properties.items():
                self.cell_properties[k] = np.array(v)[~triangles_with_nan]
    @property
    def triangle_area(self):
        """_summary_

        Returns
        -------
        _type_
            _description_


        Notes
        -----

        Area of triangle for a 3d triangle with vertices at points A, B, C is given by
        det([A-C, B-C])**.5
        """
        tri_points = self.vertices[self.triangles, :]
        mat = np.array(
            [
                [
                    tri_points[:, 0, 0] - tri_points[:, 2, 0],
                    tri_points[:, 0, 1] - tri_points[:, 2, 1],
                    tri_points[:, 0, 2] - tri_points[:, 2, 2],
                ],
                [
                    tri_points[:, 1, 0] - tri_points[:, 2, 0],
                    tri_points[:, 1, 1] - tri_points[:, 2, 1],
                    tri_points[:, 1, 2] - tri_points[:, 2, 2],
                ],
            ]
        )
        matdotmatT = np.einsum("ijm,mjk->mik", mat, mat.T)
        area = np.sqrt(np.linalg.det(matdotmatT))
        return area

    @property
    def triangle_normal(self) -> np.ndarray:
        """_summary_

        Returns
        -------
        np.ndarray
            numpy array of normals N,3 where N is the number of triangles


        Notes
        -----

        The normal of a triangle is given by the cross product of two vectors in the plane of the triangle
        """
        tri_points = self.vertices[self.triangles, :]
        normals = np.cross(
            tri_points[:, 0, :] - tri_points[:, 2, :], tri_points[:, 1, :] - tri_points[:, 2, :]
        )
        normals = normals / np.linalg.norm(normals, axis=1)[:, np.newaxis]
        return normals

    def vtk(self):
        import pyvista as pv

        surface = pv.PolyData.from_regular_faces(self.vertices, self.triangles)
        if self.values is not None:
            surface["values"] = self.values
        if self.properties is not None:
            for k, v in self.properties.items():
                surface.point_data[k] = np.array(v)
        if self.cell_properties is not None:
            for k, v in self.cell_properties.items():
                surface.cell_data[k] = np.array(v)
        return surface

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

    def to_dict(self, flatten=False):
        triangles = self.triangles
        vertices = self.vertices
        if flatten:
            vertices = self.vertices.flatten()
            triangles = (
                np.hstack([np.ones((self.triangles.shape[0], 1)) * 3, self.triangles])
                .astype(int)
                .flatten()
            )
        return {
            "vertices": vertices.tolist(),
            "triangles": triangles.tolist(),
            "normals": self.normals.tolist() if self.normals is not None else None,
            "properties": (
                {k: p.tolist() for k, p in self.properties.items()} if self.properties else None
            ),
            "cell_properties": (
                {k: p.tolist() for k, p in self.cell_properties.items()}
                if self.cell_properties
                else None
            ),
            "name": self.name,
            "values": self.values.tolist() if self.values is not None else None,
        }

    @classmethod
    def from_dict(cls, d, flatten=False):
        vertices = np.array(d['vertices'])
        triangles = np.array(d['triangles'])
        if flatten:
            vertices = vertices.reshape((-1, 3))
            triangles = triangles.reshape((-1, 4))[:, 1:]
        return cls(
            vertices,
            triangles,
            np.array(d['normals']),
            d['name'],
            np.array(d['values']),
            d.get('properties', None),
            d.get('cell_properties', None),
        )

    def save(self, filename, *, group='Loop',replace_spaces=True, ext=None):
        filename = filename.replace(' ', '_') if replace_spaces else filename
        if isinstance(filename, (io.StringIO, io.BytesIO)):
            if ext is None:
                raise ValueError('Please provide an extension for StringIO')
            ext = ext.lower()
        else:
            filename = str(filename)
            if ext is None:
                ext = filename.split('.')[-1].lower()
        if ext == 'json':
            import json

            with open(filename, 'w') as f:
                json.dump(self.to_dict(), f)
        elif ext == 'vtk':
            self.vtk().save(filename)
        elif ext == 'obj':
            import meshio

            meshio.write_points_cells(
                filename,
                self.vertices,
                [("triangle", self.triangles)],
                point_data={"normals": self.normals},
            )
        elif ext == 'ts' or ext == 'gocad':
            from LoopStructural.export.exporters import _write_feat_surfs_gocad

            _write_feat_surfs_gocad(self, filename)
        elif ext == 'geoh5':
            from LoopStructural.export.geoh5 import add_surface_to_geoh5

            add_surface_to_geoh5(filename, self, groupname=group)

        elif ext == 'pkl':
            import pickle

            with open(filename, 'wb') as f:
                pickle.dump(self, f)
        elif ext == 'csv':
            import pandas as pd

            df = pd.DataFrame(self.vertices, columns=['x', 'y', 'z'])
            if self.properties:
                for k, v in self.properties.items():
                    df[k] = v
            df.to_csv(filename, index=False)
        elif ext == 'omf':
            from LoopStructural.export.omf_wrapper import add_surface_to_omf

            add_surface_to_omf(self, filename)
        else:
            raise ValueError(f"Extension {ext} not supported")
