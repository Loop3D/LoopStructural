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
        elif ext == 'geoh5':
            from LoopStructural.export.geoh5 import add_surface_to_geoh5

            add_surface_to_geoh5(filename, self)

        elif ext == 'pkl':
            import pickle

            with open(filename, 'wb') as f:
                pickle.dump(self, f)
        else:
            raise ValueError(f"Extension {ext} not supported")
