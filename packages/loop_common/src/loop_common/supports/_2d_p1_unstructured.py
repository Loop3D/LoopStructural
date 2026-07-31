"""
Tetmesh based on cartesian grid for piecewise linear interpolation
"""
from __future__ import annotations

import logging

import numpy as np

from . import SupportType
from ._2d_base_unstructured import BaseUnstructured2d
from ._2d_structured_grid import StructuredGrid2D

logger = logging.getLogger(__name__)


class P1Unstructured2d(BaseUnstructured2d):
    """ """

    def __init__(
        self,
        elements: np.ndarray | None = None,
        vertices: np.ndarray | None = None,
        neighbours: np.ndarray | None = None,
        aabb_nsteps=None,
        origin: np.ndarray | None = None,
        step_vector: np.ndarray | None = None,
        nsteps: np.ndarray | None = None,
    ):
        if elements is None or vertices is None or neighbours is None:
            if origin is None or step_vector is None or nsteps is None:
                raise ValueError(
                    "P1Unstructured2d requires either explicit elements/vertices/"
                    "neighbours arrays, or origin/step_vector/nsteps to build a "
                    "triangular mesh over a bounding box"
                )
            vertices, elements, neighbours = self._build_from_bbox(origin, step_vector, nsteps)
        BaseUnstructured2d.__init__(self, elements, vertices, neighbours, aabb_nsteps)
        self.type = SupportType.P1Unstructured2d

    @staticmethod
    def _build_from_bbox(origin: np.ndarray, step_vector: np.ndarray, nsteps: np.ndarray):
        """Build a triangular mesh over a structured 2D grid by splitting
        every grid cell into two triangles along the (bottom-left, top-right)
        diagonal.

        Returns
        -------
        tuple of (vertices, elements, neighbours) suitable for
        BaseUnstructured2d.__init__
        """
        grid = StructuredGrid2D(origin=origin, nsteps=nsteps, step_vector=step_vector)
        vertices = grid.nodes
        quads = grid.elements

        tri_a = quads[:, [0, 1, 2]]
        tri_b = quads[:, [1, 3, 2]]
        elements = np.vstack([tri_a, tri_b])
        n_tris = elements.shape[0]

        local_edges = np.array([[0, 1], [1, 2], [2, 0]])
        edge_nodes = elements[:, local_edges]
        edge_nodes_sorted = np.sort(edge_nodes, axis=2)
        flat_edges = edge_nodes_sorted.reshape(-1, 2)
        _unique_edges, inverse = np.unique(flat_edges, axis=0, return_inverse=True)

        tri_ids = np.repeat(np.arange(n_tris), 3)
        local_edge_ids = np.tile(np.arange(3), n_tris)

        order = np.argsort(inverse, kind="stable")
        sorted_inverse = inverse[order]
        sorted_tri = tri_ids[order]
        sorted_local = local_edge_ids[order]

        same_as_next = sorted_inverse[:-1] == sorted_inverse[1:]
        pair_idx = np.where(same_as_next)[0]

        neighbours = np.full((n_tris, 3), -1, dtype=np.int64)
        neighbours[sorted_tri[pair_idx], sorted_local[pair_idx]] = sorted_tri[pair_idx + 1]
        neighbours[sorted_tri[pair_idx + 1], sorted_local[pair_idx + 1]] = sorted_tri[pair_idx]

        return vertices, elements, neighbours

    def evaluate_shape_derivatives(self, locations, elements=None):
        """
        Compute dN/ds (1st row), dN/dt(2nd row)
        """
        inside = None
        if elements is not None:
            inside = np.zeros(self.n_elements, dtype=bool)
            inside[elements] = True
        locations = np.array(locations)
        if elements is None:
            vertices, _c, tri, inside = self.get_element_for_location(locations)
        else:
            tri = elements
            M = np.ones((elements.shape[0], 3, 3))
            M[:, :, 1:] = self.vertices[self.elements[elements], :][:, :3, :]
            points_ = np.ones((locations.shape[0], 3))
            points_[:, 1:] = locations
            # minv = np.linalg.inv(M)
            # c = np.einsum("lij,li->lj", minv, points_)

        vertices = self.nodes[self.elements[tri][:, :3]]
        jac = np.zeros((tri.shape[0], 2, 2))
        jac[:, 0, 0] = vertices[:, 1, 0] - vertices[:, 0, 0]
        jac[:, 0, 1] = vertices[:, 1, 1] - vertices[:, 0, 1]
        jac[:, 1, 0] = vertices[:, 2, 0] - vertices[:, 0, 0]
        jac[:, 1, 1] = vertices[:, 2, 1] - vertices[:, 0, 1]
        # N = np.zeros((tri.shape[0], 6))

        # dN containts the derivatives of the shape functions
        dN = np.array([[-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]])

        # find the derivatives in x and y by calculating the dot product between the jacobian^-1 and the
        # derivative matrix
        #         d_n = np.einsum('ijk,ijl->ilk',np.linalg.inv(jac),dN)
        d_n = np.linalg.inv(jac)
        #         d_n = d_n.swapaxes(1,2)
        d_n = d_n @ dN
        # d_n = d_n.swapaxes(2, 1)
        # d_n = np.dot(np.linalg.inv(jac),dN)
        return d_n, tri, inside

    def evaluate_shape(self, locations):
        locations = np.array(locations)
        _vertices, c, tri, inside = self.get_element_for_location(locations, return_verts=False)
        # c = np.dot(np.array([1,x,y]),np.linalg.inv(M)) # convert to barycentric coordinates
        # order of bary coord is (1-s-t,s,t)
        N = c  # np.zeros((c.shape[0],3)) #evaluate shape functions at barycentric coordinates
        return N, tri, inside
