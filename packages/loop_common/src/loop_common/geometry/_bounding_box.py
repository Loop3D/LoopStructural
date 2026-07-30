from __future__ import annotations
from typing import Optional, Union, Dict

# from LoopStructural.utils.exceptions import LoopValueError
from loop_common.math import rng
from loop_common.supports import StructuredGrid
import numpy as np
import copy

from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)


class LoopValueError(ValueError):
    """Custom error for invalid values in LoopStructural."""

    pass


class BoundingBox:
    def __init__(
        self,
        origin: Optional[np.ndarray] = None,
        maximum: Optional[np.ndarray] = None,
        nsteps: Optional[np.ndarray] = None,
        step_vector: Optional[np.ndarray] = None,
        dimensions: Optional[int] = 3,
    ):
        """A bounding box for a model, defined by the
        origin, maximum and number of steps in each direction

        Parameters
        ----------
        dimensions : int, optional
            _description_, by default 3
        origin : Optional[np.ndarray], optional
            _description_, by default None
        maximum : Optional[np.ndarray], optional
            _description_, by default None
        nsteps : Optional[np.ndarray], optional
            _description_, by default None
        """
        self.dimensions = dimensions

        def _coerce_point(point, name):
            if point is None:
                return None
            arr = np.asarray(point, dtype=float)
            if arr.shape != (self.dimensions,):
                logger.warning(
                    f"{name} has shape {arr.shape} but bounding box has {self.dimensions} dimensions"
                )
                raise LoopValueError(f"{name} has incorrect number of dimensions")
            return arr

        origin = _coerce_point(origin, "Origin")
        maximum = _coerce_point(maximum, "Maximum")

        if nsteps is not None:
            if len(nsteps) != dimensions:
                logger.warning(f"Nsteps has {len(nsteps)} dimensions but bounding box has {dimensions}")
                raise LoopValueError("Nsteps has incorrect number of dimensions")
            if np.any(np.asarray(nsteps) <= 0):
                raise LoopValueError("Nsteps must be positive integers")

        if (
            maximum is None
            and nsteps is not None
            and step_vector is not None
            and origin is not None
        ):
            maximum = np.asarray(origin, dtype=float) + np.asarray(nsteps) * np.asarray(
                step_vector, dtype=float
            )

        if origin is not None and maximum is not None and np.any(maximum < origin):
            raise LoopValueError("Maximum must be greater than or equal to origin")

        self._origin = origin
        self._maximum = maximum

        # Local interpolation coordinate frame (world -> local affine transform).
        self._world_to_local = np.eye(4)
        self._local_to_world = np.eye(4)
        self._local_origin = np.zeros(self.dimensions, dtype=float)
        self._local_rotation = np.eye(self.dimensions, dtype=float)

        if self.valid:
            self.nelements = 10_000
        else:
            default_nsteps = np.ones(self.dimensions, dtype=int) * 50
            if self.dimensions == 3:
                default_nsteps[-1] = 25
            self.nsteps = default_nsteps
        if nsteps is not None:
            self.nsteps = np.array(nsteps)

        self.set_local_transform(local_origin=np.zeros(self.dimensions, dtype=float))

        self.name_map = {
            "xmin": (0, 0),
            "ymin": (0, 1),
            "zmin": (0, 2),
            "xmax": (1, 0),
            "ymax": (1, 1),
            "zmax": (1, 2),
            "lower": (0, 2),
            "upper": (1, 2),
            "minx": (0, 0),
            "miny": (0, 1),
            "minz": (0, 2),
            "maxx": (1, 0),
            "maxy": (1, 1),
            "maxz": (1, 2),
        }

    def set_local_transform(
        self,
        local_origin: Optional[np.ndarray] = None,
        rotation_matrix: Optional[np.ndarray] = None,
    ):
        """Set the world->local affine transform used for interpolation coordinates.

        Parameters
        ----------
        local_origin : Optional[np.ndarray]
            World-space origin of the local frame. If None, uses zeros.
        rotation_matrix : Optional[np.ndarray]
            Rotation matrix mapping world axes to local axes.
        """
        if local_origin is None:
            local_origin = np.zeros(self.dimensions, dtype=float)
        local_origin = np.asarray(local_origin, dtype=float)
        if local_origin.shape != (self.dimensions,):
            raise LoopValueError("Local origin has incorrect number of dimensions")

        if rotation_matrix is None:
            rotation_matrix = np.eye(self.dimensions, dtype=float)
        rotation_matrix = np.asarray(rotation_matrix, dtype=float)
        if rotation_matrix.shape != (self.dimensions, self.dimensions):
            raise LoopValueError(
                f"Rotation matrix must have shape ({self.dimensions}, {self.dimensions})"
            )

        self._local_origin = local_origin
        self._local_rotation = rotation_matrix

        world_to_local = np.eye(4)
        world_to_local[: self.dimensions, : self.dimensions] = rotation_matrix
        world_to_local[: self.dimensions, 3] = -rotation_matrix @ local_origin

        self._world_to_local = world_to_local
        self._local_to_world = np.linalg.inv(world_to_local)

    def _apply_affine(
        self, xyz: np.ndarray, matrix: np.ndarray, inplace: bool = False
    ) -> np.ndarray:
        arr = np.asarray(xyz, dtype=float)
        is_vector = arr.ndim == 1
        points = arr.reshape(1, -1) if is_vector else arr
        if points.shape[1] != self.dimensions:
            raise LoopValueError(
                f"locations array is {points.shape[1]}D but bounding box is {self.dimensions}"
            )

        hom = np.ones((points.shape[0], 4), dtype=float)
        hom[:, : self.dimensions] = points
        transformed = (matrix @ hom.T).T[:, : self.dimensions]

        if inplace and isinstance(xyz, np.ndarray):
            xyz[...] = transformed.reshape(arr.shape)
            return xyz
        if is_vector:
            return transformed[0]
        return transformed

    @property
    def local_origin(self):
        """World-space origin of the local interpolation frame."""
        return self._local_origin.copy()

    @property
    def local_rotation(self):
        """Rotation matrix that maps world coordinates into local coordinates."""
        return self._local_rotation.copy()

    @property
    def world_to_local_matrix(self):
        """Homogeneous 4x4 matrix for world -> local coordinates."""
        return self._world_to_local.copy()

    @property
    def local_to_world_matrix(self):
        """Homogeneous 4x4 matrix for local -> world coordinates."""
        return self._local_to_world.copy()

    @property
    def valid(self):
        """Check if the bounding box has valid origin and maximum values.

        Returns
        -------
        bool
            True if both origin and maximum are set, False otherwise
        """
        return self._origin is not None and self._maximum is not None

    @property
    def origin(self) -> np.ndarray:
        """Get the origin coordinates of the bounding box.

        Returns
        -------
        np.ndarray
            Origin coordinates

        Raises
        ------
        LoopValueError
            If the origin is not set
        """
        if self._origin is None:
            raise LoopValueError("Origin is not set")
        return self._origin

    @origin.setter
    def origin(self, origin: np.ndarray):
        """Set the origin coordinates of the bounding box.

        Parameters
        ----------
        origin : np.ndarray
            Origin coordinates
        """
        if self.dimensions != len(origin):
            logger.warning(
                f"Origin has {len(origin)} dimensions but bounding box has {self.dimensions}"
            )
        self._origin = np.asarray(origin, dtype=float)

    @property
    def maximum(self) -> np.ndarray:
        """Get the maximum coordinates of the bounding box.

        Returns
        -------
        np.ndarray
            Maximum coordinates

        Raises
        ------
        LoopValueError
            If the maximum is not set
        """
        if self._maximum is None:
            raise LoopValueError("Maximum is not set")
        return self._maximum

    @maximum.setter
    def maximum(self, maximum: np.ndarray):
        """Set the maximum coordinates of the bounding box.

        Parameters
        ----------
        maximum : np.ndarray
            Maximum coordinates
        """
        self._maximum = np.asarray(maximum, dtype=float)

    @property
    def nelements(self):
        """Get the total number of elements in the bounding box.

        Returns
        -------
        int
            Total number of elements (product of nsteps)
        """

        return self.nsteps.prod()

    @property
    def volume(self):
        """Calculate the volume of the bounding box.

        Returns
        -------
        float
            Volume of the bounding box
        """
        length = self.maximum - self.origin
        length = length.astype(float)
        return np.prod(length)

    @property
    def bb(self):
        """Get a numpy array containing origin and maximum coordinates.

        Returns
        -------
        np.ndarray
            Array with shape (2, n_dimensions) containing [origin, maximum]
        """
        return np.array([self.origin, self.maximum])

    @nelements.setter
    def nelements(self, nelements: Union[int, float]):
        """Update the number of elements in the associated grid
        This is for visualisation, not for the interpolation
        When set it will update the nsteps/step vector for cubic
        elements

        Parameters
        ----------
        nelements : int,float
            The new number of elements
        """
        box_vol = self.volume
        ele_vol = box_vol / nelements
        # calculate the step vector of a regular cube
        step_vector = np.zeros(self.dimensions)
        if self.dimensions == 2:
            step_vector[:] = ele_vol ** (1.0 / 2.0)
        elif self.dimensions == 3:
            step_vector[:] = ele_vol ** (1.0 / 3.0)
        else:
            logger.warning("Can only set nelements for 2d or 3D bounding box")
            return
        # number of steps is the length of the box / step vector
        nsteps = np.ceil((self.maximum - self.origin) / step_vector).astype(int)
        self.nsteps = nsteps

    @property
    def corners(self) -> np.ndarray:
        """Returns the corners of the bounding box in local coordinates



        Returns
        -------
        np.ndarray
            array of corners in clockwise order
        """

        return np.array(
            [
                self.origin.tolist(),
                [self.maximum[0], self.origin[1], self.origin[2]],
                [self.maximum[0], self.maximum[1], self.origin[2]],
                [self.origin[0], self.maximum[1], self.origin[2]],
                [self.origin[0], self.origin[1], self.maximum[2]],
                [self.maximum[0], self.origin[1], self.maximum[2]],
                self.maximum.tolist(),
                [self.origin[0], self.maximum[1], self.maximum[2]],
            ]
        )

    @property
    def corners_global(self) -> np.ndarray:
        """Returns the corners of the bounding box
        in the original space

        Returns
        -------
        np.ndarray
            corners of the bounding box
        """
        return self.corners

    @property
    def step_vector(self):
        if np.any(self.nsteps == 0):
            raise LoopValueError("Cannot compute step_vector: nsteps contains zero values")
        return (self.maximum - self.origin) / self.nsteps

    @property
    def length(self):
        return self.maximum - self.origin

    def fit(self, locations: np.ndarray, local_coordinate: bool = False) -> BoundingBox:
        """Initialise the bounding box from a set of points.

        Parameters
        ----------
        locations : np.ndarray
            xyz locations of the points to fit the bbox
        local_coordinate : bool, optional
            whether to set the origin to [0,0,0], by default False

        Returns
        -------
        BoundingBox
            A reference to the bounding box object, note this is not a new bounding box
            it updates the current one in place.

        Raises
        ------
        LoopValueError
            _description_
        """
        if locations.shape[1] != self.dimensions:
            raise LoopValueError(
                f"locations array is {locations.shape[1]}D but bounding box is {self.dimensions}"
            )
        origin = locations.min(axis=0)
        maximum = locations.max(axis=0)
        origin = np.array(origin)
        maximum = np.array(maximum)
        self.origin = origin
        self.maximum = maximum
        if local_coordinate:
            self.set_local_transform(local_origin=origin)
        else:
            self.set_local_transform(local_origin=np.zeros(self.dimensions, dtype=float))
        return self

    def with_buffer(self, buffer: float = 0.2) -> BoundingBox:
        """Create a new bounding box with a buffer around the existing bounding box

        Parameters
        ----------
        buffer : float, optional
            percentage to expand the dimensions by, by default 0.2

        Returns
        -------
        BoundingBox
            The new bounding box object.

        Raises
        ------
        LoopValueError
            if the current bounding box is invalid
        """
        if self.origin is None or self.maximum is None:
            raise LoopValueError("Cannot create bounding box with buffer, no origin or maximum")
        # local coordinates, rescale into the original bounding boxes global coordinates
        origin = self.origin - buffer * np.max(self.maximum - self.origin)
        maximum = self.maximum + buffer * np.max(self.maximum - self.origin)
        buffered = BoundingBox(
            origin=origin,
            maximum=maximum,
            nsteps=self.nsteps,
            dimensions=self.dimensions,
        )
        buffered.set_local_transform(
            local_origin=self.local_origin,
            rotation_matrix=self.local_rotation,
        )
        return buffered

    def __call__(self, xyz):
        xyz = np.asarray(xyz, dtype=float)
        if xyz.ndim == 1:
            xyz = xyz[None, :]
        # Calculate center and half-extents of the box
        center = (self.maximum + self.origin) / 2
        half_extents = (self.maximum - self.origin) / 2

        # Calculate the distance from point to center
        offset = np.abs(xyz - center) - half_extents

        # Inside distance: negative value based on the smallest penetration
        inside_distance = np.min(half_extents - np.abs(xyz - center), axis=1)

        # Outside distance: length of the positive components of offset
        outside_distance = np.linalg.norm(np.maximum(offset, 0), axis=1)

        # If any component of offset is positive, we're outside
        # Otherwise, we're inside and return the negative penetration distance
        distance = np.zeros(xyz.shape[0])
        mask = np.any(offset > 0, axis=1)
        distance[mask] = outside_distance
        distance[~mask] = -inside_distance[~mask]
        return distance
        # return outside_distance if np.any(offset > 0) else -inside_distance

    def get_value(self, name):
        ix, iy = self.name_map.get(name, (-1, -1))
        if ix == -1 and iy == -1:
            raise LoopValueError(f"{name} is not a valid bounding box name")
        if iy == -1:
            return self.origin[ix]

        return self.bb[ix,]

    def __getitem__(self, name):
        if isinstance(name, str):
            return self.get_value(name)
        elif isinstance(name, tuple):
            return self.origin
        return self.get_value(name)

    def is_inside(self, xyz):
        xyz = np.array(xyz)
        if len(xyz.shape) == 1:
            xyz = xyz.reshape((1, -1))
        if xyz.shape[1] != 3:
            raise LoopValueError(
                f"locations array is {xyz.shape[1]}D but bounding box is {self.dimensions}"
            )
        inside = np.ones(xyz.shape[0], dtype=bool)
        inside = np.logical_and(inside, xyz[:, 0] > self.origin[0])
        inside = np.logical_and(inside, xyz[:, 0] < self.maximum[0])
        inside = np.logical_and(inside, xyz[:, 1] > self.origin[1])
        inside = np.logical_and(inside, xyz[:, 1] < self.maximum[1])
        inside = np.logical_and(inside, xyz[:, 2] > self.origin[2])
        inside = np.logical_and(inside, xyz[:, 2] < self.maximum[2])
        return inside

    def regular_grid(
        self,
        nsteps: Optional[Union[list, np.ndarray]] = None,
        shuffle: bool = False,
        order: str = "F",
        local: bool = True,
    ) -> np.ndarray:
        """Get the grid of points from the bounding box

        Parameters
        ----------
        nsteps : Optional[Union[list, np.ndarray]], optional
            number of steps, by default None uses self.nsteps
        shuffle : bool, optional
            Whether to return points in order or random, by default False
        order : str, optional
            when flattening using numpy "C" or "F", by default "C"
        local : bool, optional
            Whether to return the points in the local coordinate system of global
            , by default True

        Returns
        -------
        np.ndarray
            numpy array N,3 of the points
        """

        if nsteps is None:
            nsteps = self.nsteps
        coordinates = [
            np.linspace(self.origin[i], self.maximum[i], nsteps[i]) for i in range(self.dimensions)
        ]
        coordinate_grid = np.meshgrid(*coordinates, indexing="ij")
        locs = np.array([coord.flatten(order=order) for coord in coordinate_grid]).T

        if local:
            locs = self.project(locs)

        if shuffle:
            # logger.info("Shuffling points")
            rng.shuffle(locs)
        return locs

    def cell_centres(self, order: str = "F") -> np.ndarray:
        """Get the cell centres of a regular grid

        Parameters
        ----------
        order : str, optional
            order of the grid, by default "C"

        Returns
        -------
        np.ndarray
            array of cell centres
        """
        locs = self.regular_grid(order=order, nsteps=self.nsteps - 1)

        return locs + 0.5 * self.step_vector

    def to_dict(self) -> dict:
        """Export the defining characteristics of the bounding
        box to a dictionary for json serialisation

        Returns
        -------
        dict
            dictionary with origin, maximum and nsteps
        """
        return {
            "origin": self.origin.tolist(),
            "maximum": self.maximum.tolist(),
            "nsteps": self.nsteps.tolist(),
            "local_origin": self.local_origin.tolist(),
            "local_rotation": self.local_rotation.tolist(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "BoundingBox":
        """Create a bounding box from a dictionary

        Parameters
        ----------
        data : dict
            dictionary with origin, maximum and nsteps

        Returns
        -------
        BoundingBox
            bounding box object
        """
        bbox = cls(
            origin=np.array(data["origin"]),
            maximum=np.array(data["maximum"]),
            nsteps=np.array(data["nsteps"]),
        )
        if "local_origin" in data or "local_rotation" in data:
            bbox.set_local_transform(
                local_origin=np.array(data.get("local_origin", np.zeros(bbox.dimensions))),
                rotation_matrix=np.array(
                    data.get("local_rotation", np.eye(bbox.dimensions).tolist())
                ),
            )
        return bbox

    def vtk(self):
        """Export the model as a pyvista RectilinearGrid

        Returns
        -------
        pv.RectilinearGrid
            a pyvista grid object

        Raises
        ------
        ImportError
            If pyvista is not installed raise import error
        """
        try:
            import pyvista as pv
        except ImportError:
            raise ImportError("pyvista is required for vtk support")
        x = np.linspace(self.origin[0], self.maximum[0], self.nsteps[0])
        y = np.linspace(self.origin[1], self.maximum[1], self.nsteps[1])
        z = np.linspace(self.origin[2], self.maximum[2], self.nsteps[2])
        return pv.RectilinearGrid(
            x,
            y,
            z,
        )

    def structured_grid(
        self,
        cell_data: Optional[Dict[str, np.ndarray]] = None,
        vertex_data: Optional[Dict] = None,
        name: str = "bounding_box",
        local_coordinates: bool = False,
    ):
        # python is passing a reference to the cell_data, vertex_data dicts so we need to
        # copy them to make sure that different instances of StructuredGrid are not sharing the same
        # underlying objects
        if cell_data is None:
            cell_data = {}
        if vertex_data is None:
            vertex_data = {}
        _cell_data = copy.deepcopy(cell_data)
        _vertex_data = copy.deepcopy(vertex_data)
        if local_coordinates:
            # Project origin/maximum directly (rather than all corners, which
            # only supports 3D) -- exact for translation-only transforms.
            local_points = self.project(np.array([self.origin, self.maximum]))
            local_origin = np.min(local_points, axis=0)
            local_maximum = np.max(local_points, axis=0)
            step_vector = (local_maximum - local_origin) / self.nsteps
            origin = local_origin
        else:
            step_vector = self.step_vector
            origin = self.origin
        return StructuredGrid(
            origin=origin,
            step_vector=step_vector,
            nsteps=self.nsteps,
            cell_properties=_cell_data,
            properties=_vertex_data,
            name=name,
        )

    def project(self, xyz, inplace=False):
        """Project a point into the bounding box

        Parameters
        ----------
        xyz : np.ndarray
            point to project
        inplace : bool, optional
            Whether to modify the input array in place, by default False

        Returns
        -------
        np.ndarray
            projected point
        """
        return self._apply_affine(xyz, self.world_to_local_matrix, inplace=inplace)

    def project_vectors(self, vectors: np.ndarray) -> np.ndarray:
        """Rotate vectors from world frame into local frame."""
        arr = np.asarray(vectors, dtype=float)
        is_vector = arr.ndim == 1
        vec = arr.reshape(1, -1) if is_vector else arr
        if vec.shape[1] != self.dimensions:
            raise LoopValueError(
                f"vector array is {vec.shape[1]}D but bounding box is {self.dimensions}"
            )
        projected = (self.local_rotation @ vec.T).T
        return projected[0] if is_vector else projected

    def scale_by_projection_factor(self, value):
        return value / np.max((self.maximum - self.origin))

    def reproject(self, xyz, inplace=False):
        """Reproject a point from the bounding box to the global space

        Parameters
        ----------
        xyz : np.ndarray
            point to reproject
        inplace : bool, optional
            Whether to modify the input array in place, by default False
        Returns
        -------
        np.ndarray
            reprojected point
        """
        return self._apply_affine(xyz, self.local_to_world_matrix, inplace=inplace)

    def reproject_vectors(self, vectors: np.ndarray) -> np.ndarray:
        """Rotate vectors from local frame back into world frame."""
        arr = np.asarray(vectors, dtype=float)
        is_vector = arr.ndim == 1
        vec = arr.reshape(1, -1) if is_vector else arr
        if vec.shape[1] != self.dimensions:
            raise LoopValueError(
                f"vector array is {vec.shape[1]}D but bounding box is {self.dimensions}"
            )
        rotation = self.local_to_world_matrix[: self.dimensions, : self.dimensions]
        reprojected = (rotation @ vec.T).T
        return reprojected[0] if is_vector else reprojected

    def __repr__(self):
        return f"BoundingBox(origin:{self.origin}, maximum:{self.maximum}, nsteps:{self.nsteps})"

    def __str__(self):
        return f"BoundingBox(origin:{self.origin}, maximum:{self.maximum}, nsteps:{self.nsteps})"

    def __eq__(self, other):
        if not isinstance(other, BoundingBox):
            return False
        return (
            np.allclose(self.origin, other.origin)
            and np.allclose(self.maximum, other.maximum)
            and np.allclose(self.nsteps, other.nsteps)
        )

    def matrix(self, normalise: bool = False) -> np.ndarray:
        """Get the world-to-local transformation matrix.

        Returns
        -------
        np.ndarray
            4x4 transformation matrix
        """
        matrix = self.world_to_local_matrix
        if normalise:
            L = np.max(self.maximum - self.origin)
            if L > 0:
                matrix[: self.dimensions, : self.dimensions] /= L
                matrix[: self.dimensions, 3] /= L
        return matrix
