from __future__ import annotations
from typing import Optional, Union, Dict
from LoopStructural.utils.exceptions import LoopValueError
from LoopStructural.utils import rng
from LoopStructural.datatypes._structured_grid import StructuredGrid
import numpy as np
import copy

from LoopStructural.utils.logging import getLogger

logger = getLogger(__name__)


class BoundingBox:
    def __init__(
        self,
        origin: Optional[np.ndarray] = None,
        maximum: Optional[np.ndarray] = None,
        global_origin: Optional[np.ndarray] = None,
        global_maximum: Optional[np.ndarray] = None,
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
        if origin is not None and len(origin) != dimensions:
            logger.warning(
                f"Origin has {len(origin)} dimensions but bounding box has {dimensions}"
            )
            raise LoopValueError("Origin has incorrect number of dimensions")
        if maximum is not None and len(maximum) != dimensions:
            logger.warning(
                f"Maximum has {len(maximum)} dimensions but bounding box has {dimensions}"
            )
            raise LoopValueError("Maximum has incorrect number of dimensions")
        if global_origin is not None and len(global_origin) != dimensions:
            logger.warning(
                f"Global origin has {len(global_origin)} dimensions but bounding box has {dimensions}"
            )
            raise LoopValueError("Global origin has incorrect number of dimensions")
        if nsteps is not None and len(nsteps) != dimensions:
            logger.warning(
                f"Nsteps has {len(nsteps)} dimensions but bounding box has {dimensions}"
            )
            raise LoopValueError("Nsteps has incorrect number of dimensions")
        # reproject relative to the global origin, if origin is not provided.
        # we want the local coordinates to start at 0
        # otherwise uses provided origin. This is useful for having multiple bounding boxes rela
        if global_origin is not None and origin is None:
            origin = np.zeros(np.array(global_origin).shape, dtype=float)
        if global_maximum is not None and global_origin is not None:
            maximum = np.array(global_maximum,dtype=float) - np.array(global_origin,dtype=float)

        if maximum is None and nsteps is not None and step_vector is not None:
            maximum = np.array(origin) + np.array(nsteps) * np.array(step_vector)
        if origin is not None and global_origin is None:
            global_origin = np.zeros(3)
        self._origin = np.array(origin, dtype=float)
        self._maximum = np.array(maximum, dtype=float)
        self.dimensions = dimensions
        if self.origin.shape:
            if self.origin.shape[0] != self.dimensions:
                logger.warning(
                    f"Origin has {self.origin.shape[0]} dimensions but bounding box has {self.dimensions}"
                )

        else:
            self.dimensions = dimensions
        self._global_origin = global_origin
        if self.origin is not None and self.maximum is not None:
            self.nelements = 10_000
        else:
            self.nsteps = np.array([50, 50, 25])
        if nsteps is not None:
            self.nsteps = np.array(nsteps)
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

    @property
    def global_origin(self):
        """Get the global origin of the bounding box.

        Returns
        -------
        np.ndarray
            The global origin coordinates
        """
        return self._global_origin

    @global_origin.setter
    def global_origin(self, global_origin):
        """Set the global origin of the bounding box.

        Parameters
        ----------
        global_origin : array_like
            The global origin coordinates
        """
        if self.dimensions != len(global_origin):
            logger.warning(
                f"Global origin has {len(global_origin)} dimensions but bounding box has {self.dimensions}"
            )
        self._global_origin = global_origin

    @property
    def global_maximum(self):
        """Get the global maximum coordinates of the bounding box.

        Returns
        -------
        np.ndarray
            The global maximum coordinates (local maximum + global origin)
        """
        return self.maximum + self.global_origin

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
        self._origin = origin

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
        self._maximum = maximum

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
        return np.array(
            [
                self.global_origin.tolist(),
                [self.global_maximum[0], self.global_origin[1], self.global_origin[2]],
                [self.global_maximum[0], self.global_maximum[1], self.global_origin[2]],
                [self.global_origin[0], self.global_maximum[1], self.global_origin[2]],
                [self.global_origin[0], self.global_origin[1], self.global_maximum[2]],
                [self.global_maximum[0], self.global_origin[1], self.global_maximum[2]],
                self.global_maximum.tolist(),
                [self.global_origin[0], self.global_maximum[1], self.global_maximum[2]],
            ]
        )

    @property
    def step_vector(self):
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
        if local_coordinate:
            self.global_origin = origin
            self.origin = np.zeros(3)
            self.maximum = maximum - origin
        else:
            self.origin = origin
            self.maximum = maximum
            self.global_origin = np.zeros(3)
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
        return BoundingBox(
            origin=origin,
            maximum=maximum,
            global_origin=self.global_origin,
            nsteps=self.nsteps,
            dimensions=self.dimensions,
        )

    # def __call__(self, xyz):
    #     xyz = np.array(xyz)
    #     if len(xyz.shape) == 1:
    #         xyz = xyz.reshape((1, -1))

    #     distances = np.maximum(0,
    #                         np.maximum(self.global_origin+self.origin - xyz,
    #                                 xyz - self.global_maximum))
    #     distance = np.linalg.norm(distances, axis=1)
    #     distance[self.is_inside(xyz)] = -1
    #     return distance

    def __call__(self, xyz):
        # Calculate center and half-extents of the box
        center = (self.maximum + self.global_origin + self.origin) / 2
        half_extents = (self.maximum - self.global_origin + self.origin) / 2

        # Calculate the distance from point to center
        offset = np.abs(xyz - center) - half_extents

        # Inside distance: negative value based on the smallest penetration
        inside_distance = np.min(half_extents - np.abs(xyz - center), axis=1)

        # Outside distance: length of the positive components of offset
        outside_distance = np.linalg.norm(np.maximum(offset, 0))

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

        if not local:
            coordinates = [
                np.linspace(self.global_origin[i]+self.origin[i], self.global_maximum[i], nsteps[i])
                for i in range(self.dimensions)
            ]
        coordinate_grid = np.meshgrid(*coordinates, indexing="ij")
        locs = np.array([coord.flatten(order=order) for coord in coordinate_grid]).T

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
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'BoundingBox':
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
        return cls(
            origin=np.array(data["origin"]),
            maximum=np.array(data["maximum"]),
            nsteps=np.array(data["nsteps"]),
        )

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
        x = np.linspace(
            self.global_origin[0] + self.origin[0], self.global_maximum[0], self.nsteps[0]
        )
        y = np.linspace(
            self.global_origin[1] + self.origin[1], self.global_maximum[1], self.nsteps[1]
        )
        z = np.linspace(
            self.global_origin[2] + self.origin[2], self.global_maximum[2], self.nsteps[2]
        )
        return pv.RectilinearGrid(
            x,
            y,
            z,
        )

    def structured_grid(
        self, cell_data: Dict[str, np.ndarray] = {}, vertex_data={}, name: str = "bounding_box"
    ):
        # python is passing a reference to the cell_data, vertex_data dicts so we need to
        # copy them to make sure that different instances of StructuredGrid are not sharing the same
        # underlying objects
        _cell_data = copy.deepcopy(cell_data)
        _vertex_data = copy.deepcopy(vertex_data)
        return StructuredGrid(
            origin=self.global_origin + self.origin,
            step_vector=self.step_vector,
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
        if inplace:
            xyz -= self.global_origin
            return xyz
        return (xyz - self.global_origin)  # np.clip(xyz, self.origin, self.maximum)

    def scale_by_projection_factor(self, value):
        return value / np.max((self.global_maximum - self.global_origin))

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
        if inplace:
            xyz += self.global_origin
            return xyz
        return xyz + self.global_origin

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
        """Get the transformation matrix from local to global coordinates

        Returns
        -------
        np.ndarray
            4x4 transformation matrix
        """
        matrix = np.eye(4)
        L = self.global_maximum - self.global_origin
        L = np.max(L)
        matrix[0, 3] = -self.global_origin[0]/L
        matrix[1, 3] = -self.global_origin[1]/L
        matrix[2, 3] = -self.global_origin[2]/L
        if normalise:
            matrix[0,0] = 1/L
            matrix[1,1] = 1/L
            matrix[2,2] = 1/L
        return matrix