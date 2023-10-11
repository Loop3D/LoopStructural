from __future__ import annotations
from typing import Optional
from LoopStructural.utils.exceptions import LoopValueError
import numpy as np


class BoundingBox:
    def __init__(
        self,
        dimensions: int = 3,
        origin: Optional[np.ndarray] = None,
        maximum: Optional[np.ndarray] = None,
        nsteps: Optional[np.ndarray] = None,
    ):
        self._origin = origin
        self._maximum = maximum
        self.dimensions = dimensions
        if nsteps is None:
            self.nsteps = np.array([50, 50, 25])
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
    def valid(self):
        return self._origin is not None and self._maximum is not None

    @property
    def origin(self) -> np.ndarray:
        if self._origin is None:
            raise LoopValueError("Origin is not set")
        return self._origin

    @origin.setter
    def origin(self, origin: np.ndarray):
        self._origin = origin

    @property
    def maximum(self) -> np.ndarray:
        if self._maximum is None:
            raise LoopValueError("Maximum is not set")
        return self._maximum

    @maximum.setter
    def maximum(self, maximum: np.ndarray):
        self._maximum = maximum

    @property
    def nelements(self):
        return self.nsteps.prod()

    @property
    def volume(self):
        return np.product(self.maximum - self.origin)

    @property
    def bb(self):
        return np.array([self.origin, self.maximum])

    @nelements.setter
    def nelements(self, nelements):
        box_vol = self.volume
        ele_vol = box_vol / nelements
        # calculate the step vector of a regular cube
        step_vector = np.zeros(3)
        step_vector[:] = ele_vol ** (1.0 / 3.0)
        # step_vector /= np.array([1,1,2])
        # number of steps is the length of the box / step vector
        nsteps = np.ceil((self.maximum - self.origin) / step_vector).astype(int)
        self.nsteps = nsteps

    @property
    def step_vector(self):
        return (self.maximum - self.origin) / self.nsteps

    def fit(self, locations: np.ndarray):
        if locations.shape[1] != self.dimensions:
            raise LoopValueError(
                f"locations array is {locations.shape[1]}D but bounding box is {self.dimensions}"
            )
        print("fitting")
        self.origin = locations.min(axis=0)
        self.maximum = locations.max(axis=0)
        return self

    def with_buffer(self, buffer: float = 0.2) -> BoundingBox:
        if self.origin is None or self.maximum is None:
            raise LoopValueError(
                "Cannot create bounding box with buffer, no origin or maximum"
            )
        origin = self.origin - buffer * (self.maximum - self.origin)
        maximum = self.maximum + buffer * (self.maximum - self.origin)
        return BoundingBox(origin=origin, maximum=maximum)

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
        inside = np.zeros(xyz.shape[0], dtype=bool)
        inside = np.logical_and(inside, xyz[:, 0] > self.origin[0])
        inside = np.logical_and(inside, xyz[:, 0] < self.maximum[0])
        inside = np.logical_and(inside, xyz[:, 1] > self.origin[1])
        inside = np.logical_and(inside, xyz[:, 1] < self.maximum[1])
        inside = np.logical_and(inside, xyz[:, 2] > self.origin[2])
        inside = np.logical_and(inside, xyz[:, 2] < self.maximum[2])
        return inside

    def regular_grid(self, nsteps=None, shuffle=False, order="C"):
        if nsteps is None:
            nsteps = self.nsteps
        x = np.linspace(self.origin[0], self.maximum[0], nsteps[0])
        y = np.linspace(self.origin[1], self.maximum[1], nsteps[1])
        z = np.linspace(self.origin[2], self.maximum[2], nsteps[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
        locs = np.array(
            [xx.flatten(order=order), yy.flatten(order=order), zz.flatten(order=order)]
        ).T
        if shuffle:
            # logger.info("Shuffling points")
            np.random.shuffle(locs)
        return locs
