import logging
import os
from typing import Any, Callable, Optional

import numpy as np

from .logging import get_logger


class LoopException(Exception):
    """Base class for LoopStructural exceptions."""


class LoopImportError(LoopException):
    def __init__(self, message, additional_information=None):
        super().__init__(message)
        self.additional_information = additional_information


class InterpolatorError(LoopException):
    pass


class LoopTypeError(LoopException):
    pass


class LoopValueError(LoopException):
    pass


class LogSink:
    def __init__(self, name: str = "sink") -> None:
        self.name = name

    def __call__(self, record: logging.LogRecord) -> None:
        return None


class StreamSink(LogSink):
    def __init__(self, stream=None, name: str = "stream") -> None:
        super().__init__(name=name)
        self.stream = stream


class FileSink(LogSink):
    def __init__(self, filename: str, name: str = "file") -> None:
        super().__init__(name=name)
        self.filename = filename


class SqliteSink(LogSink):
    def __init__(self, filename: str, name: str = "sqlite") -> None:
        super().__init__(name=name)
        self.filename = filename


_extra_sinks = []


def add_sink(sink: Callable[[logging.LogRecord], None]):
    _extra_sinks.append(sink)
    return sink


def remove_sink(sink: Callable[[logging.LogRecord], None]):
    if sink in _extra_sinks:
        _extra_sinks.remove(sink)
    return sink


def get_levels():
    return {"info": logging.INFO, "warning": logging.WARNING, "error": logging.ERROR, "debug": logging.DEBUG}


def getLogger(name: str):
    return get_logger(name)


def log_to_file(filename, overwrite=True, level="info"):
    logger = getLogger(__name__)
    if overwrite and os.path.isfile(filename):
        os.remove(filename)
    levels = get_levels()
    level_value = levels.get(level, logging.WARNING)
    handler = logging.FileHandler(filename)
    handler.setLevel(level_value)
    logger.addHandler(handler)
    logger.setLevel(level_value)
    return logger


def log_to_console(level="warning"):
    levels = get_levels()
    level_value = levels.get(level, logging.WARNING)
    logger = getLogger(__name__)
    logger.setLevel(level_value)
    return logger


def timed_stage(name: str):
    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        return func

    return decorator


def timed(func: Callable[..., Any]) -> Callable[..., Any]:
    return func


class EuclideanTransformation:
    def __init__(self, rotation=None, translation=None):
        self.rotation = rotation
        self.translation = translation


rng = np.random.default_rng()


def get_data_bounding_box_map(xyz, buffer):
    xyz = np.asarray(xyz, dtype=float)
    minx, maxx = xyz[:, 0].min(), xyz[:, 0].max()
    miny, maxy = xyz[:, 1].min(), xyz[:, 1].max()
    minz, maxz = xyz[:, 2].min(), xyz[:, 2].max()
    minx -= buffer
    maxx += buffer
    miny -= buffer
    maxy += buffer
    minz -= buffer
    maxz += buffer
    bb = np.array([[minx, miny, minz], [maxx, maxy, maxz]])

    def region(xyz):
        xyz = np.asarray(xyz, dtype=float)
        inside = np.ones(xyz.shape[0], dtype=bool)
        inside &= xyz[:, 0] > minx
        inside &= xyz[:, 0] < maxx
        inside &= xyz[:, 1] > miny
        inside &= xyz[:, 1] < maxy
        return inside

    return bb, region


def get_data_bounding_box(xyz, buffer):
    xyz = np.asarray(xyz, dtype=float)
    minx, maxx = xyz[:, 0].min(), xyz[:, 0].max()
    miny, maxy = xyz[:, 1].min(), xyz[:, 1].max()
    minz, maxz = xyz[:, 2].min(), xyz[:, 2].max()
    xlen = maxx - minx
    ylen = maxy - miny
    zlen = maxz - minz
    length = max([xlen, ylen, zlen])
    minx -= length * buffer
    maxx += length * buffer
    miny -= length * buffer
    maxy += length * buffer
    minz -= length * buffer
    maxz += length * buffer
    bb = np.array([[minx, miny, minz], [maxx, maxy, maxz]])

    def region(xyz):
        xyz = np.asarray(xyz, dtype=float)
        inside = np.ones(xyz.shape[0], dtype=bool)
        inside &= xyz[:, 0] > minx
        inside &= xyz[:, 0] < maxx
        inside &= xyz[:, 1] > miny
        inside &= xyz[:, 1] < maxy
        inside &= xyz[:, 2] > minz
        inside &= xyz[:, 2] < maxz
        return inside

    return bb, region


def create_surface(bounding_box, nstep):
    x = np.linspace(bounding_box[0, 0], bounding_box[1, 0], nstep[0])
    y = np.linspace(bounding_box[0, 1], bounding_box[1, 1], nstep[1])
    xx, yy = np.meshgrid(x, y, indexing="xy")

    def gi(i, j):
        return i + j * nstep[0]

    corners = np.array([[0, 1, 0, 1], [0, 0, 1, 1]])
    i = np.arange(0, nstep[0] - 1)
    j = np.arange(0, nstep[1] - 1)
    ii, jj = np.meshgrid(i, j, indexing="ij")
    corner_gi = gi(
        ii[:, :, None] + corners[None, None, 0, :],
        jj[:, :, None] + corners[None, None, 1, :],
    )
    corner_gi = corner_gi.reshape((nstep[0] - 1) * (nstep[1] - 1), 4)
    tri = np.vstack([corner_gi[:, :3], corner_gi[:, 1:]])
    return tri, xx.flatten(), yy.flatten()


def create_box(bounding_box, nsteps):
    from .geometry import BoundingBox

    if isinstance(bounding_box, BoundingBox):
        bounding_box = bounding_box.bb

    tri, xx, yy = create_surface(bounding_box[0:2, :], nsteps[0:2])
    zz = np.zeros(xx.shape)
    zz[:] = bounding_box[1, 2]
    tri = np.vstack([tri, tri + np.max(tri) + 1])
    xx = np.hstack([xx, xx])
    yy = np.hstack([yy, yy])
    z = np.zeros(zz.shape)
    z[:] = bounding_box[0, 2]
    zz = np.hstack([zz, z])
    t, x, z = create_surface(bounding_box[:, [0, 2]], nsteps[[0, 2]])
    tri = np.vstack([tri, t + np.max(tri) + 1])
    y = np.zeros(x.shape)
    y[:] = bounding_box[0, 1]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])
    tri = np.vstack([tri, t + np.max(tri) + 1])
    y[:] = bounding_box[1, 1]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])
    t, y, z = create_surface(bounding_box[:, [1, 2]], nsteps[[1, 2]])
    tri = np.vstack([tri, t + np.max(tri) + 1])
    x = np.zeros(y.shape)
    x[:] = bounding_box[0, 0]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])
    tri = np.vstack([tri, t + np.max(tri) + 1])
    x[:] = bounding_box[1, 0]
    xx = np.hstack([xx, x])
    zz = np.hstack([zz, z])
    yy = np.hstack([yy, y])
    points = np.zeros((len(xx), 3))
    points[:, 0] = xx
    points[:, 1] = yy
    points[:, 2] = zz
    return points, tri


__all__ = [
    "getLogger",
    "log_to_file",
    "log_to_console",
    "get_levels",
    "LogSink",
    "StreamSink",
    "FileSink",
    "SqliteSink",
    "add_sink",
    "remove_sink",
    "timed_stage",
    "timed",
    "EuclideanTransformation",
    "rng",
    "get_data_bounding_box",
    "get_data_bounding_box_map",
    "create_surface",
    "create_box",
    "LoopException",
    "LoopImportError",
    "InterpolatorError",
    "LoopTypeError",
    "LoopValueError",
]
