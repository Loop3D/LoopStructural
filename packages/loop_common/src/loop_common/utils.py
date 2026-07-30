import logging
import os

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


rng = np.random.default_rng()


__all__ = [
    "getLogger",
    "log_to_file",
    "log_to_console",
    "get_levels",
    "rng",
    "LoopException",
    "LoopImportError",
    "InterpolatorError",
    "LoopTypeError",
    "LoopValueError",
]
