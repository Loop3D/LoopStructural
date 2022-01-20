import logging

from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class LoopException(Exception):
    """
    Base loop exception
    """


class LoopImportError(LoopException):
    """ """

    pass


class InterpolatorError(LoopException):
    pass


class LoopTypeError(LoopException):
    pass


class LoopValueError(LoopException):
    pass
