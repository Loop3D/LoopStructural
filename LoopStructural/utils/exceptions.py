from ..utils import getLogger

logger = getLogger(__name__)


class LoopException(Exception):
    """
    Base loop exception
    """


class LoopImportError(LoopException):
    """ """

    def __init__(self, message, additional_information=None):
        super().__init__(message)
        self.additional_information = additional_information

    pass


class InterpolatorError(LoopException):
    pass


class LoopTypeError(LoopException):
    pass


class LoopValueError(LoopException):
    pass
