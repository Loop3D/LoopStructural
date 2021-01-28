import logging

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class LoopBaseException(Exception):
    """
    Base loop exception
    """
    # logger.error("Raising loop base exception")