import logging

logger = logging.getLogger(__name__)

class LoopBaseException(Exception):
    """
    Base loop exception
    """
    # logger.error("Raising loop base exception")