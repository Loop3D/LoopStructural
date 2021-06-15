import logging

from LoopStructural.utils import getLogger
logger = getLogger(__name__)

class LoopBaseException(Exception):
    """
    Base loop exception
    """
    
class LoopImportError(LoopBaseException):
    def __init__(self,library,additional_information=None):
        error_message = 'Could not import: {}'.format(library)
        logger.error('Import error raise message: \"{}\"'.format(error_message))
        if additional_information:
            logger.error(additional_information)
        LoopBaseException.__init__(self,error_message)
