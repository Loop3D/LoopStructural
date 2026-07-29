# Make submodules available for import

from . import geometry
from . import io
from . import logging
from . import math
from . import supports

# Expose get_logger at the package level
from .logging.logger import get_logger
