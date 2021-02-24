from LoopStructural.utils import getLogger
logger = getLogger(__name__)

try:
    import theano
except ImportError:
    logger.error("Cannot use LoopStructural.probability without theano \n"
                "Please install theano and try again")
from ._normal import normal
from ._theano_wrapper import LogLikelihood, LogLikelihoodGradient