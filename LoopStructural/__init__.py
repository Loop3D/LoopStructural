"""


LoopStructural API
=======================

.. autosummary::
    :toctree:


    interpolators
    modelling
    visualisation
    utils
"""
import logging
from logging.config import dictConfig

from .modelling.core.geological_model import GeologicalModel
from .visualisation.model_visualisation import LavaVuModelViewer
from .visualisation.map_viewer import MapView

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename='loop_structural_log.log',
                    filemode='w')
# define a Handler which writes INFO messages or higher to the sys.stderr
console = logging.StreamHandler()
console.setLevel(logging.INFO)
# add the handler to the root logger
logging.getLogger().addHandler(console)



# logging_config = dict(
#     version=1,
#     formatters={
#         'f': {'format':
#                   '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'}
#     },
#     handlers={
#         'h': {'class': 'logging.StreamHandler',
#               'formatter': 'f',
#               'level': logging.INFO
#               }
#     },
#     root={
#         'handlers': ['h'],
#         'level': logging.INFO,
#     },
# )
#
# dictConfig(logging_config)
