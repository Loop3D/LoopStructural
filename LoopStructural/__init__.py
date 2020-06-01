"""


LoopStructural API
=======================

.. autosummary::


    interpolators
    modelling
    visualisation
    utils
"""
import logging
from logging.config import dictConfig

from .modelling.core.geological_model import GeologicalModel

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
