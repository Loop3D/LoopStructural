import logging

import numpy as np

logger = logging.getLogger(__name__)


def strike_symbol(strike):
    R = np.zeros((2, 2))
    R[0, 0] = np.cos(np.deg2rad(-strike))
    R[0, 1] = -np.sin(np.deg2rad(-strike))
    R[1, 0] = np.sin(np.deg2rad(-strike))
    R[1, 1] = np.cos(np.deg2rad(-strike))
    R = np.zeros((2, 2))
    R[0, 0] = np.cos(np.deg2rad(-strike))
    R[0, 1] = -np.sin(np.deg2rad(-strike))
    R[1, 0] = np.sin(np.deg2rad(-strike))
    R[1, 1] = np.cos(np.deg2rad(-strike))

    vec = np.array([0, 1])
    rotated = R @ vec
    vec2 = np.array([-0.5, 0])
    r2 = R @ vec2
    return rotated, r2
def get_levels():
    """dict for converting to logger levels from string


    Returns
    -------
    dict
        contains all strings with corresponding logging levels.
    """
    return {'info':logging.INFO,'warning':logging.WARNING,'error':logging.ERROR,'debug':logging.DEBUG}

def log_to_file(filename,level='info'):
    """Set the logging parameters for log file


    Parameters
    ----------
    filename : string
        name of file or path to file
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    levels = get_levels()
    level = levels.get(level,logging.WARNING)
    logging.basicConfig(level=level,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    filename=filename,
                    filemode='w')

def log_to_console(level='warning'):
    """Set the level of logging to the console


    Parameters
    ----------
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    levels = get_levels()
    level = levels.get(level,logging.WARNING)

    changed_level = False
    for h in logging.getLogger().handlers:
        if type(h) is logging.StreamHandler:
            h.setLevel(level)
            changed_level = True
    if not changed_level:
        console = logging.StreamHandler()
        console.setLevel(level)
        # add the handler to the root logger
        logging.getLogger().addHandler(console)