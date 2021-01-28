import logging
import LoopStructural

def get_levels():
    """dict for converting to logger levels from string


    Returns
    -------
    dict
        contains all strings with corresponding logging levels.
    """
    return {'info':logging.INFO,'warning':logging.WARNING,'error':logging.ERROR,'debug':logging.DEBUG}

def getLogger(name):
    logger = logging.getLogger(name)
    logger.addHandler(LoopStructural.ch)
    LoopStructural.loggers[name] = logger
    return logger

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
    fh = logging.FileHandler(filename)
    fh.setFormatter(LoopStructural.formatter)
    fh.setLevel(level)
    for logger in LoopStructural.loggers.values():
        for hdlr in logger.handlers[:]:  # remove the existing file handlers
            if isinstance(hdlr,logging.FileHandler): #fixed two typos here
                logger.removeHandler(hdlr)
        logger.addHandler(fh) 
        logger.setLevel(level)
    

def log_to_console(level='warning'):
    """Set the level of logging to the console


    Parameters
    ----------
    level : str, optional
        'info', 'warning', 'error', 'debug' mapped to logging levels, by default 'info'
    """
    levels = get_levels()
    level = levels.get(level,logging.WARNING)
    for logger in LoopStructural.loggers.values():
        for hdlr in logger.handlers:
            # both stream and file are base stream, so check if not a filehandler
            if not isinstance(hdlr,logging.FileHandler):
                hdlr.setLevel(level)
