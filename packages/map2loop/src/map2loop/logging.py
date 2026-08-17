import logging

loggers = {}

ch = ch = logging.StreamHandler()
formatter = logging.Formatter("%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s")
ch.setFormatter(formatter)
ch.setLevel(logging.WARNING)
def get_levels():
    """dict for converting to logger levels from string


    Returns
    -------
    dict
        contains all strings with corresponding logging levels.
    """
    return {
        "info": logging.INFO,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "debug": logging.DEBUG,
    }


def getLogger(name: str):
    """Get a logger object with a specific name


    Parameters
    ----------
    name : str
        name of the logger object

    Returns
    -------
    logging.Logger
        logger object
    """
    if name in loggers:
        return loggers[name]
    logger = logging.getLogger(name)
    logger.addHandler(ch)
    logger.propagate = False
    loggers[name] = logger
    return logger


logger = getLogger(__name__)


def set_level(level: str):
    """Set the level of the logging object


    Parameters
    ----------
    level : str
        level of the logging object
    """
    levels = get_levels()
    level = levels.get(level, logging.WARNING)
    ch.setLevel(level)

    for name in loggers:
        logger = logging.getLogger(name)
        logger.setLevel(level)
    logger.info(f"Logging level set to {level}")


logger = getLogger(__name__)
logger.info("Imported map2loop.logging module")


def setLogging(level="info", handler=None):
    """Set the logging parameters for log file or custom handler.

    Parameters
    ----------
    level : str, optional
        Logging level to set, by default "info"
        Valid options: 'info', 'warning', 'error', 'debug'
    handler : logging.Handler, optional
        A logging handler to use instead of the default StreamHandler, by default None

    Examples
    --------
    >>> from map2loop.logging import setLogging
    >>> setLogging('debug')
    >>> setLogging('info', logging.FileHandler('loop.log'))
    """
    levels = get_levels()
    level_value = levels.get(level, logging.WARNING)

    # Create default handler if none provided
    if handler is None:
        handler = logging.StreamHandler()

    formatter = logging.Formatter(
        "%(levelname)s: %(asctime)s: %(filename)s:%(lineno)d -- %(message)s"
    )
    handler.setFormatter(formatter)
    handler.setLevel(level_value)

    # Replace handlers in all known loggers
    for name in loggers:
        logger = logging.getLogger(name)
        logger.handlers = []
        logger.addHandler(handler)
        logger.setLevel(level_value)

    # Also apply to main module logger
    main_logger = logging.getLogger(__name__)
    main_logger.handlers = []
    main_logger.addHandler(handler)
    main_logger.setLevel(level_value)

    main_logger.info(f"Set logging to {level}")
