import logging
import time
import os

_format = '%(name)s: %(message)s'

def get_logger(name):

    # Setup the logger

    logger = logging.getLogger(os.path.splitext(name)[0])

    logger.setLevel(logging.INFO)

    # Prepare the handler

    handler = logging.StreamHandler()

    formatter = logging.Formatter(_format)

    handler.setFormatter(formatter)

    if len(logger.handlers) == 0:

        logger.addHandler(handler)

    else:

        # Logger already existed, no need to add a new handler
        pass

    logger.propagate = False

    logger.info("Setup logger %s at %s" % (name, time.asctime()))

    return logger
