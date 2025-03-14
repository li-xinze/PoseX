import logging
import os
import contextlib, tempfile
from pathlib import Path
from functools import partial
from subprocess import check_output as _call
run_shell_cmd = partial(_call, shell=True)

DEBUG_MODE = bool(int(os.environ.get("DEBUG_MODE", 0)))


def create_logger(filename):
    logger_ = logging.getLogger(filename)  # type: logging.Logger

    if DEBUG_MODE:
        logger_.setLevel(logging.DEBUG)
    else:
        logger_.setLevel(logging.INFO)

    handler = logging.StreamHandler()
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger_.addHandler(handler)
    return logger_




@contextlib.contextmanager
def temprary_filename(mode='w+b', suffix=None):
    tmp_name = None
    try:
        tmp_file = tempfile.NamedTemporaryFile(mode=mode, suffix=suffix, delete=False)
        tmp_name = tmp_file.name
        tmp_file.close()
        yield tmp_name
    finally:
        Path(tmp_name).unlink(missing_ok=True)