import os
from pathlib import Path
import numpy as np


def check_rm(path):
    """
    Check and delete file

    :param path: path string of target file
    """
    if os.path.isfile(path):
        print(f"removing {path}")
        os.remove(path)


def file_exists(path):
    if os.path.isfile(path):
        return True
    else:
        raise FileNotFoundError(f"{path} not exist")


def is_empty(path):
    """
    Check if file is empty. If the file does not exist, then create one.

    :param path: path string of target file
    :returns: boolean operator
    """
    if os.path.isfile(path):
        return os.stat(path).st_size == 0
    else:
        Path(path).touch()
