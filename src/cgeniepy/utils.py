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
   """
   Check if file exists. If not, raise FileNotFoundError.

   :param path: path string of target file
   :returns: True if file exists
   """
   file_path = Path(path).expanduser()
   if file_path.is_file():
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

def efficient_log(data, replace_zero=10):
    """A shortcut to efficiently calculate log10 of data.

    :param data: data to calculate log10
    :param replace_zero: value to replace zero in data

    :returns: log10 of data
    """
    return np.where(data == 0, replace_zero, np.log10(data))
