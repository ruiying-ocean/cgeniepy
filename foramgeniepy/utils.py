import os
from pathlib import Path
import numpy as np

def check_rm(path):
    """
    Check and delete file

    :param path: path string of target file
    """
    if os.path.isfile(path):
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

def set_sns_barwidth(ax, width):
    width_list = [width] * len(ax.patches)

    for bar,newwidth in zip(ax.patches, width_list):
        # get current position
        x = bar.get_x()
        width = bar.get_width()
        centre = x+width/2.

        # update width
        bar.set_x(centre-newwidth/2.)
        bar.set_width(newwidth)

def mean_w_na(data, na_policy):
    if na_policy == "zero":
        return np.mean(np.nan_to_num(data))
    elif na_policy == "ignore":
        return np.nanmean(data)
