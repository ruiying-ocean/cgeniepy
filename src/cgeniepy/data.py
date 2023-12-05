import pathlib

from xarray import open_dataset
import numpy as np
from .grid import normalise_obs_lon


def efficient_log(data):
    "keep NA, remove zeros"
    return np.where(data == 0, -10, np.log10(data))


def foram_groups():
    """
    get a dictionary with foram abbrev (keys), pft_index and complete name (values).

    :returns: dictionary
    """
    foram_names = {
        "bn": [16, "symbiont-barren non-spinose"],
        "bs": [17, "symbiont-barren spinose"],
        "sn": [18, "symbiont-facultative non-spinose"],
        "ss": [19, "symbiont-obligate spinose"],
    }

    return foram_names
