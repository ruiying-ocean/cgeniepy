import pathlib

from xarray import open_dataset
import numpy as np
from .grid import normalise_obs_lon


def efficient_log(data):
    "keep NA, remove zeros"
    return np.where(data == 0, -10, np.log10(data))


def obs_data(source, var):
    """
    An interface to fetch prescribed observational data.
    The longitude coordinate will be assigned to fit GENIE output

    :returns: numpy 2D array
    :raise ValueError:
    """
    if source == "coretop":
        file_path = pathlib.Path(__file__).parent.parent / "data/ForCenS_regridded.nc"
    elif source == "net":
        file_path = (
            pathlib.Path(__file__).parent.parent / "data/plankton_tow_regridded_abs.nc"
        )
    elif source == "trap":
        file_path = (
            pathlib.Path(__file__).parent.parent
            / "data/sediment_trap_regridded_organic.nc"
        )
    elif source == "margo":
        file_path = pathlib.Path(__file__).parent.parent / "data/MARGO_regridded.nc"
    else:
        raise ValueError("type must be coretop, net or trap!")

    ds = open_dataset(file_path)
    long_name = foram_names()[var][1]
    obs = ds[long_name]
    modified_obs = normalise_obs_lon(obs).to_numpy()

    return modified_obs


def foram_names():
    """
    get a dictionary with foram abbrev (keys) and complete name (values).

    :returns: dictionary
    """
    foram_names = {
        "bn": [16, "symbiont-barren non-spinose"],
        "bs": [17, "symbiont-barren spinose"],
        "sn": [18, "symbiont-facultative non-spinose"],
        "ss": [19, "symbiont-obligate spinose"],
    }

    return foram_names

def filter_foramdf(dataframe, foram_group):
    """
    select data based on foraminifer group, either pass foraminifer abbrev
    or explicitly pass boolean values to `symbiosi` & `spinose`

    :param dataframe: pandas dataframe
    :param foram_group: abbreviation name of foraminifer

    :returns: a subsetted dataframe
    """
    foram_bool = {
        "bn": ["No", "No"],
        "bs": ["No", "Yes"],
        "sn": ["Yes", "No"],
        "ss": ["Yes", "Yes"],
    }

    trait = foram_bool[foram_group]

    query_string = "Symbiosis == '{}' & Spinose == '{}'".format(trait[0], trait[1])
    return dataframe.query(query_string)
