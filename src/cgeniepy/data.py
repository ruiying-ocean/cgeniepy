from typing import Union
import pathlib

from xarray import open_dataset
import numpy as np
import pandas as pd
from scipy.stats import sem

from .grid import reassign_obs
from .utils import remove_outliers


def efficient_log(data: Union[int, float]) -> float:
    "keep NA, remove zeros"
    return np.where(data == 0, -10, np.log10(data))


def obs_data(source, var, stat=None, outlier_level=None):
    """
    Quickly fetch observational data from ForamData/data directory.
    The longitude coordinate will be assigned to fit GENIE output

    :returns: numpy 2D array
    :raise ValueError:
    """
    if source == "core":
        file_path = pathlib.Path(__file__).parent.parent / "data/ForCenS_regridded.nc"
    elif source == "tow":
        file_path = (
            pathlib.Path(__file__).parent.parent / "data/plankton_tow_regridded_abs.nc"
        )
    elif source == "trap":
        file_path = (
            pathlib.Path(__file__).parent.parent
            / "data/sediment_trap_regridded_organic.nc"
        )
    elif source == "MARGO":
        file_path = pathlib.Path(__file__).parent.parent / "data/MARGO_regridded.nc"
    else:
        raise ValueError("type must be core, obs or trap!")

    ds = open_dataset(file_path)
    long_name = foram_names()[var]
    obs = ds[long_name]
    modified_obs = reassign_obs(obs).to_numpy()

    if not outlier_level:
        pass
    else:
        modified_obs = remove_outliers(modified_obs, m=outlier_level)

    if not stat:
        return modified_obs
    else:
        mean = np.nanmean(modified_obs)
        sd = np.nanstd(modified_obs)
        se = sem(modified_obs, nan_policy="omit", axis=None)
        return mean, sd, se


def obs_stat_bysource(source, *args, **kwargs) -> pd.DataFrame:
    obs = []
    foram_abbrev = foram_names().keys()
    foram_fullnames = foram_names().values()

    for i in foram_abbrev:
        tmp = obs_data(source=source, var=i, stat="Yes", *args, **kwargs)
        obs.append(tmp)

    table = pd.DataFrame(obs, index=foram_fullnames, columns=["mean", "sd", "se"])
    return table


def obs_stat_bytype(type, *args, **kwargs) -> pd.DataFrame:
    tow = obs_stat_bysource("tow", *args, **kwargs).loc[:, type]
    trap = obs_stat_bysource("trap", *args, **kwargs).loc[:, type]
    core = obs_stat_bysource("core", *args, **kwargs).loc[:, type]
    # combination
    data = pd.concat([tow, trap, core], axis=1)
    data.columns = [
        "Biomass(mmol C/m3)",
        "Carbon Export (mmol C/m3/d)",
        "Relative Abundance",
    ]

    return data

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

def filter_foramdf(dataframe, foram_group=None):
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
