from typing import Union
import pathlib

from xarray import open_dataset
import numpy as np
from scipy.stats import sem

from .grid import reassign_obs


def efficient_log(data: Union[int, float]) -> float:
    "keep NA, remove zeros"
    return np.where(data == 0, -10, np.log10(data))

def obs_data(source, var, summary=None):
    """
    Quickly fetch observational data from ForamData/data directory.
    The longitude coordinate will be assigned to fit GENIE output

    :returns: numpy 2D array
    :raise ValueError:
    """
    if source == "core":
        file_path = pathlib.Path(__file__).parent.parent / "data/ForCenS_regridded.nc"
    elif source == "tow":
        file_path = pathlib.Path(__file__).parent.parent / "data/plankton_tow_regridded_abs.nc"
    elif source == "trap":
        file_path = pathlib.Path(__file__).parent.parent /  "data/sediment_trap_regridded_organic.nc"
    elif source == "MARGO":
        file_path = pathlib.Path(__file__).parent.parent /  "data/MARGO_regridded.nc"
    else:
        raise ValueError("type must be core, obs or trap!")

    ds = open_dataset(file_path)
    long_name = foram_names()[var]
    obs = ds[long_name]
    modified_obs = reassign_obs(obs).to_numpy()

    if not summary:
        return modified_obs
    else:
        mean = np.nanmean(modified_obs)
        sd = np.nanstd(modified_obs)
        se = sem(modified_obs, nan_policy="omit", axis = None)
        return mean, sd, se


# def get_diff_xarray(filename, type, var):
#     sorted_model_data, target_obs = get_comparable_data(filename, type, var)
#     diff_xarray = sorted_model_data - target_obs
#
#     return diff_xarray

def foram_names():
    """
    get a dictionary with foram abbrev (keys) and complete name (values).

    :returns: dictionary
    """
    foram_names = {
        "bn": "symbiont-barren non-spinose",
        "bs": "symbiont-barren spinose",
        "sn": "symbiont-facultative non-spinose",
        "ss": "symbiont-obligate spinose",
    }

    return foram_names


def foram_dict():
    """
    get a dictionary with foram variable name (key) and corresponding data source name(value[0]), abbrev (value[1])

    :returns: dictionary
    """
    foram_dict = {
        "eco2D_Plankton_C_016": ["tow", "bn"],
        "eco2D_Plankton_C_017": ["tow", "bs"],
        "eco2D_Plankton_C_018": ["tow", "sn"],
        "eco2D_Plankton_C_019": ["tow", "ss"],
        "eco2D_Export_C_016": ["trap", "bn"],
        "eco2D_Export_C_017": ["trap", "bs"],
        "eco2D_Export_C_018": ["trap", "sn"],
        "eco2D_Export_C_019": ["trap", "ss"]}

    return foram_dict


def get_calcite_rate(file, time=-1):
    """
    Estimate total foraminiferal inorganic carbon flux rate (g/year)

    :param file: string, netcdf output file
    :returns: 2D array
    """
    data = open_dataset(file).isel(time=time)

    df_bn = data.eco2D_Export_C_016
    df_bs = data.eco2D_Export_C_017
    df_sn = data.eco2D_Export_C_018
    df_ss = data.eco2D_Export_C_019

    POC_total = df_bn + df_bs + df_sn + df_ss  # mmol C/m3/day
    shell_rate = POC_total * 12 * 1e3 / 0.845  # ind/m3/day
    calcite_rate = shell_rate * 0.845 * 3 * 1e-3 * 80.8  # mg/m2/day
    calcite_rate_g_yr = calcite_rate * 365 * 1e-3

    return calcite_rate_g_yr


def filter_foram(dataframe, symbiosis: Union[str, bool], spinose: Union[str, bool]):
    """
    :param dataframe: pandas dataframe
    :param symbiosis: "Yes" | "No", True | False
    :param spinose: "Yes" | "No", True | False

    :returns: a subsetted dataframe
    """
    if type(symbiosis) != type(spinose):
        raise TypeError("Argument symbiosis and spinose should be in the same type!")

    arglist = [symbiosis, spinose]
    for i, arg in enumerate(arglist):
        if arg == True:
            arglist[i] = "Yes"
        elif arg == False:
            arglist[i] = "No"

    query_string = "Symbiosis == '{}' & Spinose == '{}'".format(arglist[0], arglist[1])
    return dataframe.query(query_string)


def POC_to_PIC(POC):

    """
    convert POC (mmol C/m3/day) to PIC (g/m2/yr)
    """

    # unit: ind/m3/day
    shell_rate = POC * 12 * 1e3 / 0.845

    # unit: mg/m2/day
    calcite_rate = shell_rate * 0.845 * 3 * 1e-3 * 80.8

    # unit: g/m2/year
    calcite_rate_g_yr = calcite_rate * 365 * 1e-3

    return calcite_rate_g_yr


def groupped_obsdf(df):
    log_lst = [False, True]
    data_list = []
    for i in log_lst:
        for j in log_lst:
            data_list.append(filter_foram(df, i, j))

    return data_list
