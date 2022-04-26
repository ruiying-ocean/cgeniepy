import numpy as np
from scipy.spatial import distance
from netCDF4 import Dataset
from .data import get_obs_data

def safe_unveil(data):
    "get pure array from a numpy masked array object"
    if data.__class__ != np.ma.core.MaskedArray:
        return data
    else:
        return data.filled(np.nan)

def intersect_index(array1, array2, verbose=False):
    """
    Return the index where corresponding values are not
    nan in both input arrays. One then can filter the array
    by the output boolean array.
    """

    # If both are not NaN, return it
    array1 = safe_unveil(array1)
    array2 = safe_unveil(array2)

    indx_array = np.logical_and(~np.isnan(array1), ~np.isnan(array2))

    if verbose is True:
        num = indx_array.flatten()[indx_array.flatten() == True].shape[0]
        print("Summary: {} elements simultaneously exist.".format(num))

    return indx_array


def cal_mscore(data1: np.array, data2: np.array):
    """
    Calculate skill metric M-score. See more in the paper Watterson, I. G. (1996)

    Use 2D array as input, order causes no difference.
    """
    indx = intersect_index(data1, data2)
    sub_data1 = data1[indx]
    sub_data2 = data2[indx]

    mse = np.square(np.subtract(sub_data1, sub_data2)).mean()
    v1 = sub_data1.var()
    v2 = sub_data2.var()
    g1 = sub_data1.mean()
    g2 = sub_data2.mean()

    mscore = (2 / np.pi) * np.arcsin(1 - mse / (v1 + v2 + np.square(g1 - g2)))

    return mscore


def cal_cosine_similarity(data1, data2):
    """
    Calculate metric cosine similarity of two input arrays.

    Use 2D array as input, order causes no difference.
    """
    indx = intersect_index(data1, data2)
    sub_data1 = data1[indx]
    sub_data2 = data2[indx]

    if sub_data1.mean() != 0 and sub_data2.mean() != 0:
        cos_sim = 1 - distance.cosine(sub_data1.flatten(), sub_data2.flatten())
    else:
        cos_sim = np.nan

    return cos_sim


def cal_rmse(data1, data2):
    """
    Calculate Root Mean Sqaure Error (rmse) between two input arrays.

    Use 2D array as input, order causes no difference.
    """

    error_2d = data1 - data2
    error_1d = error_2d.ravel()[~np.isnan(error_2d.ravel())]
    rmse = np.sqrt(np.square(error_1d).mean())

    return rmse


def get_foram_prop(file_path, var):
    """
    Quick calculation of [modelled] relative abundance, based on carbon export flux
    because of little difference between biomass and export.

    :param file_path: an netcdf file with all foram-related varialbes
    :param var: foram group abbrev: bn, bs, sn, ss

    :returns: a scalar value
    """

    f = Dataset(file_path)

    bn = safe_unveil(f.variables['eco2D_Export_C_016'][-1,:,:])
    bs = safe_unveil(f.variables['eco2D_Export_C_017'][-1,:,:])
    sn = safe_unveil(f.variables['eco2D_Export_C_018'][-1,:,:])
    ss = safe_unveil(f.variables['eco2D_Export_C_019'][-1,:,:])

    total_foram = bn + bs + sn + ss

    #ignore divided by 0
    with np.errstate(divide='ignore', invalid='ignore'):                
        one_foram = locals()[var]
        proportion =np.divide(one_foram, total_foram, out=np.zeros_like(one_foram), where=total_foram!=0)

    f.close()

    return proportion

def quick_rmse(data, obs_source, var):
    "A wrapper function to calculate RMSE"
    return cal_rmse(data, get_obs_data(obs_source, var))


def quick_mscore(data, obs_source, var):
    "A wrapper function to calculate M-Score"
    return cal_mscore(data, get_obs_data(obs_source, var))


def quick_cos_sim(data, obs_source, var):
    "A wrapper function to calculate cosine cimilarity"
    return cal_cosine_similarity(data, get_obs_data(obs_source, var))


#consider convert data into getModelData(var)
