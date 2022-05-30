import os
from math import log
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

def remove_outliers(data, m=5):
    """
    remove extreme value (outliers) based on median absolute deviation (MAD) measurement around the median
    :param data: numpy array
    :param m: tolerance, the larger `m` remove less outliers
    :returns: numpy array with outlier removed
    """
    distance = np.abs(data - np.nanmedian(data))
    mdev = np.nanmedian(distance)
    s = distance / (mdev if mdev else 1.)
    np.putmask(data, s>m, np.nan)
    return data


def esd_powerlaw(N=8, esd_min=0.6, esd_max=1900):
    """
    ecosystem size structure builder
    solve size=k(N^a)
    """
    k = esd_min
    a = log(esd_max/esd_min, N)
    size_count = np.linspace(1, N, N)
    size = k * (size_count ** a)
    return size

def esd_linear(N=8, esd_min=0.6, esd_max=1900):
    """
    ecosystem size structure builder
    solve size=kN+b
    """
    k = (esd_max - esd_min)/(N-1)
    b = esd_min - k
    size_count = np.linspace(1, N, N)
    size = size_count * k + b
    return size

def esd_original(N=2, esd_min=0.6, k =1.3, esd_max=1900):
    """
    Ward et al. (2018) method
    :param N: the least size classes in each loop
    :param k: size increment
    :param esd_max: maximum, in micron
    """

    # build a smallest size group
    l = []
    for i in range(N):
        l.append(esd_min + i*k)
    esd = np.array(l)

    # growth by a exponent of 10
    while esd.max() < esd_max:
        esd = np.append(esd, esd[-N:]*10)
    return esd


def esd_ln(N=8, esd_min=0.6, esd_max=1900):
    """
    ecosystem size structure builder
    solve size= a + b ln(x)
    """
    a = esd_min
    b = (esd_max - a)/np.log(N)
    size_count = np.linspace(1, N, N)
    size = a + b * np.log(size_count)
    return size

def esd_logistic(N=8, esd_max=1900, k = 1.5):
    """
    ecosystem size structure builder following logistic model
    """
    l = esd_max
    x0 = int(N/2)
    size_count = np.linspace(1, N, N)
    size = l/(1+np.e**(-k*(size_count - x0)))
    return size
