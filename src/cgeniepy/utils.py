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


def set_sns_barwidth(ax, width):
    width_list = [width] * len(ax.patches)

    for bar, newwidth in zip(ax.patches, width_list):
        # get current position
        x = bar.get_x()
        width = bar.get_width()
        centre = x + width / 2.0

        # update width
        bar.set_x(centre - newwidth / 2.0)
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
    s = distance / (mdev if mdev else 1.0)
    np.putmask(data, s > m, np.nan)
    return data


def powlaw(x, a, b):
    return a * np.power(x, b)


def esd_pow(n=8, esd_max=1900):
    """
    ecosystem size structure builder following power law, size=k(n^a).
    build from scipy.curve_fit based on exisiting 8P8Z and 32P32Z size structure.
    Principal is more size classes (`n`) should allocate more at smaller organisms
    instead of interpolation or equal increment.

    :param n: number of size classes
    :param esd_max: maximum individual size (um)
    """

    a = 8.5 + 0.25 * (n / 8)
    k = esd_max / n**a
    x = powlaw(np.linspace(1, n, n), k, a)
    return x


def esd_ward(N=2, esd_min=0.6, k=1.3, esd_max=1900):
    """
    Ward et al. (2018) method
    :param N: the least size classes in each loop
    :param k: size increment
    :param esd_max: maximum, in micron
    """

    # build a smallest size group
    l = []
    for i in range(N):
        l.append(esd_min + i * k)
        esd = np.array(l)

    # growth by a exponent of 10
    while esd.max() < esd_max:
        esd = np.append(esd, esd[-N:] * 10)
    return esd


def distance(point1: tuple, point2: tuple):
    """
    point: a tuple (with >= 2 numbers) of coordinate values.
    Coordinate value can be float or numpy array.
    """
    if len(point1) != len(point2):
        raise ValueError("Inconsistent dimension between two points!")

    for a in point1:
        for b in point2:
            if np.shape(a) != np.shape(b):
                raise ValueError("Inconsistent coordinate data shape")

    x = np.zeros(np.shape(point1[0]))

    for i in range(len(point1)):
        x = x + np.square(point1[i] - point2[i])

    return np.sqrt(x)

def my_sum(*args):
    total = 0
    for x in args:
        total += x
    return total    
