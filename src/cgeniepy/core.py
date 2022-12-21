from os.path import join

import xarray as xr
import numpy as np
from scipy.stats import sem
import regionmask
from netCDF4 import Dataset
from .plot import GeniePlottable

from . import Q_
from .grid import (
    GENIE_grid_area,
    reassign_GENIE,
    GENIE_grid_mask,
    GENIE_grid_vol,
    mask_Arctic_Med
)
from .utils import file_exists
from .chem import rm_element
from .scores import ModelSkill
from .utils import remove_outliers

class GenieArray(GeniePlottable):

    def __init__(self, M=36, N=36):
        """
        Create an empty 2D array, default as 36x36
        """
        # set dimension
        self.M = M
        self.N = N
        # set data
        self.array = self._set_array()
        # set unit
        if hasattr(self.array, "units"):
            self.unit = self.array.units
        else:
            self.unit = ""

    def _set_array(self):
        "assign real data"
        return np.zeros((self.M, self.N))

    def __getitem__(self, item):
        "make iterable"
        return self.array()[item]

    def pure_array(self):
        "get a numpy array"
        if hasattr(self.array, "values"):
            return self.array.values
        else:
            return self.array

    def uarray(self):
        "array with unit"
        unit = rm_element(self.unit)
        uarray = Q_(self.pure_array(), unit)
        return uarray

    def sel(self, *args, **kwargs):
        "a wrapper to xarray sel method"
        self.array = self.array.sel(*args, **kwargs)
        return self

    def isel(self, *args, **kwargs):
        "a wrapper to xarray sel method"
        self.array = self.array.isel(*args, **kwargs)
        return self

    def dim(self):
        return self.pure_array().ndim

    def flip(self, *args, **kwargs):
        return np.flip(self.array,  *args, **kwargs)

    def flatten(self):
        "flatten in row-major (C-style)"
        return self.pure_array().flatten(order="C")

    def apply(self, f):
        vfunc = np.vectorize(f)
        x = GenieArray()
        x.array = vfunc(self.pure_array())
        return x

    def reassign_array(self):
        "if self.array is xarray, then reassign the coordinate"
        x = GenieArray()
        x.array = reassign_GENIE(self.array).to_numpy()
        return x

    def _run_method(self, method: str, *args, **kwargs):
        "an alias to run stat for GenieArray class"
        return getattr(self, method)(*args, **kwargs)

    def __add__(self, other):
        sum = GenieArray()
        if hasattr(other, "array"):
            sum.array = self.array + other.array
        else:
            sum.array = self.array + other
        return sum

    def __sub__(self, other):
        diff = GenieArray()
        if hasattr(other, "array"):
            diff.array = self.array - other.array
        else:
            diff.array = self.array - other
        return diff

    def __truediv__(self, other):
        quotient = GenieArray()
        if hasattr(other, "array"):
            quotient.array = np.divide(
                self.array,
                other.array,
                out=np.zeros_like(self.array),
                where=other.array != 0,
            )
        else:
            try:
                quotient.array = np.divide(self.array, other)
            except ValueError:
                print("Sorry, either number and GenieArray are accepted")

        return quotient

    def __mul__(self, other):
        product = GenieArray()
        if hasattr(other, "array"):
            product.array = self.array * other.array
        else:
            try:
                product.array = self.array * other
            except ValueError:
                print("Sorry, only number and GenieArray are accepted")

        return product

    def max(self, *args, **kwargs):
        return np.max(self.array, *args, **kwargs)

    def nanmax(self, *args, **kwargs):
        return np.nanmax(self.array, *args, **kwargs)

    def min(self, *args, **kwargs):
        return np.min(self.array, *args, **kwargs)

    def nanmin(self, *args, **kwargs):
        return np.nanmin(self.array, *args, **kwargs)

    def sum(self, *args, **kwargs):
        return np.sum(self.pure_array(), *args, **kwargs)

    def square(self, *arg, **kwargs):
        return np.square(self.pure_array())

    def sqrt(self, *arg, **kwargs):
        return np.sqrt(self.pure_array())

    def nansum(self, *args, **kwargs):
        return np.nansum(self.pure_array(), *args, **kwargs)

    def mean(self, *args, **kwargs):
        return np.mean(self.pure_array(), *args, **kwargs)

    def nanmean(self, *args, **kwargs):
        return np.nanmean(self.pure_array(), *args, **kwargs)

    def ptp(self, *args, **kwargs):
        "range of values"
        return np.ptp(self.pure_array(), *args, **kwargs)

    def sd(self, *args, **kwargs):
        return np.std(self.pure_array(), *args, **kwargs)

    def nansd(self, *args, **kwargs):
        return np.nanstd(self.pure_array(), *args, **kwargs)

    def cv(self):
        "coefficient of variance, or normalized standard deviation"
        cv = self.nansd() / self.nanmean()
        return cv

    def se(self, *args, **kwargs):
        return sem(self.array, nan_policy="omit", axis=None, *args, **kwargs)

    def select_basin(self, basin):
        ocean = regionmask.defined_regions.ar6.ocean
        index = ocean.map_keys(basin)
        mask = ocean.mask(self.reassign_array())
        regional_data = self.reassign_array().where(mask == index)

        return regional_data

    def mask_basin(self, base, basin, basin_lvl):
        # mask data
        data = self.pure_array()
        mask = GENIE_grid_mask(base=base, basin=basin, basin_lvl=basin_lvl, invert=True)

        if self.dim() > 2:
            mask = np.broadcast_to(mask, (16, 36, 36))

        mask_data = np.ma.array(data, mask=mask)
        mask_data = np.ma.masked_invalid(mask_data)

        garray = GenieArray()
        garray.array = mask_data

        return garray

    def search_grid(self, *args, **kwargs):
        return self.array.sel(*args, **kwargs, method="nearest")

    def search_range(self, lon_min=-255, lon_max=95, lat_min=0, lat_max=90):
        """
        default longitude is unassigned of cGENIE grids
        """

        if lon_min > lon_max or lat_min > lat_max:
            raise ValueError("longitude/latitude min must be less than max!")

        lon = self.array.coords["lon"]
        lat = self.array.coords["lat"]

        return self.array.loc[
            dict(
                lat=lat[(lat >= lat_min) & (lat <= lat_max)],
                lon=lon[(lon >= lon_min) & (lon <= lon_max)],
            )
        ]

    def filter(self, threshold, greater_sign=True):
        data = self.array

        if greater_sign:
            return data.where(data > threshold, drop=True)
        else:
            return data.where(data < threshold, drop=True)

    def compare_obs(self, obs, *args, **kwargs):
        return ModelSkill(model=self.pure_array(), observation=obs, *args, **kwargs)

    def remove_outliers(self, outlier_level):
        self.array = remove_outliers(self.array, m = outlier_level)
        return self

    def mask_Arctic_Med(self, *args, **kwargs):
        self.array = mask_Arctic_Med(self.array, *args, **kwargs)
        return self


class GenieModel(object):
    def __init__(self, model_path):
        self.model_path = model_path

    def _nc_path(self, gem="ecogem", dim="2d"):
        "Essential a fstring extending model_path to NETCDF file path"

        model_path = self.model_path
        nc_file = f"fields_{gem}_{dim}.nc"
        nc_path = join(model_path, gem, nc_file)
        if file_exists(nc_path):
            return nc_path

    def _auto_find_path(self, var):
        "automatically find ecogem/biogem path according to selected variable"

        for gem in ["biogem", "ecogem"]:
            for dim in ["2d", "3d"]:
                nc_path = self._nc_path(gem, dim)
                if self.has_var(var, nc_path):
                    return nc_path
        raise ValueError("Variable not found, please check the spelling!")

    def _open_nc(self, path):
        "Use xarray to open netcdf file"
        return xr.open_dataset(path)

    def _run_method(self, method: str, *args, **kwargs):
        return getattr(self, method)(*args, **kwargs)

    def has_var(self, var, nc_path):
        """
        check if variable exists
        """
        t = Dataset(nc_path, "r")
        if_exist = var in t.variables.keys()
        t.close()
        return if_exist

    def get_vars(self, var_name=None, *args, **kwargs):
        "return specified or all available variables"
        if not var_name:
            t = Dataset(self._nc_path(*args, **kwargs), "r")
            tmp = t.variables.keys()
            t.close()
            return tmp
        else:
            return xr.open_dataset(self._nc_path(*args, **kwargs))[var_name]

    def select_var(self, var: str):
        return GenieVariable(var=var, model_path=self.model_path)

    def _grid_mask(self, source="ecogem", Arctic=True, Med=True):
        """
        cGENIE continent mask array (by setting zero values),
        either calculated from existing data or using biogem.grid_mask
        """
        if source == "ecogem":
            data = self.select_var("eco2D_xGamma_T").array
            grid_mask = xr.where(~np.isnan(data), 1, 0).values
        elif source == "biogem":
            grid_mask = self.select_var("grid_mask").array
        else:
            raise ValueError("source only accept ecogem or biogem")

        if Arctic:
            grid_mask[34:36, :] = 0
        if Med:
            grid_mask[27:30, 25:30] = 0

        return grid_mask

    def _marine_area(self):
        "grid area array in km2"
        grid_mask = self.grid_mask()
        grid_area = GENIE_grid_area()
        mask_area = grid_area * grid_mask

        return mask_area

    def _marine_volume(self):
        "grid volume array in km3"
        grid_mask = self.grid_mask()
        grid_volume = GENIE_grid_vol()
        mask_volume = grid_volume * grid_mask

        return mask_volume

    def check_completeness(self):
        "tbd"
        pass

    def diff(self, model2compare, var):

        if isinstance(model2compare, GenieModel):
            B = model2compare
        else:
            B = GenieModel(model2compare)

        diff = self.select_var(var) - B.select_var(var)
        return diff

    def div(self, model2compare, var):
        B = GenieModel(model2compare)
        diff = self.select_var(var) / B.select_var(var)
        return diff


class GenieVariable(GenieArray):
    def __init__(self, model_path, var):
        self.model_path = model_path
        self.var = var
        GenieArray.__init__(self)

    def _set_array(self):
        gm = GenieModel(model_path = self.model_path)
        if isinstance(self.var, str):
            path2nc =gm._auto_find_path(var=self.var)
            array = gm._open_nc(path2nc)[self.var]
        elif isinstance(self.var, list) or isinstance(self.var, tuple):
            array = []
            for v in self.var:
                path2nc =gm._auto_find_path(v)
                array.append(gm._open_nc(path2nc)[v])

        return array
