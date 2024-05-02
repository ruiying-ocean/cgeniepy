import numpy as np
import xarray as xr
from scipy.stats import sem
import regionmask

from .grid import Interporaltor, GridOperation
from .plot import GriddedDataVis
from .chem import Chemistry
import cgeniepy.table as ct
from functools import cache

import warnings


class GriddedData:
    """
    GriddedData is a class to store and compute GENIE netcdf data.

    It stores data in xarray.DataArray format, and provides optimalised methods for GENIE model output to compute statistics.    
    """
    
    def __init__(self, array=np.nan, mutable=False, attrs={}):
        """
        Initialise an instance of GriddedData

        
        :param array: a xarray.DataArray object
        :param mutable: whether the data is mutable in place, default is False
        :param attrs: a dictionary of attributes, default is empty

        -----------
        Example
        -----------
        data = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20]}, attrs={"units": "unitless""})
        dg = GriddedData(data, mutable=False)
        """
        
        # set data
        self.data = array
        self.mutable = mutable
        self.attrs = attrs

        ## formatting the unit
        if 'units' in self.attrs:
            self.attrs['units'] = Chemistry().format_unit(self.attrs['units'])

    def __getitem__(self, item):
        "make GriddedData subscriptable like xarray.DataArray"
        return self.data[item]

    def interpolate(self, *args, **kwargs):
        """Interpolate the GriddedData to a finer grid, mostly useful for GENIE plotting

        :returns: a GriddedData object with interpolated data
        
        ------------
        Example
        ------------
        >>> model = GenieModel(path)
        >>> var = model.get_var(target_var)
        >>> var.interpolate(method='r-linear') ## regular linear interpolation        
        """
        
        coords = tuple([self.data[dim].values for dim in self.data.dims])
        values = self.data.values
        dims = self.data.dims
        output = Interporaltor(dims, coords,values, *args, **kwargs).to_xarray()
        
        if self.mutable:
            self.data = output
            self.attrs = self.attrs            
            return self            
        else:
            return GriddedData(output, mutable=False, attrs=self.attrs)
        

    def sel(self, *args, **kwargs):
        """select grid points based on the given coordinates
        a wrapper to xarray `sel` method

        :returns: a GriddedData object with selected data

        -----------
        Examples:
        -----------
        >>> model = GenieModel(path)
        >>> var = model.get_var(target_var)
        >>> var.search_grid(lon=XX, lat=XX, zt=XX, method='nearest')
        """

        if self.mutable:
            self.data = self.data.sel(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.sel(*args, **kwargs))


    def isel(self, *args, **kwargs):
        """select grid points based on the given index
        a wrapper to xarray `isel` method

        -----------
        Examples:
        -----------
        >>> model = GenieModel(path)
        >>> var = model.get_var(target_var)
        >>> var.isel(lon=XX, lat=XX, zt=XX)
        """
        if self.mutable:
            self.data = self.data.isel(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.isel(*args, **kwargs), mutable=False, attrs=self.attrs)
    
    def normalise_longitude(self, method='g2n', *args, **kwargs):
        """Normalise GENIE's longitude (eastern degree) to normal longitude (-180, 180)

        :param method: normalise method, default is 'g2n' (GENIE to normal). A full list of methods are:
            - 'g2n': GENIE to normal
            - 'n2g': normal to GENIE
            - 'e2n': eastern to normal
            - 'n2e': normal to eastern

        -----------
        Example
        -----------

        >>> Model = GenieModel("a path")
        >>> Model.get_var('abc').normalise_longitude()
        """
        gp = GridOperation()

        match method:
            case 'g2n':
                output = gp.xr_g2n(self.data, *args, **kwargs)
            case 'n2g':
                output = gp.xr_n2g(self.data, *args, **kwargs)
            case 'e2n':
                output = gp.xr_e2n(self.data, *args, **kwargs)
            case 'n2e':
                output = gp.xr_n2e(self.data, *args, **kwargs)
            case _:
                ## No change
                warnings.warn("No change applied because of invalid method")
                output = self.data

        if self.mutable:
            self.data = output
            return self        
        else:
             return GriddedData(output, mutable=False, attrs=self.attrs)
                    

    def __add__(self, other):
        """Allow GriddedData to be added by a number or another GriddedData
        """
        sum_array = GriddedData()
        if hasattr(other, "array"):
            sum_array.data = self.data + other.array
            ## check the attributes
            np.testing.assert_equal(self.attrs, other.attrs)            
            sum_array.attrs = self.attrs
        else:
            ## a scalar or xarray.DataArray
            sum_array.data = self.data + other
            sum_array.attrs = self.attrs
        return sum_array

    def __sub__(self, other):
        """
        Allow GriddedData to be subtracted by a number or another GriddedData
        """
        diff = GriddedData()
        if hasattr(other, "data"):
            
            diff.data = self.data - other.data
            ## check all the attributes
            np.testing.assert_equal(self.attrs, other.attrs)
            diff.attrs = self.attrs
        else:
            ## a scalar
            diff.data = self.data - other
        return diff

    def __rsub__(self, other):
        """
        Allow a number or another GriddedData to be subtracted from this GriddedData
        """
        diff = GriddedData()
        if hasattr(other, "data"):
            diff.data = other.data - self.data
            # check all the attributes
            np.testing.assert_equal(self.attrs, other.attrs)
            diff.attrs = self.attrs
        else:
            # a scalar
            diff.data = other - self.data

        return diff

    def __truediv__(self, other):
        """
        Allow GriddedData to be divided by a number or another GriddedData

        NA/NA -> NA
        """
        quotient_array = np.zeros_like(self.data)
        
        if isinstance(other, GriddedData):
            quotient_array = np.divide(self.data, other.data)
        else:
            # Handle division by a scalar
            quotient_array = np.divide(self.data, other)

        quotient = GriddedData()
        quotient.data = xr.DataArray(quotient_array)
        return quotient

    def __mul__(self, other):
        """
        Allow GriddedData to be multiplied by a number or another GriddedData
        """
        product = GriddedData()
        if hasattr(other, "data"):
            product.data = self.data * other.data
        else:
            try:
                product.data = self.data * other
            except ValueError:
                print("Sorry, only number and GriddedData are accepted")

        return product

    def __pow__(self, other):
        """
        Allow GriddedData to be raised to a power
        """

        if self.mutable:
            self.data = self.data ** other
            return self
        else:
            return GriddedData(self.data ** other, mutable=False, attrs=self.attrs)


    def max(self, *args, **kwargs):
        """compute the maximum value of the array
        """
        
        self.data = self.data.max(*args, **kwargs)
        return self


    def min(self, *args, **kwargs):
        """compute the minimal value of the array
        """        
        if self.mutable:
            self.data = self.data.min(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.min(*args, **kwargs), mutable=False, attrs=self.attrs)


    def sum(self, *args, **kwargs):
        """ compute the sum of the array
        """
        
        if self.mutable:
            self.data = self.data.sum(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.sum(*args, **kwargs), mutable=False, attrs=self.attrs)


    def mean(self, *args, **kwargs):
        """compute the mean of the array.
        Note this is not weighted mean, for weighted mean, use `weighted_mean` method
        """
        if self.mutable:
            self.data = self.data.mean(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.mean(*args, **kwargs), mutable=False, attrs=self.attrs)
    

    def median(self, *args, **kwargs):
        """compute the median of the array
        """
        
        if self.mutable:
            self.data = self.data.median(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.median(*args, **kwargs), mutable=False, attrs=self.attrs)

    def sd(self, *args, **kwargs):        
        "compute the standard deviation of the mean"
        if self.mutable:
            self.data = np.std(self.data, *args, **kwargs)
            return self
        else:
            return np.std(self.data, *args, **kwargs)


    def variance(self, *args, **kwargs):
        "compute the variance of the array"
        if self.mutable:
            self.data = np.var(self.data, *args, **kwargs)
            return self
        else:
            return np.var(self.data, *args, **kwargs)

    def se(self, *args, **kwargs):
        "compute the standard error of the mean"
        if self.mutable:
            self.data = sem(self.data, nan_policy="omit", axis=None, *args, **kwargs)
            return self
        else:
            return sem(self.data, nan_policy="omit", axis=None, *args, **kwargs)

    def weighted_mean(self, weights, *args, **kwargs):
        """
        compute the weighted average of the array (e.g., ocean volume)
        """
        array_ma = self.data.to_masked_array()
        return np.ma.average(array_ma, weights=weights, *args, **kwargs)        
    
    def sel_modern_basin(self, basin, norm_lon_method='g2n'):
        """
        select modern basin from regionmask.defined_regions.ar6.ocean

        The 58 defines marine regions from the sixth IPCC assessment report (AR6), Iturbide et al., (2020) ESSD
        
        https://regionmask.readthedocs.io/en/stable/_images/plotting_ar6_all.png

        :param basin: the basin index (int) or basin abbrev name (str), or a list of basin index or basin name
        :param norm_lon_method: normalise longitude method, default is 'g2n'
        
        --------------------------------
        47: North Pacific
        48: Equatorial Pacific
        49: South Pacific
        
        50: North Atlantic Ocean
        51: Equatorial Atlantic Ocean
        52: Southern Atlantic Ocean

        53: North Indian Ocean
        55: Equatorial Indian Ocean
        56: South Indian Ocean

        57: Southern Ocean
        46: Arctic Ocean
        --------------------------------

        
        Example:
        ---------        
        >>> gd.sel_modern_basin(47) ## North Pacific
        """
        ocean = regionmask.defined_regions.ar6.ocean

        target_ocean_index = ocean.map_keys(basin)
        mask = ocean.mask(self.normalise_longitude(method=norm_lon_method)) ## an masked array target regions are the index, others are NA

        if isinstance(target_ocean_index, list):
            cond_lst = []
            for i in target_ocean_index:
                cond_lst.append(mask == i)
            cond = np.logical_or.reduce(cond_lst)
        else:
            cond = mask == target_ocean_index


        if self.mutable:
            self.data = self.normalise_longitude(method=norm_lon_method).data.where(cond)                    
            return self
        else:
            return GriddedData(self.normalise_longitude(method=norm_lon_method).data.where(cond), mutable=False, attrs=self.attrs)


    def mask_basin(self, base, basin, subbasin):

        """use pre-defined grid mask to select_basin, mostly used for cGENIE model

        :param base: the base configuration name, e.g., worjh2, worlg4
        :param basin: the basin name, e.g., Atlantic, Pacific, Indian
        :param subbasin: N/S/ALL, ALL means Southern Ocean section included

        Example
        ---------
        >>> gd.mask_basin('worjh2', 'Atlantic'').mean()        
        """
        gp = GridOperation()
        data = self.data
        attr = self.attrs
        mask = gp.GENIE_grid_mask(base=base, basin=basin, subbasin=subbasin, invert=True)

        if self.data.ndim > 2:
            mask = np.broadcast_to(mask, (16, 36, 36))

        mask_data = np.ma.array(data, mask=mask)
        mask_data = np.ma.masked_invalid(mask_data)

        output = xr.DataArray(mask_data, dims=data.dims, coords=data.coords)


        if self.mutable:
            self.data = output
            self.attrs = attr        
            return self
        else:
            return GriddedData(output, mutable=False, attrs=attr)


    @cache
    def ocn_only_data(self, index=False):
        """
        remove the NA grid (i.e., land definition)
        :param index: return GeoIndex if True or real value if false
        """

        stacked_data = self.data.stack(x=self.data.dims)            
        stacked_ocn_mask =~np.isnan(stacked_data)
        ocn_only_data = stacked_data[stacked_ocn_mask]
        
        if index:
            return np.stack(ocn_only_data.x.values)
        else:
            return ocn_only_data
        

    def search_point(self, point, ignore_na=False, to_genielon=False, **kwargs):
        """
        search the nearest grid point to the given coordinates
        
        :param point: a list/tuple of coordinate values in the same order to the data dimension (self.data.dims)
        :param ignore_na: whether only check ocean data (which ignore the NA grids)
        :param to_genielon: whether convert the input longitude to genie longitude, input point must be list if True

        ----------------------
        Example
        ----------------------
        >>> lat, lon, depth = 10, 20, 30
        >>> data.search_point((lat, lon, depth))
        """
        ## ignore the first dimension (which is time)        
        ndim = self.data.ndim
        
        if len(point) != ndim:
            raise ValueError("Input point has incompatiable coordinate")
        
        if ignore_na:
            index_pool = self.ocn_only_data(index=True)
            match ndim:
                case 3:
                    ## point: (depth, lat, lon); Index: (depth, lat, lon)
                    if to_genielon: point[2] = GridOperation().lon_g2n(point[2])
                    distances = GridOperation().geo_dis3d(point, index_pool[:, 0:4])
                    idx_min = np.argmin(distances)
                case 2:
                    ## point: (lat, lon); Index: (lat, lon)
                    if to_genielon: point[1] = GridOperation().lon_g2n(point[1])
                    distances = GridOperation().geo_dis2d(point, index_pool[:, 0:3])
                    idx_min = np.argmin(distances)
            
            nearest_value = self.ocn_only_data(index=False).values[idx_min]
            return nearest_value
        else:
            ## construct a dictionary based on self.data.dims and input point, and method ("nearest")
            kwargs = dict(zip(self.data.dims, point), method='nearest')
            return self.data.sel(**kwargs).values.item()

    def to_GriddedDataVis(self):
        """
        convert GriddedData to GriddedDataVis
        """
        return GriddedDataVis(self)

    def plot(self, *args, **kwargs):
        """
        plot the data using GriddedDataVis
        """
        return self.to_GriddedDataVis().plot(*args, **kwargs)

    def to_dataframe(self):
        """
        convert GriddedData to pandas.DataFrame
        """
        return self.data.to_dataframe()

    def to_ScatterData(self):
        """
        convert GriddedData to ScatterData
        """
        return ct.ScatterData(self.to_dataframe())
