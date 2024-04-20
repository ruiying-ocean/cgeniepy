import numpy as np
import xarray as xr
from scipy.stats import sem
import regionmask

from .grid import Interporaltor, GridOperation
from .plot import GriddedDataVis
from .chem import Chemistry
from functools import cache


class GriddedData():
    """
    GriddedData is a class to store and compute GENIE netcdf data.

    It stores data in xarray.DataArray format, and provides optimalised methods for GENIE model output to compute statistics.    
    """
    
    def __init__(self, array=np.nan, mutable=False, attrs={}):
        """
        Initialise an instance of GriddedData

        Parameters
        ----------
        array : xarray.DataArray
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
        """
        Interpolate the GriddedData to a finer grid, mostly useful for GENIE plotting
        """
        
        coords = tuple([self.data[dim].values for dim in self.data.dims])
        values = self.data.values
        dims = self.data.dims
        
        if self.mutable:
            self.data = Interporaltor(dims, coords,values, *args, **kwargs).to_xarray()
            return self            
        else:
            return GriddedData(Interporaltor(dims, coords,values, *args, **kwargs).to_xarray())
        

    def sel(self, *args, **kwargs):
        """select grid points based on the given coordinates
        a wrapper to xarray `sel` method

        -----------
        Examples:
        -----------
        >>> model = GenieModel(path)
        >>> var = model.get_var(target_var)
        >>> var.search_grid(lon=XX, lat=XX, zt=XX, method='nearest')
        
        ## for loop over a data frame (df) with longitude, latitude, depth columns
        >>> modelvar = []
        >>> for i in range(len(df)):
        >>>     lat = df.latitude.iloc[i][0]
        >>>     lon = df.longitude[i].iloc[0]
        >>>     zt = df.depth[i].iloc[0]
        >>>     kwargs = {'lat': lat, 'lon':lon, 'zt': zt, 'method': 'nearest'}
        >>>     var = var.sel(**kwargs).values
        >>>     modelvar.append(var)        
        """

        if self.mutable:
            self.data = self.data.sel(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.sel(*args, **kwargs))


    def isel(self, *args, **kwargs):
        """select grid points based on the given index
        a wrapper to xarray `isel` method
        """
        if self.mutable:
            self.data = self.data.isel(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.isel(*args, **kwargs), mutable=False, attrs=self.attrs)
    
    def normalise_longitude(self):
        """Normalise GENIE's longitude (eastern degree) to normal longitude (-180, 180)

        -----------
        Example
        -----------

        >>> Model = GenieModel("a path")
        >>> Model.get_var('abc').normalise_longitude()
        """
        gp = GridOperation()
        if self.mutable:
            self.data = gp.normalise_GENIE_lon(self.data)
            return self
        else:
            return GriddedData(gp.normalise_GENIE_lon(self.data), mutable=False, attrs=self.attrs)

    def __add__(self, other):
        """Allow GriddedData to be added by a number or another GriddedData
        """
        sum = GriddedData()
        if hasattr(other, "array"):
            sum.data = self.data + other.array
        else:
            ## a scalar
            sum.data = self.data + other
        return sum

    def __sub__(self, other):
        """
        Allow GriddedData to be subtracted by a number or another GriddedData
        """
        diff = GriddedData()
        if hasattr(other, "data"):
            diff.data = self.data - other.data
        else:
            ## a scalar
            diff.data = self.data - other
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


    def max(self, *args, **kwargs):
        self.data = self.data.max(*args, **kwargs)
        return self


    def min(self, *args, **kwargs):
        if self.mutable:
            self.data = self.data.min(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.min(*args, **kwargs), mutable=False, attrs=self.attrs)


    def sum(self, *args, **kwargs):
        if self.mutable:
            self.data = self.data.sum(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.sum(*args, **kwargs), mutable=False, attrs=self.attrs)


    def mean(self, *args, **kwargs):
        if self.mutable:
            self.data = self.data.mean(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.mean(*args, **kwargs), mutable=False, attrs=self.attrs)
    

    def median(self, *args, **kwargs):
        if self.mutable:
            self.data = self.data.median(*args, **kwargs)
            return self
        else:
            return GriddedData(self.data.median(*args, **kwargs), mutable=False, attrs=self.attrs)

    def sd(self, *args, **kwargs):
        "standard deviation of the mean"
        if self.mutable:
            self.data = np.std(self.data, *args, **kwargs)
            return self
        else:
            return np.std(self.data, *args, **kwargs)


    def variance(self, *args, **kwargs):
        if self.mutable:
            self.data = np.var(self.data, *args, **kwargs)
            return self
        else:
            return np.var(self.data, *args, **kwargs)

    def se(self, *args, **kwargs):
        "standard error of the mean"
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
    
    def sel_modern_basin(self, basin):
        """
        select modern basin from regionmask.defined_regions.ar6.ocean

        The 58 defines marine regions from the sixth IPCC assessment report (AR6), Iturbide et al., (2020) ESSD
        
        https://regionmask.readthedocs.io/en/stable/_images/plotting_ar6_all.png
        
        --------------------------------
        47: North Pacific
        48: Equatorial Pacific
        49: South Pacific
        
        50: North Atlantic Ocean
        51: Equatorial Atlantic Ocean
        52: Southern Atlantic Ocean

        57: Southern Ocean
        --------------------------------        
        """
        ocean = regionmask.defined_regions.ar6.ocean
        index = ocean.map_keys(basin)
        mask = ocean.mask(self.normalise_longitude())        

        if self.mutable:
            self.data = self.normalise_longitude().data.where(mask == index)                    
            return self
        else:
            return GriddedData(self.normalise_longitude().data.where(mask == index), mutable=False, attrs=self.attrs)


    def mask_basin(self, base, basin, subbasin):

        """
        use pre-defined grid mask to select_basin, mostly used for cGENIE model
        """
        gp = GridOperation()
        data = self.data
        mask = gp.GENIE_grid_mask(base=base, basin=basin, subbasin=subbasin, invert=True)

        if self.data.ndim > 2:
            mask = np.broadcast_to(mask, (16, 36, 36))

        mask_data = np.ma.array(data, mask=mask)
        mask_data = np.ma.masked_invalid(mask_data)


        if self.mutable:
            self.data = xr.DataArray(mask_data, dims=data.dims, coords=data.coords)            
            return self
        else:
            return xr.DataArray(mask_data, dims=data.dims, coords=data.coords)


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


    def plot(self, *args, **kwargs):
        """
        plot the data using GriddedDataVis
        """
        return GriddedDataVis(self.data, self.attrs).plot(*args, **kwargs)
