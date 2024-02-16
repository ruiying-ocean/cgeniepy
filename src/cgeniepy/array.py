import numpy as np
import xarray as xr
from scipy.stats import sem
import regionmask

from .grid import Interporaltor, GridOperation
from .plot import ArrayVis

class GriddedData(ArrayVis):
    """
    GriddedData is a class to store and compute GENIE netcdf data.

    It stores data in xarray.DataArray format, and provides optimalised methods for GENIE model output to compute statistics.    
    """
    
    def __init__(self, arr=np.nan):
        """
        Initialise an instance of GriddedData

        Parameters
        ----------
        arr : xarray.DataArray
        """
        
        # set data
        self.array = arr
        
        # init plottable instance
        super().__init__(array=self.array)

    def __getitem__(self, item):
        "make GenieArray subscriptable like xarray.DataArray"
        return self.array[item]

    def interpolate(self, *args, **kwargs):
        """
        Interpolate the GENIE array to a finer grid, might be useful for plotting        
        """
        
        coords = tuple([self.array[dim].values for dim in self.array.dims])
        values = self.array.values
        dims = self.array.dims
        self.array = Interporaltor(dims, coords,values, *args, **kwargs).to_xarray()
        
        return self

    def sel(self, *args, **kwargs):
        """a wrapper to xarray `sel` method

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
        
        self.array = self.array.sel(*args, **kwargs)
        return self


    def isel(self, *args, **kwargs):
        "a wrapper to xarray `isel` method"
        try:
            self.array = self.array.isel(*args, **kwargs)
            return self
        except:
            print("This array instance does not contain xarray.DataArray")
    
    def normalise_longitude(self):
        """
        Normalise GENIE's longitude (eastern degree) to normal longitude (-180, 180)

        -----------
        Example
        -----------

        >>> Model = GenieModel("a path")
        >>> Model.get_var('abc').normalise_longitude()
        """
        gp = GridOperation()
        self.array = gp.normalise_GENIE_lon(self.array)
        return self

    def __add__(self, other):
        """
        Allow GenieArray to be added by a number or another GenieArray
        """
        sum = GriddedData()
        if hasattr(other, "array"):
            sum.array = self.array + other.array
        else:
            ## a scalar
            sum.array = self.array + other
        return sum

    def __sub__(self, other):
        """
        Allow GenieArray to be subtracted by a number or another GenieArray
        """
        diff = GriddedData()
        if hasattr(other, "array"):
            diff.array = self.array - other.array
        else:
            ## a scalar
            diff.array = self.array - other
        return diff

    def __truediv__(self, other):
        """
        Allow GenieArray to be divided by a number or another GenieArray

        NA/NA -> NA
        """
        quotient_array = np.zeros_like(self.array)
        
        if isinstance(other, GriddedData):
            quotient_array = np.divide(self.array, other.array)
        else:
            # Handle division by a scalar
            quotient_array = np.divide(self.array, other)

        quotient = GriddedData()
        quotient.array = xr.DataArray(quotient_array)
        return quotient

    def __mul__(self, other):
        """
        Allow GenieArray to be multiplied by a number or another GenieArray
        """
        product = GriddedData()
        if hasattr(other, "array"):
            product.array = self.array * other.array
        else:
            try:
                product.array = self.array * other
            except ValueError:
                print("Sorry, only number and GenieArray are accepted")

        return product


    def max(self, *args, **kwargs):
        self.array = self.array.max(*args, **kwargs)
        return self


    def min(self, *args, **kwargs):
        self.array = self.array.min(*args, **kwargs)
        return self


    def sum(self, *args, **kwargs):
        ## restore units
        units = self.array.units
        self.array = self.array.sum(*args, **kwargs)
        self.array.units = units
        return self        


    def mean(self, *args, **kwargs):
        self.array = self.array.mean(*args, **kwargs)
        return self
    

    def median(self, *args, **kwargs):
        self.array = np.median(self.array, *args, **kwargs)
        return self        


    def sd(self, *args, **kwargs):
        "standard deviation of the mean"
        self.array = np.std(self.array, *args, **kwargs)
        return self        

    def variance(self, *args, **kwargs):
        self.array = np.var(self.array, *args, **kwargs)
        return self    

    def se(self, *args, **kwargs):
        "standard error of the mean"
        self.array = sem(self.array, nan_policy="omit", axis=None, *args, **kwargs)
        return self        
    
    def select_basin(self, basin):
        """
        select (modern) basin from regionmask.defined_regions.ar6.ocean

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
        self.array = self.normalise_longitude().array.where(mask == index)

        return self


    def mask_basin(self, base, basin, subbasin):

        """
        similar to select_basin, but use pre-defined grid mask
        """
        gp = GridOperation()
        data = self.array
        mask = gp.GENIE_grid_mask(base=base, basin=basin, subbasin=subbasin, invert=True)

        if self.array.geo_ndim > 2:
            mask = np.broadcast_to(mask, (16, 36, 36))

        mask_data = np.ma.array(data, mask=mask)
        mask_data = np.ma.masked_invalid(mask_data)

        self.array = xr.DataArray(mask_data, dims=data.dims, coords=data.coords)        

        return self

    def ocn_only_data(self, index=False, center_position=None):
        """
        remove the NA grid (land in cGENIE definition)

        :param index: return GeoIndex or value
        :param center_position: a list/tuple of target position to reduce the output size
        """

        if center_position:
            lat, lon = center_position
            tolerance = 20 #20 degree
            reduced_array = self.array.sel(lon=slice(lon - 20, lon + 20),
                                           lat=slice(lat - 20, lat + 20))
            
            stacked_data = reduced_array.stack(x=self.array.dims)
        else:
            stacked_data = self.array.stack(x=self.array.dims)
            
        stacked_ocn_mask =~np.isnan(stacked_data)
        ocn_only_data = stacked_data[stacked_ocn_mask]
        
        if index:
            return ocn_only_data.x.values
        else:            
            return ocn_only_data
        

    def search_grid(self, point, ignore_na=False, to_genielon=False):
        """
        search the nearest grid point to the given coordinates
        
        :param point: a list/tuple of coordinate values
        :param ignore_na: whether only check ocean data (which ignore the NA grids)
        :param to_genielon: whether convert the input longitude to genie longitude, input point must be list if True       
        """
        ## ignore the first dimension (which is time)        
        geo_ndim = self.array.ndim - 1
        
        if len(point) != geo_ndim:
            raise ValueError("Input point has incompatiable coordinate")        
        
        if ignore_na:
            
            match geo_ndim:
                case 3:
                    ## point: (zt, lat, lon); Index: (time, zt, lat, lon)
                    if to_genielon: point[2] = GridOperation().lon_g2n(point[2])
                    index_pool = self.ocn_only_data(index=True, center_position=point[1:3])
                    distances = np.array([GridOperation().geo_dis3d(point, pt[1:4]) for pt in index_pool])
                case 2:
                    ## point: (lat, lon); Index: (time, lat, lon)                    
                    if to_genielon: point[1] = GridOperation().lon_g2n(point[1])
                    index_pool = self.ocn_only_data(index=True, center_position=point[0:2])                    
                    distances = np.array([GridOperation().geo_dis2d(point, pt[1:3]) for pt in index_pool])

            idx_min = np.argmin(distances)
            nearest_value = self.ocn_only_data(index=False).values[idx_min]
            return nearest_value
        else:
            match geo_ndim:
                case 3:
                    z, lon, lat = point
                    kwargs = {'lat': lat, 'lon':lon, 'zt': zt, 'method': 'nearest'}
                case 2:
                    lon, lat = point
                    kwargs = {'lat': lat, 'lon':lon, 'method': 'nearest'}
            return self.array.sel(**kwargs).values.item()

        
    def select_bbox(self, lon_min=-255, lon_max=95, lat_min=0, lat_max=90):
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
