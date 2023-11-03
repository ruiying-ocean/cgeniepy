import numpy as np
import xarray as xr
from scipy.stats import sem
import regionmask

from . import Q_, ureg
from .grid import normalise_GENIE_lon, GENIE_grid_mask, mask_Arctic_Med
from .chem import format_unit, pure_unit, molecular_weight
from .plot import GeniePlottable


def attr_conservation(cal_func):
    """
    A decorator to keep the attribution after array calculation

    This is because xarray will remove the units and long_name after calculation
    """
    def wrappered_func(self, *args, **kwargs):
        # Temporially store the original units
        units = self.array.units
        long_name = self.array.long_name
        
        # Call the decorated function
        result = cal_func(self, *args, **kwargs)
        
        # Assign back the unit to the array
        self.array.attrs['units'] = units
        self.array.attrs['long_name'] = long_name
        return result
    return wrappered_func


class GenieArray(GeniePlottable):
    """
    GenieArray is a class to store and compute GENIE netcdf data.

    It stores data in xarray.DataArray format, and provides a set of methods to compute statistics and plot.
    
    Particularly, it provides:

    1) ready-to-use plotting system
    2) subset grid cells (select basin, search by lat/lon, etc.)
    3) statistics (median, mean, sd, variance, se)
    """
    
    def __init__(self):
        """
        Initialise an empty array, default as NaN
        """
        
        # set data
        self.array = self._set_array()
        
        # init plottable instance
        super().__init__(array=self.array)

    def _set_array(self) -> xr.DataArray:
        """
        a function to be overwritten by subclass
        always return a xarray.DataArray
        """
        arr = np.nan
        return xr.DataArray(arr)

    def __getitem__(self, item):
        "make GenieArray subscriptable like xarray.DataArray"
        return self.array[item]

    def uarray(self):
        """convert array to pint.Quantity

        This funtion is useful when doing unit conversion
        """
        uarray = Q_(self.array.values, self.array.units)
        return uarray

    @attr_conservation
    def sel(self, *args, **kwargs):
        "a wrapper to xarray `sel` method"
        try:
            self.array = self.array.sel(*args, **kwargs)
            return self
        except:
            print("This array instance does not contain xarray.DataArray")

    @attr_conservation
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
        self.array = normalise_GENIE_lon(self.array)
        return self

    def _run_method(self, method: str, *args, **kwargs):
        "an alias to run stat for GenieArray class"
        return getattr(self, method)(*args, **kwargs)

    def __add__(self, other):
        """
        Allow GenieArray to be added by a number or another GenieArray
        """
        sum = GenieArray()
        if hasattr(other, "array"):
            sum.array = self.array + other.array
        else:
            sum.array = self.array + other
        return sum

    def __sub__(self, other):
        """
        Allow GenieArray to be subtracted by a number or another GenieArray
        """
        diff = GenieArray()
        if hasattr(other, "array"):
            diff.array = self.array - other.array
        else:
            diff.array = self.array - other
        return diff

    def __truediv__(self, other):
        """
        Allow GenieArray to be divided by a number or another GenieArray
        It will ignore the NA values (i.e., land grid points) when dividing
        """
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
                print("Only number and GenieArray are accepted")

        return quotient

    def __mul__(self, other):
        """
        Allow GenieArray to be multiplied by a number or another GenieArray
        """
        product = GenieArray()
        if hasattr(other, "array"):
            product.array = self.array * other.array
        else:
            try:
                product.array = self.array * other
            except ValueError:
                print("Sorry, only number and GenieArray are accepted")

        return product

    @attr_conservation
    def max(self, *args, **kwargs):
        self.array = self.array.max(*args, **kwargs)
        return self

    @attr_conservation
    def min(self, *args, **kwargs):
        self.array = self.array.min(*args, **kwargs)
        return self

    @attr_conservation
    def sum(self, *args, **kwargs):
        ## restore units
        units = self.array.units
        self.array = self.array.sum(*args, **kwargs)
        self.array.units = units
        return self        

    @attr_conservation
    def mean(self, *args, **kwargs):
        self.array = self.array.mean(*args, **kwargs)
        return self
    
    @attr_conservation
    def median(self, *args, **kwargs):
        self.array = np.median(self.array, *args, **kwargs)
        return self        

    @attr_conservation
    def sd(self, *args, **kwargs):
        self.array = np.std(self.array, *args, **kwargs)
        return self        

    @attr_conservation
    def variance(self, *args, **kwargs):
        self.array = np.var(self.array, *args, **kwargs)
        return self
    
    @attr_conservation
    def se(self, *args, **kwargs):
        self.array = sem(self.array, nan_policy="omit", axis=None, *args, **kwargs)
        return self        

    @attr_conservation
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

    @attr_conservation
    def mask_basin(self, base, basin, subbasin):

        """
        similar to select_basin, but use pre-defined grid mask
        """
        
        data = self.array
        mask = GENIE_grid_mask(base=base, basin=basin, subbasin=subbasin, invert=True)

        if self.array.ndim > 2:
            mask = np.broadcast_to(mask, (16, 36, 36))

        mask_data = np.ma.array(data, mask=mask)
        mask_data = np.ma.masked_invalid(mask_data)

        self.array = xr.DataArray(mask_data, dims=data.dims, coords=data.coords)        

        return self

    def search_grid(self, *args, **kwargs):
        """
        search the nearest grid point to the given coordinates
        
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
        >>>     var = var.search_grid(**kwargs).values
        >>>     modelvar.append(var)
        """
        return self.array.sel(*args, **kwargs)

    def select_box(self, lon_min=-255, lon_max=95, lat_min=0, lat_max=90):
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

    def mask_Arctic_Med(self, *args, **kwargs):
        """
        only support modern grid configuration
        """
        self.array = mask_Arctic_Med(self.array, *args, **kwargs)
        return self

            
    # def compare_obs(self, obs, *args, **kwargs):
    #     return ModelSkill(model=self.array, observation=obs, *args, **kwargs)


    ## TODO
    ## [ ] foraminifera relative abundance, calcite
    ## [ ] Plankton biomass, export
    ## [ ] Check unit when conduct calculation
