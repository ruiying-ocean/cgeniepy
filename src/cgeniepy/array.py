import numpy as np
import xarray as xr
from scipy.stats import sem
import regionmask

from cgeniepy.grid import Interpolator, GridOperation
from cgeniepy.plot import GriddedDataVis
from cgeniepy.chem import Chemistry
import cgeniepy.table as ct

import warnings

class GriddedData:

    modify_in_place = False ## modify the data in place
    keep_attrs = True ## keep the attributes during operations
    
    """
    GriddedData is a class to store and compute GENIE netcdf data.

    It stores data in xarray.DataArray format, and provides optimalised methods for GENIE model output to compute statistics.    
    """
    
    def __init__(self, array=np.nan, attrs=None):
        """
        Initialise an instance of GriddedData

        :param array: a xarray.DataArray object
        :param attrs: a dictionary of attributes, default is empty

        -----------
        Example
        -----------
        data = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20]}, attrs={"units": "unitless"})
        dg = GriddedData(data)
        """

        # set data
        self.data = array
        # Fix mutable default argument issue - use None and create new dict
        self.attrs = attrs if attrs is not None else {}

        ## formatting the unit
        if 'units' in self.attrs and self.attrs['units'] is not None:
            self.attrs['units'] = Chemistry().format_unit(self.attrs['units'])

    def __repr__(self):
        return f"GriddedData\ndata={self.data}\n"

    def __getitem__(self, item):
        "make GriddedData subscriptable like xarray.DataArray"
        return self.data[item]

    def __array__(self):
        "return the numpy array of the data"
        return self.data.values

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
        output = Interpolator(dims, coords, values, *args, **kwargs).to_xarray()

        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:            
            return GriddedData(output, attrs=output_attrs)
        

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

        output_attrs = self.attrs if GriddedData.keep_attrs else {}

        if GriddedData.modify_in_place:
            self.data = self.data.sel(*args, **kwargs)
            self.attrs = output_attrs
        else:
            return GriddedData(self.data.sel(*args, **kwargs), attrs=output_attrs)

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
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = self.data.isel(*args, **kwargs)
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(self.data.isel(*args, **kwargs), attrs=self.attrs)
    
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

        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)
        
                    

    def __add__(self, other):
        """Allow GriddedData to be added by a number or another GriddedData
        """
        sum_array = GriddedData()
        if isinstance(other, GriddedData):
            self._check_same_units(other)
            sum_array.data = self.data + other.data
            sum_array.attrs = self.attrs
        else:
            ## a scalar or xarray.DataArray
            sum_array.data = self.data + other
            sum_array.attrs = self.attrs
        return sum_array

    def __radd__(self, other):
        """
        Allow a number or another GriddedData to be added to this GriddedData
        """
        sum_array = GriddedData()
        if isinstance(other, GriddedData):
            self._check_same_units(other)
            sum_array.data = other.data + self.data
            sum_array.attrs = self.attrs
        else:
            sum_array.data = other + self.data
            sum_array.attrs = self.attrs
        return sum_array

    def __sub__(self, other):
        """
        Allow GriddedData to be subtracted by a number or another GriddedData
        """
        diff = GriddedData()
        if isinstance(other, GriddedData):
            self._check_same_units(other)
            diff.data = self.data - other.data
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
        if isinstance(other, GriddedData):
            self._check_same_units(other)
            diff.data = other.data - self.data
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

    def __rtruediv__(self, other):
        """
        Allow a number or another GriddedData to be divided by this GriddedData
        """
        quotient_array = np.zeros_like(self.data)
        if isinstance(other, GriddedData):
            quotient_array = np.divide(other.data, self.data)
        else:
            # Handle division by a scalar
            quotient_array = np.divide(other, self.data)

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

    def __rmul__(self, other):
        """
        Allow a number or another GriddedData to be multiplied by this GriddedData
        """
        product = GriddedData()
        if hasattr(other, "data"):
            product.data = other.data * self.data
        else:
            try:
                product.data = other * self.data
            except ValueError:
                print("Sorry, only number and GriddedData are accepted")

        return product

    def __pow__(self, other):
        """
        Allow GriddedData to be raised to a power
        """
        return GriddedData(self.data ** other, attrs=self.attrs)

    def _check_same_units(self, other):
        "fields can only be added or subtracted in the same units"
        units = self.attrs.get("units")
        other_units = other.attrs.get("units")
        if units != other_units:
            raise ValueError(f"Cannot combine fields in {units} and {other_units}")


    ## allow comparison
    def __lt__(self, other):
        return self.data < other

    def __le__(self, other):
        return self.data <= other

    def __gt__(self, other):
        return self.data > other
            
    def __ge__(self, other):
        return self.data >= other

    def __eq__(self, other):
        return self.data == other
            
    def __ne__(self, other):
        return not self.__eq__(other)


    def max(self, *args, **kwargs):
        """compute the maximum value of the array
        """
        
        output = self.data.max(*args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)


    def min(self, *args, **kwargs):
        """compute the minimal value of the array
        """        
        output = self.data.min(*args, **kwargs)
        
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)


    def sum(self, *args, **kwargs):
        """ compute the sum of the array
        """
        
        output = self.data.sum(*args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)


    def mean(self, *args, **kwargs):
        """compute the mean of the array.
        Note this is not weighted mean, for weighted mean, use `weighted_mean` method
        """
        output = self.data.mean(*args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)
    

    def median(self, *args, **kwargs):
        """compute the median of the array
        """
        
        output= self.data.median(*args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)

    def std(self, *args, **kwargs):
        """Compute the standard deviation using xarray's API so named dimensions can be selected."""

        if isinstance(self.data, xr.DataArray):
            output = self.data.std(*args, **kwargs)
        else:
            output = np.std(self.data, *args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)

    def sd(self, *args, **kwargs):
        """Backward-compatible alias for std."""
        return self.std(*args, **kwargs)


    def variance(self, *args, **kwargs):
        "compute the variance of the array"

        output = np.var(self.data, *args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)

    def se(self, *args, **kwargs):
        "compute the standard error of the mean"
        output = sem(self.data, nan_policy="omit", axis=None, *args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)
                
    def weighted(self, weights, *args, **kwargs):
        """
        assign weights to the data
        """
        ## using xarray's API

        output= self.data.weighted(weights, *args, **kwargs)
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)


    
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

        output = self.normalise_longitude(method=norm_lon_method).data.where(cond)                    
        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)        


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

        if isinstance(basin, str):
            basin = [basin]  # Convert to list if a single basin is provided

        # Loop and merge masks
        combined_mask = None

        for bs in basin:
            mask = gp.GENIE_grid_mask(base=base, basin=bs, subbasin=subbasin, invert=True)
            
            if combined_mask is None:
                combined_mask = mask
            else:
                combined_mask = np.logical_and(combined_mask, mask)  

        # repeat the (lat, lon) mask over any leading dimensions, e.g. time and depth
        combined_mask = np.broadcast_to(combined_mask, data.shape)

        mask_data = np.ma.array(data, mask=combined_mask)
        mask_data = np.ma.masked_invalid(mask_data)

        output = xr.DataArray(mask_data, dims=data.dims, coords=data.coords)

        output_attrs = self.attrs if GriddedData.keep_attrs else {}
        if GriddedData.modify_in_place:
            self.data = output
            self.attrs = output_attrs
            return self
        else:
            return GriddedData(output, attrs=output_attrs)
    
    def _ocn_only_data(self, index=False):
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

    def save_ocn_data(self):
        self.ocn_index = self._ocn_only_data(index=True)
        self.ocn_data = self._ocn_only_data(index=False)

    def _ensure_ocn_data(self):
        if not hasattr(self, 'ocn_index') or not hasattr(self, 'ocn_data'):
            self.save_ocn_data()        

    def _curvilinear_lat_lon(self, lat_coord=None, lon_coord=None):
        """Return 2-D latitude/longitude coordinates, when present.

        Coordinates are identified from explicit names, CF metadata, or common
        model naming conventions (in that order).
        """
        if (lat_coord is None) != (lon_coord is None):
            raise ValueError("lat_coord and lon_coord must be provided together")
        if lat_coord is not None:
            try:
                lat = self.data.coords[lat_coord]
                lon = self.data.coords[lon_coord]
            except KeyError as exc:
                raise ValueError(f"Unknown coordinate: {exc.args[0]}") from exc
            if lat.ndim != 2 or lon.ndim != 2 or lat.dims != lon.dims:
                raise ValueError(
                    "Curvilinear latitude and longitude must be 2-D and use "
                    "the same dimensions"
                )
            return lat, lon

        lat = lon = None
        for coord in self.data.coords.values():
            if coord.ndim != 2:
                continue
            name = (coord.name or "").lower()
            standard_name = coord.attrs.get("standard_name", "").lower()
            units = coord.attrs.get("units", "").lower()
            if (name in {"lat", "latitude", "tlat", "ulat", "nav_lat", "geolat"}
                    or standard_name == "latitude"
                    or units in {"degrees_north", "degree_north", "degrees_n"}):
                lat = coord
            elif (name in {"lon", "longitude", "tlong", "ulong", "nav_lon", "geolon"}
                    or standard_name == "longitude"
                    or units in {"degrees_east", "degree_east", "degrees_e"}):
                lon = coord

        if lat is not None and lon is not None and lat.dims == lon.dims:
            return lat, lon
        return None

    @staticmethod
    def _vertical_to_km(values, units):
        """Convert vertical distances to kilometres using CF-style units."""
        units = (units or "m").strip().lower()
        factors = {
            "m": 1e-3, "meter": 1e-3, "meters": 1e-3,
            "metre": 1e-3, "metres": 1e-3,
            "cm": 1e-5, "centimeter": 1e-5, "centimeters": 1e-5,
            "centimetre": 1e-5, "centimetres": 1e-5,
            "km": 1.0, "kilometer": 1.0, "kilometers": 1.0,
            "kilometre": 1.0, "kilometres": 1.0,
        }
        return np.asarray(values) * factors.get(units, 1e-3)

    def _search_curvilinear(self, point, ignore_na, lat_coord=None, lon_coord=None):
        """Search data whose geographic coordinates are two-dimensional."""
        lat, lon = self._curvilinear_lat_lon(lat_coord, lon_coord)
        horizontal_dims = lat.dims
        vertical_dims = [dim for dim in self.data.dims if dim not in horizontal_dims]

        if len(vertical_dims) > 1 or len(point) != self.data.ndim:
            raise ValueError(
                "Curvilinear search supports 2-D fields or 3-D fields with one "
                "vertical dimension; select a time index before searching"
            )

        if vertical_dims:
            vertical_dim = vertical_dims[0]
            z_target, lat_target, lon_target = point
        else:
            vertical_dim = None
            lat_target, lon_target = point

        horizontal_distance = GridOperation().haversine_distance(
            lat_target, lon_target, lat.values, lon.values
        )

        if not ignore_na:
            iy, ix = np.unravel_index(
                np.nanargmin(horizontal_distance), horizontal_distance.shape
            )
            indexers = {horizontal_dims[0]: iy, horizontal_dims[1]: ix}
            if vertical_dim is not None:
                z = self.data[vertical_dim]
                indexers[vertical_dim] = int(np.abs(z.values - z_target).argmin())
            return self.data.isel(indexers).values.item()

        values = self.data.values
        valid = ~np.isnan(values)
        if not np.any(valid):
            return np.nan

        if vertical_dim is None:
            distances = horizontal_distance
        else:
            z = self.data[vertical_dim]
            z_distance = self._vertical_to_km(
                np.abs(z.values - z_target), z.attrs.get("units")
            )
            # Transpose distance components into the data's dimension order.
            horizontal = xr.DataArray(
                horizontal_distance, dims=horizontal_dims
            ).broadcast_like(self.data)
            vertical = xr.DataArray(
                z_distance, dims=(vertical_dim,)
            ).broadcast_like(self.data)
            distances = np.sqrt(horizontal.values ** 2 + vertical.values ** 2)

        distances = np.where(valid, distances, np.inf)
        return values.flat[np.argmin(distances)]

    def search_point(self, point, ignore_na=False, to_genielon=False,
                     lat_coord=None, lon_coord=None, **kwargs):
        """
        search the nearest grid point to the given coordinates
        
        :param point: a list/tuple of coordinate values in the same order to the data dimension (self.data.dims)
        :param ignore_na: whether only check ocean data (which ignore the NA grids)
        :param to_genielon: whether convert the input longitude to genie longitude, input point must be list if True
        :param lat_coord: optional name of a 2-D curvilinear latitude coordinate
        :param lon_coord: optional name of a 2-D curvilinear longitude coordinate

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

        curvilinear = self._curvilinear_lat_lon(lat_coord, lon_coord)
        if curvilinear is not None:
            return self._search_curvilinear(
                point, ignore_na, lat_coord=lat_coord, lon_coord=lon_coord
            )
        
        if ignore_na:
            self._ensure_ocn_data()
            index_pool = self.ocn_index
            
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
            
            nearest_value = self.ocn_data.values[idx_min]
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

    def fill_poles(self):
        """
        Fill the poles with a given value, default is NaN
        """
        # Create a regular latitude grid from -90 to 90
        if hasattr(self.data, 'lat'):
            current_lats = self.data.lat.values
        else:
            raise ValueError("DataArray must have a 'lat' coordinate to fill poles.")

        n_points = len(current_lats) + 2  # or choose your desired resolution
        regular_lats = np.linspace(-90, 90, n_points)

        # If you have data associated with these latitudes, you can interpolate
        # Example assuming you have a DataArray called 'da':
        da_interp = self.data.interp(lat=regular_lats, method='linear')

        self.data = da_interp
        return self
