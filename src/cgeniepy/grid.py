from . import ureg, Q_

import numpy as np
import pandas as pd
import xarray as xr

from scipy.interpolate import RegularGridInterpolator


def lon_n2g(x):
    """
    Change parts of observational latitude [100, 180] to GENIE longitude [-260, -180]
    Note it isn't axisymmetric!
    """
    # TBD: add value range checker

    if x > 100 and x < 180:
        return x - 360
    else:
        return x


def lon_g2n(x):
    """
    Change parts of observational latitude [100, 180] to GENIE longitude [-260, -180]
    CANNOT simply use +80, or -80! It isn't axisymmetric!
    """
    # TBD: add value range checker

    if x < -180:
        return x + 360
    else:
        return x


def normalise_obs_lon(data: xr.Dataset) -> xr.Dataset:
    return data.assign_coords({"lon": list(map(lon_n2g, data.lon.values))}).sortby(
        "lon"
    )


def normalise_GENIE_lon(data: xr.Dataset) -> xr.Dataset:
    """
    Change parts of observational latitude [100, 180] to GENIE longitude [-260, -180]
    """
    ## if hasattr(data, "assign_coords"):
    return data.assign_coords({"lon": list(map(lon_g2n, data.lon.values))}).sortby(
        "lon"
    )


def mask_Arctic_Med(array, policy="na"):
    """
    mask Arctic and Meditterean Sea in cGENIE modern continent configuration
    """
    if policy == "na":
        array[34:36, :] = np.nan
        array[27:30, 25:30] = np.nan
    elif policy == "zero":
        array[34:36, :] = 0
        array[27:30, 25:30] = 0

    return array


def GENIE_grid_mask(
    base="worjh2", basin="ALL", subbasin="", mask_Arc_Med=False, invert=False
):
    """
    Get a modern GENIE 36x36 mask array from input data.
    The input array is flipped (left/right flip -> up/down flip) for easy recognition

    :continent: worjh2, worlg4, worbe2, GIteiiaa, GIteiiaa, p0055c
    :basin: Atlantic/Pacific/Indian/ALL/Tanzania
    :subbasin: N/S/ALL, ALL means Southern Ocean section included

    :returns: GENIE grid array where continent/ice cap is 0 and ocean is 1, default is 'worjh2'
    """

    file_path = (
        pathlib.Path(__file__).parent.parent
        / f"data/mask_{base}_{basin}{subbasin}.txt"
    )
    grid_mask_raw = np.loadtxt(file_path, dtype=int)
    grid_mask = np.flip(np.fliplr(grid_mask_raw))

    if mask_Arc_Med:
        grid_mask = mask_Arctic_Med(grid_mask, policy="zero")

    if invert:
        grid_mask = ~grid_mask + 2

    return grid_mask


def normal_lon(N=36, edge=False):
    """
    Normal longitude in 10 degree resolution,
    if edge is False, then return midpoint
    """
    if edge:
        lon_edge = np.linspace(-180, 180, N + 1)
        return lon_edge
    else:
        lon = np.linspace(-175, 175, N)
        return lon


## check grid area (1)
# print(marine_area.sum()*1e-8) #around 3.7

# check grid area (2)
# np.nansum(np.where(mask_array == 1, grid_area, np.nan))*1e-8 #same

# check grid area (3) plot out
# plt.pcolormesh(marine_area)
# plt.colorbar()


# check grid area (4) compare to biogem grid area
# bgem_grid_area = xr.load_dataset("~/Downloads/fields_biogem_2d.nc")
# np.nansum(np.where(mask_array == 1, bgem_grid_area, np.nan))*1e-6*1e-8


# def interp_GENIE(file_path, var):
#     """
#     Interpolate GENIE into 360x360 grids using xesmf API. Depends on Xarray.
#     Known limitations:
#         (1) unexpected NAs grids in boundaries, e.g, Arctic and 95-degree longitude;
#         (2) thus necessary offset in the map;
#         (3) ugly way to draw continent lines

#     :param file_path: string, netcdf file
#     :parm var: string, variable name
#     """

#     # ------Prepare finer grids-------------
#     coarse_data = open_dataset(file_path)
#     # get 10 times higher resolution
#     #fine_lat = np.rad2deg(np.arcsin(np.linspace(-0.97222222, 0.97222222, 360)))  # sin(x) -> radian -> degree
#     #fine_lon = np.linspace(coarse_data.lon[0], coarse_data.lon[-1], 360)
#     fine_lat = GENIE_lat(360)
#     fine_lon = GENIE_lon(360)

#     # prepare finer grid
#     fine_data = Dataset(
#         {"lat": (["lat"], fine_lat),
#          "lon": (["lon"], fine_lon)}
#     )

#     # ----- Mask coarse & fine grids ---------
#     coarse_data['mask'] = np.where(~np.isnan(coarse_data[var].isel(time=0)), 1, 0)  # NA data to 0, otherwise 1

#     # mask output grid
#     mask_array_coarse = coarse_data['mask'].values
#     mask_array_fine = np.zeros((360, 360))

#     #359 or 360?
#     for i in range(360):
#         for j in range(360):
#             if mask_array_coarse[i // 10, j // 10] != 0:
#                 mask_array_fine[i, j] = 1

#     mask_array_fine_da = DataArray(
#         data=mask_array_fine,
#         dims=["lat", "lon"],
#     )

#     fine_data['mask'] = mask_array_fine_da

#     # --------- Regrid -----------
#     regridder = xe.Regridder(coarse_data, fine_data, method="bilinear", extrap_method="inverse_dist")
#     fine_data = regridder(coarse_data[var])

#     return fine_data


def regrid_lat(x):

    """
    Transform <latitude> into cGENIE resolution to facilitate comparison between
    model and observational data.
    """
    if x >= -90 and x <= 90:
        lat_edge = np.rad2deg(np.arcsin(np.linspace(-1, 1, 37)))
        lat = GENIE_lat()

        for i in range(36):
            if x > lat_edge[i] and x <= lat_edge[i + 1]:
                x = lat[i]
    else:
        raise ValueError("Latitude must be in [-90,90]")

    return x


def regrid_lon(x):

    """
    Transform <longitude> into cGENIE resolution to facilitate comparison between
    model and observational data.
    """

    if x >= -180 and x <= 180:
        lon_edge = np.linspace(-180, 180, 37)
        for i in range(36):
            if x > lon_edge[i] and x <= lon_edge[i + 1]:
                x = (lon_edge[i] + lon_edge[i + 1]) / 2  # middle value in the bin
    else:
        raise ValueError("Longitude must be in [-180,180]")

    return x


def new_basin_mask(mask_array, basin, filename):
    """
    a dirty way to copy the basin mask from worjh2 and facilitate further manual modification
    """

    x = mask_array + GENIE_grid_mask(base="worjh2", basin=basin)
    o = np.where(x > 1, 1, 0)
    o = np.flip(np.fliplr(o))
    np.savetxt(filename, o, fmt="%i")
    print(f"array saved in {filename}")


def regrid_dataframe(
    dataframe,
    low_threshold=None,
    new_low_bound=None,
    high_threshold=None,
    new_high_bound=None,
):
    """
    Regrid a dataframe within certain format to cGENIE grids

    Input: A dataframe with three necessary columns: "Latitude", "Longitude", "Observation".
    The other columns will be unselected. Note that observation column is the data you want to process.
    Rename dataframe by using command: `df = df.rename({"old name": "new name"}, axis='columns')`

    Output:
    A 36x36 2D array.

    Optional parameters:
    low_threshold/new_low_bound: change the data that is smaller than low_threshold to new_low_bound.
    For example, covert any data < 0.03 to 0. Work same to high_threshold/new_high_bound.
    """

    df = dataframe.copy()

    # subset
    df = df[["Longitude", "Latitude", "Observation"]]

    # drop NAN
    df = df.dropna(axis="rows", how="any")

    # regrid coordinate
    df.loc[:, "Latitude"] = df.loc[:, "Latitude"].apply(regrid_lat)
    df.loc[:, "Longitude"] = df.loc[:, "Longitude"].apply(regrid_lon)

    # group and aggregate (still in long format)
    df_agg = df.groupby(["Longitude", "Latitude"]).agg("mean")

    # Create a all-nan data frame with full cGENIE grids
    lat = GENIE_lat()
    lon = normal_lon()
    data = np.zeros([36 * 36])
    index = pd.MultiIndex.from_product([lon, lat], names=["Longitude", "Latitude"])
    df_genie = pd.DataFrame(data, index=index, columns=["Observation"])
    df_genie["Observation"] = df_genie["Observation"].replace(0, np.nan)

    # transform the tuple elements in multindex to lists
    df_genie_index_list = [list(item) for item in df_genie.index.values]
    df_agg_index_list = [list(item) for item in df_agg.index.values]

    # copy aggregated dataframe to full-grid dataframe, one by one
    # if the returned dataframe is full of NA, it's very likely from this step
    for i in df_genie_index_list:

        longitude = i[0]
        latitude = i[1]

        if i in df_agg_index_list:
            df_genie.loc[longitude, latitude] = df_agg.loc[longitude, latitude]
        else:
            df_genie.loc[longitude, latitude] = np.nan

    # filter data if necessary
    if (low_threshold is not None) and (new_low_bound is not None):
        df_genie.loc[
            df_genie.Observation < low_threshold, "Observation"
        ] = new_low_bound

    if (high_threshold is not None) and (new_high_bound is not None):
        df_genie.loc[
            df_genie.Observation > high_threshold, "Observation"
        ] = new_high_bound

    # long -> wide data format
    df_genie_wide = df_genie.pivot_table(
        values="Observation", index="Latitude", columns="Longitude", dropna=False
    )

    return df_genie_wide


class regridder:
    """
    regridder class, make cGENIE to finer resolution
    """
    
    def __init__(self, array, grid_number=200, method='linear'):
        """
        initialize regridder
        
        :param array: xr data array
        :param grid_number: number of target grid points
        :param method: interpolation method, only linear is supported
        """
        self.array = array
        ## a tuple of dimensions
        self.dims = array.dims
        ## dimension values
        self.coords = tuple([array[dim].values for dim in self.dims])
        ## array values
        self.values = array.values

        ## interpolation function
        self.interp_function = self._create_interp_function(method=method)        
        ## number of grid points
        self.grid_number = grid_number        
        ## create new coordinates, meshgrid
        self.gridded_coord = self.new_coordinate(n=grid_number)      
        self.meshgrid = self.new_meshgrid(*self.gridded_coord, indexing='ij')
        self.gridded_data = self.interpolate_data(tuple(self.meshgrid))
        
    def _create_interp_function(self, method):        
        interp_function = RegularGridInterpolator(self.coords,
                                                  self.values,
                                                  method=method,
                                                  bounds_error=True, 
                                                  fill_value=None)
        return interp_function

    def new_coordinate(self, n):
        new_coords = []
        for dim_values in self.coords:
            min_val, max_val = np.min(dim_values), np.max(dim_values)
            new_values = np.linspace(min_val, max_val, n)
            new_coords.append(new_values)
        return new_coords

    def new_meshgrid(self, *args, **kwargs):
        "regrid to finer dimensions"        
        return np.meshgrid(*args, **kwargs)

    def interpolate_data(self, *args, **kwargs):        
        "Use the interpolation function to get regridded values"
        return self.interp_function(*args, **kwargs)
    
    def to_xarray(self):
        """
        convert numpy array to xr data array
        """
        return xr.DataArray(self.gridded_data, dims=self.dims, coords=self.gridded_coord)
        
    # def to_geniearray(self):
    #     """
    #     convert numpy array to cGENIE 36x36 array
    #     """
    #     return GenieArray(self.to_xarray())
