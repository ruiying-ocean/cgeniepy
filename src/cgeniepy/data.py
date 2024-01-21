import pathlib


from metpy.interpolate import natural_neighbor_to_grid, inverse_distance_to_grid
from scipy.interpolate import griddata
from xarray import open_dataset
import numpy as np
from .grid import normalise_obs_lon


def efficient_log(data):
    "keep NA, remove zeros"
    return np.where(data == 0, -10, np.log10(data))


def foram_groups():
    """
    get a dictionary with foram abbrev (keys), pft_index and complete name (values).

    :returns: dictionary
    """
    foram_names = {
        "bn": [16, "symbiont-barren non-spinose"],
        "bs": [17, "symbiont-barren spinose"],
        "sn": [18, "symbiont-facultative non-spinose"],
        "ss": [19, "symbiont-obligate spinose"],
    }

    return foram_names

class Observation:

    def __init__(self) -> None:
        pass

    def plot(self):
        pass

    def interpolate(self):
        pass
    


# def scatter_map(
#         df: pd.DataFrame,
#         var,
#         ax,
#         x="Longitude",
#         y="Latitude",
#         interpolate=None,
#         log=False,
#         land_mask=True,
#         *args,
#         **kwargs,
# ):
#     """plot map based on dataframe with latitude/longitude
#     using cartopy as engine

#     :param df: pandas dataframe
#     :param var: variable (column) in dataframe to plot
#     :param x: coordinate attribute, default "Longitude"
#     :param y: coordinate attribute, default "Latitude"
#     :param interpolate: whether interpolate scatter data

#     :returns: a map
#     """

#     if land_mask:
#         ax.set_global()
#         # plot land and coastline, zorder is the drawing order, smaller -> backer layer
#         ax.add_feature(cfeature.LAND.with_scale('110m'), zorder=2, facecolor="#B1B2B4")
#         ax.add_feature(cfeature.COASTLINE.with_scale('110m'), zorder=3)

#     if log:
#         df[var] = efficient_log(df[var])

#     if interpolate:        
#         subdf = df[[x, y, var]]
#         subdf = subdf.dropna().astype('float64').to_numpy()
#         lat = subdf[:,0]
#         lon = subdf[:,1]
#         values = subdf[:,2]
        
#         # construct meshgrid
#         min_lat=round(min(lat))
#         max_lat=round(max(lat))
#         min_lon=round(min(lon))
#         max_lon=round(max(lon))

#         # every 1x1 pixel
#         # equivalent to
#         # grid_lat, grid_lon = np.mgrid[min_lat:max_lat:nlat*1j, min_lon:max_lon:nlon*1j]
#         grid_lat, grid_lon = np.meshgrid(np.linspace(min_lat, max_lat, max_lat-min_lat),
#                                          np.linspace(min_lon, max_lon, max_lon-min_lon))
        
#         # interpolate and return data in 2D array
#         match interpolate:
#             case 'natural_neighbor':
#                 grid_values = natural_neighbor_to_grid(lat, lon, values, grid_lat, grid_lon)
#             case 'inverse_distance':
#                 grid_values = inverse_distance_to_grid(lat, lon, values, grid_lat,grid_lon,r=3, min_neighbors=0.5)
#             case 'linear':
#                 grid_values = griddata((lat, lon), values, (grid_lat, grid_lon), method='linear')
#             case 'nearest':
#                 points = subdf[:,1:3]
#                 grid_values = griddata((lat, lon), values, (grid_lat, grid_lon), method='nearest')
#             case 'cubic':
#                 points = subdf[:,1:3]
#                 grid_values = griddata((lat, lon), values, (grid_lat, grid_lon), method='cubic')
#             case _:
#                 raise ValueError(f"Interpolation method {interpolate} not supported")

#         # plot
#         p = ax.pcolormesh(grid_lat, grid_lon,
#                           grid_values,
#                           transform=ccrs.PlateCarree(),
#                           zorder=1,
#                           *args,
#                           **kwargs)
#     else:
#         p = ax.scatter(
#             x=df[x],
#             y=df[y],
#             c=df[var],
#             transform=ccrs.PlateCarree(),
#             *args,
#             **kwargs,
#         )

#     return p
    
