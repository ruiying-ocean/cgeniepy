from metpy.interpolate import natural_neighbor_to_grid, inverse_distance_to_grid
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .utils import efficient_log


class GeoTable:

    def __init__(self, path, lat_col, lon_col,  var_col,  *args, **kwargs):
        self.df = pd.read_csv(path, *args, **kwargs)
        self.lat_col = lat_col
        self.lon_col = lon_col
        self.var_col = var_col

    def scatter_map(
            self,
            ax = None,
            interpolate=None,
            log=False,
            land_mask=True,
            *args,
            **kwargs,
    ):
        """plot map based on dataframe with latitude/longitude
        using cartopy as engine

        :param df: pandas dataframe
        :param var: variable (column) in dataframe to plot
        :param x: coordinate attribute, default "Longitude"
        :param y: coordinate attribute, default "Latitude"
        :param interpolate: whether interpolate scatter data

        :returns: a map
        """
        if not ax:
            fig, ax = self._init_fig(subplot_kw={'projection': ccrs.EckertIV()})

        if land_mask:
            ax.set_global()
            # plot land and coastline, zorder is the drawing order, smaller -> backer layer
            ax.add_feature(cfeature.LAND.with_scale('110m'), zorder=2, facecolor="#B1B2B4")
            ax.add_feature(cfeature.COASTLINE.with_scale('110m'), zorder=3)

        if log:
            self.df[self.var_col] = efficient_log(self.df[self.var_col])

        if interpolate:        
            subdf = self.df[[self.lon_col, self.lat_col, self.var_col]]
            subdf = subdf.dropna().astype('float64').to_numpy()
            lat = subdf[:,0]
            lon = subdf[:,1]
            values = subdf[:,2]
            
            # construct meshgrid
            min_lat=round(min(lat))
            max_lat=round(max(lat))
            min_lon=round(min(lon))
            max_lon=round(max(lon))

            # every 1x1 pixel
            # equivalent to
            # grid_lat, grid_lon = np.mgrid[min_lat:max_lat:nlat*1j, min_lon:max_lon:nlon*1j]
            grid_lat, grid_lon = np.meshgrid(np.linspace(min_lat, max_lat, max_lat-min_lat),
                                            np.linspace(min_lon, max_lon, max_lon-min_lon))
            
            # interpolate and return data in 2D array
            match interpolate:
                case 'natural_neighbor':
                    grid_values = natural_neighbor_to_grid(lat, lon, values, grid_lat, grid_lon)
                case 'inverse_distance':
                    grid_values = inverse_distance_to_grid(lat, lon, values, grid_lat,grid_lon,r=3, min_neighbors=0.5)
                case 'linear':
                    grid_values = griddata((lat, lon), values, (grid_lat, grid_lon), method='linear')
                case 'nearest':
                    points = subdf[:,1:3]
                    grid_values = griddata((lat, lon), values, (grid_lat, grid_lon), method='nearest')
                case 'cubic':
                    points = subdf[:,1:3]
                    grid_values = griddata((lat, lon), values, (grid_lat, grid_lon), method='cubic')
                case _:
                    raise ValueError(f"Interpolation method {interpolate} not supported")

            # plot
            p = ax.pcolormesh(grid_lat, grid_lon,
                            grid_values,
                            transform=ccrs.PlateCarree(),
                            zorder=1,
                            *args,
                            **kwargs)
        else:
            p = ax.scatter(
                x=self.df[self.lon_col],
                y=self.df[self.lat_col],
                c=self.df[self.var_col],
                transform=ccrs.PlateCarree(),
                *args,
                **kwargs,
            )

        return p
    
    def _init_fig(self, *args, **kwargs):
        return plt.subplots(dpi=120, *args, **kwargs)
