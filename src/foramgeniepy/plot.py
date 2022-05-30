import numpy as np
import cartopy.crs as ccrs
from cartopy.feature import LAND
from matplotlib import cm

from .data import efficient_log
from .grid import GENIE_lat, GENIE_lon

## TODO [] add contour layer

def plot_GENIE(ax, data, log=False, grid_line = False, continent_outline=True, contour_layer=False, cmap=cm.Spectral_r, *args, **kwargs):
    """
    plot map for 2D GENIE time-slice data

    :param data: 2D numpy array
    :param ax: matplotlib axes object with cartopy projection. e.g., projection=ccrs.EckertIV()
    :param vmin/vmax: float, colorbar limit
    :param log: logic, whether take log for data
    :param grid_line: logic, whthere plot grid line
    :param continent_outline: logic whethere plot continent
    :param cmap: color map object

    :returns: a mappable plot object
    """

    data_crs = ccrs.PlateCarree()  # should not change

    # -------------------Box line------------------------
    ax.spines['geo'].set_edgecolor('black')
    ax.spines['geo'].set_linewidth(1)

    # -------------------Facecolor-----------------------
    ax.patch.set_color("silver")

    # -------------------Plot-----------------------
    if log: data = efficient_log(data)
    lon_edge = GENIE_lon(edge=True)
    lat_edge = GENIE_lat(edge=True)

    #cartopy transform seems to help reassign the GENIE longitude to normal
    p = ax.pcolormesh(lon_edge, lat_edge, data, cmap=cmap, transform=data_crs, shading="flat", *args, **kwargs)

    if contour_layer:
        lat = GENIE_lat(edge=False)
        lon = GENIE_lon(edge=False)
        cs = ax.contour(lon, lat ,data, transform=ccrs.PlateCarree(), cmap=cmap)
        #label every three levels
        ax.clabel(cs, cs.levels[::3], colors=['black'], fontsize=8, inline=False)

    # -------------------Grid lines-----------------------
    if grid_line:
        ax.gridlines(crs=data_crs, draw_labels=False, linewidth=0.5, color='gray', alpha=0.5, linestyle='-')

    # -------------------Continent lines-----------------------
    if continent_outline:
        outline_color = "black"
        outline_width = 1
        mask_array = np.where(~np.isnan(data), 1, 0)

        Nlat_edge  = len(lat_edge)
        Nlon_edge  = len(lon_edge)
        # dimension of array
        Nlat_array = Nlat_edge - 1
        Nlon_array = Nlon_edge - 1
        # index of array
        Nlat_index = Nlat_array - 1
        Nlon_index = Nlon_array - 1

        # i stands for latitude index, j stands for longitude index
        # e.g., mask_array[0, ] is Antarctic ice cap
        # i,j are for mask array and will be converted to lon/lat edge lines
        for i in range(Nlat_array):
            for j in range(Nlon_array):
                # compare with the right grid, and plot vertical line if different
                if j < Nlon_index and mask_array[i, j] != mask_array[i, j+1]:
                    ax.vlines(
                        lon_edge[j+1],
                        lat_edge[i],
                        lat_edge[i+1],
                        color=outline_color, linewidth=outline_width, transform=data_crs)

                # connect the circular longitude axis
                if j == Nlon_index and mask_array[i, j] != mask_array[i, 0]:
                    ax.vlines(
                        lon_edge[j+1],
                        lat_edge[i],
                        lat_edge[i+1],
                        color=outline_color, linewidth=outline_width, transform=data_crs)

                # compare with the above grid, and plot horizontal line if different
                if i < Nlat_index and mask_array[i, j] != mask_array[i+1, j]:
                    ax.hlines(
                        lat_edge[i+1],
                        lon_edge[j],
                        lon_edge[j+1],
                        colors=outline_color,linewidth=outline_width, transform=data_crs)

    return p


def scatter_map(ax, df, cmap=cm.Spectral_r, *args, **kwargs):
    ax.set_global()
    ax.coastlines()
    ax.add_feature(LAND, zorder=0, facecolor="#B1B2B4", edgecolor="white")  # zorder is drawing sequence

    p = ax.scatter(y=df.Latitude,
                   x=df.Longitude,
                   s=23,
                   transform=ccrs.PlateCarree(),
                   cmap=cmap,
                   *args, **kwargs)

    return p

def scatter_on_GENIE(ax, df, var, cmap=cm.Spectral_r, log=False, *args, **kwargs):

    if 'Latitude' not in df.columns or 'Longitude' not in df.columns:
        raise ValueError("Input data lack Latitude/Longitude column")

    if log:
        df[var] = efficient_log(df[var])

    p = ax.scatter(y=df.Latitude, x=df.Longitude, c=df[var], linewidths=.5,
                   cmap=cmap, edgecolors='black',
                   transform = ccrs.PlateCarree(), *args, **kwargs)

    return p
