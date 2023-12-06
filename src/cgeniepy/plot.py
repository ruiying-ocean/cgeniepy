import pathlib
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from metpy.interpolate import natural_neighbor_to_grid, inverse_distance_to_grid
from scipy.interpolate import griddata


import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import ListedColormap, to_rgb as hex_to_rgb
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf

from .data import efficient_log
from .grid import GENIE_lat, GENIE_lon, GENIE_depth


def scatter_map(
    df: pd.DataFrame,
    var,
    ax,
    x="Longitude",
    y="Latitude",
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

    if land_mask:
        ax.set_global()
        # plot land and coastline, zorder is the drawing order, smaller -> backer layer
        ax.add_feature(cfeature.LAND.with_scale('110m'), zorder=2, facecolor="#B1B2B4")
        ax.add_feature(cfeature.COASTLINE.with_scale('110m'), zorder=3)

    if log:
        df[var] = efficient_log(df[var])

    if interpolate:        
        subdf = df[[x, y, var]]
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
            x=df[x],
            y=df[y],
            c=df[var],
            linewidths=0.5,
            edgecolors="black",
            transform=ccrs.PlateCarree(),
            *args,
            **kwargs,
        )

    return p

def community_palette(cmap_name, N=256, reverse=False, alpha=None):
    """
    community-driven colormaps with multiple sources

    :param cmap_name: colormap name, can be found in avail_palette()
    :type cmap_name: str

    :returns: colormap

    XML data: https://sciviscolor.org/colormaps/
    txt data: from original packages
    """

    if cmap_name not in avail_palette():
        raise ValueError(f"{cmap_name} not found, accepted values are {avail_palette()}")
    
    data_dir = pathlib.Path(__file__).parent.parent

    file_path = None
    file_ext = None
    colors = []
    c = None

    # Search for the directory based on list comprehension
    file_path = [f for f in data_dir.glob("**/*") if cmap_name in str(f)]
    if len(file_path) > 1:
        raise ValueError("Multiple colormaps found")
    else:
        file_path = file_path[0]
        file_ext = file_path.suffix    

    if file_ext == ".txt":
        colors = pd.read_csv(file_path, header=None).values.tolist()
        colors = [colors[i][0] for i in range(len(colors))]
        c = ListedColormap(colors, name=cmap_name)
    elif file_ext == ".xml":
        tree = ET.parse(file_path)
        root = tree.getroot()
        positions = []

        # Extract color points and their positions from the XML
        for point in root.findall('.//Point'):
            r = float(point.get('r'))
            g = float(point.get('g'))
            b = float(point.get('b'))
            colors.append((r, g, b))
            positions.append(float(point.get('x')))

        c = ListedColormap(colors, name=cmap_name)

    # Optionally add transparency
    if alpha is not None:
        if file_ext == ".txt":
            rgb_array = np.array([hex_to_rgb(i) for i in colors])
            rgb_array[:, -1] = alpha
            c = ListedColormap(rgb_array, N=N)
        elif file_ext == ".xml":
            # Assuming alpha is to be applied uniformly to all colors in the XML-based colormap
            rgba_colors = [(color[0], color[1], color[2], alpha) for color in colors]
            c = ListedColormap(rgba_colors, name=cmap_name)

    if reverse and c is not None:
        return c.reversed()

    ## if _r exsit in the string, reverse the colormap
    if "_r" in cmap_name and c is not None:
        return c.reversed()
    
    if c is None:
        raise ValueError("Colormap could not be created")

    return c

def avail_palette():
    """return a list of colormap names"""
    data_dir = pathlib.Path(__file__).parent.parent

    return [f.stem for f in data_dir.glob("data/colormaps/*") if f.suffix in [".txt", ".xml"]]

def cbar_wrapper(plotting_func):
    """a decorator to add color bar

    :param plotting_func: function returning a mappable object
    :returns: plotting function with colorbar
    """

    def wrappered_func(*args, **kwargs):
        p = plotting_func(*args, **kwargs)
        cbar = plt.colorbar(p, fraction=0.046, pad=0.04, orientation="horizontal")
        cbar.ax.tick_params(color="k", direction="in")
        cbar.outline.set_edgecolor('black')
        cbar.minorticks_on()
        
    wrappered_func.__unwrapped__ =  plotting_func

    return wrappered_func

def cond_cbar(func):
    def funct_w_cond(colorbar):
        return func if colorbar else func.__unwrappered__
    return funct_w_cond

## TODO: add more layers: quiver,

class GeniePlottable:

    transform_crs = ccrs.PlateCarree()  # do not change

    grid_dict = {
    "lon": GENIE_lon(edge=False),
    "lat": GENIE_lat(edge=False),
    "zt":  GENIE_depth(edge=False)/1000,

    "lon_edge": GENIE_lon(edge=True),
    "lat_edge": GENIE_lat(edge=True),
    "zt_edge": GENIE_depth(edge=True)/1000,
    }

    
    def __init__(self, array):
        self.array = array

    def plot(self, *args, **kwargs):
        if self.array.ndim == 1:
            return self._plot_1d(*args, **kwargs)
        elif self.array.ndim == 2:
            return self._plot_2d(*args, **kwargs)
        elif self.array.ndim == 3:
            return self._plot_3d(*args, **kwargs)
        else:
            raise ValueError(f"{self.ndim} dimensions not supported")

    def _plot_1d(self, *args, **kwargs):
        """
        plot 1D data, e.g., zonal_average, time series
        dim: time/lon/lat
        """
        if 'ax' not in kwargs:
            fig, local_ax = self._init_fig()
        else:
            local_ax = kwargs.pop('ax')

        ## get the only dimension as x
        dim = self.array[self.array.dims[0]]

        if not 'swap_xy' in kwargs:            
            local_ax.set_xlabel(dim.name)
            local_ax.set_ylabel(f"{self.array.long_name} ({self.array.units})")
            p = local_ax.plot(dim, self.array, *args, **kwargs)
        else:
            local_ax.set_xlabel(f"{self.array.long_name} ({self.array.units})")
            local_ax.set_ylabel(dim.name)
            p = local_ax.plot(self.array, dim, *args, **kwargs)

        return p

    def _plot_2d(self, *args, **kwargs):
        ## if lon, lat then plot map
        ## if lon, zt then plot cross section
        ## if lat, zt then plot cross section
        dims = self.array.dims
        if 'lon' in dims and 'lat' in dims:
            return self._plot_map(*args, **kwargs)
        elif 'zt' in dims and 'lon' in dims:
            return self._plot_cross_section(*args, **kwargs)
        elif 'zt' in dims and 'lat' in dims:
            return self._plot_cross_section(*args, **kwargs)
        else:
            raise ValueError(f"{dims} not supported")

    def _plot_map(self, x_edge="lon_edge", y_edge="lat_edge", contour=False, colorbar=True, *args, **kwargs):

        if 'ax' not in kwargs:
            fig, local_ax = self._init_fig(subplot_kw={'projection': ccrs.EckertIV()})
        else:
            local_ax = kwargs.pop('ax')

        x_edge_arr = self.grid_dict.get(x_edge)
        y_edge_arr = self.grid_dict.get(y_edge)

        self._set_facecolor(local_ax)
        self._set_borderline(local_ax)
        self._add_outline(local_ax, x_edge=x_edge_arr, y_edge=y_edge_arr,  transform=self.transform_crs)
        p = self._add_pcolormesh(local_ax, x_edge=x_edge_arr, y_edge=y_edge_arr, transform=self.transform_crs, *args, **kwargs)

        if colorbar:
            self._add_colorbar(p)
        
        if contour:
            x_arr = self.grid_dict.get(x_edge[:3:1])
            y_arr = self.grid_dict.get(y_edge[:3:1])
            ## if cmap is specified, also do it in contour
            # if "cmap" in kwargs:
            #     p = self._add_contour(ax, x=x_arr, y=y_arr, transform=self.transform_crs, *args, **kwargs)
            # else:
            #     p = self._add_contour(ax, x=x_arr, y=y_arr, transform=self.transform_crs)
            p = self._add_contour(local_ax, x=x_arr, y=y_arr, transform=self.transform_crs, colors="k", linewidths=0.5, zorder=1)      

        return p

    def plot_polar(self, ax=None, hemisphere="South", x_edge="lon_edge", y_edge="lat_edge", contour=False, colorbar=True, *args, **kwargs):

        if not ax:
            match hemisphere:
                case "North":
                    fig, ax = self._init_fig(subplot_kw={'projection': ccrs.Orthographic(0, 90)})
                case "South":
                    fig, ax = self._init_fig(subplot_kw={'projection': ccrs.Orthographic(180, -90)})

        x_edge_arr = self.grid_dict.get(x_edge)
        y_edge_arr = self.grid_dict.get(y_edge)

        self._set_facecolor(ax)
        self._set_borderline(ax)
        p = self._add_pcolormesh(ax, x_edge=x_edge_arr, y_edge=y_edge_arr, transform=self.transform_crs, *args, **kwargs)
        self._add_outline(ax, x_edge=x_edge_arr, y_edge=y_edge_arr, transform=self.transform_crs)
        self._add_gridline(ax,
            transform=self.transform_crs,
            draw_labels=False,
            linewidth=0.5,
            color="gray",
            alpha=0.5,
            linestyle="-",
        )

        if colorbar:
            self._add_colorbar(p)

        
        if contour:
            x_arr = self.grid_dict.get(x_edge[:3:1])
            y_arr = self.grid_dict.get(y_edge[:3:1])
            p = self._add_contour(ax, x=x_arr, y=y_arr, transform=self.transform_crs)

            
        return p

    
    def _plot_cross_section(self, ax=None, x_edge="lat_edge", y_edge="zt_edge", contour=False, colorbar=True, *args, **kwargs):
        """
        Examples

        import matplotlib.pyplot as plt

        model = EcoModel("path_to_model")
        
        fig, axs=plt.subplots(nrows=1, ncols=3, figsize=(15, 3))
        
        basins = ['Atlantic', 'Pacific', 'Indian']
        
        for i in range(3):
            model.get_var('ocn_PO4').isel(time=-1).mask_basin(base='worjh2',basin=basins[i], subbasin='').mean(dim='lon').plot(ax=axs[i])
            axs[i].title.set_text(basins[i])
        """
        if not ax:
            fig, ax = self._init_fig(figsize=(5, 2.5))

        x_edge_arr = self.grid_dict.get(x_edge)
        y_edge_arr = self.grid_dict.get(y_edge)

        
        self._set_facecolor(ax)
        self._set_borderline(ax, geo=False)

        p = self._add_pcolormesh(ax, x_edge=x_edge_arr, y_edge=y_edge_arr, *args, **kwargs)
        self._add_outline(ax, x_edge=x_edge_arr, y_edge=y_edge_arr)
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel(x_edge[:3:1])
        ax.set_ylabel("Depth (km)")

        if colorbar:
            self._add_colorbar(p, location="vertical")
        
        if contour:
            ## e.g., lat_edge[:3:1] -> "lat"
            ##       zt_edge[:2:1] -> "zt"
            x_arr = self.grid_dict.get(x_edge[:3:1])
            y_arr = self.grid_dict.get(y_edge[:2:1])
            
            p = self._add_contour(ax, x=x_arr, y=y_arr, linewidths=0.6, colors="black", linestyles="solid")
            ax.clabel(p, p.levels[::1], colors=["black"], fontsize=8.5, inline=False)

        return p


    ## ------- Below is implementations -------------------------

    def _init_fig(self, *args, **kwargs):
        # self._init_style()
        return plt.subplots(dpi=120, *args, **kwargs)

    def _init_style(self):
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams['image.cmap'] = 'viridis'
        plt.rcParams['grid.linestyle'] = '--'
        plt.rcParams['grid.width'] = 0.3
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')

    def _add_pcolormesh(self, ax, x_edge, y_edge, *args, **kwargs):
        return ax.pcolormesh(
            x_edge,
            y_edge,
            self.array,
            *args,
            **kwargs
        )

    def _add_contour(self, ax, x, y, *args, **kwargs):
        cs = ax.contour(x, y, self.array, *args, **kwargs)
        # label every three levels
        ax.clabel(cs, cs.levels[::3], colors=["black"], fontsize=8, inline=False)
        return cs

    def _add_gridline(self, ax, *args, **kwargs):
        ax.gridlines(*args, **kwargs)

    def _set_borderline(self, ax, geo=True, width=1):
        if geo:
            ax.spines["geo"].set_edgecolor("black")
            ax.spines["geo"].set_linewidth(width)
        else:
            for direction in ['top','bottom','left','right']:
                ax.spines[direction].set_linewidth(width)
                ax.spines[direction].set_edgecolor("black")
                # vertical layer order
                ax.spines[direction].set_zorder(0)

    def _set_facecolor(self, ax):
        ax.patch.set_color("silver")

    def _add_outline(self, ax, x_edge, y_edge, *args, **kwargs):
        outline_color = "black"
        outline_width = 1
        mask_array = np.where(~np.isnan(self.array), 1, 0)

        Ny_edge = len(y_edge)
        Nx_edge = len(x_edge)
        # dimension of array
        Ny_array = Ny_edge - 1
        Nx_array = Nx_edge - 1
        # index of array
        Ny_index = Ny_array - 1
        Nx_index = Nx_array - 1

        # i stands for yitude index, j stands for xgitude index
        # e.g., mask_array[0, ] is Antarctic ice cap
        # i,j are for mask array and will be converted to x/y edge lines
        for i in range(Ny_array):
            for j in range(Nx_array):
                # compare with the right grid, and plot vertical line if different
                if j < Nx_index and mask_array[i, j] != mask_array[i, j + 1]:
                    ax.vlines(
                        x_edge[j + 1],
                        y_edge[i],
                        y_edge[i + 1],
                        color=outline_color,
                        linewidth=outline_width,
                        *args, **kwargs
                    )

                # connect the circular xgitude axis
                if j == Nx_index and mask_array[i, j] != mask_array[i, 0]:
                    ax.vlines(
                        x_edge[j + 1],
                        y_edge[i],
                        y_edge[i + 1],
                        color=outline_color,
                        linewidth=outline_width,
                        *args, **kwargs
                    )

                # compare with the above grid, and plot horizontal line if different
                if i < Ny_index and mask_array[i, j] != mask_array[i + 1, j]:
                    ax.hlines(
                        y_edge[i + 1],
                        x_edge[j],
                        x_edge[j + 1],
                        colors=outline_color,
                        linewidth=outline_width,
                        *args, **kwargs
                    )

    def _add_colorbar(self, mappable_object, location="horizontal"):
        cbar = plt.colorbar(mappable_object, fraction=0.05, pad=0.04, orientation=location)
        cbar.ax.tick_params(color="k", direction="in")
        cbar.outline.set_edgecolor('black')
        cbar.minorticks_on()
        ## set colorbar label
        cbar.set_label(f"{self.array.long_name}\n{self.array.units}", size=10, labelpad=10)
        
    def plot_quiver(x,y):
        pass

    
class TaylorDiagram(object):
    """
    Taylor diagram.
    Plot model standard deviation and correlation to reference data in a single-quadrant polar plot,
    with r=stddev and theta=arccos(correlation).

    modified from Yannick Copin, https://gist.github.com/ycopin/3342888
    reference: https://matplotlib.org/stable/gallery/axisartist/demo_floating_axes.html
    """

    def __init__(
        self,
        fig=None,
        figscale=1,
        subplot=111,
        xmax=None,
        tmax=np.pi / 2,
        ylabel="Standard Deviation",
        rotation=None,
    ):

        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.

        Parameters:

        * fig: input Figure or None
        * subplot: subplot definition
        * xmax: the length of radius, xmax can be 1.5* reference std
        """

        # --------------- tickers --------------------------
        # Correlation labels (if half round)
        cor_label = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1])

        # add the negative ticks if more than half round
        excess_theta = tmax - np.pi / 2
        if excess_theta > 0:
            cor_label = np.concatenate((-cor_label[:0:-1], cor_label))

        # convert to radian
        rad = np.arccos(cor_label)
        # tick location
        gl = gf.FixedLocator(rad)
        # tick formatting: bind radian and correlation coefficient
        tf = gf.DictFormatter(dict(zip(rad, map(str, cor_label))))

        # --------------- coordinate -----------------------
        # Standard deviation axis extent (in units of reference stddev)
        # xmin must be 0, which is the centre of round

        self.xmin = 0
        self.xmax = xmax
        self.tmax = tmax

        # ------- curvilinear coordinate definition -------
        # use built-in polar transformation (i.e., from theta and r to x and y)
        tr = PolarAxes.PolarTransform()
        ghelper = fa.GridHelperCurveLinear(
            tr,
            extremes=(0, self.tmax, self.xmin, self.xmax),
            grid_locator1=gl,
            tick_formatter1=tf,
        )

        # ------- create floating axis -------
        if fig is None:
            fig_height = 4.5 * figscale
            fig_width = fig_height * (1 + np.sin(excess_theta))
            fig = plt.figure(figsize=(fig_width, fig_height), dpi=100)

        ax = fa.FloatingSubplot(fig, subplot, grid_helper=ghelper)
        fig.add_subplot(ax)

        # Adjust axes
        # Angle axis
        ax.axis["top"].label.set_text("Correlation")
        ax.axis["top"].toggle(ticklabels=True, label=True)
        # inverse the direction
        ax.axis["top"].set_axis_direction("bottom")
        ax.axis["top"].major_ticklabels.set_axis_direction("top")
        ax.axis["top"].label.set_axis_direction("top")

        # X axis
        ax.axis["left"].set_axis_direction("bottom")

        # Y axis direction & label
        ax.axis["right"].toggle(all=True)
        ax.axis["right"].label.set_text(ylabel)
        ax.axis["right"].set_axis_direction("top")
        # ticklabel direction
        ax.axis["right"].major_ticklabels.set_axis_direction("left")

        ax.axis["bottom"].set_visible(False)

        # ------- Set instance attribute ----------
        self.fig = fig
        # Graphical axes
        self._ax = ax
        # grid line
        self._ax.grid(True, zorder=0, linestyle="--")
        # aspect ratio
        self._ax.set_aspect(1)
        # A parasite axes for further plotting data
        self.ax = ax.get_aux_axes(tr)
        # Collect sample points for latter use (e.g. legend)
        self.samplePoints = []

    def add_ref(self, refstd, reflabel="Observation", linestyle="-", color="k"):
        """add a reference point"""
        self.refstd = refstd
        # Add reference point
        # slightly higher than 0 so star can be fully seen
        l = self.ax.plot(0.01, self.refstd, "k*", ls="", ms=10)
        # xy for the point, xytext for the text (the coordinates are
        # defined in xycoords and textcoords, respectively)
        self.ax.annotate(
            reflabel,
            xy=(0.01, self.refstd),
            xycoords="data",
            xytext=(-25, -30),
            textcoords="offset points",
        )
        # add stddev contour
        t = np.linspace(0, self.tmax)
        r = np.zeros_like(t) + self.refstd
        self.ax.plot(t, r, linestyle=linestyle, color=color)
        self.samplePoints.append(l)

    def add_scatter(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """

        l = self.ax.scatter(
            np.arccos(corrcoef), stddev, *args, **kwargs
        )  # (theta, radius)
        self.samplePoints.append(l)

        return l

    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """

        rs, ts = np.meshgrid(
            np.linspace(self.xmin, self.xmax), np.linspace(0, self.tmax)
        )
        # Compute centered RMS difference
        crmse = np.sqrt(self.refstd**2 + rs**2 - 2 * self.refstd * rs * np.cos(ts))
        contours = self.ax.contour(ts, rs, crmse, levels, linestyles="--", **kwargs)
        self.ax.clabel(contours, contours.levels[::1], inline=False)

        return contours

    def add_legend(self, *args, **kwargs):
        return self.ax.legend(*args, **kwargs)

    def add_annotation(self, *args, **kwargs):
        return self.ax.annotation(*args, **kwargs)

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)
