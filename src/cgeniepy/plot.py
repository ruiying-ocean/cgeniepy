import pathlib

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
from cartopy.feature import LAND
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf

from .data import efficient_log
from .grid import GENIE_lat, GENIE_lon, GENIE_depth

# [] decorator for boundary line in pcolormesh object


def plot_genie(
    ax,
    data,
    log=False,
    grid_line=False,
    continent_outline=True,
    contour_layer=False,
    *args,
    **kwargs,
):
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
    ax.spines["geo"].set_edgecolor("black")
    ax.spines["geo"].set_linewidth(1)

    # -------------------Facecolor-----------------------
    ax.patch.set_color("silver")

    # -------------------Plot-----------------------
    if log:
        data = efficient_log(data)

    lon_edge = GENIE_lon(edge=True)
    lat_edge = GENIE_lat(edge=True)

    # cartopy transform seems to help reassign the GENIE longitude to normal
    p = ax.pcolormesh(
        lon_edge,
        lat_edge,
        data,
        transform=data_crs,
        shading="flat",
        *args,
        **kwargs,
    )

    if contour_layer:
        lat = GENIE_lat(edge=False)
        lon = GENIE_lon(edge=False)
        cs = ax.contour(lon, lat, data, transform=ccrs.PlateCarree(), cmap=cmap)
        # label every three levels
        ax.clabel(cs, cs.levels[::3], colors=["black"], fontsize=8, inline=False)

    # -------------------Grid lines-----------------------
    if grid_line:
        ax.gridlines(
            crs=data_crs,
            draw_labels=False,
            linewidth=0.5,
            color="gray",
            alpha=0.5,
            linestyle="-",
        )

    # -------------------Continent lines-----------------------
    if continent_outline:
        outline_color = "black"
        outline_width = 1
        mask_array = np.where(~np.isnan(data), 1, 0)

        Nlat_edge = len(lat_edge)
        Nlon_edge = len(lon_edge)
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
                if j < Nlon_index and mask_array[i, j] != mask_array[i, j + 1]:
                    ax.vlines(
                        lon_edge[j + 1],
                        lat_edge[i],
                        lat_edge[i + 1],
                        color=outline_color,
                        linewidth=outline_width,
                        transform=data_crs,
                    )

                # connect the circular longitude axis
                if j == Nlon_index and mask_array[i, j] != mask_array[i, 0]:
                    ax.vlines(
                        lon_edge[j + 1],
                        lat_edge[i],
                        lat_edge[i + 1],
                        color=outline_color,
                        linewidth=outline_width,
                        transform=data_crs,
                    )

                # compare with the above grid, and plot horizontal line if different
                if i < Nlat_index and mask_array[i, j] != mask_array[i + 1, j]:
                    ax.hlines(
                        lat_edge[i + 1],
                        lon_edge[j],
                        lon_edge[j + 1],
                        colors=outline_color,
                        linewidth=outline_width,
                        transform=data_crs,
                    )

    return p


def scatter_map(
    ax,
    df,
    var,
    x="Longitude",
    y="Latitude",
    log=False,
    *args,
    **kwargs,
):

    ax.set_global()
    ax.coastlines()
    ax.add_feature(
        LAND, zorder=0, facecolor="#B1B2B4", edgecolor="white"
    )  # zorder is drawing sequence

    if "Latitude" not in df.columns or "Longitude" not in df.columns:
        raise ValueError("Input data lack Latitude/Longitude column")

    if log:
        df[var] = efficient_log(df[var])

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


def genie_cmap(cmap_name, N=256, reverse=False):
    """
    :param cmap_name: Zissou1, FantasticFox, Rushmore, Darjeeling, ODV
    :return: colormap
    """
    file_name = f"data/{cmap_name}.txt"
    file_path = pathlib.Path(__file__).parent.parent / file_name

    colors = pd.read_csv(file_path, header=None).values.tolist()
    colors = [colors[i][0] for i in range(len(colors))]

    c = ListedColormap(colors, N=N)

    if reverse:
        return c.reversed()
    return c


def cbar_wrapper(func, *args, **kwarags):
    def wrappered_func():
        p = func(*args, **kwarags)
        cbar = plt.colorbar(p, fraction=0.05, pad=0.04, orientation="vertical")
        cbar.ax.tick_params(color="k", direction="in")

    return wrappered_func


class GeniePlottable(object):
    def plot_line(self, dim, *args, **kwargs):
        """plot 1D data, e.g., zonal_average"""

        if self.dim() > 1:
            raise ValueError("Data should be of 1 Dimension")
        if dim == "lon":
            x = GENIE_lat()
        if dim == "lat":
            x = GENIE_lon()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.minorticks_on()

        p = ax.plot(x, self.array, "k", *args, **kwargs)

        return p

    def plot_map(self, ax=None, cbar=True, *args, **kwargs):

        """plot lat-lon 2D array"""

        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = ["Arial"]

        if not ax:
            fig = plt.figure(dpi=75)
            ax = fig.add_subplot(111, projection=ccrs.EckertIV())

        if np.any(self.array < 0) and np.any(self.array > 0):
            # x = max(abs(self.nanmin()), self.nanmax())
            # vmax = x/2
            # vmin = vmax * -1
            p = plot_genie(ax=ax, data=self.array, cmap="RdBu_r", *args, **kwargs)
        else:
            p = plot_genie(ax=ax, data=self.array, *args, **kwargs)

        if cbar:
            cax = fig.add_axes([0.15, 0.1, 0.73, 0.07])  # xmin, ymin, dx, dy
            cbar = fig.colorbar(p, cax=cax, orientation="horizontal", pad=0.04)
            cbar.minorticks_on()
            if hasattr(self, "unit"):
                cbar.set_label(self.unit, size=12)

        return p

    def cross_section(self, ax=None, cmap="YlGnBu_r", *args, **kwargs):
        """plot 3D array, tracers and stream function"""

        data = self.array

        # visualise
        if not ax:
            fig, ax = plt.subplots(figsize=(8, 4))

        lat_edge = GENIE_lat(edge=True)
        z_edge = GENIE_depth(edge=True) / 1000
        lat = GENIE_lat(edge=False)
        z = GENIE_depth(edge=False) / 1000

        # pcolormesh
        p = ax.pcolormesh(lat_edge, z_edge, data, cmap=cmap, *args, **kwargs)

        # contourf
        # contourf_max = np.max(moc)//10 * 10
        # controuf_N = int(contourf_max*2/10 + 1)
        # contourf_level = np.linspace(contourf_max*-1, contourf_max, controuf_N)
        # cs = ax.contourf(lat, z, data, corner_mask=False, levels=contourf_level, cmap=contourf_cmap)

        # contour
        cs = ax.contour(
            lat, z, data, linewidths=0.6, colors="black", linestyles="solid"
        )
        ax.clabel(cs, cs.levels[::1], colors=["black"], fontsize=8.5, inline=False)

        # colorbar, ticks & labels
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = plt.colorbar(p, fraction=0.05, pad=0.04, orientation="vertical")
        cbar.ax.tick_params(color="k", direction="in")
        # cbar.set_label("Stream function (Sv)")

        plt.rcParams["font.family"] = "Arial"
        ax.patch.set_color("#999DA0")
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel(r"Latitude ($\degree$N)", fontsize=13)
        ax.set_ylabel("Depth (km)", fontsize=12)
        ax.minorticks_on()
        ax.tick_params("both", length=4, width=1, which="major")

        return p


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
