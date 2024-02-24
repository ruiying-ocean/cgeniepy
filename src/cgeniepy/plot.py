import pathlib
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import ListedColormap, to_rgb as hex_to_rgb
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf

from .utils import efficient_log
from .chem import Chemistry

class ArrayVis:

    transform_crs = ccrs.PlateCarree()  # do not change

    def __init__(self, array):

        self.array = array

        if hasattr(self.array, "units"):
            self.units = Chemistry().format_unit(self.array.units)
        else:
            self.units = ""

        if hasattr(self.array, "long_name"):
            self.long_name = self.array.long_name
        else:
            self.long_name = ""

        ## aesthetic parameters
        self.aes_dict = {
            ## x,ylabel
            "general_kwargs": {"cmap": "viridis", "font": "Helvetica", "fontsize": 10},
            "facecolor_kwargs": {"c": "silver"},
            "borderline_kwargs": {"c": "black", "linewidth": 0.5},
            "outline_kwargs": {"colors": "black", "linewidth": 0.5},
            "gridline_kwargs": {"colors": "gray", "linewidth": 0.5},
            "pcolormesh_kwargs": {"shading": "auto", "cmap": plt.get_cmap("viridis")},
            "contour_kwargs": {
                "linewidths": 0.6,
                "colors": "black",
                "linestyles": "solid",
                "zorder": 10,
            },
            "contour_label_kwargs": {
                "colors": ["black"],
                "fontsize": 8,
                "inline": False,
            },
            "contourf_kwargs": {"levels": 20},
            "colorbar_label_kwargs": {
                "label": f"{self.long_name}\n({self.units})",
                "size": 10,
                "labelpad": 10,
            },
            "colorbar_kwargs": {"fraction": 0.046, "pad": 0.04},
        }

    def plot(self, *args, **kwargs):
        """
        visualise the data based on the dimension of the array

        *args and **kwargs are passed to plotting functions
        in practical, turn on/off the plotting elements
        e.g., plot(x, y, pcolormesh=False, contour=True)
        each plotting function will seek for the corresponding kwargs
        in self.aes_dict
        """

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
        if "ax" not in kwargs:
            fig, local_ax = self._init_fig()
        else:
            local_ax = kwargs.pop("ax")

        ## get the only dimension as x
        dim = self.array[self.array.dims[0]]

        # local_ax.set_xlabel(dim.name)
        # local_ax.set_ylabel(f"{self.array.long_name} ({self.array.units})")
        p = local_ax.plot(dim, self.array, *args, **kwargs)

        return p

    def _plot_2d(self, *arg, **kwargs):
        "plot 2D data, e.g., map, trasect"
        ## if lon, lat then plot map
        ## if lon, zt then plot transect
        ## if lat, zt then plot transect
        dims = self.array.dims
        if "lon" in dims and "lat" in dims:
            return self._plot_map(*arg, **kwargs)
        elif "zt" in dims and "lon" in dims:
            return self._plot_transect(*arg, **kwargs)
        elif "zt" in dims and "lat" in dims:
            return self._plot_transect(*arg, **kwargs)
        else:
            raise ValueError(f"{dims} not supported")

    def _plot_3d(self, *args, **kwargs):
        "plot 3D data = plot mutiple 2D plots"
        print("3D plot not supported yet")
        pass

    def _plot_map(
        self,
        pcolormesh=True,
        contour=False,
        colorbar=False,
        contourf=False,
        outline=False,
        facecolor=True,
        borderline=True,
        gridline=False,
        *args,
        **kwargs,
    ):

        x_name = self.array.dims[1]  ## lon
        y_name = self.array.dims[0]  ## lat

        x_arr = self.array[x_name]
        y_arr = self.array[y_name]

        x_edge = np.linspace(-260, 100, x_arr.size + 1)
        y_edge = np.rad2deg(np.arcsin(np.linspace(-1, 1, y_arr.size + 1)))

        if "ax" not in kwargs:
            fig, local_ax = self._init_fig(subplot_kw={"projection": ccrs.EckertIV()})
        else:
            local_ax = kwargs.pop("ax")

        if facecolor:
            self._set_facecolor(local_ax, **self.aes_dict["facecolor_kwargs"])

        if borderline:
            self._set_borderline(
                local_ax, geo=True, **self.aes_dict["borderline_kwargs"]
            )

        if outline:
            ## outline uses edge coordinates
            ## need to be transformed to PlateCarree
            self._add_outline(
                local_ax,
                x=x_edge,
                y=y_edge,
                transform=self.transform_crs,
                **self.aes_dict["outline_kwargs"],
            )

        if gridline:
            self._add_gridline(
                local_ax,
                transform=self.transform_crs,
                **self.aes_dict["gridline_kwargs"],
            )

        if pcolormesh:
            ## pcolormesh uses edge coordinates
            ## need to be transformed to PlateCarree
            p_pcolormesh = self._add_pcolormesh(
                local_ax,
                x=x_edge,
                y=y_edge,
                transform=self.transform_crs,
                *args,
                **self.aes_dict["pcolormesh_kwargs"],
            )
            if colorbar:
                cbar = self._add_colorbar(p_pcolormesh, orientation="horizontal")
                self._add_colorbar_label(cbar, **self.aes_dict["colorbar_label_kwargs"])

        if contour:
            ## contour uses center coordinates
            ## need to be transformed to PlateCarree
            p_contour = self._add_contour(
                local_ax,
                x_arr,
                y_arr,
                transform=self.transform_crs,
                **self.aes_dict["contour_kwargs"],
            )
            self._add_contour_label(
                local_ax, p_contour, **self.aes_dict["contour_label_kwargs"]
            )
            ## contour will not be used to plot colorbar because it's set to black

        if contourf:
            p_contourf = self._add_contourf(
                local_ax,
                x_arr,
                y_arr,
                transform=self.transform_crs,
                **self.aes_dict["contourf_kwargs"],
            )

        if pcolormesh:
            return p_pcolormesh
        elif contour:
            return p_contour
        elif contourf:
            return p_contourf
        else:
            return local_ax

    def _plot_transect(
        self,
        x="lat_edge",
        y="zt_edge",
        pcolormesh=True,
        contour=False,
        colorbar=True,
        contourf=False,
        outline=False,
        facecolor=True,
        borderline=True,
        *args,
        **kwargs,
    ):
        """
        Examples

        import matplotlib.pyplot as plt

        model = EcoModel("path_to_model")

        fig, axs=plt.subplots(nrows=1, ncols=3, figsize=(15, 3))

        basins = ['Atlantic', 'Pacific', 'Indian']

        for i in range(3):
            model.get_var('ocn_PO4').isel(time=-1).mask_basin(base='worjh2',basin=basins[i], subbasin='').mean(dim='lon').plot(ax=axs[i])
            axs[i].title.set_text(basins[i])


        available kwargs:
        outline_kwargs = {'outline_color': 'red', 'outline_width': .5}
        """

        x_name = self.array.dims[1]  ## lat
        y_name = self.array.dims[0]  ## zt

        x_arr = self.array[x_name]
        y_arr = self.array[y_name]

        if "ax" not in kwargs:
            fig, local_ax = self._init_fig(figsize=(6, 3))
        else:
            local_ax = kwargs.pop("ax")

        if facecolor:
            self._set_facecolor(local_ax, **self.aes_dict["facecolor_kwargs"])

        if borderline:
            self._set_borderline(
                local_ax, geo=False, **self.aes_dict["borderline_kwargs"]
            )

        if outline:
            ## outline uses edge coordinates
            self._add_outline(
                local_ax, x=x_arr, y=y_arr, **self.aes_dict["outline_kwargs"]
            )

        if pcolormesh:
            ## pcolormesh uses edge coordinates
            p_pcolormesh = self._add_pcolormesh(
                local_ax, x=x_arr, y=y_arr, *args, **self.aes_dict["pcolormesh_kwargs"]
            )
            if colorbar:
                cbar = self._add_colorbar(p_pcolormesh, orientation="vertical")
                self._add_colorbar_label(cbar, **self.aes_dict["colorbar_label_kwargs"])

        if contour:
            ## contour uses center coordinates
            ## need to be transformed to PlateCarree
            p_contour = self._add_contour(
                local_ax, x_arr, y_arr, **self.aes_dict["contour_kwargs"]
            )
            self._add_contour_label(
                local_ax, p_contour, **self.aes_dict["contour_label_kwargs"]
            )
            ## contour will not be used to plot colorbar because it's set to black

        if contourf:
            p_contourf = self._add_contourf(local_ax, x_arr, y_arr)

        ## reverse y axis
        local_ax.set_ylim(local_ax.get_ylim()[::-1])
        local_ax.set_ylabel(
            "Depth (m)",
            fontsize=self.aes_dict["general_kwargs"]["fontsize"],
            font=self.aes_dict["general_kwargs"]["font"],
        )
        local_ax.set_xlabel(
            x.split("_")[0],
            fontsize=self.aes_dict["general_kwargs"]["fontsize"],
            font=self.aes_dict["general_kwargs"]["font"],
        )

        ## x/y tick label
        local_ax.tick_params(
            axis="both",
            which="major",
            labelsize=self.aes_dict["general_kwargs"]["fontsize"],
            labelfontfamily=self.aes_dict["general_kwargs"]["font"],
        )

        if pcolormesh:
            return p_pcolormesh
        elif contour:
            return p_contour
        elif contourf:
            return p_contourf
        else:
            return local_ax


    ## ------- Below is implementations -------------------------

    def _init_fig(self, *args, **kwargs):
        return plt.subplots(dpi=120, *args, **kwargs)

    def _add_pcolormesh(self, ax, x, y, *args, **kwargs):
        return ax.pcolormesh(x, y, self.array, *args, **kwargs)

    def _add_contour(self, ax, x, y, *args, **kwargs):
        return ax.contour(x, y, self.array, *args, **kwargs)

    def _add_contourf(self, ax, x, y, *args, **kwargs):
        return ax.contourf(x, y, self.array, *args, **kwargs)

    def _add_contour_label(self, ax, cs, *args, **kwargs):
        ax.clabel(cs, cs.levels[::2], *args, **kwargs)

    def _add_gridline(self, ax, *args, **kwargs):
        ax.gridlines(*args, **kwargs)

    def _set_borderline(self, ax, geo=True, **kwargs):
        if geo:
            ax.spines["geo"].set_edgecolor(kwargs["c"])
            ax.spines["geo"].set_linewidth(kwargs["linewidth"])
        else:
            for direction in ["top", "bottom", "left", "right"]:
                ax.spines[direction].set_edgecolor(kwargs["c"])
                ax.spines[direction].set_linewidth(kwargs["linewidth"])
                # vertical layer order
                ax.spines[direction].set_zorder(0)

    def _set_facecolor(self, ax, **kwargs):
        ax.patch.set_color(**kwargs)

    def _add_outline(self, ax, x, y, **kwargs):
        """
        draw outlines for the NA grids
        """
        mask_array = np.where(~np.isnan(self.array), 1, 0)

        Ny_edge = len(y)
        Nx_edge = len(x)
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
                    ax.vlines(x[j + 1], y[i], y[i + 1], **kwargs)

                # connect the circular xgitude axis
                if j == Nx_index and mask_array[i, j] != mask_array[i, 0]:
                    ax.vlines(x[j + 1], y[i], y[i + 1], **kwargs)

                # compare with the above grid, and plot horizontal line if different
                if i < Ny_index and mask_array[i, j] != mask_array[i + 1, j]:
                    ax.hlines(y[i + 1], x[j], x[j + 1], **kwargs)

    def _add_colorbar(self, mappable_object, *args, **kwargs):
        cbar = plt.colorbar(mappable_object, *args, **kwargs)
        return cbar

    def _add_colorbar_label(self, cbar, *args, **kwargs):
        ## set colorbar label
        cbar.set_label(*args, **kwargs)
        cbar.ax.tick_params(
            axis="both",
            which="major",
            labelsize=self.aes_dict["general_kwargs"]["fontsize"],
            labelfontfamily=self.aes_dict["general_kwargs"]["font"],
        )
        cbar.outline.set_edgecolor("black")
        cbar.outline.set_linewidth(0.5)


class ScatterVis:
    "Visualisation object based on dataframe"

    def __init__(self, df):
        self.df = df

    def _init_fig(self, *args, **kwargs):
        return plt.subplots(dpi=120, *args, **kwargs)

    def plot_map(
        self,
        ax=None,
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
            fig, ax = self._init_fig(subplot_kw={"projection": ccrs.EckertIV()})

        if land_mask:
            ax.set_global()
            # plot land and coastline, zorder is the drawing order, smaller -> backer layer
            ax.add_feature(
                cfeature.LAND.with_scale("110m"), zorder=2, facecolor="#B1B2B4"
            )
            ax.add_feature(cfeature.COASTLINE.with_scale("110m"), zorder=3)

        if log:
            self.df[self.var] = efficient_log(self.df[self.var])

        if self.df[self.var].dtype != float:
            self.df[self.var] = self.df[self.var].astype(float)

        p = ax.scatter(
            x=self.df[self.lon],
            y=self.df[self.lat],
            c=self.df[self.var],
            transform=ccrs.PlateCarree(),
            *args,
            **kwargs,
        )

        return p

    def plot_transect(
        self,
        ax=None,
        bathy_lon=None,
        *args,
        **kwargs,
    ):
        if not ax:
            fig, ax = self._init_fig()

        if self.df[self.var].dtype != float:
            self.df[self.var] = self.df[self.var].astype(float)
            
        p = ax.scatter(
            x=self.df[self.lat],
            y=self.df[self.depth],
            c=self.df[self.var],
            *args,
            **kwargs,
        )        

        if bathy_lon:
            data_dir = pathlib.Path(__file__).parent.parent
            bathy = xr.open_dataset(data_dir / "data/GEBCO2002_bathy.nc")

            ## get the lat range
            min_lat = self.df[self.lat].min()
            max_lat = self.df[self.lat].max()
            min_depth = self.df[self.depth].min()

            up = bathy.z.sel(lon=slice(*bathy_lon)).mean(dim="lon")
            # up = up.sel(lat=slice(min_lat, max_lat))
            bottom = np.ones(len(up)) * -5500

            ax.fill_between(up.lat, up, bottom, color="black")

        return p


class CommunityPalette:

    def get_palette(self, cmap_name, N=256, reverse=False, alpha=None):
        """
        community-driven colormaps with multiple sources

        :param cmap_name: colormap name, can be found in avail_palette()
        :type cmap_name: str

        :returns: colormap

        XML data: https://sciviscolor.org/colormaps/
        txt data: from original packages
        """

        if cmap_name not in self.avail_palettes():
            raise ValueError(
                f"{cmap_name} not found, accepted values are {self.avail_palette()}"
            )

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
            for point in root.findall(".//Point"):
                r = float(point.get("r"))
                g = float(point.get("g"))
                b = float(point.get("b"))
                colors.append((r, g, b))
                positions.append(float(point.get("x")))

            c = ListedColormap(colors, name=cmap_name)

        # Optionally add transparency
        if alpha is not None:
            if file_ext == ".txt":
                rgb_array = np.array([hex_to_rgb(i) for i in colors])
                rgb_array[:, -1] = alpha
                c = ListedColormap(rgb_array, N=N)
            elif file_ext == ".xml":
                # Assuming alpha is to be applied uniformly to all colors in the XML-based colormap
                rgba_colors = [
                    (color[0], color[1], color[2], alpha) for color in colors
                ]
                c = ListedColormap(rgba_colors, name=cmap_name)

        if reverse and c is not None:
            return c.reversed()

        ## if _r exsit in the string, reverse the colormap
        if "_r" in cmap_name and c is not None:
            return c.reversed()

        if c is None:
            raise ValueError("Colormap could not be created")

        return c

    def avail_palettes(self):
        """return a list of colormap names"""
        data_dir = pathlib.Path(__file__).parent.parent

        return [
            f.stem
            for f in data_dir.glob("data/colormaps/*")
            if f.suffix in [".txt", ".xml"]
        ]
