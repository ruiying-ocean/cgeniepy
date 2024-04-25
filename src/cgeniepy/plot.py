import xml.etree.ElementTree as ET
from importlib.resources import files
import itertools

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb as hex_to_rgb
import matplotlib.patheffects as pe
from matplotlib.patheffects import Stroke, Normal
import cartopy.mpl.geoaxes
from cgeniepy.grid import GridOperation
from .utils import efficient_log


class GriddedDataVis:

    """A class to visualise the GriddedData object"""

    transform_crs = ccrs.PlateCarree()  # do not change

    def __init__(self, gd):
        """
        initialise the GriddedDataVis object

        :param gd: GriddedData object
        """
        self.data = gd.data
        self.attrs = gd.attrs
        
        ## aesthetic parameters
        self.aes_dict = {
            ## x,ylabel
            "general_kwargs": {"cmap": "viridis", "font": "Helvetica", "fontsize": 10},
            "facecolor_kwargs": {"c": "silver"},
            "borderline_kwargs": {"c": "black", "linewidth": 0.5},
            "outline_kwargs": {"colors": "black", "linewidth": 0.5},
            "gridline_kwargs": {"color": "black", "linewidth": 0.5, "linestyle": "dashed", "draw_labels": False},
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
                "label": f"{self.attrs['long_name']}\n({self.attrs['units']})",
                "size": 10,
                "labelpad": 10,
            },
            "colorbar_kwargs": {"fraction": 0.046, "pad": 0.03},
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
        plt.rcParams['font.family'] = 'sans-serif'

        if self.data.ndim == 1:
            return self._plot_1d(*args, **kwargs)
        elif self.data.ndim == 2:
            return self._plot_2d(*args, **kwargs)
        elif self.data.ndim == 3:
            return self._plot_3d(*args, **kwargs)
        else:
            raise ValueError(f"{self.data.ndim} dimensions not supported")

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
        dim = self.data[self.data.dims[0]]

        # local_ax.set_xlabel(dim.name)
        # local_ax.set_ylabel(f"{self.array.long_name} ({self.array.units})")
        p = local_ax.plot(dim, self.data, *args, **kwargs)

        return p

    def _plot_2d(self, *arg, **kwargs):
        "plot 2D data, e.g., map, trasect"
        ## if lon, lat then plot map
        ## if zt or depth then plot transec   

        has_lat, has_lon, has_depth, has_time = GridOperation().check_dimension(self.data.dims)
        if has_lat and has_lon:
            return self._plot_map(*arg, **kwargs)
        elif has_depth:
            return self._plot_transect(*arg, **kwargs)
        else:
            raise ValueError("2D plotting not supported")

    def _plot_3d(self, *args, **kwargs):
        "plot 3D data = plot mutiple 2D plots"
        has_lat, has_lon, has_depth, has_time = GridOperation().check_dimension(self.data.dims)
        
        if has_time:
            time_order = GridOperation().dim_order(self.data.dims)[0]
            time_name = self.data.dims[time_order]
            time_arr = self.data[time_name]

            # Determine the optimal number of columns and rows for subplots
            num_plots = len(time_arr)
            ncols = min(num_plots, 4)  # Maximum of 4 columns
            nrows = (num_plots + ncols - 1) // ncols  # Calculate the required number of rows

            # Create a figure and axis objects
            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6 * ncols, 3 * nrows), squeeze=False)

            # Flatten the axes array for easier iteration
            axs = axs.flatten()

            # Plot the data for each time step
            for i, ax in enumerate(axs):
                if i < num_plots:
                    im = self.data.isel(time=i).plot(ax=ax,add_colorbar=False)

                # Hide unused axes
                if i >= num_plots:
                    ax.set_visible(False)

            # Add a common colorbar to the figure
            plt.colorbar(im, ax=axs.tolist(), orientation='horizontal', label=f"{self.attrs['long_name']} ({self.attrs['units']})", fraction=0.046, pad=0.05)

            # Adjust spacing between subplots
            fig.subplots_adjust(hspace=0.5, wspace=0.5)
        else:
            raise ValueError("Not support 3D plot iterating over other dimension than time")
    

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
        cfeature=None,
        zebra_frame=False,
        *args,
        **kwargs,
    ):

        dim_order = GridOperation().dim_order(self.data.dims)
        lat_order = dim_order[0] ## in the case of 2D, lat is the first dimension
        lon_order = dim_order[1] ## in the case of 2D, lon is the second dimension

        x_name = self.data.dims[lon_order]  ## lon
        y_name = self.data.dims[lat_order]  ## lat

        x_arr = self.data[x_name]
        y_arr = self.data[y_name]

        x_min = x_arr.min()
        x_max = x_arr.max()
        x_res = x_arr[1] - x_arr[0]
        x_edge = np.linspace(x_min-x_res/2, x_max+x_res/2, x_arr.size + 1)

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

        if cfeature:
            self._add_cfeature(local_ax, cfeature)
        
        if zebra_frame:
            self._add_zebra_frame(local_ax)

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

        dim_order = GridOperation().dim_order(self.data.dims)        
        zt_order = dim_order[0] ## in the case of 2D, zt is the first dimension
        lat_order = dim_order[1] ## in the case of 2D, lat is the second dimension
        
        x_name = self.data.dims[lat_order]  ## lat
        y_name = self.data.dims[zt_order]  ## zt

        x_arr = self.data[x_name]
        y_arr = self.data[y_name]

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
            y_name,
            fontsize=self.aes_dict["general_kwargs"]["fontsize"],
            font=self.aes_dict["general_kwargs"]["font"],
        )
        local_ax.set_xlabel(
            x_name,
            fontsize=self.aes_dict["general_kwargs"]["fontsize"],
            font=self.aes_dict["general_kwargs"]["font"],
        )

        ## x/y tick label
        local_ax.tick_params(
            axis="both",
            which="major",
            direction="in",
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
        return ax.pcolormesh(x, y, self.data, *args, **kwargs)

    def _add_contour(self, ax, x, y, *args, **kwargs):
        return ax.contour(x, y, self.data, *args, **kwargs)

    def _add_contourf(self, ax, x, y, *args, **kwargs):
        return ax.contourf(x, y, self.data, *args, **kwargs)

    def _add_contour_label(self, ax, cs, *args, **kwargs):
        ax.clabel(cs, cs.levels[::2], *args, **kwargs)

    def _add_gridline(self, ax, *args, **kwargs):
        gl = ax.gridlines(*args, **kwargs)     

        if ax.projection == ccrs.PlateCarree():
            gl.xlines = False  # removing gridlines
            gl.ylines = False            
            ## add ticks
            xticks = np.linspace(-180, 180, 7)
            yticks = np.linspace(-90, 90, 7)
            ax.set_xticks(xticks, crs=ccrs.PlateCarree())
            ax.set_yticks(yticks, crs=ccrs.PlateCarree())
            ax.set_yticklabels("")
            ax.set_xticklabels("")
            bbox_pixels = ax.get_window_extent().get_points()
            width_pixels = bbox_pixels[1, 0] - bbox_pixels[0, 0]

            # Set the tick length as a fraction of the width
            tick_length = 0.015 * width_pixels            
            ax.tick_params(axis='x', length=tick_length)
            ax.tick_params(axis='y', length=tick_length)

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
        mask_array = np.where(~np.isnan(self.data), 1, 0)

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

    def _add_cfeature(self, ax, feature_name):
        """add natural earth mask"""

        feature = cfeature.NaturalEarthFeature('physical', feature_name, '110m',
                                        edgecolor='black', facecolor='silver')
        ## use ax's transform
        ax.add_feature(feature, zorder=10)


    def _add_zebra_frame(self,ax, lw=1.2):
        """
        copied from https://github.com/SciTools/cartopy/issues/1830
        """
        ax.spines["geo"].set_visible(False)
        left, right, bot, top = ax.get_extent()
        
        # Alternate black and white line segments
        bws = itertools.cycle(["k", "white"])

        xticks = sorted([left, *ax.get_xticks(), right])
        xticks = np.unique(np.array(xticks))
        yticks = sorted([bot, *ax.get_yticks(), top])
        yticks = np.unique(np.array(yticks))
        for ticks, which in zip([xticks, yticks], ["lon", "lat"]):
            for idx, (start, end) in enumerate(zip(ticks, ticks[1:])):
                bw = next(bws)
                if which == "lon":
                    xs = [[start, end], [start, end]]
                    ys = [[bot, bot], [top, top]]
                else:
                    xs = [[left, left], [right, right]]
                    ys = [[start, end], [start, end]]

                # For first and lastlines, used the "projecting" effect
                capstyle = "butt" if idx not in (0, len(ticks) - 2) else "projecting"
                for (xx, yy) in zip(xs, ys):
                    ax.plot(
                        xx,
                        yy,
                        color=bw,
                        linewidth=lw,
                        clip_on=False,
                        transform=self.transform_crs,
                        solid_capstyle=capstyle,
                        # Add a black border to accentuate white segments
                        path_effects=[
                            pe.Stroke(linewidth=lw + 1, foreground="black"),
                            pe.Normal(),
                        ],
                    )




class ScatterDataVis:
    "Visualisation object based on ScatterData object"

    def __init__(self, sd):
        self.data = sd.data.reset_index(inplace=False)
        self.index = sd.index

        gp = GridOperation()
        gp.set_coordinates(obj=self, index=self.index)

    def _init_fig(self, *args, **kwargs):
        return plt.subplots(*args, **kwargs)

    def plot(self, var, *args, **kwargs):
        """
        visualise the data based on the dimension of the data
        """
        plt.rcParams['font.family'] = 'sans-serif'

        if isinstance(self.index, str):
            self.index = [self.index]

        match len(self.index):
            case 1:
                return self._plot_1d(var=var, *args, **kwargs)
            case 2:
                if hasattr(self, 'lat') and hasattr(self, 'lon'):
                    return self._plot_map(var=var, *args, **kwargs)
                elif hasattr(self, 'depth') and hasattr(self, 'lat'):
                    return self._plot_transect(var=var,*args, **kwargs)
                else:
                    raise ValueError("self.dim must be lat/lon or depth not real column name")
            case _:
                raise ValueError("Not supported 3D and higher dimensions")

    def _plot_1d(self,var, *args, **kwargs):
        """
        plot 1D data, e.g., zonal_average, time series
        dim: time/lon/lat
        """
        if "ax" not in kwargs:
            fig, local_ax = self._init_fig()
        else:
            local_ax = kwargs.pop("ax")

        ## get the only dimension as x
        x = self.data[self.index[0]]
        y = self.data[var]

        p = local_ax.scatter(x,y, *args, **kwargs)
        local_ax.set_xlabel(self.index[0])
        local_ax.set_ylabel(var)

        return p


    def _plot_map(
        self,
        var,
        ax=None,
        log=False,
        land_mask=True,
        colorbar=True,
        *args,
        **kwargs,
    ):
        """plot map based on dataframe with latitude/longitude
        using cartopy as engine

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
            self.data[var] = efficient_log(self.data[var])

        if self.data[var].dtype != float:
           self.data[var] = self.data[var].astype(float)

        p = ax.scatter(
            self.data[self.lon],
            self.data[self.lat],
            c=self.data[var],
            transform=ccrs.PlateCarree(),
            *args,
            **kwargs,
        )

        if colorbar:
            cbar = plt.colorbar(p, ax=ax, orientation="horizontal")
            cbar.ax.tick_params(axis="both", which="major", labelsize=8)
            cbar.set_label(f"{var}")

        return p    

    def _plot_transect(
        self,
        var,
        ax=None,
        bathy_lon=None,
        *args,
        **kwargs,
    ):
        """plot transect based on dataframe with depth/latitude
        bathy_lon: longitude range for the bathymetry plot
        """
        if not ax:
            fig, ax = self._init_fig()

        if self.data[var].dtype != float:
            self.data[var] = self.data[var].astype(float)
            
        p = ax.scatter(
            x=self.data[self.lat],
            y=self.data[self.depth],
            c=self.data[var],
            *args,
            **kwargs,
        )        

        if bathy_lon:
            data_path = str(files('data').joinpath('GEBCO2002_bathy.nc'))
            bathy = xr.open_dataset(data_path)

            ## get the lat range
            min_lat = self.data[self.lat].min()
            max_lat = self.data[self.lat].max()
            min_depth = self.data[self.depth].min()

            up = bathy.z.sel(lon=slice(*bathy_lon)).mean(dim="lon")
            # up = up.sel(lat=slice(min_lat, max_lat))
            bottom = np.ones(len(up)) * -5500

            ax.fill_between(up.lat, up, bottom, color="black")

        return p


class CommunityPalette:

    """A class to handle community-driven colormaps"""    

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
                f"{cmap_name} not found, accepted values are {self.avail_palettes()}"
            )

        data_dir = files("data").joinpath("colormaps")

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
        data_dir = files("data").joinpath("colormaps")

        return [
            f.stem
            for f in data_dir.glob("*")
            if f.suffix in [".txt", ".xml"]
        ]
