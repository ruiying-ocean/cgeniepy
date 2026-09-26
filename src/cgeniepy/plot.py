import xml.etree.ElementTree as ET
from importlib.resources import files
import itertools

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, to_rgb as hex_to_rgb, LinearSegmentedColormap, rgb2hex
import matplotlib.patheffects as pe
import matplotlib as mpl
from cgeniepy.grid import GridOperation
from .utils import efficient_log
import warnings


# MOC and psi latitudes and depths are nodes that include the domain limits.
_NODE_DIMS = ("lat_moc", "lat_psi", "lat_opsi", "zt_moc")
# GOLDSTEIN depth levels are evenly spaced in log(depth + ez0 * dsc).
_GENIE_DEPTH_OFFSET = 500.0  # m; ez0 = 0.1, dsc = 5000 m


def _evenly_spaced(values, rtol):
    steps = np.diff(values)
    return np.allclose(steps, steps[0], rtol=rtol)


def _goldstein_depth_edges(values):
    """Edges of GOLDSTEIN depth levels recovered from the depths BIOGEM writes as zt, or None.

    BIOGEM puts the top zt halfway down the top layer and spaces the rest like the
    GOLDSTEIN tracer points, so zt + a = B * g**(k + 1/2) and the edges are
    B * g**k - B, with B = ez0 * dsc and a = B * (1 - (sqrt(g) - 1)**2 / 2).
    Fitting a, g and B to the values finds the edges for any depth scale and
    number of levels.
    """
    if values.size < 3:
        return None
    m = (values.size - 1) // 2
    z0, z1, z2 = values[0], values[m], values[2 * m]
    curvature = z0 + z2 - 2 * z1
    if curvature <= 0:
        return None
    a = (z1**2 - z0 * z2) / curvature
    u = values + a
    g = (u[m] / u[0]) ** (1 / m)
    if not np.allclose(u, u[0] * g ** np.arange(values.size), rtol=1e-6):
        return None
    B = a / (1 - (np.sqrt(g) - 1) ** 2 / 2)
    return np.concatenate([[u[0] / np.sqrt(g)], u * np.sqrt(g)]) - B


def _cell_edges(coordinate):
    """Edges of the cells drawn around each value of a coordinate.

    Depths written by BIOGEM get the edges of their GOLDSTEIN levels. Otherwise
    edges sit halfway between neighbouring values, measured in the space where
    GENIE spaces its grid evenly: sine of latitude, and log(depth + 500 m).
    Other coordinates are bisected linearly. The outer cells of a node grid
    stop at the first and last node.
    """
    values = coordinate.values.astype(float)
    if values.size < 2:
        raise ValueError(f"Cannot derive cell edges from one {coordinate.name} value")
    name = coordinate.name.lower()
    is_lat, _, is_depth, _ = GridOperation.check_dimension((name,))
    is_depth = is_depth and values.min() >= 0
    if is_depth and name not in _NODE_DIMS:
        edges = _goldstein_depth_edges(values)
        if edges is not None:
            return np.clip(edges, 0, None)
    space, kind = values, "linear"
    if is_lat and _evenly_spaced(np.sin(np.deg2rad(values)), rtol=1e-5):
        space, kind = np.sin(np.deg2rad(values)), "sine"
    elif is_depth and _evenly_spaced(np.log(values + _GENIE_DEPTH_OFFSET), rtol=1e-2):
        space, kind = np.log(values + _GENIE_DEPTH_OFFSET), "log"

    middle = (space[:-1] + space[1:]) / 2
    if name in _NODE_DIMS:
        first, last = space[0], space[-1]
    else:
        first, last = 2 * space[0] - middle[0], 2 * space[-1] - middle[-1]
    edges = np.concatenate([[first], middle, [last]])
    if kind == "sine":
        edges = np.rad2deg(np.arcsin(np.clip(edges, -1, 1)))
    elif kind == "log":
        edges = np.exp(edges) - _GENIE_DEPTH_OFFSET
    if is_lat:
        return np.clip(edges, -90, 90)
    if is_depth:
        return np.clip(edges, 0, None)
    return edges


class GriddedDataVis:

    """A class to visualise the GriddedData object"""

    transform_crs = ccrs.PlateCarree()  # do not change
    _MAP_COLORBAR = {"fraction": 0.046, "pad": 0.03, "orientation": "horizontal"}
    _TRANSECT_COLORBAR = {"orientation": "vertical"}

    def __init__(self, gd):
        """
        initialise the GriddedDataVis object

        :param gd: GriddedData object
        """
        self.data = gd.data
        self.attrs = gd.attrs

        ## if long name and units are not provided
        if 'long_name' not in self.attrs:
            self.plot_name = ''
        else:
            self.plot_name = self.attrs['long_name']
            
        if 'units' not in self.attrs:
            self.plot_units = ''
        else:
            self.plot_units = self.attrs['units']

        colourbar_label = f"{self.plot_name}\n{self.plot_units}"

        self.has_negative_and_positive = (self.data.min() < 0 and self.data.max() > 0)

        pal = CommunityPalette(name='parula').colormap
        self._default_cmap = pal  # store default to check if user changed it

        self.aes_dict = {
            ## font None follows Matplotlib's rcParams
            "general_kwargs": {"font": None, "fontsize": 10},
            "facecolor_kwargs": {"c": "silver"}, #silver
            "borderline_kwargs": {"c": "black", "linewidth": 1.0},
            "outline_kwargs": {"colors": "black", "linewidth": 1.0},
            "gridline_kwargs": {"color": "black", "linewidth": 0.25, "linestyle": "--", "draw_labels": False},
            "pcolormesh_kwargs": {"shading": "auto", "cmap": pal},
            "contour_kwargs": {
                "linewidths": 1.35,
                "colors": "black",
                "linestyles": "solid",
                "zorder": 10,
                "levels": 15,
            },
            "contour_label_kwargs": {
                "colors": ["black"],
                "fontsize": 10,
                "inline": False,
            },
            "contourf_kwargs": {"levels": 15},
            "colorbar_label_kwargs": {
                "label": colourbar_label,
                "size": 10,
                "labelpad": 10,
            },
            ## applied over _MAP_COLORBAR or _TRANSECT_COLORBAR
            "colorbar_kwargs": {},
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
        plt.rcParams['lines.antialiased'] = True


        def validate_and_get_cmap(cmap_name):
            "even if our cmap is not registered successfully, it will still be used"
            try:
                if isinstance(cmap_name, mpl.colors.Colormap):
                    return cmap_name
                if cmap_name is None:
                    cmap_name = mpl.rcParams["image.cmap"]
                return mpl.colormaps[cmap_name]
            except (KeyError, ValueError):
                try:
                    return CommunityPalette(name=cmap_name).colormap
                except (FileNotFoundError, ValueError) as exc:
                    raise ValueError(
                        f"Colormap '{cmap_name}' not found in matplotlib or "
                        "CommunityPalette"
                    ) from exc

            
        ## if value has both negative and positive values, set cmap to diverging
        if kwargs.get('pcolormesh', True):

            ## only override if cmap not in kwargs AND user hasn't changed aes_dict cmap from default
            cmap_is_default = self.aes_dict['pcolormesh_kwargs']['cmap'] is self._default_cmap
            if 'cmap' not in kwargs and self.has_negative_and_positive and cmap_is_default:
                div_pal = CommunityPalette(name='cspace_BlRd').colormap
                self.aes_dict['pcolormesh_kwargs']['cmap'] = div_pal

            ## try to update pcolor_kwargs
            if 'vmin' in kwargs:
                    self.aes_dict['pcolormesh_kwargs']['vmin'] = kwargs['vmin']
            if 'vmax' in kwargs:
                self.aes_dict['pcolormesh_kwargs']['vmax'] = kwargs['vmax']
            if 'cmap' in kwargs:
                validated_cmap = validate_and_get_cmap(kwargs['cmap'])
                self.aes_dict['pcolormesh_kwargs']['cmap'] = validated_cmap


        ## for contourf plots, do the same
        if kwargs.get('contourf', True):
            ## only override if cmap not in kwargs AND user hasn't changed aes_dict cmap from default
            contourf_cmap = self.aes_dict['contourf_kwargs'].get('cmap')
            cmap_is_default = contourf_cmap is None or contourf_cmap is self._default_cmap
            if 'cmap' not in kwargs and self.has_negative_and_positive and cmap_is_default:
                div_pal = CommunityPalette(name='cspace_BlRd').colormap
                self.aes_dict['contourf_kwargs']['cmap'] = div_pal

            ## try to update contourf_kwargs
            if 'vmin' in kwargs:
                    self.aes_dict['contourf_kwargs']['vmin'] = kwargs['vmin']
            if 'vmax' in kwargs:
                    self.aes_dict['contourf_kwargs']['vmax'] = kwargs['vmax']
            if 'cmap' in kwargs:
                validated_cmap = validate_and_get_cmap(kwargs['cmap'])
                self.aes_dict['contourf_kwargs']['cmap'] = validated_cmap
        

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
            fig, local_ax = self._init_fig(figsize=(5, 3.5))
        else:
            local_ax = kwargs.pop("ax")

        ## get the only dimension as x
        dim = self.data[self.data.dims[0]]
        
        p = local_ax.plot(dim, self.data, *args, **kwargs)
            
        
        local_ax.grid(True, linestyle='--', alpha=0.5)
        local_ax.minorticks_on()
        local_ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)        
        local_ax.set_xlabel(dim.name)
        local_ax.set_ylabel(f"{self.plot_name} ({self.plot_units})")

        return p

    def _plot_2d(self, *arg, **kwargs):
        "plot 2D data, e.g., map, trasect"
        ## if lon, lat then plot map
        ## if zt or depth then plot transec   

        has_lat, has_lon, has_depth, _ = GridOperation().check_dimension(self.data.dims)
        if has_lat and has_lon:
            return self._plot_map(*arg, **kwargs)
        
        elif has_depth:
            return self._plot_transect(*arg, **kwargs)
        else:
            raise ValueError("2D plotting not supported")

    def _plot_3d(self, *args, **kwargs):
        "plot 3D data = plot mutiple 2D plots"
        _, _, has_depth, has_time = GridOperation().check_dimension(self.data.dims)
        
        if has_time or has_depth:
            target_order = GridOperation().dim_order(self.data.dims)[0]
            target_name = self.data.dims[target_order]
            target_arr = self.data[target_name]

            # Determine the optimal number of columns and rows for subplots
            num_plots = len(target_arr)
            ncols = min(num_plots, 3)
            nrows = (num_plots + ncols - 1) // ncols  # Calculate the required number of rows

            fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3 * ncols, 2.2 * nrows),
                                    sharex=True, sharey=True, squeeze=False)

            # Flatten the axes array for easier iteration
            axs = axs.flatten()

            # Plot the data for each time step                
            for i, ax in enumerate(axs):
                if i < num_plots:
                    im = self.data.isel({target_name: i}).plot(ax=ax, add_colorbar=False)
                    title = f"{target_name} = {target_arr[i].values:.2f} {self.data[target_name].attrs['units']}"
                    ax.set_title(title, fontsize=10)
                    ax.set_xlabel("")
                    ax.set_ylabel("")
                else:
                    ax.set_visible(False)

            # Add a common colorbar to the figure
            fig.colorbar(im, ax=axs.tolist(), orientation='horizontal', label=f"{self.plot_name} ({self.plot_units})",
                         fraction=0.046, pad=0.05)            
        else:
            raise ValueError("Not support 3D plot iterating over other dimension than time and depth")

    def _map_coordinates(self):
        lat_order, lon_order = GridOperation().dim_order(self.data.dims)
        lon = self.data[self.data.dims[lon_order]]
        lat = self.data[self.data.dims[lat_order]]
        return lon.values, lat.values, _cell_edges(lon), _cell_edges(lat)

    def _split_at_map_edge(self, ax, x_edge):
        """Split cells that straddle the map's edge, which cartopy fails to wrap (#4).

        Both halves keep the cell's value, so nothing moves on the map.
        """
        data = np.asarray(self.data)
        projection = getattr(ax, "projection", None)
        if projection is None or not np.all(np.diff(x_edge) > 0):
            return x_edge, data
        lon_axis = GridOperation().dim_order(self.data.dims)[1]
        map_edge = projection.proj4_params.get("lon_0", 0.0) + 180
        turns = np.arange(
            np.ceil((x_edge[0] - map_edge) / 360),
            np.floor((x_edge[-1] - map_edge) / 360) + 1,
        )
        for boundary in map_edge + 360 * turns:
            index = np.searchsorted(x_edge, boundary)
            if (
                0 < index < x_edge.size
                and not np.isclose(x_edge[index - 1 : index + 1], boundary).any()
            ):
                x_edge = np.insert(x_edge, index, boundary)
                data = np.insert(
                    data, index - 1, data.take(index - 1, axis=lon_axis), axis=lon_axis
                )
        return x_edge, data

    def _close_longitude(self, x_arr, x_edge):
        """Repeat the first column after the last so contours close on a global grid."""
        data = np.asarray(self.data)
        if not np.isclose(x_edge[-1] - x_edge[0], 360):
            return x_arr, data
        lon_axis = GridOperation().dim_order(self.data.dims)[1]
        data = np.concatenate([data, data.take([0], axis=lon_axis)], axis=lon_axis)
        return np.append(x_arr, x_arr[0] + 360), data

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
        engine='cartopy',
        *args,
        contour_label=True,
        **kwargs,
    ):

        if engine == 'pygmt':
            import pygmt
            fig = pygmt.Figure()
            ## plot xarray data
            fig.grdimage(self.data, projection="Q12c", cmap="turbo", frame=True)
            return fig

        x_arr, y_arr, x_edge, y_edge = self._map_coordinates()


        if "ax" not in kwargs:
            fig, local_ax = self._init_fig(subplot_kw={"projection": ccrs.EqualEarth()})
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


        plots = {}
        if pcolormesh:
            ## pcolormesh uses edge coordinates
            ## need to be transformed to PlateCarree
            mesh_x_edge, mesh_data = self._split_at_map_edge(local_ax, x_edge)
            plots["pcolormesh"] = self._add_pcolormesh(
                local_ax,
                x=mesh_x_edge,
                y=y_edge,
                data=mesh_data,
                transform=self.transform_crs,
                *args,
                **self.aes_dict["pcolormesh_kwargs"],
            )

        ## contours use center coordinates, closed across the seam
        contour_x, contour_data = self._close_longitude(x_arr, x_edge)
        if contour:
            ## need to be transformed to PlateCarree
            plots["contour"] = self._add_contour(
                local_ax,
                contour_x,
                y_arr,
                data=contour_data,
                transform=self.transform_crs,
                **self.aes_dict["contour_kwargs"],
            )
            if contour_label:
                self._add_contour_label(
                    local_ax, plots["contour"], **self.aes_dict["contour_label_kwargs"]
                )

        if contourf:
            plots["contourf"] = self._add_contourf(
                local_ax,
                contour_x,
                y_arr,
                data=contour_data,
                transform=self.transform_crs,
                **self.aes_dict["contourf_kwargs"],
            )

        if colorbar:
            self._add_labelled_colorbar(local_ax, plots, self._MAP_COLORBAR)

        if cfeature:
            self._add_cfeature(local_ax, cfeature)

        if zebra_frame:
            add_zebra_frame(local_ax)

        return self._main_plot(plots, local_ax)

    def _transect_coordinates(self):
        ## the horizontal axis is latitude, or longitude for zonal sections
        zt_order, x_order = GridOperation().dim_order(self.data.dims)
        x = self.data[self.data.dims[x_order]]
        zt = self.data[self.data.dims[zt_order]]
        return x.name, zt.name, x.values, zt.values, _cell_edges(x), _cell_edges(zt)

    def _plot_transect(
        self,
        pcolormesh=True,
        contour=False,
        colorbar=True,
        contourf=False,
        outline=False,
        facecolor=True,
        borderline=True,
        *args,
        contour_label=True,
        gridline=True,
        **kwargs,
    ):
        """
        Examples
        ---------------

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
        x_name, y_name, x_arr, y_arr, x_edge, y_edge = self._transect_coordinates()

        if "ax" not in kwargs:
            fig, local_ax = self._init_fig(figsize=(5, 3))
        else:
            local_ax = kwargs.pop("ax")

        if gridline:
            local_ax.grid(which='major', linestyle='--', linewidth=0.5, alpha=0.7)


        if facecolor:
            self._set_facecolor(local_ax, **self.aes_dict["facecolor_kwargs"])

        if borderline:
            self._set_borderline(
                local_ax, geo=False, **self.aes_dict["borderline_kwargs"]
            )

        if outline:
            ## outline uses edge coordinates
            self._add_outline(
                local_ax, x=x_edge, y=y_edge, **self.aes_dict["outline_kwargs"]
            )

        plots = {}
        if pcolormesh:
            ## pcolormesh uses edge coordinates
            plots["pcolormesh"] = self._add_pcolormesh(
                local_ax, x=x_edge, y=y_edge, *args, **self.aes_dict["pcolormesh_kwargs"],

            )

        if contour:
            ## contour uses center coordinates
            plots["contour"] = self._add_contour(
                local_ax, x_arr, y_arr, **self.aes_dict["contour_kwargs"]
            )
            if contour_label:
                self._add_contour_label(
                    local_ax, plots["contour"], **self.aes_dict["contour_label_kwargs"]
                )


        if contourf:
            plots["contourf"] = self._add_contourf(local_ax, x_arr, y_arr, **self.aes_dict["contourf_kwargs"])

        if colorbar:
            self._add_labelled_colorbar(local_ax, plots, self._TRANSECT_COLORBAR)

        ## reverse y axis
        local_ax.set_ylim(local_ax.get_ylim()[::-1])
        fontsize = self.aes_dict["general_kwargs"]["fontsize"]
        family = self._font_family()
        local_ax.set_ylabel(y_name, fontsize=fontsize, **family)
        local_ax.set_xlabel(x_name, fontsize=fontsize, **family)

        ## x/y tick label
        local_ax.tick_params(
            axis="both",
            which="major",
            direction="out",
            labelsize=fontsize,
            **{f"label{key}": value for key, value in family.items()},
        )
        # local_ax.minorticks_on()

        return self._main_plot(plots, local_ax)

    def _font_family(self):
        font = self.aes_dict["general_kwargs"]["font"]
        return {} if font is None else {"fontfamily": font}

    @staticmethod
    def _main_plot(plots, ax):
        """Return the filled plot first, so it can feed a colorbar."""
        for plot_type in ("pcolormesh", "contourf", "contour"):
            if plot_type in plots:
                return plots[plot_type]
        return ax

    def _add_labelled_colorbar(self, ax, plots, defaults):
        ## contour lines are black, so only filled plots get a colorbar
        mappable = plots.get("pcolormesh", plots.get("contourf"))
        if mappable is None:
            return None
        cbar = self._add_colorbar(
            mappable, **{"ax": ax, **defaults, **self.aes_dict["colorbar_kwargs"]}
        )
        self._add_colorbar_label(cbar, **self.aes_dict["colorbar_label_kwargs"])
        return cbar


    ## ------- Below is implementations -------------------------

    def _init_fig(self,dpi=120, *args, **kwargs):
        return plt.subplots(dpi=dpi, *args, **kwargs)

    def _init_pygmt_fig(self, *args, **kwargs):
        import pygmt
        fig = pygmt.Figure()
        fig.basemap(*args, **kwargs)
        return fig

    def _add_pcolormesh(self, ax, x, y, *args, data=None, **kwargs):
        return ax.pcolormesh(x, y, self.data if data is None else data, *args, **kwargs)

    def _add_contour(self, ax, x, y, *args, data=None, **kwargs):
        return ax.contour(x, y, self.data if data is None else data, *args, **kwargs)

    def _add_contourf(self, ax, x, y, *args, data=None, **kwargs):
        return ax.contourf(x, y, self.data if data is None else data, *args, **kwargs)

    def _add_contour_label(self, ax, cs, *args, levels=None, **kwargs):
        """Label ``levels``, or every other contour level by default."""
        if levels is None:
            levels = cs.levels[::2]
        ax.clabel(cs, levels, *args, **kwargs)

    def _add_gridline(self, ax, *args, **kwargs):
        gl = ax.gridlines(*args, **kwargs)        

        if ax.projection == ccrs.PlateCarree():
            gl.xlines = False  # removing gridlines
            gl.ylines = False
            add_border_ticks(ax, 0.015)


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
            **{f"label{key}": value for key, value in self._font_family().items()},
        )
        cbar.outline.set_edgecolor("black")
        cbar.outline.set_linewidth(0.5)

    def _add_cfeature(self, ax, feature_name):
        """add natural earth mask"""

        feature = cfeature.NaturalEarthFeature('physical', feature_name, '110m',
                                        edgecolor='black', facecolor='silver')
        ## use ax's transform
        ax.add_feature(feature, zorder=10)



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
                    ## drop kind if it exsit in kwargs
                    del kwargs['kind']
                    warnings.warn("kind is not used for map plotting")
                    return self._plot_map(var=var, *args, **kwargs)
                elif hasattr(self, 'depth') and hasattr(self, 'lat'):
                    del kwargs['kind']
                    warnings.warn("kind is not used for transect plotting")
                    return self._plot_transect(var=var,*args, **kwargs)
                else:
                    raise ValueError("self.dim must be lat/lon or depth not real column name")
            case _:
                raise ValueError("Not supported 3D and higher dimensions")

    def _plot_1d(self, var, kind='scatter', *args, **kwargs):
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

        if kind == 'line':
            p = local_ax.plot(x, y, *args, **kwargs)
        elif kind == 'scatter':
            p = local_ax.scatter(x, y, *args, **kwargs)
        else:
            raise ValueError("Invalid plot type. Choose 'line' or 'scatter'.")
        
        local_ax.grid(True, linestyle='--', alpha=0.5)
        local_ax.minorticks_on()
        local_ax.tick_params(axis='both', which='both', direction='in', top=True, right=True)

        return p


    def _plot_map(
        self,
        var,
        ax=None,
        log=False,
        land_mask=True,
        colorbar=True,
        gridline=True,
        zebra_frame=False,
        mask_age = None,
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
            fig, ax = self._init_fig(subplot_kw={"projection": ccrs.PlateCarree()})

        if land_mask:
            if not mask_age:
                ax.set_global()
                ## plot land and coastline, zorder is the drawing order, smaller -> backer layer
                # ax.stock_img()
                ax.add_feature(cfeature.LAND.with_scale("110m"), zorder=2, facecolor="lightgrey")
                ax.add_feature(cfeature.COASTLINE.with_scale("110m"), linewidth=1.4)
                ax.add_feature(cfeature.LAKES.with_scale('110m'), zorder=3, facecolor='black')
            else:
                from gwspy import PlateModel

                plate_model = PlateModel("Muller2022")
                coastlines_shapely = plate_model.get_coastlines(
                    time=mask_age, 
                    format="shapely"
                )

                ax.add_geometries(
                    coastlines_shapely,
                    crs=ccrs.PlateCarree(),
                    facecolor="grey",
                    edgecolor=None,
                    alpha=0.5,
                )

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
            cbar = plt.colorbar(p, ax=ax, orientation="horizontal",
                                label=f"{var}", pad=0.1)
            cbar.ax.tick_params(axis="both", which="major", direction="out", length=5)
            cbar.ax.minorticks_on()

        if zebra_frame:
            add_zebra_frame(ax)

        if gridline:
            gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=.5,
                              color='black', alpha=0.5,
                              linestyle='--', zorder=1)
            gl.top_labels = False
            gl.right_labels = False
            if isinstance(ax.projection, ccrs.PlateCarree):
                add_border_ticks(ax, 0.01)

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
            data_path = files(__package__ + '.data').joinpath('GEBCO2002_bathy.nc')
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

    def __init__(self, name=None, *args, **kwargs):
        self.name = name
        
        if self.name:
            self.colormap = self.get_palette(name, *args, **kwargs)
        else:
            self.colormap=None

    def get_palette(self, cmap_name, N=256, reverse=False, alpha=None):
        """
        community-driven colormaps with multiple sources

        :param cmap_name: colormap name, can be found in avail_palette()
        :type cmap_name: str

        :returns: colormap

        XML data: https://sciviscolor.org/colormaps/
        txt data: from original packages
        """

        ## if _r exsit in the end of the string
        if cmap_name[-2:] == "_r":
            cmap_name = cmap_name[:-2]
            reverse = True

        if cmap_name not in self.avail_palettes():
            raise ValueError(
                f"{cmap_name} not found, accepted values are {self.avail_palettes()}"
            )


        data_dir = files(__package__ + '.data').joinpath('colormaps')

        file_path = None
        file_ext = None
        colors = []
        c = None

        file_path = [f for f in data_dir.glob("**/*") if f.stem == cmap_name]        
        file_path = file_path[0]
        file_ext = file_path.suffix


        match file_ext:
            case ".txt":
                colors = self._parse_txt_data(file_path)

                if alpha is not None:
                    rgba_colors = [(*hex_to_rgb(color), alpha) for color in colors]
                    c = ListedColormap(rgba_colors, name=cmap_name)
                else:
                    c = ListedColormap(colors, name=cmap_name)
                    
            case ".xml":
                colors = self._parse_xml_data(file_path)
                if alpha is not None:
                    rgba_colors = [(color[0], color[1], color[2], alpha) for color in colors]                    
                    c = ListedColormap(rgba_colors, name=cmap_name)
                else:
                    c = ListedColormap(colors, name=cmap_name)
            case ".spk":
                data = self._parse_spk_data(file_path)
                pos = data[:,0]
                ## normalise pos to 0-1
                pos_min = np.min(pos)
                pos_max = np.max(pos)
                scaled_pos = [(val - pos_min) / (pos_max - pos_min) for val in pos]

                rgb = data[:,1:]/100

                ## if alpha is not None, add alpha to the rgb
                if alpha is not None:
                    rgba = np.zeros((len(rgb), 4))
                    rgba[:,0:3] = rgb
                    rgba[:,3] = alpha
                    c = self.create_colormap(scaled_pos, rgba)
                else:
                    c = self.create_colormap(scaled_pos, rgb)
            case _:
                raise ValueError("File extension not supported")
            
        interval = np.linspace(0, 1, N)
        c = ListedColormap(c(interval), name=cmap_name)
        
        if reverse and c is not None:
            return c.reversed()

        if c is None:
            raise ValueError("Colormap could not be created")

        return c

    def avail_palettes(self, show_ferret_data=True):
        """return a list of colormap names"""
        data_dir = files(__package__ + '.data').joinpath('colormaps')

        if show_ferret_data:
            return [
                f.stem
                for f in data_dir.glob("*")
                if f.suffix in [".txt", ".xml", ".spk"]
            ]
        else:
            return [
                f.stem
                for f in data_dir.glob("*")
                if f.suffix in [".txt", ".xml"]
            ]

    def _parse_spk_data(self, filename):
        data = []
        with open(filename, 'r') as file:
            for line in file:
                # Skip empty lines and comment lines starting with "!"
                if "RGB_Mapping" in line or not line.strip() or line.strip().startswith('!'):
                    continue
                # Split the line by whitespace and convert to floats
                if "!" in line:
                    line = line.split("!")[0]
                values = list(map(float, line.split()))
                data.append(values)
        return np.array(data)

    def _parse_txt_data(self, filename):
        colors = pd.read_csv(filename, header=None).values.tolist()
        colors = [colors[i][0] for i in range(len(colors))]
        return colors

    def _parse_xml_data(self, filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        colors = []
        for point in root.findall(".//Point"):
            r = float(point.get("r"))
            g = float(point.get("g"))
            b = float(point.get("b"))
            colors.append((r, g, b))
        return colors
        

    def to_rgb(self):
        return self.colormap.colors

    def to_hex(self, unique=True):
        df = pd.DataFrame(self.colormap.colors)
        ## get unique rows
        if unique:
            df = df.drop_duplicates()
            
        x_unique = df.to_numpy()

        hex_color = [rgb2hex(i) for i in x_unique]
        return hex_color

    
    def create_colormap(self, positions, colors):
        """
        Create a colormap with specified positions and colors.

        Args:
            positions (list): List of floats indicating the positions of colors in the colormap.
            colors (list): List of RGB tuples representing the colors.

        Returns:
            LinearSegmentedColormap: The created colormap.
        """
        segments = {'red': [], 'green': [], 'blue': []}
        for pos, color in zip(positions, colors):
            for i, key in enumerate(['red', 'green', 'blue']):
                segments[key].append((pos, color[i], color[i]))

        cmap = LinearSegmentedColormap('custom_colormap', segments)
        return cmap


    def __repr__(self):
        x = np.linspace(0, 1, 256).reshape(1,-1)
        plt.imshow(x, cmap=self.colormap, aspect=20)
        plt.axis('off')
        plt.show()
        ## title

        return f"Colormap: {self.name}"


def add_zebra_frame(ax, lw=1.2):
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
                    transform=ccrs.PlateCarree(),
                    solid_capstyle=capstyle,
                    # Add a black border to accentuate white segments
                    path_effects=[
                        pe.Stroke(linewidth=lw + 1, foreground="black"),
                        pe.Normal(),
                    ],
                )

def add_border_ticks(ax, tick_len_scale=0.015):
    """
    copied from https://github.com/SciTools/cartopy/issues/2003
    """
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
    tick_length = tick_len_scale * width_pixels            
    ax.tick_params(axis='x', length=tick_length)
    ax.tick_params(axis='y', length=tick_length)

    
    
