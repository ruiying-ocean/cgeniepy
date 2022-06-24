from functools import reduce
from os.path import join
import string

import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import cartopy.crs as ccrs
import numpy as np
from scipy.stats import sem
import regionmask
from netCDF4 import Dataset
from pandas import DataFrame, read_fwf
import seaborn as sns

from . import ureg, Q_
from .plot import plot_GENIE
from .grid import GENIE_grid_area, reassign_GENIE, GENIE_grid_mask, GENIE_grid_vol, GENIE_lat, GENIE_lon, GENIE_depth, normal_lon
from .data import foram_dict, foram_names, obs_stat_bytype, obs_stat_bysource
from .scores import quick_mscore, quick_rmse, quick_cos_sim
from .utils import file_exists, set_sns_barwidth
from .chem import rm_element, molecular_weight
# from .GENIE_GRID import interp_GENIE


# DONE [X] reassign_array
# DONE [X] basin
# DONE [X] M-score
# DONE [X] search grid
# DONE [X] basic add/subtract
# DONE [X] add streamfunction plot
# DONE [X] redefine calcite class
# DONE [X] add totol POC/calcite/zonal average to GenieModel
# DONE [X] include biogem/ecogem, 2d/3d
# DONE [X] add quiver plot
# DONE [X] add plot statisitcs
# TODO [] add functional diversity
# TODO [] add unit

class GeniePlottable(object):

    def plot_line(self, dim, *args, **kwargs):
        """plot 1D data, e.g., zonal_average
        """

        if self.dim() > 1:
            raise ValueError("Data should be of 1 Dimension")
        if dim=="lon":
            x = GENIE_lat()
        if dim=="lat":
            x = GENIE_lon()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.minorticks_on()

        p = ax.plot(x, self.array, 'k', *args, **kwargs)

        return p

    def plot_map(self, ax=None, cbar=True, *args, **kwargs):

        """plot lat-lon 2D array"""

        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial']

        if not ax:
            fig = plt.figure(dpi=75)
            ax = fig.add_subplot(111, projection=ccrs.EckertIV())

        p = plot_GENIE(ax=ax, data=self.array, *args, **kwargs)

        if cbar:
            cax = fig.add_axes([0.15, 0.1, 0.73, 0.07]) #xmin, ymin, dx, dy
            cbar = fig.colorbar(p, cax = cax, orientation = 'horizontal', pad=0.04)
            cbar.minorticks_on()
            if hasattr(self, "unit"):
                cbar.set_label(self.unit, size=12)

        return p

    def cross_section(self, ax=None, cmap='YlGnBu_r', *args, **kwargs):
        """plot 3D array, tracers and stream function
        """

        data = self.array

        # visualise
        if not ax:
            fig, ax = plt.subplots(figsize=(8, 4))

        lat_edge = GENIE_lat(edge=True)
        z_edge = GENIE_depth(edge=True)/1000
        lat = GENIE_lat(edge=False)
        z = GENIE_depth(edge=False)/1000

        # pcolormesh
        p = ax.pcolormesh(lat_edge, z_edge, data, cmap=cmap, *args, **kwargs)

        # contourf
        # contourf_max = np.max(moc)//10 * 10
        # controuf_N = int(contourf_max*2/10 + 1)
        # contourf_level = np.linspace(contourf_max*-1, contourf_max, controuf_N)
        # cs = ax.contourf(lat, z, data, corner_mask=False, levels=contourf_level, cmap=contourf_cmap)

        # contour
        cs = ax.contour(lat, z, data, linewidths=0.6, colors='black', linestyles='solid')
        ax.clabel(cs, cs.levels[::1], colors=['black'], fontsize=8.5, inline=False)

        # colorbar, ticks & labels
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = plt.colorbar(p, fraction=0.05, pad=0.04, orientation = 'vertical')
        cbar.ax.tick_params(color='k', direction='in')
        # cbar.set_label("Stream function (Sv)")

        plt.rcParams["font.family"] = "Arial"
        ax.patch.set_color("#999DA0")
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel(r"Latitude ($\degree$N)", fontsize=13)
        ax.set_ylabel("Depth (km)", fontsize=12)
        ax.minorticks_on()
        ax.tick_params('both', length=4, width=1, which='major')

        return p


class GenieArray(GeniePlottable):

    def __init__(self, M=36, N=36):
        """
        Create a 36x36 empty 2D array
        """
        # set dimension
        self.M = M
        self.N = N
        # set data
        self.array = self._set_array()
        # set unit
        if hasattr(self.array, "unit"):
            self.unit = self.array.units
        else:
            self.unit = ""

    def _set_array(self):
        return np.zeros((self.M, self.N))

    def pure_array(self):
        "get a numpy array"
        if hasattr(self.array, 'values'):
            return self.array.values
        else:
            return self.array

    def uarray(self):
        "array with unit"
        unit = rm_element(self.unit)
        uarray = Q_(self.pure_array(), unit)
        return uarray

    def dim(self):
        return self.pure_array().ndim

    def flip(self, axis=0):
        return np.flip(self.array, axis=axis)

    def flatten(self):
        "flatten in row-major (C-style)"
        return self.pure_array().flatten(order="C")

    def _value_sign(self, x):
        if np.isnan(x) or x == 0:
            return x
        elif x > 0:
            return 1
        elif x < 0:
            return -1

    def value_sign(self):
        vfunc = np.vectorize(self._value_sign)
        x = GenieArray()
        x.array = vfunc(self.pure_array())
        return x

    def reassign_array(self):
        return reassign_GENIE(self.array)

    def _to_genie_array(self):
        """
        designed for sub-classes to remove all attributes
        """
        empty = GenieArray()
        empty.array = self.array
        return empty

    def __add__(self, other):
        sum = GenieArray()
        sum.array = self.array + other.array
        return sum

    def __sub__(self, other):
        diff = GenieArray()
        diff.array = self.array - other.array
        return diff

    def __truediv__(self, other):
        quotient = GenieArray()
        quotient.array = np.divide(self.array, other.array,
                               out=np.zeros_like(self.array),
                               where=other.array!=0)
        return quotient

    def __mul__(self, other):
        product = GenieArray()
        product.array = self.array * other.array
        return product

    def max(self, *args, **kwargs):
        return np.max(self.array, *args, **kwargs)

    def min(self, *args, **kwargs):
        return np.min(self.array, *args, **kwargs)

    def sum(self, *args, **kwargs):
        return np.sum(self.pure_array(), *args, **kwargs)

    def nansum(self, *args, **kwargs):
        return np.nansum(self.pure_array(), *args, **kwargs)

    def mean(self, *args, **kwargs):
        return np.mean(self.pure_array(), *args, **kwargs)

    def nanmean(self, *args, **kwargs):
        return np.nanmean(self.pure_array(), *args, **kwargs)

    def ptp(self, *args, **kwargs):
        return np.ptp(self.pure_array(), *args, **kwargs)

    def sd(self, *args, **kwargs):
        return np.std(self.pure_array(), *args, **kwargs)

    def nansd(self, *args, **kwargs):
        return np.nanstd(self.pure_array(), *args, **kwargs)

    def se(self, *args, **kwargs):
        return sem(self.array, nan_policy="omit", axis=None, *args, **kwargs)

    def select_basin(self, basin):
        ocean = regionmask.defined_regions.ar6.ocean
        index = ocean.map_keys(basin)
        mask = ocean.mask(self.reassign_array())
        regional_data = self.reassign_array().where(mask == index)

        return regional_data

    def mask_basin(self, base, basin, basin_lvl):
        # mask data
        data = self.pure_array()
        mask = GENIE_grid_mask(base=base, basin=basin, basin_lvl=basin_lvl, invert=True)

        if self.dim() > 2:
            mask = np.broadcast_to(mask, (16, 36, 36))

        mask_data = np.ma.array(data, mask=mask)
        mask_data = np.ma.masked_invalid(mask_data)

        garray = GenieArray()
        garray.array = mask_data

        return garray

    def mean_along(self, dim):
        # time, zt, lat, lon
        if self.dim() == 2:
            if dim=='lat':
                ax = 0
            elif dim=='lon':
                ax = 1
        elif self.dim() == 3:
            if dim=='z':
                ax = 0
            if dim=='lat':
                ax = 1
            elif dim=='lon':
                ax = 2

        garray = GenieArray()
        garray.array = self.nanmean(axis=ax)

        return garray


    def search_grid(self, *args, **kwargs):
        return self.array.sel(*args, **kwargs, method="nearest")

    def search_range(self, lon_min=-255, lon_max=95, lat_min=0, lat_max=90):
        """
        default longitude is unassigned of cGENIE grids
        """

        if lon_min > lon_max or lat_min > lat_max:
            raise ValueError("longitude/latitude min must be less than max!")

        lon = self.array.coords['lon']
        lat = self.array.coords['lat']

        return self.array.loc[
            dict(lat=lat[(lat >= lat_min) & (lat <= lat_max)],
                 lon=lon[(lon >= lon_min) & (lon <= lon_max)])
        ]

    def filter(self, threshold, greater_sign=True):
        data = self.array

        if greater_sign:
            return data.where(data > threshold, drop=True)
        else:
            return data.where(data < threshold, drop=True)


class GenieModel(object):

    def __init__(self, model_path):
        self.model_path = model_path
        self._time = -1

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        if not isinstance(value, int):
            time_array = self.vars(var_name="time")
            raise ValueError(f'Please select [time index] (integer) from {time_array}')
        self._time = value

    def nc_path(self, gem="ecogem", dim="2d"):
        "find netcdf file model_path, default is ecosystem model output"
        model_path = self.model_path
        nc_file = f"fields_{gem}_{dim}.nc"
        nc_path = join(model_path, gem, nc_file)
        if file_exists(nc_path):
            return nc_path

    def GENIE_lon(self, *args, **kwargs):
        return GENIE_lon(*args, **kwargs)

    def GENIE_lat(self, *args, **kwargs):
        return GENIE_lat(*args, **kwargs)

    def normal_lon(self, *args, **kwargs):
        return normal_lon(*args, **kwargs)

    def open_nc(self, path):
        "Note the time default is -1, and will depend on the instance, e.g., ForamVariable"
        return xr.open_dataset(path).isel(time=self.time)

    def auto_find_path(self, var):
        for gem in ["biogem", "ecogem"]:
            for dim in ["2d", "3d"]:
                nc_path = self.nc_path(gem, dim)
                if self.has_var(var, nc_path):
                    return nc_path
        raise ValueError("Variable not found, please check the spelling!")

    def has_var(self, var, nc_path):
        """
        check if variable exists
        """
        t = Dataset(nc_path, "r")
        if_exist = var in t.variables.keys()
        t.close()
        return if_exist

    def get_vars(self, var_name=None, *args, **kwargs):
        "return specified or all available variables"
        if not var_name:
            t = Dataset(self.nc_path(*args, **kwargs), "r")
            tmp = t.variables.keys()
            t.close()
            return tmp
        else:
            return xr.open_dataset(self.nc_path(*args, **kwargs))[var_name]

    def select_foram(self, foram_name:str):
        "a optimised version of select_var()"
        return ForamVariable(foram_name = foram_name, model_path = self.model_path)

    def select_var(self, var:str):
        return GenieVariable(var = var, model_path = self.model_path)

    def grid_mask(self, source="ecogem", Arctic=True, Med=True):
        """
        cGENIE continent mask array (by setting zero values),
        either calculated from existing data or using biogem.grid_mask
        """
        if source=="ecogem":
            data = self.select_var("eco2D_xGamma_T").array
            grid_mask = xr.where(~np.isnan(data), 1, 0).values
        elif source=="biogem":
            grid_mask = self.select_var("grid_mask").array
        else:
            raise ValueError("source only accept ecogem or biogem")

        if Arctic:
            grid_mask[34:36,:]=0
        if Med:
            grid_mask[27:30,25:30] = 0

        return grid_mask

    def marine_area(self):
        "grid area array in km2"
        grid_mask = self.grid_mask()
        grid_area = GENIE_grid_area()
        mask_area = grid_area * grid_mask

        return mask_area

    def marine_volume(self):
        "grid volume array in km3"
        grid_mask = self.grid_mask()
        grid_volume = GENIE_grid_vol()
        mask_volume = grid_volume * grid_mask

        return mask_volume

    def mscore_table(self, table_styler=True, *args, **kwargs):
        "summarised model M-score compared to modern observations"

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())
        df = {
            "Biomass": [ForamVariable(i, self.model_path).carbon_biomass().m_score(*args, **kwargs) for i in foram_abbrev],
            "Carbon Export": [ForamVariable(i, self.model_path).POC_export().m_score(*args, **kwargs) for i in foram_abbrev],
            "Relative Abundance": [ForamVariable(i, self.model_path).POC_export().proportion().m_score(*args, **kwargs) for i in foram_abbrev],
        }

        df = DataFrame(df, index=foram_fullname)
        df['Column Total'] = df.sum(axis=1)
        df.loc['Row Total',:]= df.sum(axis=0)

        if table_styler:
            cm = plt.get_cmap("RdBu_r")
            df = (df.style.set_caption("M-score across foraminifer groups and variables compared to modern observation").
                text_gradient(cmap=cm, subset=(df.index[0:4], df.columns[0:3])))

        return df


    def rmse_table(self, table_styler=True, *args, **kwargs):
        "summarised model M-score compared to modern observations"

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())
        df = {
            "Biomass": [ForamVariable(i, self.model_path).carbon_biomass().rmse(*args, **kwargs) for i in foram_abbrev],
            "Carbon Export": [ForamVariable(i, self.model_path).POC_export().rmse(*args, **kwargs) for i in foram_abbrev],
            "Relative Abundance": [ForamVariable(i, self.model_path).POC_export().proportion().rmse(*args, **kwargs) for i in foram_abbrev],
        }

        df = DataFrame(df, index=foram_fullname)
        df['Column Total'] = df.sum(axis=1)
        df.loc['Row Total',:]= df.sum(axis=0)

        if table_styler:
            cm = plt.get_cmap("RdBu_r")
            df = (df.style.set_caption("RMSE across foraminifer groups and variables compared to modern observation").
                text_gradient(cmap=cm, subset=(df.index[0:4], df.columns[0:3])))

        return df

    def summarise(self, method="mean", diff=False,
            diff_method="percentage", table_styler=True, *args, **kwargs):
        """
        summarise basic statistics of foraminifer groups
        """

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())

        if method=="mean":
            dic = {
                "Biomass(mmol C/m3)": [ForamVariable(i, self.model_path).carbon_biomass().nanmean() for i in foram_abbrev],
                "Carbon Export (mmol C/m3/d)": [ForamVariable(i, self.model_path).POC_export().nanmean() for i in foram_abbrev],
                "Relative Abundance": [ForamVariable(i, self.model_path).POC_export().proportion().nanmean() for i in foram_abbrev],
            }
        elif method=="sd":
            dic = {
                "Biomass(mmol C/m3)": [ForamVariable(i, self.model_path).carbon_biomass().nansd() for i in foram_abbrev],
                "Carbon Export (mmol C/m3/d)": [ForamVariable(i, self.model_path).POC_export().nansd() for i in foram_abbrev],
                "Relative Abundance": [ForamVariable(i, self.model_path).POC_export().proportion().nansd() for i in foram_abbrev],
            }
        elif method=="se":
            dic = {
                "Biomass(mmol C/m3)": [ForamVariable(i, self.model_path).carbon_biomass().se().values for i in foram_abbrev],
                "Carbon Export (mmol C/m3/d)": [ForamVariable(i, self.model_path).POC_export().se().values for i in foram_abbrev],
                "Relative Abundance": [ForamVariable(i, self.model_path).POC_export().proportion().se().values for i in foram_abbrev],
            }

        df = DataFrame(dic, index=foram_fullname)

        if diff:
            obs = obs_stat_bytype(type=method, *args, **kwargs)
            diff = df - obs

            if diff_method == "percentage":
                diff_perc = diff/obs
                return diff_perc
            elif diff_method == "absolute":
                return diff

        if table_styler:
            cm = plt.get_cmap("RdBu_r")
            df = df.style.set_caption("Mean values across foraminifer groups").text_gradient(cmap=cm)

        return df

    def check_completeness(self):
        pass

    def eco_pars(self):
        """
        return ecophysiological parameter table
        """
        path = f"{self.model_path}/ecogem/Plankton_params.txt"
        df = read_fwf(path)
        return df

    def foram_POC(self):
        "Estimate total foraminiferal organic carbon flux rate"

        foram_poc = GenieArray()

        for foram in ["bn", "bs", "sn", "ss"]:
            foram_poc += self.select_foram(foram).POC_export()

        poc_total = ForamCarbonFlux(self.model_path, "ALL_FORAM")
        poc_total.array = foram_poc.array

        return poc_total

    def foram_biomass(self):
        "Estimate total foraminiferal biomass"
        foram_biomass = GenieArray()

        for foram in ["bn", "bs", "sn", "ss"]:
            foram_biomass += self.select_foram(foram).carbon_biomass()

        biomass_total = ForamBiomass(self.model_path, "ALL_FORAM")
        biomass_total.array = foram_biomass.array

        return biomass_total

    def foram_calcite(self):
        "Estimate total foraminiferal inorganic carbon flux rate"
        return self.foram_POC().to_calcite()

    def foram_ldg(self, legend=True):
        "plot ldg: latitudinal diversity gradient"
        ## TODO: change colormap; change line width; add minor ticks
        ## add background grids
        foram_array_list = []
        foram_abbrev = list(foram_names().keys())
        foram_fullname = list(foram_names().values())

        for foram in foram_abbrev:
            array_1d = self.select_foram(foram).POC_export().proportion().mean(axis=1)
            foram_array_list.append(array_1d)

        lat = GENIE_lat()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)

        for n in range(4):
            ax.plot(lat, foram_array_list[n], label = foram_fullname[n])

        if legend:
            ax.legend(loc = 'lower center', ncol = 2, bbox_to_anchor=(0.5, -0.35), edgecolor="black")

        ax.set_xlabel("Latitude")
        ax.set_ylabel("Relative abundance")
        ax.tick_params(axis="y",direction="in")
        ax.tick_params(axis="x",direction="in")

        return ax

    def plot_currents(self, z_level=0):
        """
        plot velocity field of horizontal currents, default is the toppest layer
        Not able to add upon other field map yet.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

        v = self.select_var("phys_v").reassign_array().isel(zt=z_level)
        u = self.select_var("phys_u").reassign_array().isel(zt=z_level)
        m = np.hypot(u, v)

        lat = v.lat
        lon = v.lon
        lon2d, lat2d = np.meshgrid(lon, lat)

        ax.set_global()
        ax.stock_img()
        ax.coastlines()
        ax.quiver(lon2d, lat2d, m, u, v, transform=ccrs.PlateCarree())

        return ax

    def plot_biomass(self, *args, **kwargs):
        """
        quick wrapper function to plot carbon biomass for 4 foram groups
        """

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,7),
                         subplot_kw=dict(projection=ccrs.EckertIV()))

        varlst = ["eco2D_Plankton_C_01" + str(i) for i in range(6,10)]
        most_max = max([self.select_var(varlst[i]).max() for i in range(4)])
        foram_fullnames = list(foram_names().values())

        for i,ax in enumerate(axes.flat):
                vdata = self.select_var(varlst[i]).array
                mean = np.nanmean(vdata)
                p = plot_GENIE(ax=ax, data=vdata, vmin=0, vmax=most_max, *args, **kwargs)
                ax.set_title(f"({string.ascii_lowercase[i]}) {foram_fullnames[i]} {mean:.2E}", pad=10)

        cbar = fig.colorbar(p, ax=axes, orientation="horizontal", pad=0.05, shrink=0.7)
        cbar.minorticks_on()
        cbar.set_label(r'carbon biomass $mmol C m$^{-3}$)', size=12)

        return p

    def plot_export(self):
        """
        quick wrapper function to plot carbon export for 4 foram groups
        """

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,7),
                         subplot_kw=dict(projection=ccrs.EckertIV()))

        varlst = ["eco2D_Export_C_01" + str(i) for i in range(6,10)]
        most_max = max([self.select_var(varlst[i]).max() for i in range(4)])
        foram_fullnames = list(foram_names().values())

        for i,ax in enumerate(axes.flat):
                vdata = self.select_var(varlst[i]).array
                mean = np.nanmean(vdata)
                p = plot_GENIE(ax=ax, data=vdata, vmin=0, vmax=most_max)
                ax.set_title(f"({string.ascii_lowercase[i]}) {foram_fullnames[i]} {mean:.2E}", pad=10)

        cbar = fig.colorbar(p, ax=axes, orientation="horizontal", pad=0.05, shrink=0.7)
        cbar.minorticks_on()
        cbar.set_label(r'POC export at 80.8 m (mmol C m$^{-3}$ d$^{-1})$', size=12)

        return p

    def plot_abundance(self):
        """
        quick wrapper function to plot relative abundance for 4 foram groups
        """

        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,7),
                         subplot_kw=dict(projection=ccrs.EckertIV()))

        varlst = list(foram_names().keys())
        foram_fullnames = list(foram_names().values())

        for i, ax in enumerate(axes.flat):
            vdata = self.select_foram(varlst[i]).POC_export().proportion().array
            mean = np.nanmean(vdata) * 100
            p = plot_GENIE(ax = ax, data = vdata, vmin=0, vmax=1)
            ax.set_title(f"({string.ascii_lowercase[i]}) {foram_fullnames[i]} {mean:.2f}%", pad=10)

        cbar = fig.colorbar(p, ax=axes, orientation="horizontal", pad=0.05, shrink=0.7)
        cbar.minorticks_on()
        cbar.set_label('Relative abundance', size=12)

        return p

    def barplot_comparison(self, *args, **kwargs) -> plt.axes:

        """
        Overview barplot of biomass and POC export compared to observed data
        :returns: matplotlib axes object
        """

        fname = list(foram_names().values())
        foram_abbrev = foram_names().keys()

        model_biomass_mean = [self.select_foram(i).carbon_biomass().nanmean() for i in foram_abbrev]
        model_export_mean = [self.select_foram(i).POC_export().nanmean() for i in foram_abbrev]
        model_biomass_se = [self.select_foram(i).carbon_biomass().se() for i in foram_abbrev]
        model_export_se = [self.select_foram(i).POC_export().se() for i in foram_abbrev]
        model_export_sum = [self.select_foram(i).POC_export().sum() for i in foram_abbrev]
        model_biomass_sum = [self.select_foram(i).carbon_biomass().sum() for i in foram_abbrev]
        obs_biomass_mean = obs_stat_bysource("tow", *args, **kwargs).loc[:, "mean"]
        obs_biomass_se = obs_stat_bysource("tow", *args, **kwargs).loc[:, "se"]
        obs_export_mean = obs_stat_bysource("trap", *args, **kwargs).loc[:, "mean"]
        obs_export_se = obs_stat_bysource("trap", *args, **kwargs).loc[:, "se"]

        model_biomass_sum = DataFrame({"group": fname, "value": model_biomass_sum})
        model_export_sum = DataFrame({"group": fname, "value": model_export_sum})

        data_to_plot = [[[model_biomass_mean, model_biomass_se,obs_biomass_mean, obs_biomass_se],
                        [model_export_mean, model_export_se,obs_export_mean, obs_export_se]],
                        [model_biomass_sum, model_export_sum]]
        ## Plot starts
        fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)
        bar_width = 0.3

        # The x position of bars
        x = np.arange(4)
        xlabels = [w.replace(" ", "\n") for w in fname]

        for i in range(2):
            for j in range(2):
                axes[i,j].yaxis.set_minor_locator(AutoMinorLocator(4))
                axes[i,j].set_axisbelow(True)
                axes[i,j].yaxis.grid(color='gray', linestyle='dashed')
                if i == 0:
                    axes[i,j].bar(x - bar_width/2,
                    data_to_plot[i][j][0],
                    width = bar_width,
                    color = sns.color_palette("Set1")[0],
                    edgecolor = 'black',
                    yerr=data_to_plot[i][j][1],
                    capsize=7,
                    label='model')

                    axes[i,j].bar(x + bar_width/2,
                    data_to_plot[i][j][2],
                    width = bar_width,
                    color = sns.color_palette("Set1")[1],
                    edgecolor = 'black',
                    yerr=data_to_plot[i][j][3],
                    capsize=7,
                    label='obs')

                    axes[i,j].set_xticks(x)
                    axes[i,j].set_xticklabels(xlabels, rotation = 45, ha="right")
                    axes[i,j].legend()
                else:
                    sns.barplot(data=data_to_plot[i][j], x="group", y="value", ax=axes[i,j], edgecolor="black", palette="deep")
                    axes[i,j].set_xticklabels(xlabels, rotation=45, ha="right")
                    axes[i,j].set_xlabel("")
                    axes[i,j].set_xticks(x)
                    set_sns_barwidth(axes[i,j], bar_width)


        axes[0,0].set_ylabel(r"mmol C m$^{-3}$")
        axes[0,1].set_ylabel(r"mmol C m$^{-3}$ d$^{-1}$")

        axes[0,0].set_title("(a)    global biomass mean/se",loc="left")
        axes[0,1].set_title("(b)    global POC flux mean/se", loc="left")
        axes[1,0].set_title("(c)    total biomass production",loc="left")
        axes[1,1].set_title("(d)    total POC production rate", loc="left")

        axes[1,0].set_ylabel("Tg C")
        axes[1,1].set_ylabel(r"Tg C yr$^{-1}$")

        fig.tight_layout()

        return axes

    def diff(self, model2compare, var):
        B =  GenieModel(model2compare)
        diff = self.select_var(var) - B.select_var(var)
        return diff

    def div(self, model2compare, var):
        B =  GenieModel(model2compare)
        diff = self.select_var(var) / B.select_var(var)
        return diff

    def _ptf_presence(self, x, tol=1e-8):
        """
        to determine whethere a functional group present or not
        :param tol: threshold of biomass (mmol C/m3)
        """
        if np.isnan(x):
            return x
        elif x >= tol:
            return 1
        else:
            return 0

    def pft_n(self):
        full_lst = list(self.get_vars())
        name_lst = [x for x in full_lst if 'eco2D_Export_C' in x]
        n = len(name_lst)
        return n

    def pft_richness(self):
        """
        plankton functional group richness, note it is different from species richness,
        because there should be more species in low size classes.
        See `body size-species richness relationship` for more details.
        """
        vfunc = np.vectorize(self._ptf_presence)
        total_sp = np.zeros((36, 36))
        n = self.pft_n()

        for i in range(n):
            name  = f'eco2D_Plankton_C_0{i+1:02d}'
            #select
            arr = self.select_var(name).pure_array()
            #conditional mask
            sp = vfunc(arr)
            #sum
            total_sp += sp

        x = GenieArray()
        x.array = total_sp
        return x


    def select_pft(self, start, end):
        total = np.zeros((36, 36))
        end += 1
        for i in range(start, end):
            name  = f'eco2D_Plankton_C_0{i:02d}'
            arr = self.select_var(name).pure_array()
            total += arr

        x = GenieArray()
        x.array = total
        return x


class GenieVariable(GenieModel, GenieArray):

    def __init__(self, var, model_path):
        # initialise super class to inherite attributes
        self._var = var
        GenieModel.__init__(self, model_path=model_path)
        GenieArray.__init__(self)

    def _set_array(self):
        nc_path = self.auto_find_path(self._var)
        source_data = self.open_nc(nc_path)
        return source_data[self._var]


class ForamVariable(GenieModel):

    def __init__(self, foram_name, model_path):
        self.foram_name = foram_name
        GenieModel.__init__(self, model_path = model_path)

    def carbon_biomass(self):
        return ForamBiomass(model_path = self.model_path, foram_name = self.foram_name)

    def POC_export(self):
        return ForamCarbonFlux(model_path = self.model_path, foram_name = self.foram_name)

    def sum_mscore(self):
        "compared to modern observations"
        return (
            self.POC_export().m_score() +
            self.carbon_biomass().m_score() +
            self.POC_export().proportion().m_score()
        )


class ForamBiogeochem(GenieModel, GenieArray):

    def __init__(self, model_path:str, foram_name:str, biogeo_var:str):
        self.biogeo_var = biogeo_var
        self.foram_name = foram_name
        GenieModel.__init__(self, model_path)
        self.ecogem_path = self.nc_path("ecogem", "2d")
        GenieArray.__init__(self)

    def _set_array(self):
        source_data = self.open_nc(self.ecogem_path)
        foramdict = foram_dict()
        for key, value in foramdict.items():
            if value[1] == self.foram_name and value[0] == self.biogeo_var:
                return source_data[key]

    def proportion(self):
        return ForamProportion(self.model_path, self.foram_name, self.biogeo_var)

    def m_score(self, *args, **kwargs):
        return quick_mscore(self.pure_array(), self.biogeo_var, self.foram_name, *args, **kwargs)

    def rmse(self, *args, **kwargs):
        return quick_rmse(self.pure_array(), self.biogeo_var, self.foram_name, *args, **kwargs)

    def cos_sim(self, *args, **kwargs):
        return quick_cos_sim(self.pure_array(), self.biogeo_var, self.foram_name, *args, **kwargs)


class ForamBiomass(ForamBiogeochem):

    def __init__(self, model_path, foram_name):
        self.biogeo_var = "tow"
        super(ForamBiomass, self).__init__(model_path, foram_name, self.biogeo_var)
        self.unit = "mmol m$^-3$"

    def sum(self):
        C_ = molecular_weight("C")
        c = self.uarray().to_base_units()
        v = self.marine_volume().to_base_units()
        s = c * v
        s = s.to("mol").to('g', 'chemistry', mw = C_ * ureg('g/mole')).to('Tg')

        return np.nansum(s)


class ForamCarbonFlux(ForamBiogeochem):

    def __init__(self, model_path, foram_name):
        self.biogeo_var = "trap"
        super(ForamCarbonFlux, self).__init__(model_path, foram_name, self.biogeo_var)
        self.unit = "mmol m$^-3$ d$^-1$"

    def to_calcite(self):
        calcite = ForamCalcite(model_path = self.model_path)
        calcite.array = self.array * 0.36

        return calcite

    @ureg.with_context('bgc')
    def sum(self):
        C_ = molecular_weight("C")
        c = self.uarray().to_base_units()
        v = self.marine_volume().to_base_units()
        s = c * v
        s = s.to("mol d^-1").to('g d^-1','bgc', mw = C_ * ureg('g/mol')).to("Tg yr^-1")

        return np.nansum(s)

class ForamProportion(ForamBiogeochem):

    def __init__(self, model_path, foram_name, biogeo_var):
        self.biogeo_var = biogeo_var
        ForamBiogeochem.__init__(self, model_path, foram_name, biogeo_var)

    def _set_array(self):
        source_data = self.open_nc(self.ecogem_path)
        foramdict = foram_dict()
        all_foram_vars = []

        for key, value in foramdict.items():
            if value[0] == self.biogeo_var:
                all_foram_vars.append(source_data[key])

        total_foram = reduce(np.add, all_foram_vars)

        # ignore divided by 0
        # and set grid with total_foram == 0 to 0 instead of NA
        with np.errstate(divide='ignore', invalid='ignore'):
            one_foram = super()._set_array()
            proportion =np.divide(one_foram, total_foram, out=np.zeros_like(one_foram), where= total_foram!=0)

        return proportion

    def m_score(self, observation="core", *args, **kwargs):
        return quick_mscore(self.pure_array(), observation, self.foram_name, *args, **kwargs)

    def rmse(self, observation="core", *args, **kwargs):
        return quick_rmse(self.pure_array(), observation, self.foram_name, *args, **kwargs)

    def cos_sim(self, observation="core", *args, **kwargs):
        return quick_cos_sim(self.pure_array(), observation, self.foram_name, *args, **kwargs)


class ForamCalcite(GenieArray, GenieModel):

    def __init__(self, model_path):
        "default empty array, need to assign the array manually, such as ForamCarbonFlux.to_PIC()"
        GenieModel.__init__(self, model_path=model_path)
        self.unit = "mmol m$^-3$ d$^-1$"

    @ureg.with_context('bgc')
    def sum(self):
        CaCO3 = molecular_weight("CaCO3")
        c = self.uarray().to_base_units()
        v = self.marine_volume().to_base_units()
        s = c * v
        s = s.to("mol d^-1").to('g d^-1','bgc', mw = CaCO3 * ureg('g/mol')).to("Pg yr^-1")

        return np.nansum(s)
