from functools import reduce
from os.path import join

import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import cartopy.crs as ccrs
import numpy as np
from scipy.stats import sem
import regionmask
from netCDF4 import Dataset
from pandas import DataFrame, IndexSlice

from .plot import plot_GENIE
from .grid import get_grid_area, reassign_GENIE, get_grid_volume, sum_grids, get_GENIE_lat, get_GENIE_lon, get_normal_lon
from .data import get_foramdict, POC_to_PIC, get_foram_longname
from .scores import quick_mscore, quick_rmse, quick_cos_sim
from .utils import file_exists
#from .GENIE_GRID import interp_GENIE


## DONE [X] reassign_array
## DONE [X] basin
## DONE [X] M-score
## DONE [X] search grid
## DONE [X] basic add/subtract
## DONE [X] add streamfunction plot
## DONE [X] redefine calcite class
## DONE [X] add totol POC/calcite/zonal average to ForamModel
## DONE [X] include biogem/ecogem, 2d/3d
## DONE [X] add quiver plot
## TODO [] add functional diversity
## TODO [] interpolate
## TODO [] add plot statisitcs

class ForamArray(object):

    def __init__(self, M=36, N=36):
        """
        Create a 36x36 empty 2D array
        """
        self.M = M
        self.N = N
        self.array = self._set_array()

    def _set_array(self):
        return np.zeros((self.M, self.N))

    def flip(self, axis=0):
        return np.flip(self.array, axis=axis)

    def flatten(self):
        "flatten in row-major (C-style)"
        return self.pure_array().flatten(order="C")

    def pure_array(self):
        if hasattr(self.array, 'values'):
            return self.array.values
        else:
            return self.array

    def _to_foram_array(self):
        """
        designed for sub-classes to remove all attributes
        """
        empty = ForamArray()
        empty.array = self.array
        return empty

    def plot_zonal_average(self, *args, **kwargs):
        array_1d = self.mean(axis=1)
        lat = get_GENIE_lat()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        p = ax.plot(lat, array_1d, 'k', *args, **kwargs)
        return p

    def plot_map(self, ax=None, cbar=True, *args, **kwargs):

        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Helvetica']

        if not ax:
            fig = plt.figure(dpi=75)
            ax = fig.add_subplot(111, projection=ccrs.EckertIV())

        p = plot_GENIE(self.array, ax, *args, **kwargs)

        if cbar:
            cax = fig.add_axes([0.15, 0.1, 0.73, 0.07]) #xmin, ymin, dx, dy
            cbar = fig.colorbar(p, cax = cax, orientation = 'horizontal', pad=0.04)
            cbar.minorticks_on()
            if hasattr(self, "unit"):
                cbar.set_label(self.unit, size=12)

        return p

    def plot_overturning(self):
        "only for phys_opsi variable"

        if not hasattr(self.array, "lat_moc"):
            raise AttributeError("Not availble MOC array such as 'phys_opsi'!")

        moc = self.array
        lat = self.array.lat_moc
        zt = self.array.zt_moc/1000
        X, Y = np.meshgrid(lat, zt)

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)
        ax.patch.set_color("grey")
        fig.gca().invert_yaxis()

        #contour/controuf
        contourf_cmap = plt.get_cmap("RdBu_r")
        contourf_max = np.max(moc)//10 * 10
        controuf_N = int(contourf_max*2/10 + 1)
        contourf_level = np.linspace(contourf_max*-1, contourf_max, controuf_N)

        CS = ax.contourf(X, Y, moc, corner_mask=False, levels=contourf_level, cmap=contourf_cmap)
        ax.clabel(CS, inline=False, fontsize=8, colors='k')

        ax.contour(X, Y, moc, colors=('k',), linewidths=(1,))

        # x/y labels
        ax.set_ylabel("Depth (km)")
        ax.set_xlabel("Latitude (North)")

        # colorbar
        cax = fig.add_axes([0.92, 0.12, 0.05, 0.75]) #xmin, ymin, dx, dy
        cbar = fig.colorbar(CS, cax = cax, orientation = 'vertical', extend="both")
        cbar.set_label("Stream function (Sv)")

        return ax

    def __add__(self, other):
        sum = ForamArray()
        sum.array = self.array + other.array
        return sum

    def __sub__(self, other):
        diff = ForamArray()
        diff.array = self.array - other.array
        return diff

    def __truediv__(self, other):
        quotient = ForamArray()
        quotient.array = np.divide(self.array, other.array,
                               out=np.zeros_like(self.array),
                               where=other.array!=0)
        return quotient

    def __mul__(self, other):
        product = ForamArray()
        product.array = self.array * other.array
        return product

    def sum(self, *args, **kwargs):
        return np.sum(self.array, *args, **kwargs)

    def nansum(self, *args, **kwargs):
        return np.nansum(self.array, *args, **kwargs)

    def mean(self, *args, **kwargs):
        return np.mean(self.array, *args, **kwargs)

    def sd(self, *args, **kwargs):
        return np.std(self.array, *args, **kwargs)

    def se(self, *args, **kwargs):
        return sem(self.array, nan_policy="omit", axis=None, *args, **kwargs)

    def reassign_array(self):
        return reassign_GENIE(self.array)

    def select_basin(self, basin_name):
        ocean = regionmask.defined_regions.ar6.ocean
        index = ocean.map_keys(basin_name)
        mask_array = ocean.mask(self.reassign_array())
        regional_data = self.reassign_array().where(mask_array == index)

        return regional_data

    def search_grid(self, lat, lon):
        return self.array.sel(lat=lat, lon=lon, method="nearest")

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

# #def interpolate(self):
# #    interp_GENIE()
## return a new foramarray instance


class ForamModel(object):

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
        return get_GENIE_lon(*args, **kwargs)

    def GENIE_lat(self, *args, **kwargs):
        return get_GENIE_lat(*args, **kwargs)

    def normal_lon(self, *args, **kwargs):
        return get_normal_lon(*args, **kwargs)

    def open_nc(self, path):
        "Note the time default is -1, and will depend on the instance, e.g., ForamVariable"
        return xr.open_dataset(path).isel(time=self.time)

    def auto_find_path(self, var):
        for gem in ["biogem", "ecogem"]:
            for dim in ["2d", "3d"]:
                nc_path = self.nc_path(gem, dim)
                if self.search_var(var, nc_path):
                    return nc_path
        raise ValueError("Variable not found in both ecogem and biogem, please check the spelling!")

    def search_var(self, var, nc_path):
        t = Dataset(nc_path, "r")
        if_exist = var in t.variables.keys()
        t.close()
        return if_exist

    def get_vars(self, var_name=None, *args, **kwargs):
        "return specified or all available variables"
        if not var_name:
            t = Dataset(self.nc_path(*args, **kwargs), "r")
            print(t.variables.keys())
            t.close()
            return None
        else:
            return xr.open_dataset(self.nc_path(*args, **kwargs))[var_name]

    def select_foram(self, foram_name:str):
        "a optimised version of select_var()"
        return ForamVariable(foram_name = foram_name, model_path = self.model_path)

    def select_var(self, var:str):
        return GenieVariable(var = var, model_path = self.model_path)

    def mask_array(self, source="ecogem"):
        "cGENIE mask array, either calculated from existing data or use biogem.grid_mask"
        if source=="ecogem":
            data = self.select_var("eco2D_xGamma_T").array
            return xr.where(~np.isnan(data), 1, 0).values
        elif source=="biogem":
            return self.select_var("grid_mask").array
        else:
            raise ValueError("source only accept ecogem or biogem")

    def marine_area(self):
        "grid area array in km2"
        mask_array = self.mask_array()
        grid_area = get_grid_area()
        return mask_array * grid_area

    def marine_volume(self):
        "grid volume array in km3"
        mask_array = self.mask_array()
        grid_volume = get_grid_volume()
        return mask_array * grid_volume

    def mscore_table(self):
        "summarised model M-score compared to modern observations"

        foram_abbrev = list(get_foram_longname().keys())
        foram_longname = tuple(get_foram_longname().values())
        df = {
            "Biomass": [ForamVariable(i, self.model_path).carbon_biomass().m_score() for i in foram_abbrev],
            "Carbon Export": [ForamVariable(i, self.model_path).POC_export().m_score() for i in foram_abbrev],
            "Relative Abundance": [ForamVariable(i, self.model_path).POC_export().proportion().m_score() for i in foram_abbrev],
        }

        df = DataFrame(df, index=foram_longname)
        df['Column Total'] = df.sum(axis=1)
        df.loc['Row Total',:]= df.sum(axis=0)

        cm = plt.get_cmap("RdBu_r")
        df = (df.style.set_caption("M-score across foraminifer groups and variables compared to modern observation").
              text_gradient(cmap=cm, subset=(df.index[0:4], df.columns[0:3])))

        return df

    def foram_POC(self):
        "Estimate total foraminiferal organic carbon flux rate"

        foram_poc = ForamArray()

        for foram in ["bn", "bs", "sn", "ss"]:
            foram_poc += self.select_foram(foram).POC_export()

        poc_total = ForamCarbonFlux(self.model_path, "ALL_FORAM")
        poc_total.array = foram_poc.array

        return poc_total

    def foram_biomass(self):
        "Estimate total foraminiferal biomass"
        foram_biomass = ForamArray()

        for foram in ["bn", "bs", "sn", "ss"]:
            foram_biomass += self.select_foram(foram).carbon_biomass()

        biomass_total = ForamBiomass(self.model_path, "ALL_FORAM")
        biomass_total.array = foram_biomass.array

        return biomass_total

    def foram_PIC(self):
        "Estimate total foraminiferal inorganic carbon flux rate"
        return self.foram_POC().to_PIC()

    def foram_ldg(self, legend=True):
        "plot ldg: latitudinal diversity gradient"
        ## TODO: change colormap; change line width; add minor ticks
        ## add background grids
        foram_array_list = []
        foram_abbrev = list(get_foram_longname().keys())
        foram_longname = list(get_foram_longname().values())

        for foram in foram_abbrev:
            array_1d = self.select_foram(foram).POC_export().proportion().mean(axis=1)
            foram_array_list.append(array_1d)

        lat = get_GENIE_lat()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)

        for n in range(4):
            ax.plot(lat, foram_array_list[n], label = foram_longname[n])

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


class GenieVariable(ForamModel, ForamArray):

    def __init__(self, var, model_path):
        #initialise super class to inherite attributes
        self._var = var
        ForamModel.__init__(self, model_path=model_path)
        ForamArray.__init__(self)

    def _set_array(self):
        nc_path = self.auto_find_path(self._var)
        source_data = self.open_nc(nc_path)
        return source_data[self._var]


class ForamVariable(ForamModel):

    def __init__(self, foram_name, model_path):
        self.foram_name = foram_name
        ForamModel.__init__(self, model_path = model_path)

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


class ForamBiogeochem(ForamModel, ForamArray):

    def __init__(self, model_path:str, foram_name:str, biogeo_var:str):
        self.biogeo_var = biogeo_var
        self.foram_name = foram_name
        ForamModel.__init__(self, model_path)
        self.ecogem_path = self.nc_path("ecogem", "2d")
        ForamArray.__init__(self)

    def _set_array(self):
        source_data = self.open_nc(self.ecogem_path)
        foramdict = get_foramdict()
        for key, value in foramdict.items():
            if value[1] == self.foram_name and value[0] == self.biogeo_var:
                return source_data[key]

    def proportion(self):
        return ForamProportion(self.model_path, self.foram_name, self.biogeo_var)

    def m_score(self):
        return quick_mscore(self.pure_array(), self.biogeo_var, self.foram_name)

    def rmse(self):
        return quick_rmse(self.pure_array(), self.biogeo_var, self.foram_name)

    def cos_sim(self):
        return quick_cos_sim(self.pure_array(), self.biogeo_var, self.foram_name)


class ForamBiomass(ForamBiogeochem):

    def __init__(self, model_path, foram_name):
        self.biogeo_var = "tow"
        self.unit= 'mmol C m$^{-3}$'
        super(ForamBiomass, self).__init__(model_path, foram_name, self.biogeo_var)

    def sum(self):
        """
        A sum of foram biomass across grids
        mmol C/m3 -> Tg C
        """

        print("Unit: TgC")
        return sum_grids(self.array)


class ForamCarbonFlux(ForamBiogeochem):

    def __init__(self, model_path, foram_name):
        self.biogeo_var = "trap"
        self.unit = 'mmol C m$^{-3}$ d$^{-1}$'
        super(ForamCarbonFlux, self).__init__(model_path, foram_name, self.biogeo_var)

    def to_PIC(self):
        calcite = ForamCalcite(model_path = self.model_path)
        calcite.array = POC_to_PIC(self.array)
        return calcite

    def sum(self):
        """
        A sum of organic carbon export across grids
        mmol C/m3/d -> Tg C/d
        """

        print("Unit: TgC/d")
        return sum_grids(self.array)


class ForamProportion(ForamBiogeochem):

    def __init__(self, model_path, foram_name, biogeo_var):
        self.biogeo_var = biogeo_var
        ForamBiogeochem.__init__(self, model_path, foram_name, biogeo_var)

    def _set_array(self):
        source_data = self.open_nc(self.ecogem_path)
        foramdict = get_foramdict()
        all_foram_vars = []

        for key, value in foramdict.items():
            if value[0] == self.biogeo_var:
                all_foram_vars.append(source_data[key])

        total_foram = reduce(np.add, all_foram_vars)

        #ignore divided by 0
        #and set grid with total_foram == 0 to 0 instead of NA
        with np.errstate(divide='ignore', invalid='ignore'):
            one_foram = super()._set_array()
            proportion =np.divide(one_foram, total_foram, out=np.zeros_like(one_foram), where=total_foram!=0)

        return proportion

    def m_score(self, observation="core"):
        return quick_mscore(self.pure_array(), observation, self.foram_name)

    def rmse(self, observation="core"):
        return quick_rmse(self.pure_array(), observation, self.foram_name)

    def cos_sim(self, observation="core"):
        return quick_cos_sim(self.pure_array(), observation, self.foram_name)


class ForamCalcite(ForamArray, ForamModel):

    def __init__(self, model_path):
        "default empty array, need to assign the array manually, such as ForamCarbonFlux.to_PIC()"
        ForamModel.__init__(self, model_path=model_path)

    def sum(self):
        """
        A sum of calcite export across grids
        g/m2/year -> Gt/m2/year -> Gt/year
        """

        area_m2 = self.marine_area() * 1e6
        calcite_g_m2_yr = self.array
        calcite_Gt_m2_yr = calcite_g_m2_yr * 1e-15
        calcite_Gt_yr = calcite_Gt_m2_yr * area_m2

        print("Unit: Gt/yr")
        return calcite_Gt_yr.sum().values
