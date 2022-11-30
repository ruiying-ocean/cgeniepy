import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

import cartopy.crs as ccrs
import string
import numpy as np
from functools import reduce
from pandas import DataFrame
import seaborn as sns

from .core import GenieArray, GenieVariable
from .ecology import GenieModel, PlanktonType, PlanktonBiomass, PlanktonExport
from . import ureg
from .plot import plot_genie
from .data import foram_names, obs_stat_bytype, obs_stat_bysource
from .scores import quick_mscore, quick_rmse, quick_cos_sim, quick_corr
from .utils import set_sns_barwidth, distance
from .chem import molecular_weight
from .grid import GENIE_lat, GENIE_grid_area


class ForamModel(GenieModel):
    """A further customized GenieModel subclass
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def select_foram(self, foram_type):
        "a optimised version of select_var()"
        return ForamType(foram_type = foram_type, model_path=self.model_path)

    def skill_score(self, cost_function="m_score", table_styler=True, *args, **kwargs):
        "summarised model skill score compared to modern observations"

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())
        df = {
            "Biomass": [
                ForamType(i, self.model_path)
                .biomass_c()
                ._run_method(method=cost_function, *args, **kwargs)
                for i in foram_abbrev
            ],
            "Carbon Export": [
                ForamType(i, self.model_path)
                .export_c()
                ._run_method(method=cost_function, *args, **kwargs)
                for i in foram_abbrev
            ],
            "Relative Abundance": [
                ForamType(i, self.model_path)
                .export_c()
                .proportion()
                ._run_method(method=cost_function, *args, **kwargs)
                for i in foram_abbrev
            ],
        }

        df = DataFrame(df, index=foram_fullname)

        df["Column Total"] = df.sum(axis=1)
        df.loc["Row Total", :] = df.sum(axis=0)

        if table_styler:
            df = df.style.set_caption(
                f"{cost_function} across foraminifer groups and variables compared to modern observation"
            ).text_gradient(cmap="viridis", subset=(df.index[0:4], df.columns[0:3]))

        return df

    def summarise(
        self,
        stat="nanmean",
        diff=False,
        diff_method="percentage",
        table_styler=True,
        *args,
        **kwargs,
    ):
        """
        summarise basic statistics of foraminifer groups
        """

        foram_abbrev = list(foram_names().keys())
        foram_fullname = tuple(foram_names().values())

        dic = {
            "Biomass": [
                ForamType(i, self.model_path).biomass_c()._run_method(method=stat)
                for i in foram_abbrev
            ],
            "Carbon Export": [
                ForamType(i, self.model_path).export_c()._run_method(method=stat)
                for i in foram_abbrev
            ],
            "Relative Abundance": [
                ForamType(i, self.model_path)
                .export_c()
                .proportion()
                ._run_method(method=stat)
                for i in foram_abbrev
            ],
        }
        df = DataFrame(dic, index=foram_fullname)

        df["Column Total"] = df.sum(axis=1)
        df.loc["Row Total", :] = df.sum(axis=0)

        if diff:
            obs = obs_stat_bytype(type=stat, *args, **kwargs)
            diff = df - obs

            if diff_method == "percentage":
                diff_perc = diff / obs
                return diff_perc
            elif diff_method == "absolute":
                return diff

        if table_styler:
            df = df.style.set_caption(
                f"{stat} across foraminifer groups"
            ).text_gradient(cmap="viridis")

        return df

    def foram_poc(self):
        "Estimate integraded foraminiferal organic carbon flux rate"

        foram_poc = GenieArray()

        for foram in ["bn", "bs", "sn", "ss"]:
            foram_poc += self.select_foram(foram).export_c()

        poc_total = ForamPOC(self.model_path, "ALL_FORAM")
        poc_total.array = foram_poc.array

        return poc_total

    def foram_biomass(self):
        "Estimate total foraminiferal biomass"
        foram_biomass = GenieArray()

        for foram in ["bn", "bs", "sn", "ss"]:
            foram_biomass += self.select_foram(foram).biomass_c()

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
            array_1d = self.select_foram(foram).export_c().proportion().mean(axis=1)
            foram_array_list.append(array_1d)

        lat = GENIE_lat()

        fig = plt.figure(figsize=(6, 4))
        ax = fig.add_subplot(111)

        for n in range(4):
            ax.plot(lat, foram_array_list[n], label=foram_fullname[n])

        if legend:
            ax.legend(
                loc="lower center",
                ncol=2,
                bbox_to_anchor=(0.5, -0.35),
                edgecolor="black",
            )

        ax.set_xlabel("Latitude")
        ax.set_ylabel("Relative abundance")
        ax.tick_params(axis="y", direction="in")
        ax.tick_params(axis="x", direction="in")

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

        fig, axes = plt.subplots(
            nrows=2,
            ncols=2,
            figsize=(10, 7),
            subplot_kw=dict(projection=ccrs.EckertIV()),
        )

        varlst = ["eco2D_Plankton_C_01" + str(i) for i in range(6, 10)]
        most_max = max([self.select_var(varlst[i]).max() for i in range(4)])
        foram_fullnames = list(foram_names().values())

        for i, ax in enumerate(axes.flat):
            vdata = self.select_var(varlst[i]).array
            mean = np.nanmean(vdata)
            p = plot_genie(ax=ax, data=vdata, vmin=0, vmax=most_max, *args, **kwargs)
            ax.set_title(
                f"({string.ascii_lowercase[i]}) {foram_fullnames[i]} {mean:.2E}", pad=10
            )

        cbar = fig.colorbar(p, ax=axes, orientation="horizontal", pad=0.05, shrink=0.7)
        cbar.minorticks_on()
        cbar.set_label(r"carbon biomass $mmol C m$^{-3}$)", size=12)

        return p

    def plot_export(self):
        """
        quick wrapper function to plot carbon export for 4 foram groups
        """

        fig, axes = plt.subplots(
            nrows=2,
            ncols=2,
            figsize=(10, 7),
            subplot_kw=dict(projection=ccrs.EckertIV()),
        )

        varlst = ["eco2D_Export_C_01" + str(i) for i in range(6, 10)]
        most_max = max([self.select_var(varlst[i]).max() for i in range(4)])
        foram_fullnames = list(foram_names().values())

        for i, ax in enumerate(axes.flat):
            vdata = self.select_var(varlst[i]).array
            mean = np.nanmean(vdata)
            p = plot_genie(ax=ax, data=vdata, vmin=0, vmax=most_max)
            ax.set_title(
                f"({string.ascii_lowercase[i]}) {foram_fullnames[i]} {mean:.2E}", pad=10
            )

        cbar = fig.colorbar(p, ax=axes, orientation="horizontal", pad=0.05, shrink=0.7)
        cbar.minorticks_on()
        cbar.set_label(r"POC export at 80.8 m (mmol C m$^{-3}$ d$^{-1})$", size=12)

        return p

    def plot_abundance(self):
        """
        quick wrapper function to plot relative abundance for 4 foram groups
        """

        fig, axes = plt.subplots(
            nrows=2,
            ncols=2,
            figsize=(10, 7),
            subplot_kw=dict(projection=ccrs.EckertIV()),
        )

        varlst = list(foram_names().keys())
        foram_fullnames = list(foram_names().values())

        for i, ax in enumerate(axes.flat):
            vdata = self.select_foram(varlst[i]).export_c().proportion().array
            mean = np.nanmean(vdata) * 100
            p = plot_genie(ax=ax, data=vdata, vmin=0, vmax=1)
            ax.set_title(
                f"({string.ascii_lowercase[i]}) {foram_fullnames[i]} {mean:.2f}%",
                pad=10,
            )

        cbar = fig.colorbar(p, ax=axes, orientation="horizontal", pad=0.05, shrink=0.7)
        cbar.minorticks_on()
        cbar.set_label("Relative abundance", size=12)

        return p

    def barplot_comparison(self, *args, **kwargs) -> plt.axes:

        """
        Overview barplot of biomass and POC export compared to observed data
        :returns: matplotlib axes object
        """

        fname = list(foram_names().values())
        foram_abbrev = foram_names().keys()

        model_biomass_mean = [
            self.select_foram(i).biomass_c().nanmean() for i in foram_abbrev
        ]
        model_export_mean = [
            self.select_foram(i).export_c().nanmean() for i in foram_abbrev
        ]
        model_biomass_se = [self.select_foram(i).biomass_c().se() for i in foram_abbrev]
        model_export_se = [self.select_foram(i).export_c().se() for i in foram_abbrev]
        model_export_sum = [
            self.select_foram(i).export_c().sum().magnitude for i in foram_abbrev
        ]
        model_biomass_sum = [
            self.select_foram(i).biomass_c().sum().magnitude for i in foram_abbrev
        ]
        obs_biomass_mean = obs_stat_bysource("tow", *args, **kwargs).loc[:, "mean"]
        obs_biomass_se = obs_stat_bysource("tow", *args, **kwargs).loc[:, "se"]
        obs_export_mean = obs_stat_bysource("trap", *args, **kwargs).loc[:, "mean"]
        obs_export_se = obs_stat_bysource("trap", *args, **kwargs).loc[:, "se"]

        model_biomass_sum = DataFrame({"group": fname, "value": model_biomass_sum})
        model_export_sum = DataFrame({"group": fname, "value": model_export_sum})

        data_to_plot = [
            [
                [
                    model_biomass_mean,
                    model_biomass_se,
                    obs_biomass_mean,
                    obs_biomass_se,
                ],
                [model_export_mean, model_export_se, obs_export_mean, obs_export_se],
            ],
            [model_biomass_sum, model_export_sum],
        ]
        ## Plot starts
        fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)
        bar_width = 0.3

        # The x position of bars
        x = np.arange(4)
        xlabels = [w.replace(" ", "\n") for w in fname]

        for i in range(2):
            for j in range(2):
                axes[i, j].yaxis.set_minor_locator(AutoMinorLocator(4))
                axes[i, j].set_axisbelow(True)
                axes[i, j].yaxis.grid(color="gray", linestyle="dashed")
                if i == 0:
                    axes[i, j].bar(
                        x - bar_width / 2,
                        data_to_plot[i][j][0],
                        width=bar_width,
                        color=sns.color_palette("Set1")[0],
                        edgecolor="black",
                        yerr=data_to_plot[i][j][1],
                        capsize=7,
                        label="model",
                    )

                    axes[i, j].bar(
                        x + bar_width / 2,
                        data_to_plot[i][j][2],
                        width=bar_width,
                        color=sns.color_palette("Set1")[1],
                        edgecolor="black",
                        yerr=data_to_plot[i][j][3],
                        capsize=7,
                        label="obs",
                    )

                    axes[i, j].set_xticks(x)
                    axes[i, j].set_xticklabels(xlabels, rotation=45, ha="right")
                    axes[i, j].legend()
                else:
                    sns.barplot(
                        data=data_to_plot[i][j],
                        x="group",
                        y="value",
                        ax=axes[i, j],
                        edgecolor="black",
                        palette="deep",
                    )
                    axes[i, j].set_xticklabels(xlabels, rotation=45, ha="right")
                    axes[i, j].set_xlabel("")
                    axes[i, j].set_xticks(x)
                    set_sns_barwidth(axes[i, j], bar_width)

        axes[0, 0].set_ylabel(r"mmol C m$^{-3}$")
        axes[0, 1].set_ylabel(r"mmol C m$^{-2}$ d$^{-1}$")

        axes[0, 0].set_title("(a)    global biomass mean/se", loc="left")
        axes[0, 1].set_title("(b)    global POC export mean/se", loc="left")
        axes[1, 0].set_title("(c)    globally integrated biomass", loc="left")
        axes[1, 1].set_title("(d)    globally integrated EP", loc="left")

        axes[1, 0].set_ylabel("Gt C")
        axes[1, 1].set_ylabel(r"Gt C yr$^{-1}$")

        fig.tight_layout()

        return axes

    def cwm(self, trait: str, abundance="relative"):
        """community weighted mean, also representing position of community centroid (i,j),
        usually weighted in relative abundance
        """

        fg_dict = {"symbiont": ["sn", "ss"], "spine": ["bs", "ss"]}
        trait_dict = {"symbiont": [0.5, 1], "spine": [1, 1]}
        fg = fg_dict[trait]
        tr = trait_dict[trait]

        if abundance == "absolute":
            weighted_trait = (
                self.select_foram(fg[0]).biomass_c() * tr[0]
                + self.select_foram(fg[1]).biomass_c() * tr[1]
            )
            sum_abundance = (
                self.select_foram(fg[0]).biomass_c()
                + self.select_foram(fg[1]).biomass_c()
            )
            cwm = weighted_trait / sum_abundance
        elif abundance == "relative":
            cwm = (
                self.select_foram(fg[0]).export_c().proportion() * tr[0]
                + self.select_foram(fg[1]).export_c().proportion() * tr[1]
            )

        return cwm

    def fspec(self, *args, **kwargs):
        "Functional specialisation, distance from the origin (non-symbiont non-spinose)"
        fspec = distance(
            self.cwm(trait="spine", *args, **kwargs),
            self.cwm(trait="symbiont", *args, **kwargs),
        )

        x = GenieArray()
        x.array = fspec
        return x

    def fored(self):
        "functional over redundancy (only used for modern), unit: redundant species"
        fred = (
            self.select_foram("bn").export_c().proportion() * 5.25 / 47
            + self.select_foram("ss").export_c().proportion() * 11.25 / 47
        )

        return fred

    # def fdis(self, *args, **kwargs):
    # "functional dispersion"
    # weight (relative abundance/biomass)
    # weight = self.select_foram("bn").biomass_c()
    # each group's distance to centroid
    # dist =
    # weighted sum


class ForamType(PlanktonType):

    def __init__(self, foram_type, model_path):
        self.model_path = model_path
        self.foram_type = foram_type
        self.pft_n = foram_names()[self.foram_type][0]

    def biomass(self, element="C"):
        return ForamBiomass(pft_n=self.pft_n, element = element, model_path=self.model_path)

    def export(self, element="C"):
        return ForamExport(pft_n=self.pft_n, element = element, model_path=self.model_path)

    def calcite(self):
        export_C_var = self.export(element="C").full_varstr
        return ForamCalcite(var=export_C_var, model_path=self.model_path)

    def relative_abundance(self, element="C"):
        export_var = self.export(element=element).full_varstr
        return ForamAbundance(var=export_var, model_path=self.model_path)


class ForamBiomass(PlanktonBiomass):
    obs = "net"
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ForamExport(PlanktonExport):
    obs = "trap"
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ForamCalcite(GenieVariable):

    unit = "mmol m$^-2$ d$^-1$"    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _set_array(self):
        """
        convert POC to Calcite (in mmol m-2 d-1) given POC:PIC:CaCO3 mol ratio = 100:36:36 (mass ratio = 100:36:300)
        """
        array = super()._set_array() * 0.36
        return array

    @ureg.with_context("bgc")
    def sum(self):
        CaCO3 = molecular_weight("CaCO3")
        c = self.uarray().to_base_units()
        v = GENIE_grid_area().to_base_units()
        s = c * v
        s = (
            s.to("mol d^-1")
            .to("g d^-1", "bgc", mw=CaCO3 * ureg("g/mol"))
            .to("Gt yr^-1")
        )

        return np.nansum(s)

class ForamAbundance(GenieVariable):
    obs = "coretop"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _total_foram(self):
        "if total foram is ptf No.16 to 19"
        # get the source data
        gm = GenieModel(model_path = self.model_path)
        path2nc =gm._auto_find_path(var=self.var)
        src_data = gm._open_nc(path2nc)

        # get all the foram variables
        variable_template = self.var[:-2]
        foram_variables = [variable_template + str(i) for i in range(16, 20)]

        # use foram variables to get data
        target_data = []
        for i in foram_variables:
            target_data.append(src_data[i])

        total_foram = reduce(np.add, target_data)

        return total_foram

    def _set_array(self):
        # one foram
        one_foram = super()._set_array()

        # total foram
        total_foram = self._total_foram()

        # ignore divided by 0
        # and set grid with total_foram == 0 to 0 instead of NA
        with np.errstate(divide="ignore", invalid="ignore"):
            proportion = np.divide(
                one_foram,
                total_foram,
                out=np.zeros_like(one_foram),
                where=total_foram != 0,
            )

        return proportion
