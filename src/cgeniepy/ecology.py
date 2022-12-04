from pandas import read_fwf
import numpy as np
from .core import GenieModel, GenieVariable
from .chem import molecular_weight
from .grid import GENIE_grid_vol, GENIE_grid_area
from . import ureg


class EcoModel(GenieModel):

    "EcoModel is an subclass of GenieModel, as EcoGENIE to cGENIE"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def eco_pars(self):
        """
        return ecophysiological parameter table
        """
        path = f"{self.model_path}/ecogem/Plankton_params.txt"
        df = read_fwf(path)
        return df

    def pft_num(self):
        "the total number of plankton functional type"
        full_lst = list(self.get_vars())
        name_lst = [x for x in full_lst if "eco2D_Export_C" in x]
        n = len(name_lst)
        return n

    def select_pft(self, pft_n):
        "pft can be an integer or a index list"
        return PlanktonType(pft_n=pft_n, model_path=self.model_path)


class PlanktonType:
    def __init__(self, pft_n, model_path):
        self.pft_n = pft_n
        self.model_path = model_path

    def biomass(self, element="C"):
        return PlanktonBiomass(
            pft_n=self.pft_n, element=element, model_path=self.model_path
        )

    def export(self, element="C"):
        return PlanktonExport(
            pft_n=self.pft_n, element=element, model_path=self.model_path
        )

    # def presence(self, x, tol=1e-8):
    #     """
    #     to determine whethere a functional group present or not
    #     :param tol: threshold of biomass (mmol C/m3)
    #     """
    #     if np.isnan(x):
    #         return x
    #     elif x >= tol:
    #         return 1
    #     else:
    #         return 0

    # def pft_richness(self):
    #     """
    #     plankton functional group richness, note it is different from species richness,
    #     because there are more species in low size classes (i.e., body size-species richness relationship)
    #     """
    #     vfunc = np.vectorize(self.presence)
    #     total_sp = np.zeros((36, 36))
    #     n = self.pft_num()

    #     for i in range(n):
    #         name = f"eco2D_Plankton_C_0{i+1:02d}"
    #         # select
    #         arr = self.select_var(name).pure_array()
    #         # conditional mask
    #         sp = vfunc(arr)
    #         # sum
    #         total_sp += sp

    #     x = GenieArray()
    #     x.array = total_sp
    #     return x


class PlanktonBiomass(GenieVariable):

    bgc_prefix = "Plankton"
    unit = "mmol m$^-3$"

    def __init__(self, model_path, pft_n, element):
        self.pft_n = pft_n
        self.element = element
        self.full_varstr = f"eco2D_{self.bgc_prefix}_{self.element}_{self.pft_n:03}"
        GenieVariable.__init__(self, model_path=model_path, var=self.full_varstr)

    def sum(self):
        "print in Tg, depending on the element"
        X_ = molecular_weight(self.element)
        c = self.uarray().to_base_units()
        v = GENIE_grid_vol().to_base_units()
        s = c * v
        s = s.to("mol").to("g", "chemistry", mw=X_ * ureg("g/mole")).to("Gt")
        return np.nansum(s)


class PlanktonExport(GenieVariable):

    unit = "mmol m$^-2$ d$^-1$"
    bgc_prefix = "Export"

    def __init__(self, model_path, pft_n, element):
        self.pft_n = pft_n
        self.element = element
        self.full_varstr = f"eco2D_{self.bgc_prefix}_{self.element}_{self.pft_n:03}"
        GenieVariable.__init__(self, model_path=model_path, var=self.full_varstr)

    def _set_array(self):
        # tested
        return super()._set_array() * 80.8

    @ureg.with_context("bgc")
    def sum(self):
        # concentration data
        c = self.uarray().to_base_units()
        # make volume in pint type
        v = GENIE_grid_area().to_base_units()
        # globall integrated value
        s = c * v

        # unit conversion
        C_ = molecular_weight(self.element)
        s = s.to("mol d^-1").to("g d^-1", "bgc", mw=C_ * ureg("g/mol")).to("Gt yr^-1")

        return np.nansum(s)
