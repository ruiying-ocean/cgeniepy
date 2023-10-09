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

    def num_of_pft(self):
        "the total number of plankton functional type"
        full_lst = list(self.get_vars())
        name_lst = [x for x in full_lst if "eco2D_Export_C" in x]
        n = len(name_lst)
        return n

    def select_pft(self, pft_index):
        "pft can be an integer or a index list"
        return PlanktonType(pft_index=pft_index, model_path=self.model_path)


class PlanktonType:
    def __init__(self, pft_index, model_path):
        self.pft_index = pft_index
        self.model_path = model_path

    def biomass(self, element="C", *args, **kwargs):
        return PlanktonBiomass(pft_index=self.pft_index, element=element, model_path=self.model_path, *args, **kwargs)

    def export(self, element="C", *args, **kwargs):
        return PlanktonExport(pft_index=self.pft_index, element=element, model_path=self.model_path, *args, **kwargs)

    
class PlanktonBiomass(GenieVariable):

    bgc_prefix = "Plankton"
    unit = "mmol m$^-3$"

    def __init__(self, model_path, pft_index, element, *args, **kwargs):
        self.pft_index = pft_index
        self.element = element

        # single index
        if isinstance(self.pft_index, int) or isinstance(self.pft_index, str):
            self.full_varstr = f"eco2D_{self.bgc_prefix}_{self.element}_{self.pft_index:03}"
        # multiple indices
        elif isinstance(self.pft_index, list) or isinstance(self.pft_index, tuple):
            ## create a list of variable names
            self.full_varstr = []
            for i in pft_index:
                self.full_varstr.append(f"eco2D_{self.bgc_prefix}_{self.element}_{i:03}")

        GenieVariable.__init__(self, model_path=model_path, var=self.full_varstr, *args, **kwargs)

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

    def __init__(self, model_path, pft_index, element, *args, **kwargs):
        self.pft_index = pft_index
        self.element = element

        # single index
        if isinstance(self.pft_index, int) or isinstance(self.pft_index, str):
            self.full_varstr = f"eco2D_{self.bgc_prefix}_{self.element}_{self.pft_index:03}"

        # multiple indices
        elif isinstance(self.pft_index, list) or isinstance(self.pft_index, tuple):
            self.full_varstr = []
            for i in pft_index:
                self.full_varstr.append(f"eco2D_{self.bgc_prefix}_{self.element}_{i:03}")

        GenieVariable.__init__(self, model_path=model_path, var=self.full_varstr, *args, **kwargs)

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
