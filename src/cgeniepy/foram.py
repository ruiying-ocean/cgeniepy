# import matplotlib.pyplot as plt
# from matplotlib.ticker import AutoMinorLocator
# import string
# from pandas import DataFrame

import numpy as np
from functools import reduce
from . import ureg
from .chem import molecular_weight

from .ecology import GenieModel, PlanktonType, PlanktonBiomass, PlanktonExport
from .core import GenieVariable
from .data import foram_names, obs_data
from .scores import ModelSkill

# from .plot import plot_genie

# from .scores import quick_mscore, quick_rmse, quick_cos_sim, quick_corr
# from .utils import set_sns_barwidth, distance
from .grid import GENIE_grid_area

# observation_data, functional diversity, more plot options

class ForamModel(GenieModel):
    """A further customized GenieModel subclass"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def select_foram(self, foram_type, *args, **kwargs):
        "a optimised version of select_var, can be int or a list/tuple"
        return ForamType(foram_type=foram_type, model_path=self.model_path,  *args, **kwargs)


class ForamType(PlanktonType):

    def __init__(self, foram_type, model_path):
        self.model_path = model_path
        self.foram_type = foram_type

    def biomass(self, element="C",  *args, **kwargs):
        return ForamBiomass(
            foram_type=self.foram_type,
            element=element,
            model_path=self.model_path,
            *args, **kwargs
        )

    def export(self, element="C",  *args, **kwargs):
        return ForamExport(
            foram_type=self.foram_type,
            element=element,
            model_path=self.model_path,
             *args, **kwargs
        )

    def calcite(self,  *args, **kwargs):
        export_C_var = self.export(element="C").full_varstr
        return ForamCalcite(foram_type = self.foram_type, var=export_C_var, model_path=self.model_path,  *args, **kwargs)

    def relative_abundance(self, element="C",  *args, **kwargs):
        export_var = self.export(element=element).full_varstr
        return ForamAbundance(foram_type = self.foram_type, var=export_var, model_path=self.model_path,  *args, **kwargs)


class ForamBiomass(PlanktonBiomass):
    obs = "net"

    def __init__(self, foram_type, element, model_path,  *args, **kwargs):
        self.foram_type = foram_type
        # convert foram type to pft index
        if isinstance(self.foram_type, list) or isinstance(self.foram_type, tuple):
            self.pft_index = [foram_names().get(foram)[0] for foram in self.foram_type]
        else:
            self.pft_index = foram_names()[self.foram_type][0]
        # pass pft_index to father class
        super().__init__(pft_index = self.pft_index,
                         element=element,
                         model_path=model_path,  *args, **kwargs)

    def compare_obs(self, **kwargs):
        if "obs" in kwargs:
            data = obs_data(source = kwargs["obs"], var=self.foram_type)
        else:
            data = obs_data(source = self.obs, var=self.foram_type)

        if "mask_MedArc" in kwargs:
            return ModelSkill(model=self.pure_array(), observation=data, mask_MedArc=kwargs["mask_MedArc"])
        else:
            return ModelSkill(model=self.pure_array(), observation=data)


class ForamExport(PlanktonExport):
    obs = "trap"

    def __init__(self, foram_type, element, model_path,  *args, **kwargs):
        self.foram_type = foram_type
        # convert foram type to pft index
        if isinstance(self.foram_type, list) or isinstance(self.foram_type, tuple):
            self.pft_index = [foram_names().get(foram)[0] for foram in self.foram_type]
        else:
            self.pft_index = foram_names()[self.foram_type][0]
        # pass pft_index to father class
        super().__init__(pft_index = self.pft_index,
                         element=element,
                         model_path=model_path,  *args, **kwargs)

    def compare_obs(self, **kwargs):
        if "obs" in kwargs:
            data = obs_data(source = kwargs["obs"], var=self.foram_type)
        else:
            data = obs_data(source = self.obs, var=self.foram_type)

        if "mask_MedArc" in kwargs:
            return ModelSkill(model=self.pure_array(), observation=data, mask_MedArc=kwargs["mask_MedArc"])
        else:
            return ModelSkill(model=self.pure_array(), observation=data)


class ForamCalcite(GenieVariable):

    unit = "mmol m$^-2$ d$^-1$"

    def __init__(self, foram_type, *args, **kwargs):
        self.foram_type = foram_type
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

    def __init__(self, foram_type, *args, **kwargs):
        self.foram_type = foram_type
        super().__init__(*args, **kwargs)

    def _total_foram(self):
        "if total foram is ptf No.16 to 19"
        # get the source data
        gm = GenieModel(model_path=self.model_path)
        path2nc = gm._auto_find_path(var=self.var)
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

    def compare_obs(self, **kwargs):
        if "obs" in kwargs:
            data = obs_data(source = kwargs["obs"], var=self.foram_type)
        else:
            data = obs_data(source = self.obs, var=self.foram_type)

        if "mask_MedArc" in kwargs:
            return ModelSkill(model=self.pure_array(), observation=data, mask_MedArc=kwargs["mask_MedArc"])
        else:
            return ModelSkill(model=self.pure_array(), observation=data)
