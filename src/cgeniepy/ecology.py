from pandas import read_fwf
import numpy as np
import re

from .model import GenieModel
from .chem import molecular_weight
from .grid import GENIE_grid_vol, GENIE_grid_area
from . import ureg


class EcoModel(GenieModel):
    """
    EcoModel is an subclass of GenieModel, as EcoGENIE to cGENIE

    It facilitates the access of ecophysiological parameters and plankton variables
    """

    def __init__(self, *args, **kwargs):
        
        super().__init__(*args, **kwargs)
        
        ## get the list of plankton variables
        self.eco_varlist = self.ncvar_list()[self._get_ncpath('ecogem', '2d')]

        plankton_pattern = r'eco2D_Plankton_C_(\d{3})'
        self.plank_indices = [var.split("_")[3] for var in self.eco_varlist if re.search(plankton_pattern, var)]
        self.plank_n = len(self.plank_indices)

        phyto_pattern = r'eco2D_Plankton_Chl_(\d{3})'
        self.phyto_indices = [var.split("_")[3] for var in self.eco_varlist if re.search(phyto_pattern, var)]
        
        ## exclude phytoplankton to get zooplankton indices
        self.zoo_indices = list(set(self.plank_indices) - set(self.phyto_indices))

        ## not support mixotrophy yet

    def eco_pars(self):
        """
        return ecophysiological parameter table
        """
        path = f"{self.model_path}/ecogem/Plankton_params.txt"
        df = read_fwf(path)
        return df

    def select_pft(self, pft_index, bgc_prefix="Plankton", element="C"):
        """
        a variant of get_var, to select plankton functional type

        :param pft_index: the index of plankton functional type
        :param bgc_prefix: 'Plankton' or 'Export'
        :param element: 'C', 'Fe', 'P', 'Si', 'N'

        :return: a GenieArray object

        ----------------------------------------
        Example:

        from cgeniepy.ecology import EcoModel

        ### initialise a model object
        model = EcoModel("path_to_GENIE_output")
        
        ### get PFT-1 biomass data
        model.select_pft(1, "Plankton", "C")

        ### get a list of PFT biomass data
        model.select_pft([1, 2, 3], "Plankton", "C")

        ### get all phytoplankton biomass data
        model.select_pft('phyto', "Plankton", "C")

        ### get all zooplankton biomass data
        model.select_pft('zoo', "Plankton", "C")

        ### get all plankton biomass data
        model.select_pft('all', "Plankton", "C")
        """

        ### >>> preprocess the bgc_prefix
        ## if bgc_prefix is "Biomass", replace it with "Plankton"
        if bgc_prefix == "Biomass": bgc_prefix = "Plankton"        
        ## always capitalize the bgc_prefix
        bgc_prefix = bgc_prefix.capitalize()

        ### >>> process the pft_index
        ## if pft_index is "phyto", replace it with phyto_indices
        if pft_index == "phyto": pft_index = self.phyto_indices
        if pft_index == "zoo": pft_index = self.zoo_indices

        ## >>> construct variable name using pft_index, bgc_prefix and element
        if isinstance(pft_index, (int, str)):
            fullstring = f"eco2D_{bgc_prefix}_{element}_{pft_index:03}"
        elif isinstance(pft_index, (list, tuple)):
            fullstring = []
            for i in pft_index:
                fullstring.append(f"eco2D_{bgc_prefix}_{element}_{i:03}")
        else:
            raise ValueError("pft_index must be an integer or a list/tuple of integers")

        ## modify the fullstring if pft_index is "all"
        ## eco2D_Plankton_C_Total/eco2D_Plankton_Chl_Total
        if pft_index == "all":
            fullstring = f"eco2D_{bgc_prefix}_{element}_Total"
        
        return self.get_var(fullstring)

    # def select_foram(self, foram_type):
    #     "a more specific version of select_pft, to select foram"

    #     ## convert foram type to pft index
    #     if isinstance(foram_type, (list, tuple)):
    #         self.pft_index = [foram_names().get(foram)[0] for foram in self.foram_type]
    #     else:
    #         self.pft_index = foram_names()[self.foram_type][0]
