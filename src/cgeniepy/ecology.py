import pandas as pd
import numpy as np
import re

from .model import GenieModel
from .data import foram_groups

class EcoModel(GenieModel):
    """
    EcoModel is an subclass of GenieModel, as EcoGENIE to cGENIE

    It facilitates the access of ecophysiological parameters and plankton variables
    """

    def __init__(self, *args, **kwargs):
        
        super().__init__(*args, **kwargs)        

        ## get the list of plankton variables
        ecogem2d_path = self._model_ncpath('ecogem', '2d')            
        self.eco_varlist = self.ncvar_dict()[ecogem2d_path]

        photo_par = self.eco_pars()['vmax_C'].to_numpy()
        graz_par = self.eco_pars()['max_graz_C'].to_numpy()

        autotrophy = (photo_par != 0)
        heterotrophy = (graz_par != 0)
        mixotrophy = np.logical_and(autotrophy, heterotrophy)

        phyto_idx = np.where(np.logical_and(autotrophy, ~mixotrophy))[0]
        zoo_idx = np.where(np.logical_and(heterotrophy, ~mixotrophy))[0]
        mixo_idx = np.where(mixotrophy)[0]

        ## convert to ECOGENIE format
        self.phyto_indices = [f'{i+1:03d}' for i in phyto_idx]
        self.zoo_indices = [f'{i+1:03d}' for i in zoo_idx]
        self.mixo_indices = [f'{i+1:03d}' for i in mixo_idx]
        
        ## some basic information (number of PFTs)
        self.phyto_n = len(self.phyto_indices)
        self.zoo_n = len(self.zoo_indices)
        self.mixo_n = len(self.mixo_indices)
        self.plank_n = np.sum([self.phyto_n, self.zoo_n, self.mixo_n])

    def eco_pars(self):
        """
        return ecophysiological parameter table
        """
        if self.is_ensemble:
            df_list = []
            for path in self.model_path:
                path = f"{path}/ecogem/Plankton_params.txt"
                df = pd.read_fwf(path)
                df_list.append(df)
            df_all = pd.concat(df_list)
            return df_all
        else:
            path = f"{self.model_path}/ecogem/Plankton_params.txt"
            df = pd.read_fwf(path)
            return df
            

    def get_pft(self, pft_index, bgc_prefix="Plankton", element="C"):
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
        model.get_pft(1, "Plankton", "C")

        ### get a list of PFT biomass data
        model.get_pft([1, 2, 3], "Plankton", "C")

        ### get all phytoplankton biomass data
        model.get_pft('phyto', "Plankton", "C")

        ### get all zooplankton biomass data
        model.get_pft('zoo', "Plankton", "C")

        ### get all plankton biomass data
        model.get_pft('all', "Plankton", "C")
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

    def get_foram(self, foram_type, *args, **kwargs):

        "a more specific version of get_pft, to select foram"
        
        ## convert foram type to pft index
        if isinstance(foram_type, (list, tuple)):
            self.pft_index = [foram_groups().get(foram)[0] for foram in foram_type]
        else:
            self.pft_index = foram_groups()[foram_type][0]

        return self.get_pft(self.pft_index, *args, **kwargs)
