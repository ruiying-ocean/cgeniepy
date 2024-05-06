import pandas as pd
import numpy as np
import re

from cgeniepy.model import GenieModel

class EcoModel(GenieModel):
    """
    EcoModel is an subclass of GenieModel, as EcoGENIE to cGENIE

    It facilitates the access of ecophysiological parameters and plankton variables

    Initialise a EcoModel object with a path to cGENIE output directory


    Example
    -----------
    
    >>> from cgeniepy.ecology import EcoModel
    >>> model = EcoModel("path_to_GENIE_output")
    """

    def __init__(self, *args, **kwargs):
        
        super().__init__(*args, **kwargs)        

        ## get the list of plankton variables
        ecogem2d_path = self._model_ncpath('ecogem', '2d')            
        self.eco_varlist = self.ncvar_dict[ecogem2d_path]
        
        pattern = r'eco2D_Plankton_C_\d+'
        plank_vars = [var for var in self.eco_varlist if re.match(pattern, var)]
        self.plank_n = len(plank_vars)

        if self.is_ensemble:
            sub_ecopars = self.eco_pars().iloc[range(self.plank_n),]
        else:
            sub_ecopars = self.eco_pars()

        photo_par = sub_ecopars['vmax_C'].to_numpy()
        graz_par = sub_ecopars['max_graz_C'].to_numpy()

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

        ## check everything is correct
        if (self.phyto_n + self.zoo_n + self.mixo_n != self.plank_n):
            raise ValueError("The number of PFTs does not match the number of plankton variables")


    def eco_pars(self):
        """
        Get all the ecophysiological parameters used in the model

        :return: a pandas DataFrame object
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
            

    def get_pft(self, pft_index, prefix="Plankton", element="C"):
        """
        a variant of GenieModel's `get_var`, to select plankton functional type

        :param pft_index: the index of plankton functional type
        :param prefix: 'Plankton' or 'Export'
        :param element: 'C', 'Fe', 'P', 'Si', 'N'

        :return: a GenieArray object


        Example
        ----------------------------------------

        >>> from cgeniepy.ecology import EcoModel

        >>> ### initialise a model object
        >>> model = EcoModel("path_to_GENIE_output")
        
        >>> ### get PFT-1 carbon biomass data
        >>> model.get_pft(1, "Plankton", "C")

        >>> ### get a list of PFT carbon biomass data
        >>> model.get_pft([1, 2, 3], "Plankton", "C")

        >>> ### get all phytoplankton carbon biomass data
        >>> model.get_pft('phyto', "Plankton", "C")

        >>> ### get all zooplankton carbon biomass data
        >>> model.get_pft('zoo', "Plankton", "C")

        >>> ### get all plankton carbon biomass data
        >>> model.get_pft('all', "Plankton", "C")
        """

        ### >>> preprocess the prefix
        ## if prefix is "Biomass", replace it with "Plankton"
        if prefix == "Biomass": prefix = "Plankton"        
        ## always capitalize the prefix
        ## prefix = prefix.capitalize()

        ### >>> process the pft_index
        ## if pft_index is "phyto", replace it with phyto_indices
        if pft_index == "phyto": pft_index = self.phyto_indices
        if pft_index == "zoo": pft_index = self.zoo_indices
        if pft_index == "mixo": pft_index = self.mixo_indices

        ## >>> construct variable name using pft_index, prefix and element
        if isinstance(pft_index, (int, str)):
            fullstring = f"eco2D_{prefix}_{element}_{pft_index:03}"
        elif isinstance(pft_index, (list, tuple)):
            fullstring = []
            for i in pft_index:
                fullstring.append(f"eco2D_{prefix}_{element}_{i:03}")
        else:
            raise ValueError("pft_index must be an integer or a list/tuple of integers")

        ## modify the fullstring if pft_index is "all"
        ## eco2D_Plankton_C_Total/eco2D_Plankton_Chl_Total
        if pft_index == "all":
            fullstring = f"eco2D_{prefix}_{element}_Total"
        
        return self.get_var(fullstring)
