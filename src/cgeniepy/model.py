from os.path import join
from typing import Union, List, Tuple
from pathlib import Path

import pandas as pd
from netCDF4 import Dataset
import xarray as xr

from .grid import (
    GENIE_grid_area,
    GENIE_grid_vol,
)
from .utils import file_exists
from .chem import format_unit
from .array import GenieArray


# Add method to read in time series and parameter table

class GenieModel(object):
    """
    GenieModel is the fundamental class for users to access cGENIE output

    Initialise a GenieModel object with a path to cGENIE output directory
    -------
    Example
    >>> from cgeniepy.model import GenieModel
    >>> model = GenieModel("path_to_GENIE_output")
    """
    
    def __init__(self, model_path: Union[str, List, Tuple], gemflag=None):
        
        ## check if model_path is a valid directory
        if isinstance(model_path, (list, tuple)):
            self.is_ensemble = True
            for path in model_path:
                if not Path(path).is_dir():
                    raise ValueError(f"{path} is not a valid directory")
        elif isinstance(model_path, str):
            self.is_ensemble = False
            if not Path(model_path).is_dir():
                raise ValueError(f"{model_path} is not a valid directory")


        if not gemflag:
            print("No gemflag is provided, assuming the model includes biogem and ecogem")
            self.gemflag = ["biogem", "ecogem"]
        else:
            self.gemflag = gemflag

        
        self.model_path = model_path

    def _model_ncpath(self, gem="ecogem", dim="2d"):
        """
        get the path of model NetCDF output of target gem and dimension
        if the model is an ensemble, return a list of paths
        """

        nc_file = f"fields_{gem}_{dim}.nc"

        if not self.is_ensemble:
            model_path = self.model_path
            nc_path = join(model_path, gem, nc_file)
            if file_exists(nc_path):
                return nc_path
        else:
            nc_paths = []
            for path in self.model_path:
                model_path = path
                nc_path = join(model_path, gem, nc_file)
                if file_exists(nc_path):
                    nc_paths.append(nc_path)
            ## make `nc_paths` hashable
            nc_paths = tuple(nc_paths)
            return nc_paths            
            

    def ncvar_dict(self) -> dict:
        """
        Return all available variables and related biogem/ecogem (2d or 3d)
        NetCDF path for each model        
        """

        ## initialise a dictionary, key: model_path, value: ncvar_list
        var_path = {}

            
        for gem in self.gemflag:
            for dim in ["2d", "3d"]:
                ## ignore ecogem 3d which is not available now (Oct 2023)
                if gem == "ecogem" and dim == "3d": continue

                if self.is_ensemble:
                    ## use one model to get the nc_path
                    nc_path = self._model_ncpath(gem, dim)[0]
                    all_nc_path = self._model_ncpath(gem, dim)
                else:
                    nc_path = self._model_ncpath(gem, dim)

                nc = Dataset(nc_path, "r")
                ncvar_list = list(nc.variables.keys())
                nc.close()

                ## add key, value to the dictionary
                if self.is_ensemble:
                    var_path[all_nc_path] = ncvar_list
                else:
                    var_path[nc_path] = ncvar_list
                
        return var_path

    def _lookup_ncpath(self, var):
        """
        find the netcdf path of the target variable

        if more than one model satisfies the condition,
        return all model's NetCDF paths
        """
        for path, value_list in self.ncvar_dict().items():
            if var in value_list:
                return path
        raise ValueError(f"variable {var} not found in the ncvar_dict")

    def _open_nc(self, nc_path):
        "Use xarray to open netcdf file"
        ## if path is a list of paths
        if not isinstance(nc_path, (list, tuple)):
            ## `open_dataset` to lazy load the data
            return xr.open_dataset(nc_path)
        else:            
            datasets = [xr.open_dataset(file) for file in nc_path]
            combined_ds  = xr.concat(datasets, "model")
            return combined_ds
    
    def get_var(self, var: Union[str, List, Tuple]):
        """
        Get the data of target variable. 
        A list of variables is supported as well.
        """
        self.target_var = var

        ## if varstr is a string
        if isinstance(self.target_var, str):
            ## find the path to the netcdf file
            path2nc =self._lookup_ncpath(var=self.target_var)
            ## open the netcdf file
            array = self._open_nc(path2nc)[self.target_var]
        ## if varstr is a list/tuple of strings
        ## concat the data along the "variable" dimension
        elif isinstance(self.target_var, (list, tuple)):
            array_container = []
            for v in self.target_var:
                path2nc =self._lookup_ncpath(var=v)
                single_array = self._open_nc(path2nc)[v]
                array_container.append(single_array)            
            
            array = xr.concat(array_container, "variable")
            array.name = "ensemble_variable"

        ## initialise GenieArray
        target_data = GenieArray()
        target_data.array = array
        target_data.array.attrs['units'] = format_unit(target_data.array.attrs['units'] )
                
        return target_data
    
    def tsvar_list(self):
        """list all files biogem timeseries files

        If the model is an ensemble, it assumes that all models share the same 
        file/data structure and only return the first model's timeseries list
        """
        
        if self.is_ensemble:
            biogem_path = join(self.model_path[0],"biogem")
        else:
            biogem_path = join(self.model_path, "biogem")

        ## return all files in the directory
        all_bg_files = Path(biogem_path)

        ## find out files with .res extension
        tsvar_list = [f for f in all_bg_files.glob("*.res")]
                
        return tsvar_list

    def get_ts(self, filename: str):
        """
        read in time series output of GENIE

        :param filename: the name of the time series file
        :return: a pandas DataFrame
        """
        
        if not self.is_ensemble:
            f = join(self.model_path, "biogem", filename)
            if not file_exists(f): raise ValueError(f"{f} does not exist")

            # Read the text file into a DataFrame
            with open(f, 'r') as file:
                lines = file.readlines()

            # Remove " %" from the header and split it using "/"
            header = [col.strip().lstrip('% ') for col in lines[0].split(' / ')]
            data = [line.split() for line in lines[1:]]
            df= pd.DataFrame(data, columns=header)

            ## convert to numeric
            df = df.apply(pd.to_numeric, errors='coerce')
            
            df['model'] = self.model_path.split('/')[-1]
            return df
        else:
            ## concatenate all models' data into one data frame
            df_list = []
            for path in self.model_path:
                f = join(path, "biogem", filename)
                if not file_exists(f): raise ValueError(f"{f} does not exist")
                with open(f, 'r') as file:
                    lines = file.readlines()
                header = [col.strip().lstrip('%') for col in lines[0].split('/')]
                data = [line.split() for line in lines[1:]]
                df= pd.DataFrame(data, columns=header)
                df = df.apply(pd.to_numeric, errors='coerce')
                df['model'] = path.split('/')[-1]
                df_list.append(df)
            ## concatenate by row
            all_df = pd.concat(df_list, axis=0)
 
            return all_df

    def _grid_mask(self, Arctic=True, Med=True):
        """
        cGENIE continent mask array (by setting zero values),
        either calculated from existing data or using biogem.grid_mask
        """
        try:
            grid_mask = self.get_var("grid_mask").array
            if Arctic: grid_mask[34:36, :] = 0
            if Med: grid_mask[27:30, 25:30] = 0
            return grid_mask
        except ValueError:
            print("grid_mask not found!")

    def _marine_area(self):
        "grid area array in km2"
        grid_mask = self._grid_mask()
        grid_area = GENIE_grid_area()
        mask_area = grid_area * grid_mask

        return mask_area

    def _marine_volume(self):
        "grid volume array in km3"
        grid_mask = self._grid_mask()
        grid_volume = GENIE_grid_vol()
        mask_volume = grid_volume * grid_mask

        return mask_volume

    def diff(self, model2compare, var):

        if isinstance(model2compare, GenieModel):
            B = model2compare
        else:
            B = GenieModel(model2compare)

        diff = self.get_var(var) - B.get_var(var)
        return diff

    def div(self, model2compare, var):
        B = GenieModel(model2compare)
        diff = self.get_var(var) / B.get_var(var)
        return diff
        
    
    def _run_method(self, method: str, *args, **kwargs):
        return getattr(self, method)(*args, **kwargs)
