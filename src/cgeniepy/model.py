from os.path import join
from typing import Union, List, Tuple
from pathlib import Path


import numpy as np
import pandas as pd
from scipy.stats import sem, t
import regionmask
from netCDF4 import Dataset
import xarray as xr

from . import Q_
from .grid import (
    GENIE_grid_area,
    normalise_GENIE_lon,
    GENIE_grid_mask,
    GENIE_grid_vol,
    mask_Arctic_Med
)
from .utils import file_exists
from .chem import pure_unit, format_unit
from .scores import ModelSkill
from .utils import remove_outliers
from .array import GenieArray


# Add method to read in time series and parameter table

class GenieModel(object):
    
    def __init__(self, model_path):
        self.model_path = model_path
        
        if Path(self.model_path).is_dir():
            pass
        else:
            raise ValueError(f"{self.model_path} is not a valid directory")

    def _get_ncpath(self, gem="ecogem", dim="2d"):
        "Extending model_path to the NETCDF file path"

        model_path = self.model_path
        nc_file = f"fields_{gem}_{dim}.nc"
        nc_path = join(model_path, gem, nc_file)
        if file_exists(nc_path):
            return nc_path

    def ncvar_list(self):
        "return all available variables in biogem/ecogem (2d or 3d)"

        ## initialise a dictionary, key: model_path, value: ncvar_list
        var_path = {}
        
        for gem in ["biogem", "ecogem"]:
            for dim in ["2d", "3d"]:
                ## ignore ecogem 3d
                if gem == "ecogem" and dim == "3d": continue
                nc_path = self._get_ncpath(gem, dim)
                nc = Dataset(nc_path, "r")
                ncvar_list = list(nc.variables.keys())
                nc.close()
                ## add key, value to the dictionary
                var_path[nc_path] = ncvar_list
                
        return var_path

    def tsvar_list(self):
        "return all available variables in biogem timeseries"

        ## biogem/xxx.res
        biogem_path = join(self.model_path, "biogem")

        ## return all files in the directory
        biogem_files = Path(biogem_path)

        ## find out files with .res extension
        tsvar_list = [f for f in biogem_files.glob("*.res")]
                
        return tsvar_list
        
    def _lookup_ncpath(self, var):
        "get netcdf path according to the given variable"
        for path, value_list in self.ncvar_list().items():
            if var in value_list:
                return path        
        raise ValueError(f"variable {var} not found in the model")

    def _open_nc(self, path):
        "Use xarray to open netcdf file"
        ## if path is a list of paths
        if isinstance(path, list):
            return xr.open_mfdataset(path, parallel=True)
        else:
            return xr.open_dataset(path)

    def get_ts(self, filename, col_name=None):

        ## check if file exists
        filename = join(self.model_path, "biogem", filename)
        
        if not file_exists(filename):
            raise ValueError(f"{filename} does not exist")

        # Read the text file into a DataFrame
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Remove "%" from the header and split it using "/"
        header = [col.strip().lstrip('%') for col in lines[0].split('/')]
        data = [line.split() for line in lines[1:]]
        
        df= pd.DataFrame(data, columns=header)
        ## convert to numeric
        df = df.apply(pd.to_numeric, errors='coerce')
        
        return df
    
    def get_var(self, var: Union[str, List, Tuple]):
        """
        get the data of target variable(s)

        variable name and model path will be used to find
        the corresponding netcdf file
        """
        self.target_var = var

        ## if varstr is a string
        if isinstance(self.target_var, str):
            path2nc =self._lookup_ncpath(var=self.target_var)
            array = self._open_nc(path2nc)[self.target_var]
        ## if varstr is a list/tuple of strings
        elif isinstance(self.target_var, (list, tuple)):
            array_container = []
            for v in self.target_var:
                path2nc =self._lookup_ncpath(var=v)
                single_array = self._open_nc(path2nc)[v]
                array_container.append(single_array)            
            
            array = xr.concat(array_container, "variable")
            array.name = "ensemble"

        ## initialise GenieArray
        target_data = GenieArray()
        target_data.array = array
        target_data.array.attrs['units'] = format_unit(target_data.array.attrs['units'] )
        
        
        return target_data

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
