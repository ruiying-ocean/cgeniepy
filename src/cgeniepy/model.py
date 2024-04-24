from os.path import join
from typing import Union, List, Tuple
import re
import warnings
from pathlib import Path

import pandas as pd
from netCDF4 import Dataset
import xarray as xr
from scipy.ndimage import binary_dilation
import numpy as np

from .utils import file_exists
from .array import GriddedData


class GenieModel(object):
    """
    GenieModel is the interface to cGENIE output
    """

    def __init__(self, model_path: Union[str, List, Tuple], gemflag=None):
        """
        Initialise a GenieModel object with a path to cGENIE output directory

        :param model_path: the path to the cGENIE output directory
        :param gemflag: the cgenie component, default is ["biogem"]
        -------
        Example
        >>> from cgeniepy.model import GenieModel
        >>> model = GenieModel("path_to_GENIE_output")

        ## specify cgenie component
        >>> model = GenieModel("path_to_GENIE_output", gemflag=["biogem"])
        """

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
            warnings.warn("No gemflag is provided, use default gemflags: [biogem]")
            self.gemflag = ["biogem"]
        else:
            self.gemflag = gemflag
            ## if gemflag is a string, convert it to a list
            if isinstance(gemflag, str):
                self.gemflag = [gemflag]

        self.model_path = model_path
        self.ncvar_dict = self._ncvar_dict()
        self.tsvar_list = self.tsvar_list()

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

    def _ncvar_dict(self) -> dict:
        """
        Return all available variables and related biogem/ecogem (2d or 3d)
        NetCDF path for each model
        """

        ## initialise a dictionary, key: model_path, value: ncvar_list
        var_path = {}

        for gem in self.gemflag:
            for dim in ["2d", "3d"]:
                ## ignore ecogem 3d which is not available now (Oct 2023)
                if gem == "ecogem" and dim == "3d":
                    continue

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

        :param var: the name of the target variable
        """
        for path, value_list in self.ncvar_dict.items():
            if var in value_list:
                return path
        raise ValueError(f"variable {var} not found in the ncvar_dict")

    def _open_nc(self, nc_path):
        """Use xarray to open netcdf file

        :param nc_path: the path to the netcdf file
        """
        ## if path is a list of paths
        if not isinstance(nc_path, (list, tuple)):
            ## `open_dataset` to lazy load the data
            return xr.open_dataset(nc_path)
        else:
            datasets = [xr.open_dataset(file) for file in nc_path]
            combined_ds = xr.concat(datasets, "model")
            combined_ds.coords["model"] = self.model_path
            return combined_ds

    def get_var(self, var: Union[str, List, Tuple], attrs=None, mutable=True):
        """
        Get the data of target variable.
        A list of variables is supported as well.

        :param var: the name of target variable
        :param unit: the unit of target variable, usually provided in the model output

        Example
        ----------
        >>> from cgeniepy.model import GenieModel
        >>> model = GenieModel("path_to_GENIE_output")
        >>> po4 = model.get_var("ocn_PO4")
        """
        self.target_var = var

        ## if varstr is a string
        if isinstance(self.target_var, str):
            ## find the path to the netcdf file
            path2nc = self._lookup_ncpath(var=self.target_var)
            ## open the netcdf file
            array = self._open_nc(path2nc)[self.target_var]
        ## if varstr is a list/tuple of strings
        ## concat the data along the "variable" dimension
        elif isinstance(self.target_var, (list, tuple)):
            array_container = []
            for v in self.target_var:
                path2nc = self._lookup_ncpath(var=v)
                single_array = self._open_nc(path2nc)[v]
                array_container.append(single_array)

            array = xr.concat(array_container, "variable")
            array.name = "ensemble_variable"

        ## initialise GriddedData object
        if not attrs:
            return GriddedData(array, mutable=mutable, attrs=array.attrs)
        else:
            return GriddedData(array, mutable=mutable, attrs=attrs)
            

    def tsvar_list(self):
        """list all files biogem timeseries files

        If the model is an ensemble, it assumes that all models share the same
        file/data structure and only return the first model's timeseries list
        """

        if self.is_ensemble:
            biogem_path = join(self.model_path[0], "biogem")
        else:
            biogem_path = join(self.model_path, "biogem")

        ## return all files in the directory
        all_bg_files = Path(biogem_path)

        ## find out files with .res extension
        tsvar_list = [f for f in all_bg_files.glob("*.res")]

        return tsvar_list

    def get_ts(self, var: str):
        """
        read in time series output of GENIE

        :param var: the name of the target variable
        :return: a pandas DataFrame
        """

        if not self.is_ensemble:
            filename = f"biogem_series_{var}.res"            
            f = join(self.model_path, "biogem", filename)
            if not file_exists(f):
                raise ValueError(f"{f} does not exist")

            # Read the text file into a DataFrame
            with open(f, "r") as file:
                lines = file.readlines()

            # Remove " %" from the header and split it using " / "
            # not '/' which is also included in the isotope unit (e.g. 'd13C o/oo')
            header = [col.strip().lstrip("% ") for col in lines[0].split(" / ")]
            data = [line.split() for line in lines[1:]]
            df = pd.DataFrame(data, columns=header)

            ## convert to numeric
            df = df.apply(pd.to_numeric, errors="coerce")

            ## add a column of model name
            df["model"] = self.model_path.split("/")[-1]
            return df
        else:
            ## concatenate all models' data into one data frame
            df_list = []

            for path in self.model_path:
                ## do the same thing as above
                filename = f"biogem_series_{var}.res"
                f = join(path, "biogem", filename)
                if not file_exists(f):
                    raise ValueError(f"{f} does not exist")

                ## read text file into a DataFrame
                with open(f, "r") as file:
                    lines = file.readlines()

                ## Remove " %" from the header and split it using " / "
                header = [col.strip().lstrip("%") for col in lines[0].split(" / ")]
                data = [line.split() for line in lines[1:]]

                df = pd.DataFrame(data, columns=header)

                ## convert to numeric
                df = df.apply(pd.to_numeric, errors="coerce")

                df["model"] = path.split("/")[-1]
                df_list.append(df)

            ## concatenate by row
            all_df = pd.concat(df_list, axis=0)

            return all_df
        
    def get_diag_avg(self, filename):
        """print out the summary of a diagnostic file

        :param filename: the name of the diagnostic file

        Example
        ----------
        >>> from cgeniepy.model import GenieModel
        >>> model = GenieModel("path_to_GENIE_output")
        >>> model.get_diag_avg("biogem_year_09999_500_diag_GLOBAL_AVERAGE.res")
        """
        f = join(self.model_path, "biogem", filename)
        with open(f, 'r') as f:
            text = f.read()

        # Split data by lines, remove empty lines
        lines = [line.strip() for line in text.splitlines() if line.strip()]

        # Initialize an empty dictionary to store all data
        variable = []
        value1 = []
        value2 = []

        # Loop through each line
        for line in lines:
        # Split line by separator (colon)
            try:
                key, value = line.split(':', 1)
            except ValueError:  # Ignore lines without a colon separator
                continue

            # Remove leading/trailing whitespaces from key and value
            key = key.strip()
            value = value.strip()

            if '<->' in value:
                v1 = value.split('<->')[0].strip()
                v2 = value.split('<->')[1].strip()
            else:
                v1 = value
                v2 = None

            key = re.sub(r"\s*\.+$", "", key)

            variable.append(key)
            value1.append(v1)
            value2.append(v2)
    
        
        ## get the model year and make it a separate column
        year = value1[0]

        ## convert data to data frame
        df = pd.DataFrame({'year': year,
                           'variable': variable, 
                           'value1': value1, 
                           'value2': value2})
        
        ## remove the duplicated first row
        df = df.iloc[1:]

        return df

    def grid_mask(self):
        """
        cGENIE continent mask array (by setting zero values),
        either calculated from existing data or using biogem.grid_mask
        """
        try:
            grid_mask = self.get_var("grid_mask").data
            return grid_mask
        except ValueError:
            print("grid_mask not found!")

    def grid_category(self):
        """an alogirthm to define surface grid catogories depending on the land-sea mask
        0: coastal sea
        1: land
        2: open ocean
        """

        is_sea = np.isnan(self.grid_mask())
        land_around = binary_dilation(is_sea)
        grid_catogories = xr.where(land_around, np.logical_and(land_around, is_sea), 2)

        return grid_catogories

    def grid_area(self):
        "return the grid area array in used in this model experiment, unit: m2"
        ## read grid_area from biogem
        print("grid area returned in the unit of 'm2'")
        return self.get_var("grid_area")
    
    def grid_mask_3d(self):
        "return the 3d mask of the grid used in this model experiment"
        return self.get_var("grid_mask_3d")
    
    def grid_topo(self):
        "return the topography of the grid used in this model experiment"
        return self.get_var("grid_topo")
    
    def grid_zt_edges(self):
        "return the depth edges of the grid used in this model experiment"
        return self.get_var("zt_edges")
    
    def grid_lat_edges(self):
        "return the latitude edges of the grid used in  this model experiment"
        return self.get_var("lat_edges")
    
    def grid_lon_edges(self):
        "return the longitude edges of the grid used in  this model experiment"
        return self.get_var("lon_edges")
    
    def grid_zt_depths(self):
        "return the depth of the grid used in this model experiment"
        return self.get_var("grid_dD")
    
    def grid_volume(self):
        "return the grid volume array (3d) in m3"
        grid_area = self.grid_area()  ## m2
        depth = self.grid_zt_depths() ## m
        ocn_mask = self.grid_mask_3d()
        try:
            grid_volume = depth * grid_area * ocn_mask 
            grid_volume.data.attrs["units"] = "m$^{3}$"
            grid_volume.data.attrs["long_name"] = "grid volume"
            print("grid volume calculated in the unit of 'm3'")
        except ValueError:
            print(
                "Depth array not found! Please ensure 3d data is exported in the model!"
            )

        return grid_volume
