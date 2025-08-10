from typing import Union, List, Tuple
import re
import warnings
from pathlib import Path
import itertools

import pandas as pd
from netCDF4 import Dataset
import xarray as xr
from scipy.ndimage import binary_dilation
import numpy as np

from cgeniepy.utils import file_exists
from cgeniepy.array import GriddedData
from cgeniepy.table import ScatterData


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
            expanded_paths = []
            for path in model_path:
                expanded_path = Path(path).expanduser()
                if not expanded_path.is_dir():
                    raise ValueError(f"{path} is not a valid directory")
                expanded_paths.append(str(expanded_path))
            self.model_path = expanded_paths
        elif isinstance(model_path, str):
            self.is_ensemble = False
            expanded_path = Path(model_path).expanduser()
            if not expanded_path.is_dir():
                raise ValueError(f"{model_path} is not a valid directory")
            self.model_path = str(expanded_path)
        else:
            self.is_ensemble = False

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
        self.tsvar_list = self._tsvar_list()

    def __repr__(self):
        if self.is_ensemble:
            return f"GenieModel(ensemble of {len(self.model_path)} models)"
        else:
            poc_export = self.get_diag_avg(-1, is_index=True).query("variable == 'Total POC export'").iloc[:, 2].to_string(index=False)
            ocn_temp = self.get_diag_avg(-1, is_index=True).query("variable == 'Ocean temp'").iloc[:, 2].to_string(index=False)
            atm_co2 = self.get_diag_avg(-1, is_index=True).query("variable == 'Atmospheric pCO2'").iloc[:, 2].to_string(index=False)
            repr_str1 = f"GenieModel({self.model_path})\nLast Year Diagnostics:\n"
            repr_str2 = f"Ocean temp: {ocn_temp}\nTotal POC export: {poc_export}\natm pCO2: {atm_co2}"
            return repr_str1+repr_str2

    def _model_ncpath(self, gem="ecogem", dim="2d"):
        """
        get the path of model NetCDF output of target gem and dimension
        if the model is an ensemble, return a list of paths
        """
        nc_file = f"fields_{gem}_{dim}.nc"
        if not self.is_ensemble:
            model_path = Path(self.model_path).expanduser()
            nc_path = model_path / gem / nc_file
            if nc_path.exists():
                return str(nc_path)
        else:
            nc_paths = []
            for path in self.model_path:
                model_path = Path(path).expanduser()
                nc_path = model_path / gem / nc_file
                if nc_path.exists():
                    nc_paths.append(str(nc_path))
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

    def get_var(self, var: Union[str, List, Tuple], attrs=None):
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
            return GriddedData(array, attrs=array.attrs)
        else:
            return GriddedData(array, attrs=attrs)
            

    def _tsvar_list(self):
        """list all files biogem timeseries files

        If the model is an ensemble, it assumes that all models share the same
        file/data structure and only return the first model's timeseries list
        """

        if self.is_ensemble:
            biogem_paths = []
            for path in self.model_path:
                biogem_paths.append(Path(path).expanduser() / "biogem")

            ## return all files in all directories
            tsvar_list = [list(Path(path).glob("*.res")) for path in biogem_paths]
            ## unnest the list using itertools.chain
            tsvar_list = list(itertools.chain(*tsvar_list))
            return tsvar_list
        else:
            biogem_path = Path(self.model_path).expanduser() / "biogem"
            ## return all files in the directory
            all_bg_files = Path(biogem_path)

            ## find out files with .res extension
            tsvar_list = [f for f in all_bg_files.glob("*.res")]         
            return tsvar_list

    def get_config(self, config_name='BIOGEM'):
        file = Path(self.model_path).expanduser() / f"data_{config_name}"
        configdata = self._parse_namelist(file)[f"INI_{config_name}_NML"]
        return configdata

    def _parse_namelist(self, file_path):
        namelist = {}
        current_group = None

        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if not line or line.startswith('!'):
                    continue

                # Check for the start of a namelist group
                if line.startswith('&'):
                    current_group = line[1:].strip()  # Remove the '&' and get the group name
                    namelist[current_group] = {}
                # Check for the end of a namelist group
                elif line.startswith('/'):
                    current_group = None
                elif '=' in line:
                    # Parse the key-value pairs within the group
                    key, value = self._split_key_value(line)

                    # Convert Fortran logicals (.true., .false.) and numeric values
                    if value.lower() in ['.true.', 't']:
                        value = True
                    elif value.lower() in ['.false.', 'f']:
                        value = False
                    else:
                        try:
                            value = float(value) if '.' in value or 'e' in value.lower() else int(value)
                        except ValueError:
                            # Remove trailing comma if present
                            value = value.rstrip(',')
                            # Remove surrounding quotes if present
                            value = value.strip("'\"")

                    if current_group:
                        namelist[current_group][key] = value

        return namelist

    def _split_key_value(self, line):
        # Split the line into key and value, respecting quotes
        parts = line.split('=', 1)
        key = parts[0].strip()
        value = parts[1].strip() if len(parts) > 1 else ''

        # Handle quoted values with potential commas
        if value.startswith(("'", '"')) and value.endswith(("'", '"')):
            return key, value.strip("'\"")
        
        # For unquoted values, remove trailing comma if present
        return key, value.rstrip(',')

    def get_ts(self, var: str, to_ScatterData=False):
        """
        read in time series output of GENIE

        :param var: the name of the target variable
        :return: a pandas DataFrame
        """
        if not self.is_ensemble:
            filename = f"biogem_series_{var}.res"                        
            f = Path(self.model_path).expanduser() / "biogem" / filename

            if not f.is_file():
                raise ValueError(f"{f} does not exist")

            with f.open("r") as file:
                lines = file.readlines()

            # Remove " %" from the header and split it using " / "
            # not '/' which is also included in the isotope unit (e.g. 'd13C o/oo')
            header = [col.strip().lstrip("% ") for col in lines[0].split(" / ")]
            data = [line.split() for line in lines[1:]]
            df = pd.DataFrame(data, columns=header)

            ## convert to numeric
            df = df.apply(pd.to_numeric, errors="coerce")

            ## add a column of model name                        
            if self.model_path.__class__.__name__ == "PosixPath":
                df["model"] = self.model_path.name
            elif self.model_path.__class__.__name__ == "str":
                df["model"] = self.model_path.split("/")[-1]
            else:
                # do nothing
                pass

            if to_ScatterData:
                df = ScatterData(df)
                df.set_index('time (yr)')
            
            return df
        else:
            ## concatenate all models' data into one data frame
            df_list = []

            for path in self.model_path:
                ## do the same thing as above
                filename = f"biogem_series_{var}.res"
                f = Path(self.model_path).expanduser() / "biogem" / filename

                if not f.is_file():
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

            if to_ScatterData:
                all_df = ScatterData(df)
                all_df.set_index('time (yr)')

            return all_df

    def _render_diag_avg(self, f):
        """read the diagnostic file of cGENIE

        :param f: the name of the diagnostic file
        """
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
        
        
    def get_diag_avg(self, target_year, is_index=False, pattern_year=r"\d+",pattern_year_digit=r"\d{3}"):
        """
        read the diagnostic file of cGENIE

        :param target_year: the target year of the diagnostic file
        :param is_index: if True, target_year is the index of the sorted years
        :param pattern: the pattern of the diagnostic file name
        
        :return: a pandas DataFrame

        Example
        ----------
        >>> from cgeniepy.model import GenieModel
        >>> model = GenieModel("path_to_GENIE_output")
        >>> model.get_diag_avg(9999)
        """
        pattern = r"biogem_year_(?P<year>" + pattern_year + ")_" + pattern_year_digit + "_diag_GLOBAL_AVERAGE"

        all_years, diag_files = [], []
        for path in self.tsvar_list:
          match = re.search(pattern, path.name)
          if match:
            # Extract the year using the named capture group 'year'
            year = match.group('year')
            all_years.append(int(year))
            diag_files.append(path)

        if self.is_ensemble:
            assert is_index == False, "is_index is not supported for ensemble model"
            output = []
            ## get all files matching target years
            for year, file in zip(all_years, diag_files):
                if year == target_year:
                    single_df = self._render_diag_avg(file)
                    ## add label
                    single_df["model"] = file.parent.parent.name
                    output.append(single_df)
            return pd.concat(output)
        else:
            ## A single model        
            if is_index:
                sorted_diagfiles = dict(sorted(zip(all_years, diag_files)))        
                sorted_years = sorted(all_years)
                try:
                    target_year = sorted_years[target_year]
                except IndexError:
                    print(f"Index {target_year} is out of range. Available indices are from 0 to {len(sorted_years) - 1}.")
                return self._render_diag_avg(sorted_diagfiles[target_year])
            else:
                # merge two lists into a dictionary (and sort by year)
                diagfiles_dict = dict(zip(all_years, diag_files))            
                if target_year not in all_years:
                    raise ValueError(f"{target_year} not found in the diagnostic files. Available years are {all_years}")    
                return self._render_diag_avg(diagfiles_dict[target_year])
        

    def grid_mask(self):
        """
        cGENIE continent mask array (by setting zero values),
        either calculated from existing data or using biogem.grid_mask
        """
        try:
            grid_mask = self.get_var("grid_mask")
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
        ##print("grid area returned in the unit of 'm2'")
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
            ## print("grid volume calculated in the unit of 'm3'")
        except ValueError:
            print(
                "Depth array not found! Please ensure 3d data is exported in the model!"
            )

        return grid_volume
