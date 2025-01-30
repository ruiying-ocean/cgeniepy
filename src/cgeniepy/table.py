import numpy as np
from cgeniepy.skill import DFComparison
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd

from io import StringIO
import re
from typing import Union

from cgeniepy.grid import Interpolator, GridOperation

from cgeniepy.plot import ScatterDataVis
from cgeniepy.grid import GridOperation
import cgeniepy.array as ca
from importlib.resources import files


class ScatterData:

    modify_in_place = True
    
    """ScatterData is a class to store non-gridded data with columns of coordinates."""

    def __init__(self, data: Union[pd.DataFrame, int, str], **kwargs):
        """
        Initialize a ScatterData object.

        Parameters:
        data: The path to the file, the data, or a PanDataSet ID.
        **kwargs: Additional keyword arguments for pandas read functions.
        """
        self.data = self._process_data(data, **kwargs)

        # if the index is already set in the data, then set the coordinates
        if not isinstance(self.data.index, pd.core.indexes.range.RangeIndex):
            self.index= list(self.data.index.names)
            GridOperation().set_coordinates(obj=self, index=self.index)

    def _process_data(self, data: Union[pd.DataFrame, int, str], **kwargs) -> pd.DataFrame:
        if isinstance(data, pd.DataFrame):
            return data
        if isinstance(data, int):
            try:
                from pangaeapy.pandataset import PanDataSet
                return PanDataSet(data).data
            except ImportError:
                print("Unable to import PanDataSet from pangaeapy. Please make sure the package is installed.")                            
        if isinstance(data, str):
            return self._process_string_data(data, **kwargs)
        raise ValueError("Unsupported data type. Expected DataFrame, int, or str.")

    def _process_string_data(self, data: str, **kwargs) -> pd.DataFrame:
        if data.endswith(".tab"):
            return pd.read_csv(StringIO(self._parse_tab_file(data)), **kwargs)
        if data.endswith(".xlsx"):
            return pd.read_excel(data, **kwargs)
        if "PANGAEA" in data:
            try:
                from pangaeapy.pandataset import PanDataSet
            except ImportError:
                print("Unable to import PanDataSet from pangaeapy. Please make sure the package is installed.")   

            doi = self._extract_doi(data)
            return PanDataSet(doi).data
        return pd.read_csv(data, **kwargs)

    @staticmethod
    def _extract_doi(url: str) -> str:
        doi_pattern = r"10\.\d{4,9}/[-._;()/:A-Z0-9]+"
        doi_match = re.search(doi_pattern, url, re.IGNORECASE)
        if doi_match:
            return doi_match.group(0)
        raise ValueError("No valid DOI found in the URL.")

    def __repr__(self):
        prefix = "ScatterData\n"
        columns = f"Columns: {self.data.columns}\n"
        rows = f"Rows: {len(self.data)}\n"
        if hasattr(self, "index"):
            index = f"Index: {self.index}\n"
        else:
            index = ""

        return prefix + columns + index + rows

    def _parse_tab_file(self, filename, begin_cmt = '/*', end_cmt = '*/'):
        """
        Read a tab-delimited file and return a pandas that is optimised for pangea-format data.

        :param filename: The name of the file to read.
        :param begin_cmt: The string that marks the beginning of a comment.
        :param end_cmt: The string that marks the end of a comment.

        :return: A pandas dataframe.
        """
        lines = []
        in_comment = False

        with open(filename) as f:
            for line in f:
                if line.strip().startswith(begin_cmt):
                    in_comment = True
                    continue

                if line.strip().startswith(end_cmt):
                    in_comment = False
                    continue

                if not in_comment:
                    lines.append(line.rstrip('\n'))

        data = '\n'.join(lines)
        return data    

    def __getitem__(self, item):
        return self.data[item]
    
    def _check_cols(self, cols):
        """
        Check if the columns are present in the dataframe.

        :param cols: A list of columns to check.

        :raises: ValueError if a column is not found.
        """
        for col in cols:
            if col not in self.data.columns:
                raise ValueError(f"{col} not found in the dataframe")


    def set_index(self, index):
        """Tell the object which columns are the coordinates.        
        """
        self.data.set_index(index, inplace=True)
        self.index= index
        GridOperation().set_coordinates(obj=self, index=self.index)


    def reset_index(self):
        if ScatterData.modify_in_place:
            self.data = self.data.reset_index()
            return self
        else:
            return self.data.reset_index()
        
    def detect_basin(self):
        """use point-in-polygon strategy to detect modern ocean basin according to lon/lat column

        Example
        ----------
        >>> data = ScatterData(data)
        >>> data.set_index(['lon', 'lat'])
        >>> data.detect_basin()
        """
        
        def sub_detect_basin(lon, lat):
            p = Point(lon, lat)
            ocean_name = oceans[oceans.contains(p)].Oceans.values
            if ocean_name.size > 0:
                return ocean_name[0]
            else:
                return ""
            
        file_path = str(files("data").joinpath("oceans/oceans.shp"))
        oceans = gpd.read_file(file_path)
        
        if len(self.index) != 2:
            raise ValueError("The index must have two columns: lon and lat")
        
        assert self.lon in self.index, f"{self.lon} not found in the index"
        assert self.lat in self.index, f"{self.lat} not found in the index"

        result = self.data.reset_index().apply(lambda row: sub_detect_basin(row[self.lon], row[self.lat]), axis=1)        
        
        ## add to self.data
        if ScatterData.modify_in_place:            
            self.data['basin'] = list(result)
            return self
        else:        
            tmp_data = self.data.reset_index()
            tmp_data['basin'] = result
            tmp_index = self.index
            y = ScatterData(tmp_data)
            y.set_index(tmp_index)
            return y


    def to_xarray(self):
        "convert to xarray dataset"
        ## set the coordinate using the data from self.coordinates         
        return self.data.to_xarray()

    def to_GriddedData(self, var):
        "convert to gridded data"
        output = ca.GriddedData(self.to_xarray()[var])
        return output

    def interpolate(self, var):
        ## a tuple of coordinate arrays
        coords =  tuple([self.data.reset_index(inplace=False)[dim].values for dim in self.index])
        ## array values
        values = self.data[var].values
        dims = self.index

        ## if 1D
        if len(dims) == 1:            
            output= Interpolator(dims, coords, values, 200, '1d').to_dataframe()
            output.rename(columns={"interpolated_values": var}, inplace=True)
            output.set_index(dims, inplace=True)
        else:
            output= Interpolator(dims, coords, values, 200, 'ir-linear')

        if ScatterData.modify_in_place:
            self.data = output
            return self
        else:
            return output


    def to_geniebin(
        self,
        var,
        agg_method="mean"
    ):
        """
        Regrid a dataframe within certain format to cGENIE grids.
        This method does not consider the land-sea mask. So if convert from a higher resolution data to GENIE grids,
        the land-sea mask could be different (Not a problem in comparison though).

        :param var: The variable to regrid.
        :param agg_method: The aggregation method to use when regridding.
        :return: An indexed data frame
        """

        go= GridOperation()

        src_df = self.data.reset_index(inplace=False)

        # subset
        select_cols = self.index + [var]
        src_df = src_df[select_cols]

        # drop NAN
        src_df = src_df.dropna(axis="rows", how="any")

        # regrid coordinate: genie lat x normal lon
        ## lat
        if hasattr(self, "lat"):
            src_df.loc[:, self.lat] = src_df.loc[:, self.lat].apply(go.geniebin_lat)
        if hasattr(self, "lon"):
            src_df.loc[:, self.lon] = src_df.loc[:, self.lon].apply(go.normbin_lon)
        if hasattr(self, "depth"):
            src_df.loc[:, self.depth] = src_df.loc[:, self.depth].apply(go.geniebin_depth)

        # aggregated source data (in long format)
        src_df_agg = src_df.groupby(self.index).agg(agg_method)
        return src_df_agg

    def drop_na(self, *args, **kwargs):
        "drop rows with NA values"
        self.data = self.data.dropna(*args, **kwargs)
        return self

    def to_ScatterDataVis(self):
        "convert to ScatterDataVis"
        return ScatterDataVis(self)

    def plot(self, var, kind='scatter',*args, **kwargs):
        "plot the data"
        return self.to_ScatterDataVis().plot(var, kind=kind,*args, **kwargs)

    def compare(self, var1, var2, model_name=None, obs_name=None):
        "compare two ScatterData objects"
        if not model_name:
            model_name = var1            
        if not obs_name:
            obs_name = var2
            
        return DFComparison(self.data, var1, var2, model_name=model_name, obs_name=obs_name)

    def rolling(self, window, *args, **kwargs):
        "apply rolling to the data"
        if ScatterData.modify_in_place:
            self.data = self.data.rolling(window, *args, **kwargs)
            return self
        else:
            return self.data.rolling(window, *args, **kwargs)

