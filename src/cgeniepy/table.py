import numpy as np
from cgeniepy.skill import DFComparison
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd

from io import StringIO

from cgeniepy.grid import Interporaltor, GridOperation

from .plot import ScatterDataVis
from .grid import GridOperation
import cgeniepy.array as ca
from importlib.resources import files
    

class ScatterData:
    """ScatterData is a class to store non-gridded data with columns of coordinates.
    """

    def __init__(self, data, *args, **kwargs):
        """
        Initialize a ScatterData object.

        Parameters:
        data: The path to the file or the data.
        coord_cols (dict): A dictionary specifying the coordinate columns.
        """
        ## if already a dataframe
        if isinstance(data, pd.DataFrame):
            self.data = data

        ## if a file path then read the file into a dataframe
        if isinstance(data, str):
            if data.endswith(".tab"):
                data = self.parse_tab_file(data)
                self.data= pd.read_csv(StringIO(data), *args, **kwargs)
            elif data.endswith("xlsx"):
                self.data = pd.read_excel(data, *args, **kwargs)                
            else:
                self.data = pd.read_csv(data, *args, **kwargs)

        ## if the index is already set in the data, then set the coordinates
        if not isinstance(self.data.index, pd.core.indexes.range.RangeIndex):
            self.index= list(self.data.index.names)
            GridOperation().set_coordinates(obj=self, index=self.index)

    def set_index(self, index):
        """Tell the object which columns are the coordinates.        
        """
        self.data.set_index(index, inplace=True)
        self.index= index
        GridOperation().set_coordinates(obj=self, index=self.index)


    def reset_index(self, inplace=True):
        if inplace:
            self.data = self.data.reset_index()
        else:
            return self.data.reset_index()

    def parse_tab_file(self, filename, begin_cmt = '/*', end_cmt = '*/'):
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
        output= Interporaltor(dims, coords, values, 200, 'ir-linear')
        return output

    def to_geniebin(
        self,
        var,
    ):
        """
        Regrid a dataframe within certain format to cGENIE grids

        Output:
        A 36x36 2D array.
        """

        go= GridOperation()
        
        src_df = self.data.reset_index(inplace=False)

        # subset
        src_df = src_df[[self.lon, self.lat, var]]

        # drop NAN
        src_df = src_df.dropna(axis="rows", how="any")

        # regrid coordinate: genie lat x normal lon
        ## lat
        src_df.loc[:, self.lat] = src_df.loc[:, self.lat].apply(go.geniebin_lat)
        ## lon 
        src_df.loc[:, self.lon] = src_df.loc[:, self.lon].apply(go.normbin_lon)

        # aggregated source data (in long format)
        src_df_agg = src_df.groupby([self.lat,self.lon]).agg("mean")
        return src_df_agg

    def drop_na(self, *args, **kwargs):
        "drop rows with NA values"
        self.data = self.data.dropna(*args, **kwargs)

    def to_ScatterDataVis(self):
        "convert to ScatterDataVis"
        return ScatterDataVis(self)

    def plot(self, var, *args, **kwargs):
        "plot the data"
        return self.to_ScatterDataVis().plot(var, *args, **kwargs)

    def compare(self, var1, var2):
        "compare two ScatterData objects"
        return DFComparison(self.data, var1, var2)

