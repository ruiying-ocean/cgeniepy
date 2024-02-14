import pathlib
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd

from io import StringIO

from cgeniepy.grid import Interporaltor

from .skill import DataFrameComp
from .plot import ScatterVis
from .array import GriddedData

class ScatterData(ScatterVis):
    """ScatterData is a class to store non-gridded data, often seen in the observations
    """

    def __init__(self, path: str, *args, **kwargs):
        """
        Initialize a ScatterData object.

        Parameters:
        path (str): The path to the file or the data.
        coord_cols (dict): A dictionary specifying the coordinate columns.

        *args: Additional positional arguments to be passed to pd.read_csv().
        **kwargs: Additional keyword arguments to be passed to pd.read_csv().
        """
        if path.endswith(".tab"):
            data = self.parse_tab_file(path)
            self.df= pd.read_csv(StringIO(data), *args, **kwargs)
        elif path.endswith("xlsx"):
            self.df = pd.read_excel(path, *args, **kwargs)                
        else:
            self.df = pd.read_csv(path, *args, **kwargs)

        self.dims = []


    def parse_tab_file(self, filename, begin_cmt = '/*', end_cmt = '*/'):
        """
        Read a tab-delimited file and return a pandas DataFrame

        This function is optimised for pangea-format data.
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

    def specify_cols(self, dim_cols: dict, var_col: str):
        "assign coordinate and variable columns"
        self._check_cols(list(dim_cols.values()) + [var_col])
        
        ## if all the cols are given in name
        if 'lon' in dim_cols:
            self.lon = dim_cols['lon']
            self.dims.append(self.lon)
        if 'lat' in dim_cols:
            self.lat = dim_cols['lat']
            self.dims.append(self.lat)
        if 'age' in dim_cols:
            self.age = dim_cols['age']
            self.dims.append(self.age)
        if 'depth' in dim_cols:
            self.depth = dim_cols['depth']
            self.dims.append(self.depth)
            
        self.var = var_col

        ## check the datatype of the coordinate columns
        for col in self.dims:
            if self.df[col].dtype != 'float64':
                self.df[col] = self.df[col].astype('float64')

    def __getitem__(self, item):
        return self.df[item]
    
    def _check_cols(self, cols):
        """
        Check if the columns are present in the dataframe.

        Args:
            cols (list): List of column names to check.

        Raises:
            ValueError: If any of the columns are not found in the dataframe.
        """
        for col in cols:
            if col not in self.df.columns:
                raise ValueError(f"{col} not found in the dataframe")

    def detect_basin(self):
        "use point-in-polygon strategy to detect modern ocean basin according to lon/lat column"
        
        def detect_basin(lon, lat):
            p = Point(lon, lat)
            ocean_name = oceans[oceans.contains(p)].Oceans.values
            if ocean_name.size > 0:
                return ocean_name[0]
            else:
                return ""
        file_path = pathlib.Path(__file__).parent.parent / "data/oceans/oceans.shp"
        oceans = gpd.read_file(file_path)
        
        self.df['basin'] = self.df.apply(lambda row: detect_basin(row[self.lon], row[self.lat]), axis=1)
        print("basin column added to the dataframe!")

    def lookup_model(self, gridded_data, new_col='model_var'):
        ## find the nearest model value given the coordinate columns
        modelvar = []
        for i in range(len(self.df)):
            lat = self.df[self.lat].iloc[i]
            lon = self.df[self.lon].iloc[i]
            kwargs = {'lat': lat, 'lon':lon, 'method': 'nearest'}
            var = gridded_data.search_grid(**kwargs).values
            modelvar.append(var)
        self.df[new_col] = modelvar
        print("model column added to the dataframe!")

        ## create a model-data comparison object
        comp = DataFrameComp(self.df, new_col, self.var)
        return comp

    def to_xarray(self):
        "convert to xarray dataset"
        ## set the coordinate using the data from self.coordinates         
        return self.df.set_index(self.dims, inplace=False).to_xarray()

    def to_gridded(self):
        "convert to gridded data"
        output = GriddedData(self.to_xarray()[self.var])
        return output
    
    def to_genie(self, example):
        "convert to GENIE grid"
        ## to be done
        pass

    def interpolate(self):
        ## a tuple of coordinate arrays
        coords =  tuple([self.df[dim].values for dim in self.dims])
        ## array values
        values = self.df[self.var].values
        dims = self.dims
        return Interporaltor(dims, coords, values, 200, 'ir-linear')
    
