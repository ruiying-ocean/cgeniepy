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

    def interpolate(self):
        ## a tuple of coordinate arrays
        coords =  tuple([self.df[dim].values for dim in self.dims])
        ## array values
        values = self.df[self.var].values
        dims = self.dims
        return Interporaltor(dims, coords, values, 200, 'ir-linear')


    def to_genie(
        self,
        low_threshold=None,
        new_low_bound=None,
        high_threshold=None,
        new_high_bound=None,
    ):
        """
        Regrid a dataframe within certain format to cGENIE grids

        Input: A dataframe with three necessary columns: "Latitude", "Longitude", "Observation".
        The other columns will be unselected. Note that observation column is the data you want to process.
        Rename dataframe by using command: `df = df.rename({"old name": "new name"}, axis='columns')`

        Output:
        A 36x36 2D array.

        Optional parameters:
        low_threshold/new_low_bound: change the data that is smaller than low_threshold to new_low_bound.
        For example, covert any data < 0.03 to 0. Work same to high_threshold/new_high_bound.
        """

        df = self.df.copy()

        # subset
        df = df[["Longitude", "Latitude", "Observation"]]

        # drop NAN
        df = df.dropna(axis="rows", how="any")

        # regrid coordinate
        df.loc[:, "Latitude"] = df.loc[:, "Latitude"].apply(regrid_lat)
        df.loc[:, "Longitude"] = df.loc[:, "Longitude"].apply(regrid_lon)

        # group and aggregate (still in long format)
        df_agg = df.groupby(["Longitude", "Latitude"]).agg("mean")

        # Create a all-nan data frame with full cGENIE grids
        lat = GENIE_lat()
        lon = normal_lon()
        data = np.zeros([36 * 36])
        index = pd.MultiIndex.from_product([lon, lat], names=["Longitude", "Latitude"])
        df_genie = pd.DataFrame(data, index=index, columns=["Observation"])
        df_genie["Observation"] = df_genie["Observation"].replace(0, np.nan)

        # transform the tuple elements in multindex to lists
        df_genie_index_list = [list(item) for item in df_genie.index.values]
        df_agg_index_list = [list(item) for item in df_agg.index.values]

        # copy aggregated dataframe to full-grid dataframe, one by one
        # if the returned dataframe is full of NA, it's very likely from this step
        for i in df_genie_index_list:

            longitude = i[0]
            latitude = i[1]

            if i in df_agg_index_list:
                df_genie.loc[longitude, latitude] = df_agg.loc[longitude, latitude]
            else:
                df_genie.loc[longitude, latitude] = np.nan

        # filter data if necessary
        if (low_threshold is not None) and (new_low_bound is not None):
            df_genie.loc[
                df_genie.Observation < low_threshold, "Observation"
            ] = new_low_bound

        if (high_threshold is not None) and (new_high_bound is not None):
            df_genie.loc[
                df_genie.Observation > high_threshold, "Observation"
            ] = new_high_bound

        # long -> wide data format
        df_genie_wide = df_genie.pivot_table(
            values="Observation", index="Latitude", columns="Longitude", dropna=False
        )

        return df_genie_wide
        
    def drop_na(self):
        "drop rows with NA values"
        self.df = self.df.dropna()
