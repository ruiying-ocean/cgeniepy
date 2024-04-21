import pathlib
from cgeniepy.skill import DFComparison
import geopandas as gpd
from shapely.geometry import Point
import pandas as pd

from io import StringIO

from cgeniepy.grid import Interporaltor

from .plot import ScatterDataVis
from .grid import GridOperation
import cgeniepy.array as ca


class ScatterData:
    """ScatterData is a class to store non-gridded data, often seen in the observations
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
        "set the index of the dataframe"
        self.data.set_index(index, inplace=True)
        self.index= index

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

    def __getitem__(self, item):
        return self.data[item]
    
    def _check_cols(self, cols):
        """
        Check if the columns are present in the dataframe.

        Args:
            cols (list): List of column names to check.

        Raises:
            ValueError: If any of the columns are not found in the dataframe.
        """
        for col in cols:
            if col not in self.data.columns:
                raise ValueError(f"{col} not found in the dataframe")

    def detect_basin(self):
        "use point-in-polygon strategy to detect modern ocean basin according to lon/lat column"
        
        def sub_detect_basin(lon, lat):
            p = Point(lon, lat)
            ocean_name = oceans[oceans.contains(p)].Oceans.values
            if ocean_name.size > 0:
                return ocean_name[0]
            else:
                return ""
        file_path = pathlib.Path(__file__).parent.parent / "data/oceans/oceans.shp"
        oceans = gpd.read_file(file_path)
        
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

        df = self.data.copy()

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

