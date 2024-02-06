import xarray as xr
import pandas as pd
from .array import Array

class Observation:

    def __init__(self, file):
        self.file = file
        
        
    def load_data(self):
        if self.file.endswith(".nc"):
            data =  Array()
            data.array = xr.open_dataset(self.file)
            return data
        elif self.file.endswith(".csv"):
            return pd.read_csv(self.file)
        else:
            raise ValueError("Only csv and netcdf are supported")

    def plot(self):
        pass

    def interpolate(self):
        pass

    def regrid(self):
        "regrid to cGENIE grid"
        pass
