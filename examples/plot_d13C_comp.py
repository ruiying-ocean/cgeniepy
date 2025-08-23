"""
=======================================================
Search the nearest grid point for a given location
=======================================================

This example is a combination of ScatterData and GriddedData to search the nearest grid point for a given location. I use CESM model output(You can download them from https://zenodo.org/records/13786013) and LGM d13C data from Peterson et al. 2014 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013PA002552) as an example.
"""

from cgeniepy.table import ScatterData
from cgeniepy.array import GriddedData
import xarray as xr
from cgeniepy.utils import download_zenodo_file

record_id = "13786013"
filename = "CESM_LGM_var_regrid.nc"
local_file_path = download_zenodo_file(record_id, filename)

## read in the data and construct GriddedData object
cesm_lgm = xr.load_dataset(local_file_path)
cesm_13C = GriddedData(cesm_lgm['CISO_DIC_d13C'], attrs=cesm_lgm['CISO_DIC_d13C'].attrs)
cesm_13C_last = cesm_13C.isel(time=-1)

## read in the proxy data and construct ScatterData object
filename = "LGM_d13c_CLEAN.txt"
local_file_path = download_zenodo_file(record_id, filename)

proxy_d13C = ScatterData(local_file_path, delimiter="\t", header=None)
proxy_d13C.data.columns = ["Longitude", "Latitude", "Depth", "observational d13C","Event"]
proxy_d13C.set_index(["Latitude", "Longitude", "Depth"])

model_data = []

for i in proxy_d13C.data.index:
    lat, lon, depth = i
    depth = depth * 100 ## not necessary for cGENIE whose depth is in m
    pos = (depth, lat, lon)
    
    data = cesm_13C_last.search_point(pos, ignore_na=True)

    model_data.append(data)

## add the model data to the dataframe
proxy_d13C.data["CESM_d13C"] = model_data

## plot the comparison
## by default, model data is in the col, and observational col is in the second
proxy_d13C.compare("CESM_d13C","observational d13C").plot()
