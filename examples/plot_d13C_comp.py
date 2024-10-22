"""
=======================================================
Search the nearest grid point for a given location
=======================================================

This example is a combination of ScatterData and GriddedData to search the nearest grid point for a given location. I use CESM model output(You can download them from https://zenodo.org/records/13786013) and LGM d13C data from Peterson et al. 2014 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013PA002552) as an example.
"""

from cgeniepy.table import ScatterData
from cgeniepy.array import GriddedData
import xarray as xr
import subprocess

## Download necessary files, you need to install zenodo_get first
## by `pip install zenodo_get`, or, just download it from the link above
subprocess.call(["zenodo_get", "10.5281/zenodo.13786013", "-o", "~/Downloads/"])


## read in the data and construct GriddedData object
cesm_lgm = xr.load_dataset("~/Downloads/CESM_LGM_var_regrid.nc")
cesm_13C = GriddedData(cesm_lgm['CISO_DIC_d13C'], attrs=cesm_lgm['CISO_DIC_d13C'].attrs)
cesm_13C_last = cesm_13C.isel(time=-1)

## read in the proxy data and construct ScatterData object
proxy_d13C = ScatterData("~/Downloads/LGM_d13c_CLEAN.txt", delimiter="\t", header=None)
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
