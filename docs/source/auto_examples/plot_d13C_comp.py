"""
=======================================================
Search the nearest grid point for a given location
=======================================================

This example is a combination of ScatterData and GriddedData to search the nearest grid point for a given location. I use CESM model output and LGM d13C data from Peterson et al. 2014 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013PA002552) as an example.
"""

from cgeniepy.table import ScatterData
from cgeniepy.array import GriddedData
import xarray as xr

cesm_lgm = xr.load_dataset("/Users/yingrui/cgeniepy/examples/CESM_LGM_var_regrid.nc")
proxy_d13C = ScatterData("~/Science/lgm_bcp/data/LGM_d13c_CLEAN.txt", delimiter="\t", header=None)
proxy_d13C.data.columns = ["Longitude", "Latitude", "Depth", "d13C","Event"]
proxy_d13C.set_index(["Latitude", "Longitude", "Depth"])

model_data = []
cesm_13C = GriddedData(cesm_lgm['CISO_DIC_d13C'], mutable=False, attrs=cesm_lgm['CISO_DIC_d13C'].attrs)
cesm_13C_last = cesm_13C.isel(time=-1)

for i in proxy_d13C.data.index:
    lat, lon, depth = i
    depth = depth * 100
    pos = (depth, lat, lon)
    
    data = cesm_13C_last.search_point(pos, ignore_na=True)

    model_data.append(data)

## add the model data to the dataframe
proxy_d13C.data["model_d13C"] = model_data

print(proxy_d13C.compare("d13C", "model_d13C").pearson_r())
proxy_d13C.compare("d13C", "model_d13C").plot()
