"""
==========================================================
Regrid observational data and compare with cGENIE output
==========================================================

This example shows how to regrid observational data to cGENIE grids and compare with cGENIE output.

The following features are used in this example:

#.

"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from cgeniepy.table import ScatterData
from cgeniepy.grid import GridOperation
from cgeniepy.array import GriddedData
from cgeniepy.model import GenieModel
from cgeniepy.skill import ArrComparison

glodap_temp = xr.load_dataset("/Users/yingrui/Science/GLODAPv2.2016b_MappedClimatologies/GLODAPv2.2016b.temperature.nc")['temperature']
glodap_temp = GridOperation().xr_e2n(glodap_temp)
glodap_temp = ScatterData(glodap_temp.isel(depth_surface=0).to_dataframe())
glodap_temp = glodap_temp.to_geniebin(var='temperature').to_xarray()['temperature']

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
worlg4_sst = model.get_var('ocn_sur_temp').isel(time=-1).normalise_longitude(method='g2n')

## calculate the skill score
print("M-score",ArrComparison(glodap_temp.values, worlg4_sst.data.values).mscore())

## add mask to the glodap data (just for plot)
worlg4_mask = np.isnan(model.get_var('grid_mask').normalise_longitude(method='g2n').data)
masked_glodap_temp = glodap_temp.where(~worlg4_mask)
masked_glodap_temp = GriddedData(masked_glodap_temp, mutable=False, attrs=glodap_temp.attrs)
masked_glodap_temp.attrs['long_name'] = 'GLODAPv2.2016b temperature'
masked_glodap_temp.attrs['units'] = 'deg C'


fig, axs = plt.subplots(1,2,subplot_kw={"projection": ccrs.PlateCarree()})
masked_glodap_temp.plot(ax=axs[0], outline=True, colorbar=True)
worlg4_sst.plot(ax=axs[1], outline=True, colorbar=True)
