"""
=============================
Plot Other Models' output
=============================

This example shows how to use cgeniepy to plot gridded data from other models (CESM and HadCM3 here).
You can download them from https://zenodo.org/records/13786014.
"""
from cgeniepy.array import GriddedData
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

## read in the data
cesm_lgm = xr.load_dataset("/Users/yingrui/cgeniepy/examples/CESM_LGM_var_regrid.nc")
## construct GriddedData object
cesm_temp = GriddedData(cesm_lgm['TEMP'], mutable=False, attrs=cesm_lgm['TEMP'].attrs)

## same for HadCM3L
hadcm3_lgm=  xr.load_dataset("/Users/yingrui/cgeniepy/examples/teitu_020_o.pgclann.nc", decode_times=False)
hadcm3_temp = GriddedData(hadcm3_lgm['temp_ym_dpth'], mutable=False, attrs=hadcm3_lgm['temp_ym_dpth'].attrs)

fig, axs = plt.subplots(1, 2, figsize=(10, 5), subplot_kw={'projection': ccrs.Robinson()})

p = cesm_temp.isel(time=0,z_t=0).to_GriddedDataVis()
p.aes_dict['pcolormesh_kwargs']['vmax'] = 30
p.plot(ax=axs[0], outline=True, colorbar=False)
axs[0].set_title('CESM LGM')

p2 = hadcm3_temp.isel(t=0,depth_1=0).to_GriddedDataVis()
p2.aes_dict['pcolormesh_kwargs']['vmax'] = 30
im = p2.plot(ax=axs[1], outline=True, colorbar=False)
axs[1].set_title('HadCM3L LGM')

fig.colorbar(im, ax=axs, orientation='horizontal', label='Sea surface temperature (°C)', fraction=0.05, pad=0.07)
