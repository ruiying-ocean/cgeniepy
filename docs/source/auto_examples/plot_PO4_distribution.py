"""
=========================================
Plot 2D transect of tracers in each basin
=========================================

This example shows how to plot the PO4 distribution in each basin.
"""

from cgeniepy.model import GenieModel
import matplotlib.pyplot as plt

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
ocn_po4 = model.get_var("ocn_PO4").isel(time=-1)

fig, axs=plt.subplots(nrows=1, ncols=3, figsize=(15, 3), tight_layout=True)        

basins = ['Atlantic', 'Pacific', 'Indian']

for i in range(3):
	basin_data = model.get_var('ocn_PO4').isel(time=-1).mask_basin(base='worjh2',basin=basins[i], subbasin='')
	basin_data.array.values = basin_data.array.values * 1E6
	basin_data.mean(dim='lon').interpolate().plot(ax=axs[i], contour=True)
	axs[i].title.set_text(basins[i])

plt.show()        
