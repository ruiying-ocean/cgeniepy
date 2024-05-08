"""
=========================================
Plot 2D transect of tracers in each basin
=========================================

This example plots the modelled oxygen distribution in cGENIE.

The following features in the package are used:

#. Access data through `cgeniepy.model` module

#. A basin-mask operation

#. A linear interpolation

#. Get pretty color palette

#. Customise the plotting details
"""

from cgeniepy.model import GenieModel
from cgeniepy.plot import CommunityPalette
import matplotlib.pyplot as plt

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
ocn_po4 = model.get_var("ocn_PO4").isel(time=-1)

fig, axs=plt.subplots(nrows=3, ncols=1, figsize=(6,9), tight_layout=True)

basins = ['Atlantic', 'Pacific', 'Indian']

cmap = CommunityPalette('tol_rainbow').colormap

for i in range(3):
    basin_data = model.get_var('ocn_O2').isel(time=-1).mask_basin(base='worjh2',basin=basins[i], subbasin='')
    basin_data_interp = basin_data.mean(dim='lon').interpolate(grid_number=50).to_GriddedDataVis()
    basin_data_interp.aes_dict['pcolormesh_kwargs']['cmap'] = cmap

    basin_data_interp.plot(ax=axs[i], contour=False, outline=True)
    axs[i].title.set_text(basins[i])

plt.show()
