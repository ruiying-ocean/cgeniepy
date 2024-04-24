"""
=========================================
Plot 2D transect of tracers in each basin
=========================================

This example plots the modelled PO4 distribution in cGENIE.

The following features in the package are used:

#. Access data through `cgeniepy.model` module

#. A basin-mask operation

#. A linear interpolation

#. Get additional color palette (mirrors the one in ODV)

#. Customise the plotting details
"""

from cgeniepy.model import GenieModel
from cgeniepy.plot import CommunityPalette
import matplotlib.pyplot as plt

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
ocn_po4 = model.get_var("ocn_PO4").isel(time=-1)

fig, axs=plt.subplots(nrows=1, ncols=3, figsize=(15, 3), tight_layout=True)

basins = ['Atlantic', 'Pacific', 'Indian']

odv_cmap = CommunityPalette().get_palette('ODV', reverse=True)

for i in range(3):
    basin_data = model.get_var('ocn_PO4').isel(time=-1).mask_basin(base='worjh2',basin=basins[i], subbasin='')
    basin_data_interp = basin_data.mean(dim='lon').interpolate().to_GriddedDataVis()
    basin_data_interp.aes_dict['pcolormesh_kwargs']['cmap'] = odv_cmap
    basin_data_interp.plot(ax=axs[i], contour=True)
    axs[i].title.set_text(basins[i])

plt.show()
