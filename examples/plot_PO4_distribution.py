"""
=========================================
Plot 2D transect of tracers in each basin
=========================================

This example plots the modelled oxygen distribution in cGENIE.

The following features in the package are used:

#. Access data through `cgeniepy.model` module

#. A basin-mask operation

#. Filled contours that reach the model's coastline and sea floor

#. Get pretty color palette

#. Customise the plotting details
"""

import cgeniepy
from cgeniepy.plot import CommunityPalette
import matplotlib.pyplot as plt

model = cgeniepy.sample_model()

fig, axs=plt.subplots(nrows=3, ncols=1, figsize=(6,9), tight_layout=True)

basins = ['Atlantic', 'Pacific', 'Indian']

cmap = CommunityPalette('tol_rainbow').colormap

for i in range(3):
    basin_data = model.get_var('ocn_O2').isel(time=-1).mask_basin(base='worjh2',basin=basins[i], subbasin='')
    basin_transect = basin_data.mean(dim='lon').to_GriddedDataVis()
    basin_transect.aes_dict['contourf_kwargs']['cmap'] = cmap

    basin_transect.plot(ax=axs[i], pcolormesh=False, contourf=True, outline=True)
    axs[i].title.set_text(basins[i])

plt.show()
