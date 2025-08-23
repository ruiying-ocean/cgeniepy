"""
=======================================================
Plot multiple variables
=======================================================

This example shows how to plot multiple variables from a GenieModel object. I use a modern model run as an example and plot the surface temperature, PO4, Fe, and O2.

The following features are particularly demonstrated:

#. Customizing the plot projection with cartopy

#. Plotting multiple variables in one figure

#. Controlling the elements in cgeniepy.GriddedDataVis object

#. Using CommunityPalette to get a pretty color palette
"""


import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cgeniepy
from cgeniepy.plot import CommunityPalette

fig, axs = plt.subplots(2, 2, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

pi_model = cgeniepy.sample_model()

variable = ['ocn_sur_temp', 'bio_fexport_POC', 'ocn_sur_PO4','ocn_sur_O2']
cmap = ['ocean_temp', 'tol_rainbow', 'parula','ODV']

for i in range(4):
    data = pi_model.get_var(variable[i]).isel(time=-1).to_GriddedDataVis()
    data.aes_dict['pcolormesh_kwargs']['cmap'] = CommunityPalette(cmap[i]).colormap
    data.aes_dict['contourf_kwargs']['cmap'] = CommunityPalette(cmap[i]).colormap    

    data.aes_dict['facecolor_kwargs']['c'] = 'white'
    data.plot(ax=axs.flatten()[i],colorbar=True, outline=True, gridline=True, contourf=True)
