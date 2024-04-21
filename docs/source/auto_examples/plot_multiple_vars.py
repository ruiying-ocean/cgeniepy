"""
=======================================================
Plot multiple variables
=======================================================

This example shows how to plot multiple variables from a GenieModel object. I use a modern model run as an example and plot the surface temperature, PO4, Fe, and O2.
"""

## many variables
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cgeniepy.model import GenieModel

fig, axs = plt.subplots(2, 2, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})

pi_model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN", gemflag='biogem')

variable = ['ocn_sur_temp', 'ocn_sur_PO4', 'ocn_sur_TDFe', 'ocn_sur_O2']

for i in range(4):
    pi_model.get_var(variable[i]).isel(time=-1).plot(ax=axs.flatten()[i],colorbar=True, outline=True, contourf=True)
