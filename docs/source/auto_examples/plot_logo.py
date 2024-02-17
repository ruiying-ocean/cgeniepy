"""
================================
Customise the 2D map projection
================================

This example shows how to customise the 2D map including the projection, the color map, which is used as the logo of this package.
"""
from cgeniepy.model import GenieModel
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

## Read in the model
model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
sst = model.get_var("ocn_sur_temp").isel(time=-1)

## use the Orthographic projection
## a full list of projections can be found at
## https://scitools.org.uk/cartopy/docs/latest/reference/projections.html#cartopy-projections
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Orthographic()})

## set the color map
sst.aes_dict['pcolormesh_kwargs']['cmap'] = plt.cm.inferno
sst.plot(ax=ax, outline=True)
plt.show()
