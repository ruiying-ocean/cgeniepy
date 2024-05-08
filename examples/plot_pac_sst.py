"""
=========================================
Plot basin-specific data
=========================================

Plot the model variable for a spficic basin, here we use sea surface temperature in the Atlantic Ocean as an example.
"""

from cgeniepy.model import GenieModel

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
sst = model.get_var('ocn_sur_temp').isel(time=-1)
## plot the sea surface temperature in the Atlantic Ocean only
sst.sel_modern_basin(['NAO','EAO','SAO']).plot(outline=True, colorbar=True)
