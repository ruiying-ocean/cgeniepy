"""
=========================================
Plot basin-specific data
=========================================

Plot the model variable for a spficic basin, here we use sea surface temperature in the Atlantic Ocean as an example.
"""

import cgeniepy

model = cgeniepy.sample_model()
sst = model.get_var('ocn_sur_temp').isel(time=-1)
## plot the sea surface temperature in the Atlantic Ocean only
sst.sel_modern_basin(['NAO','EAO','SAO']).plot(outline=True, colorbar=True)
