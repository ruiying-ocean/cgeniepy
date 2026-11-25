"""
=======================================
Map styles
=======================================

A map's style is set by two choices: its projection (``projection=``) and its frame
(``gridline=`` and ``zebra_frame=``). Each call below draws one style of the same field.
"""

import cgeniepy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

model = cgeniepy.sample_model()
sst = model.get_var("ocn_sur_temp").isel(time=-1)

###############################################################################
# The default: an Equal Earth map, with the model's land outlined.

sst.plot(colorbar=True, outline=True)

###############################################################################
# Any cartopy projection, here Robinson centred on the Pacific.

sst.plot(projection=ccrs.Robinson(central_longitude=180), colorbar=True, outline=True)

###############################################################################
# GENIE's own grid, as genieplot draws it: latitude spaced by its sine, so every
# row of the equal-area grid is equally tall, centred on the grid's longitudes.

sst.plot(projection="sine", colorbar=True, gridline=True, outline=True)

###############################################################################
# A rectangular map with the frame NCL draws: outward ticks on all four sides,
# labelled in degrees, under filled contours.

sst.plot(projection=ccrs.PlateCarree(), contourf=True, colorbar=True, gridline=True, outline=True)

###############################################################################
# GMT's black-and-white frame on the same ticks. On a sine axis the latitude
# labels crowd near the poles, so the ticks stop at 60°.

p = sst.to_GriddedDataVis()
p.aes_dict["gridline_kwargs"]["ylocs"] = [-60, -30, 0, 30, 60]
p.plot(projection="sine", colorbar=True, gridline=True, zebra_frame=True)

plt.show()
