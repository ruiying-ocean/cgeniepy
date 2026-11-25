"""
=======================================
Maps on GENIE's own grid
=======================================

GENIE's 36 x 36 grid is equal-area: its rows are spaced evenly in the sine of latitude,
so a polar row spans far more degrees than a tropical one. ``projection="sine"`` draws
the map in that space, as genieplot does: every row is equally tall, every cell takes the
same area on the page, and the map is centred on the grid's own longitudes (-260 to 100 °E).
"""

import cgeniepy
import matplotlib.pyplot as plt

model = cgeniepy.sample_model()
sst = model.get_var("ocn_sur_temp").isel(time=-1)

## GENIE's grid view, with the frame labelled in degrees
sst.plot(projection="sine", colorbar=True, gridline=True, outline=True)

###############################################################################
# On a sine axis the latitude labels crowd together near the poles (60° sits at
# 0.87 of the height), so leaving out the poles reads better.

p = sst.to_GriddedDataVis()
p.aes_dict["gridline_kwargs"]["ylocs"] = [-60, -30, 0, 30, 60]
p.plot(projection="sine", contourf=True, colorbar=True, gridline=True, zebra_frame=True)

###############################################################################
# The default Equal Earth map, for comparison: the same cells, with the polar
# rows bent to fit the globe's outline.

sst.plot(colorbar=True, outline=True)

plt.show()
