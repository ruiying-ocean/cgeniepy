"""
=======================================
Map frames in the style of NCL and GMT
=======================================

On a rectangular projection, ``gridline=True`` draws the frame NCL draws: tick marks
pointing out on all four sides, labelled in degrees, with minor ticks between them.
``zebra_frame=True`` adds GMT's fancy black-and-white frame on top of those ticks.
"""

import cgeniepy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

model = cgeniepy.sample_model()
sst = model.get_var("ocn_sur_temp").isel(time=-1)

fig, axes = plt.subplots(2, 1, figsize=(6, 8), layout="constrained", subplot_kw={"projection": ccrs.PlateCarree()})
## NCL: filled contours with the labelbar below the frame
sst.plot(ax=axes[0], contourf=True, colorbar=True, gridline=True, outline=True)
## GMT: the fancy frame
sst.plot(ax=axes[1], colorbar=True, gridline=True, zebra_frame=True, outline=True)

plt.show()
