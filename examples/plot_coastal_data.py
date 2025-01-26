"""
=========================================
Plot coastal only data
=========================================

This example shows how to plot the PO4 distribution in each basin.

This example is independent from GENIE's mask
"""

import xarray as xr
import cgeniepy
import numpy as np

model = cgeniepy.sample_model()

sst = model.get_var('ocn_sur_temp').isel(time=-1)

## grid category is a pre-defined 2D data
## 0: coastal
## 1: land
## 2: open ocean
gc = model.grid_category()

## coastal region
xr.where(gc == 0, sst.data, np.nan).plot()
