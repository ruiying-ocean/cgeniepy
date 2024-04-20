"""
=========================================
Plot coastal only data
=========================================

This example shows how to plot the PO4 distribution in each basin.
"""

import xarray as xr
from cgeniepy.model import GenieModel
import numpy as np

pi_model = GenieModel("/Users/yingrui/Science/lgm_bcp/model/muffin.CB.worlg4.BASESFeTDTL.SPIN", gemflag='biogem')

sst = pi_model.get_var('ocn_sur_temp').isel(time=-1)

## grid category is a pre-defined 2D data
## 0: coastal
## 1: land
## 2: open ocean
gc = pi_model.grid_category()

## coastal region
xr.where(gc == 0, sst.data.values, np.nan).plot()
