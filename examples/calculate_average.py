"""
==================
Calculate average
==================

This example shows how to calculate weighted or unweighted average of a variable
"""
from cgeniepy.model import GenieModel
import numpy as np

pi_model = GenieModel("/Users/yingrui/Science/lgm_bcp/model/muffin.CB.worlg4.BASESFeTDTL.SPIN", gemflag='biogem')

o2 = pi_model.get_var('ocn_O2').isel(time=-1) ##mol/kg
ocn_vol = pi_model.grid_volume().isel(time=-1) ##m3

print("average of o2 weighted by ocean grid volume", o2.weighted_average(ocn_vol.array.values))

## unweighted average of o2
print("unweighted average of o2", o2.array.mean().values)
