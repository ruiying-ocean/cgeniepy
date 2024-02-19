"""
================
Plot Time Series
================

This example shows how to read in and plot time series data
"""
from cgeniepy.model import GenieModel
import matplotlib.pyplot as plt

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.historical")
fig, ax = plt.subplots(1,2,figsize=(10,5))

model.get_ts("biogem_series_ocn_temp.res").plot(x='time (yr)',y='temperature (C)',ax=ax[0])
model.get_ts("biogem_series_ocn_O2.res").plot(x='time (yr)',y='global total O2 (mol)',ax=ax[1])

ax[0].set_title("Ocean mean temperature")
ax[1].set_title("Global total ocean oxygen")
plt.show()
