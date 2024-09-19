"""
==============================================
Extract and Plot cGENIE Time Series data
==============================================

This example shows how to read in and plot cGENIE time series data.
"""
import cgeniepy
import matplotlib.pyplot as plt
from cgeniepy.table import ScatterData


## read in the model
model = cgeniepy.sample_model()
temp = model.get_ts("ocn_temp")
o2 = model.get_ts("ocn_O2")

## merge both and convert to ScatterData format
ts_data = ScatterData(temp.merge(o2, on="time (yr)"))
ts_data.set_index("time (yr)")
fig, axs = plt.subplots(2, 1)
ts_data.plot(var="temperature (C)", ax=axs[0], kind='line')
ts_data.plot(var="surface O2 (mol kg-1)", ax=axs[1], kind='line')
