"""
==============================================
Extract and Plot cGENIE Time Series data
==============================================

This example shows how to read in and plot cGENIE time series data. The example historical run (1750-2022) is available from https://zenodo.org/records/10575295.
"""
from cgeniepy.model import GenieModel
import seaborn as sns
import matplotlib.pyplot as plt
from cgeniepy.table import ScatterData

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.historical")
fig, axs = plt.subplots(1,2,figsize=(8, 3), tight_layout=True)

temp = model.get_ts("ocn_temp")
o2 = model.get_ts("ocn_O2")

## merge both and convert to ScatterData format
ts_data = ScatterData(temp.merge(o2, on="time (yr)"))
ts_data.set_index("time (yr)")
ts_data.plot(var="temperature (C)", ax=axs[0], kind='line')
ts_data.plot(var="surface O2 (mol kg-1)", ax=axs[1], kind='line')
