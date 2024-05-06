"""
====================
Plot Time Series
====================

This example shows how to read in and plot time series data
"""
from cgeniepy.model import GenieModel
import seaborn as sns
import matplotlib.pyplot as plt

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.historical")
fig, axs = plt.subplots(1,2,figsize=(8, 3), tight_layout=True)

temp = model.get_ts("ocn_temp")
o2 = model.get_ts("ocn_O2")

## merge both
merged = temp.merge(o2, on="time (yr)")

## use seaborn to plot
sns.set_theme(context='notebook', style='ticks', palette='deep')
sns.lineplot(data=merged, x="time (yr)", y="temperature (C)", ax=axs[0])
sns.lineplot(data=merged, x="time (yr)", y="global mean O2 (mol kg-1)", ax=axs[1])
sns.despine()
