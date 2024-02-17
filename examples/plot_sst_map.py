"""
====================================
Plot 2D map from the model output
====================================

I plot the sea surface temperature
"""

from cgeniepy.model import GenieModel
import matplotlib.pyplot as plt

model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.worlg4.BASESFeTDTL.SPIN")
sst = model.get_var("ocn_sur_temp").isel(time=-1)
sst.plot()
plt.show()

