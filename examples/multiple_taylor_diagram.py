"""
===========================
Customise Taylor Diagram
===========================

This examples shows hot to plot multiple model results in one figure

For simplicity, we separate different variables of the sample model into two lists and we pretend
that they are two models.
"""
import matplotlib.pyplot as plt

import cgeniepy
from cgeniepy.skill import ArrComparison, TaylorDiagram

sample_model = cgeniepy.sample_model()
bgc_data = cgeniepy.load_obs('worjh2')

genie_var  = ['ocn_O2', 'ocn_ALK','ocn_PO4', 'ocn_temp', 'ocn_sal', 'ocn_SiO2', 'ocn_DIC']
obs_var = ['o2', 'alk', 'po4', 'temp', 'sal', 'si', 'dic']

arr_comps = []
for x,y in zip(genie_var, obs_var):
    obs = bgc_data[y]
    model = sample_model.get_var(x).isel(time=-1)
    tmp = ArrComparison(model.data.to_numpy(), obs.to_numpy(), label=y)
    arr_comps.append(tmp)

# you can set the figure size and dpi here
fig = plt.figure(dpi=100,figsize=(4,4))

# Assume this your first model
td = TaylorDiagram(arr_comps[:3])
td.setup_ax(crmse_contour=True,fig=fig)
td.plot(s=100, cmap=plt.get_cmap('tab10',3))

# Assume this your second model
td.ac = arr_comps[3:6]
td.extract_data() # refresh the data
td.plot(s=100, cmap=plt.get_cmap('tab10',3),linestyle='--')

plt.show()
