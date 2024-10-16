"""
======================================================
Evaluation of 3D ocean biogeochemistry
======================================================

This example uses WOA13/GLODAP data to evaluate the performance of cGENIE model.

The regridded data is downloaded from https://www.seao2.info/mymuffin.html
"""

import cgeniepy
from cgeniepy.skill import ArrComparison
from cgeniepy.skill import TaylorDiagram
import matplotlib.pyplot as plt

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

tay_diagram = TaylorDiagram(arr_comps)
tay_diagram.setup_ax(crmse_contour=True)
tay_diagram.plot(s=100)
plt.show()
