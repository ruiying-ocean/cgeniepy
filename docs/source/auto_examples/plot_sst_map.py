"""
=========================================
Plot 2D global map from the model output
=========================================

Here use the sea surface temperature as an example to plot the 2D global map from the model output.
"""

import cgeniepy
import matplotlib.pyplot as plt

model = cgeniepy.sample_model()
sst = model.get_var("ocn_sur_temp").isel(time=-1)

## start to plot with customised cmap
sst_plotter = sst.to_GriddedDataVis()
sst_plotter.aes_dict['pcolormesh_kwargs']['cmap'] = plt.get_cmap("Spectral_r", 15)
sst_plotter.plot(colorbar=True, outline=True)

plt.show()
