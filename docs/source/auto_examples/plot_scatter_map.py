"""
=========================================
Plot the 2D ScatterData
=========================================

This example shows how to plot the 2D ScatterData object. I use a LGM d13C data from Peterson et al. 2014 (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013PA002552) as an example.
"""

from cgeniepy.table import ScatterData
from cgeniepy.plot import CommunityPalette
import matplotlib.pyplot as plt
import subprocess

## Download necessary files, you need to install zenodo_get first
## by `pip install zenodo_get`, or, just download it from the link above
subprocess.call(["zenodo_get", "10.5281/zenodo.13786013", "-o", "~/Downloads/"])

proxy_d13C = ScatterData("~/Downloads/LGM_d13c_CLEAN.txt", delimiter="\t", header=None)
proxy_d13C.data.columns = ["Longitude", "Latitude", "Depth", "d13C","Event"]
proxy_d13C.set_index(["Latitude", "Longitude"])
cmap = CommunityPalette("BuDaRd18").colormap
proxy_d13C.plot(var='d13C', edgecolor='k', cmap=cmap)

plt.show()
