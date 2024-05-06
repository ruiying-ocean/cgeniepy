"""
=========================================
Plot the 1D ScatterData
=========================================

This example shows how to plot the 1D ScatterData object. I use a data file from the EDC ice core (https://doi.pangaea.de/10.1594/PANGAEA.472488) as an example.
"""

from cgeniepy.table import ScatterData
data= ScatterData("/Users/yingrui/cgeniepy/src/data/EDC_CO2.tab", sep='\t')
data.set_index(['Age [ka BP]'])
data.to_ScatterDataVis().plot(var='CO2 [ppmv]')
