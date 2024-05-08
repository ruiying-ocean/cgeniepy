"""
=========================================
Plot the 1D ScatterData
=========================================

This example shows how to plot the 1D ScatterData object. I use a CO2 data file from the Antarctic EDC ice core (https://doi.pangaea.de/10.1594/PANGAEA.472488) as an example.

The following features are used:

#. Load the data (which just download from Pangaea without any modification)

#. Plot the raw data

#. Plot the interpolated data (based on cubic spline interpolation)

#. Plot the rolling averaged data
"""

from cgeniepy.table import ScatterData
import matplotlib.pyplot as plt

# Load data
edc_co2 = ScatterData("/Users/yingrui/cgeniepy/src/data/EDC_CO2.tab", sep='\t')
edc_co2.set_index(['Age [ka BP]'])

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111)

# Plot the raw data
edc_co2.plot(var='CO2 [ppmv]', ax=ax, 
        label='Raw Data', kind='scatter',
        edgecolor='black', facecolor='none', marker='o')

# # Plot the interpolated data (based on cubic spline interpolation)
interpolated_data = edc_co2.interpolate(var='CO2 [ppmv]')
interpolated_data = ScatterData(interpolated_data)
interpolated_data.plot(var='CO2 [ppmv]', ax=ax, label='Interpolated', kind='line')

# # Use rolling mean to smooth the data
smoothed_data = edc_co2.rolling(window=2).mean()
smoothed_data = ScatterData(smoothed_data)
smoothed_data.plot(var='CO2 [ppmv]', ax=ax, label='Smoothed', kind='line')

ax.legend()
plt.show()
