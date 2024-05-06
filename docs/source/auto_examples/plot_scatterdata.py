"""
=========================================
Plot the 1D ScatterData
=========================================

This example shows how to plot the 1D ScatterData object. I use a CO2 data file from the Antarctic EDC ice core (https://doi.pangaea.de/10.1594/PANGAEA.472488) as an example.
"""

from cgeniepy.table import ScatterData
import matplotlib.pyplot as plt

# Load data
edc_co2 = ScatterData("/Users/yingrui/cgeniepy/src/data/EDC_CO2.tab", sep='\t')
edc_co2.set_index(['Age [ka BP]'])

# Create subplots
fig, ax = plt.subplots()

# Plot the raw data
edc_co2.to_ScatterDataVis().plot(var='CO2 [ppmv]', ax=ax, label='Raw Data', edgecolor='black', facecolor='none', marker='o')
# Plot the interpolated data (based on cubic spline interpolation)
interpolated_data = edc_co2.interpolate(var='CO2 [ppmv]').to_dataframe()
interpolated_data.plot(ax=ax, x='Age [ka BP]', y='interpolated_values', label='Interpolated', linewidth=2)

# Use rolling mean to smooth the data
smoothed_data = edc_co2.rolling(window=2).mean()
smoothed_data['CO2 [ppmv]'].plot(ax=ax, label='2pt Rolling Mean', linewidth=2)

ax.legend()
ax.grid(True)
ax.minorticks_on()

plt.show()

