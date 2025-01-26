"""
==============================================
Add paleogeography mask in scatter data plot
==============================================

This example shows how to add paleogeography mask according to input age
"""
from cgeniepy.table import ScatterData
import pandas as pd
import numpy as np

# Generate random data
np.random.seed(1239124)  # Set seed for reproducibility

# Create random latitude, longitude, and variable values
random_data = pd.DataFrame({
    'lat': np.random.uniform(-90, 90, 100),
    'lon': np.random.uniform(-180, 180, 100),
    'dummy_var': np.random.random(100) * 100 
})

# Convert to ScatterData object
random_data = ScatterData(random_data)

# Set MultiIndex for the data
random_data.set_index(['lat', 'lon'])

# Plot with specified parameters
random_data.plot(var='dummy_var', edgecolor='k', mask_age=100)
