"""
=========================================
Add observational data on model output
=========================================

"""

import cgeniepy
from cgeniepy.table import ScatterData
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

## prepare model data, replace sample model with your own model
model = cgeniepy.sample_model()
sst = model.get_var("ocn_sur_temp").isel(time=-1)

# Generate fake "observational" data
np.random.seed(319)  # Set seed for reproducibility

# Create random latitude, longitude, and variable values
random_data = pd.DataFrame({
    'lat': np.random.uniform(-90, 90, 100),  # Random latitudes between -90 and 90
    'lon': np.random.uniform(-180, 180, 100),  # Random longitudes between -180 and 180
    'temp': np.random.random(100) * 30  # Random variable values between 0 and 30
})

random_data = ScatterData(random_data)
random_data.set_index(['lat', 'lon'])


## start plotting
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

sst.plot(ax=ax,vmin=0, vmax=35, colorbar=True, outline=True)
random_data.plot(var='temp', edgecolor='k', ax=ax, vmin=0, vmax=35,
                 gridline=False, colorbar=False, land_mask=False)

plt.show()
