"""
=========================================
Detect and visualize ocean basins
=========================================

This example demonstrates how to:
1. List all available ocean basins from the IPCC AR6 regionmask
2. Detect and select specific ocean basins
3. Plot basin-specific data
"""

import cgeniepy
import matplotlib.pyplot as plt
import regionmask
import cartopy.crs as ccrs

# %%
# List all available ocean basins
# ================================
# The basin detection feature uses the IPCC AR6 marine regions defined in regionmask.
# Let's first see what basins are available.

ocean = regionmask.defined_regions.ar6.ocean

print("Available ocean basins in IPCC AR6:")
print("=" * 60)
for idx, (number, name, abbrev) in enumerate(zip(ocean.numbers, ocean.names, ocean.abbrevs)):
    print(f"{number:2d}: {abbrev:5s} - {name}")

# %%
# Key basin indices
# =================
# The most commonly used basin indices are:
#
# - 46: Arctic Ocean (AO)
# - 47: North Pacific Ocean (NPO)
# - 48: Equatorial Pacific Ocean (EPO)
# - 49: South Pacific Ocean (SPO)
# - 50: North Atlantic Ocean (NAO)
# - 51: Equatorial Atlantic Ocean (EAO)
# - 52: Southern Atlantic Ocean (SAO)
# - 53: North Indian Ocean (NIO)
# - 55: Equatorial Indian Ocean (EIO)
# - 56: South Indian Ocean (SIO)
# - 57: Southern Ocean (SO)

# %%
# Visualize a specific basin
# ===========================
# Now let's plot data for a specific basin using the basin index

model = cgeniepy.sample_model()
sst = model.get_var('ocn_sur_temp').isel(time=-1)

# Plot North Pacific Ocean only (basin index 47)
fig, ax = plt.subplots(1, 1, figsize=(10, 5), subplot_kw={"projection": ccrs.PlateCarree()})
sst.sel_modern_basin(47).plot(ax=ax, outline=True, colorbar=True)
ax.set_title('Sea Surface Temperature - North Pacific Ocean')
ax.coastlines()

# %%
# Use basin abbreviation instead of index
# ========================================
# You can also use basin abbreviations (string) instead of indices

fig, ax = plt.subplots(1, 1, figsize=(10, 5), subplot_kw={"projection": ccrs.PlateCarree()})
sst.sel_modern_basin('NAO').plot(ax=ax, outline=True, colorbar=True)
ax.set_title('Sea Surface Temperature - North Atlantic Ocean')
ax.coastlines()

# %%
# Select multiple basins
# ======================
# You can select multiple basins by passing a list of basin indices or abbreviations

# Plot entire Atlantic Ocean (North + Equatorial + Southern)
fig, ax = plt.subplots(1, 1, figsize=(10, 5), subplot_kw={"projection": ccrs.PlateCarree()})
atlantic_basins = ['NAO', 'EAO', 'SAO']
sst.sel_modern_basin(atlantic_basins).plot(ax=ax, outline=True, colorbar=True)
ax.set_title('Sea Surface Temperature - Atlantic Ocean')
ax.coastlines()

# %%
# Combine basin selection with calculations
# ==========================================
# You can chain basin selection with other operations

# Calculate mean SST for North Pacific
npo_mean_sst = sst.sel_modern_basin(47).mean()
print(f"Mean SST in North Pacific Ocean: {npo_mean_sst.data.values:.2f} °C")

# Calculate mean SST for each major basin
basins_to_compare = {
    'North Pacific': 47,
    'North Atlantic': 50,
    'North Indian': 53
}

print("\nMean SST by basin:")
print("=" * 40)
for basin_name, basin_idx in basins_to_compare.items():
    mean_sst = sst.sel_modern_basin(basin_idx).mean()
    print(f"{basin_name:20s}: {mean_sst.data.values:.2f} °C")

# %%
# Plot multiple basins side by side
# ==================================

fig, axes = plt.subplots(2, 2, figsize=(15, 8),
                         subplot_kw={"projection": ccrs.Robinson()})

basins = [
    (47, 'North Pacific', 'NPO'),
    (50, 'North Atlantic', 'NAO'),
    (49, 'South Pacific', 'SPO'),
    (52, 'South Atlantic', 'SAO')
]

for ax, (basin_idx, basin_name, basin_abbr) in zip(axes.flat, basins):
    sst.sel_modern_basin(basin_idx).plot(ax=ax, outline=True, colorbar=False, cmap='coolwarm')
    ax.set_title(f'{basin_name} ({basin_abbr})')
    ax.coastlines()

plt.tight_layout()
plt.show()
