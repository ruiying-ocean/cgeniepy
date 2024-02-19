"""
================
Additional cmap
================

This example shows the available color palettes in the `cgeniepy` package
"""

import numpy as np
import matplotlib.pyplot as plt
from cgeniepy.plot import CommunityPalette

def plot_colormaps(cmaps):
    ncols = 4
    nrows = int(np.ceil(len(cmaps) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, nrows))

    for i, cmap_name in enumerate(cmaps):
        row = i // ncols
        col = i % ncols
        ax = axes[row, col] if nrows > 1 else axes[col]

        # Create a gradient image using the colormap
        gradient = np.linspace(0, 1, 256).reshape(1, -1)
        ax.imshow(gradient, aspect='auto', cmap=CommunityPalette().get_palette(cmap_name))
        ax.set_title(cmap_name, fontsize=14, fontweight='bold')
        ax.axis('off')

    ## remove the unused axes
    for i in range(len(cmaps), ncols * nrows):
        row = i // ncols
        col = i % ncols
        fig.delaxes(axes[row, col])
        
    plt.tight_layout()

# List of colormaps from cgeniepy
cmaps_list = CommunityPalette().avail_palettes()
plot_colormaps(cmaps_list)
