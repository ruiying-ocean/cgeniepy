"""
=========================================
Plot Comparison Between Arrays
=========================================

This example shows how to build the ArrComparison class plot the comparison.
"""

from cgeniepy.skill import ArrComparison
import numpy as np

np.random.seed(2024)
x = np.random.rand(100)
y = x + np.random.normal(0, 0.1, 100)

## calculate skill score
ac = ArrComparison(x, y)
ac.plot(marker='x')
