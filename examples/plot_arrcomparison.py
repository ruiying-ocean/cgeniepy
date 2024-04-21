"""
=========================================
Plot Comparison Between Arrays
=========================================

This example shows how to use the ArrComparison class to compare two arrays and plot the comparison.
"""

from cgeniepy.skill import ArrComparison
import numpy as np

np.random.seed(1923)
x = np.random.rand(36, 36)
y = x + np.random.rand(36, 36)

## calculate skill score
ac = ArrComparison(x, y)
print("Pearson coefficient:", ac.pearson_r())
ac.plot()
