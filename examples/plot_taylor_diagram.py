"""
===========================
Plot Taylor Diagram
===========================

This example uses fake data and fitted model prediction to plot taylor diagram

The R package open-air has provided a nice explanation on what taylor diagram is:
     https://bookdown.org/david_carslaw/openair/sections/model-evaluation/taylor-diagram.html (Figure 20.2)
"""
import numpy as np
from scipy.optimize import curve_fit
from cgeniepy.skill import ArrComparison, TaylorDiagram

def generate_data(x, a, b, c, noise=0.5):
    y = a * np.exp(-b * x) + c  # Exponential function
    np.random.seed(90148)
    noise = np.random.normal(0, noise, size=len(x))
    return y + noise

def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c

def linear_func(x, a, b):
    return a * x + b

# Generate random data
x = np.linspace(0, 10, 50)
y = generate_data(x, 5, 0.3, 2)

## linear model
popt, pcov = curve_fit(linear_func, x, y)
fit1 = linear_func(x, *popt)

## exp model
popt2, pcov2 = curve_fit(exp_func, x, y)
fit2 = exp_func(x, *popt2)

## Create Comparison instance
ac1 = ArrComparison(y, fit1, 'linear')
ac2 = ArrComparison(y, fit2, 'exponential')

## Create TaylorDiagram instance
diagram = TaylorDiagram([ac1, ac2])
diagram.setup_ax(crmse_contour=True)
diagram.plot()
