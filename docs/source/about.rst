About
=====

What is this package (cgeniepy)?
--------------------------------

This package provides a handful of tools in Python for dealing with the cGENIE model outputs, similar to the existing `MATLAB scripts <https://github.com/derpycode/muffinplot>`_.

What is cGENIE?
----------------
cGENIE is an Earth System Model that simulate the climate and element cycle of Earth. It is for scientists to study the different component of the Earth System, such as the ocean, atmosphere, and land, and their complex interactions, particularly in the long-term scale. The model is used and developed by scientists from around the world. Therefore, it is important to have a package that can deal with the cGENIE model outputs in a consistent and reproducible way.


What does this package provide?
--------------------------------
Here are the core features of this package:

- **Easy to access**: It is based on Python ecosystem and included in the Python's package index pool (PyPI). You can download this software for free in any operating system, remotely or locally. You have a online documentation website with examples for you to use.
  
- **cGENIE analysis**: This package aims to analyse the cGENIE model outputs easily. This includes time slice, time series, and ensemble model results.

- **Observational data**: You can compare the cGENIE outputs with observational data using this package, including finding the nearest valid grid point, bin the observational data (gridded or scatter data), calculating the skill metric, plotting Taylor diagram.
  
- **Independent functionality**: The package can be used without cGENIE knowledge. This includes the `GriddedData` and `ScatterData` classes.


How was this package created?
--------------------------------
I used cGENIE in my PhD and I needed to efficiently analyse the model results as my wish.


Citation
------------------
You can cite this peer-reviewed paper:

Ying, R. (2024). cgeniepy: A Python package for analysing cGENIE Earth System Model output. Journal of Open Source Software, 9(101), 6762. https://doi.org/10.21105/joss.06762
