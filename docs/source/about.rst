About
=====

What is this package (cgeniepy)?
--------------------------------

This package provides a handful of tools in Python for dealing with the cGENIE model outputs, similar to the existing `MATLAB scripts <https://github.com/derpycode/muffinplot>`_.

.. note::

   This project is under active development.

What is cGENIE?
----------------
cGENIE is an Earth System Model that simulate the climate and element cycle of Earth. It is for scientists to study the different component of the Earth System, such as the ocean, atmosphere, and land, and their complex interactions, particularlly in the long-term scale. The model is used and developed by scientists from around the world. Therefore, it is important to have a package that can deal with the cGENIE model outputs in a consistent and reproducible way.


How was this package created?
------------------------------
I use cGENIE for my PhD research and I kept writing analysing codes in Python to achieve ideas in my mind. For example, I want to plot the data in a different projection. I want to achieve some clean syntax using method chaining (e.g., :code:`data.do_a().do_b().do_c()`). I also want to make all my codes sustianble, reproducible, and performs quickly. Finally, I feel it good to package my codes to the cGENIE community.


What does this package provide?
--------------------------------
Here are the core features of this package:

- **Python-based**: It is based on Python and many other packages that are open-source. It is easy and free to access compared to other commercial software (e.g., MATLAB). For example, you are reading this documentation online, which is powered by GitHub, Sphinx, and Read the doc.
- **Essentile functionalities**: You can read, compute, and plot the cGENIE model outputs using this package.
- **Observational data support**: For wet lab scientists, you can use this package to deal with observational data (e.g., core data downloaded from pangea.de).
- **Model-data comparison**: You can compare the cGENIE outputs with observational data using this package, including finding the nearest valid grid point, calculating the skill metric.


