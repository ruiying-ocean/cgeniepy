About
=====

What is this package (cgeniepy)?
--------------------------------

This package provides a handful of tools in Python for dealing with the cGENIE model outputs, similar to the existing `MATLAB scripts <https://github.com/derpycode/muffinplot>`_.

.. note::

   This project is under continuous development.

What is cGENIE?
----------------
cGENIE is an Earth System Model that simulate the climate and element cycle of Earth. It is for scientists to study the different component of the Earth System, such as the ocean, atmosphere, and land, and their complex interactions, particularly in the long-term scale. The model is used and developed by scientists from around the world. Therefore, it is important to have a package that can deal with the cGENIE model outputs in a consistent and reproducible way.


What does this package provide?
--------------------------------
Here are the core features of this package:

- **Python-based**: It is based on Python and many other packages that are open-source. It is easy and free to access compared to other commercial software (e.g., MATLAB). For example, you are reading this documentation online, which is powered by GitHub, Sphinx, and Read the doc.
- **Essential functionalities**: You can read, compute, and plot the cGENIE model outputs using this package.
- **Observational data support**: For wet lab scientists, you can use this package to deal with observational data (e.g., core data downloaded from pangea.de).
- **Model-data comparison**: You can compare the cGENIE outputs with observational data using this package, including finding the nearest valid grid point, bin the observational data (gridded or scatter data), calculating the skill metric, plotting Taylor diagram.


How was this package created?
------------------------------
I use cGENIE for my PhD research and I kept writing analysing codes in Python to implement ideas in my mind. I also want to make all my codes sustainable, reproducible, without losing performance. Finally, these codes become a package and I think it fun to share with the cGENIE community.


Citation
------------------
You can cite this preprint (https://www.researchsquare.com/article/rs-3967633/v1).
