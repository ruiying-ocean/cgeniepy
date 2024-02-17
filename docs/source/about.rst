About
=====

What is cGENIE?
--------------
cGENIE is an Earth System Model. It is a numerical model that simulates the biogeochemical cycles of carbon, nitrogen, phosphorus, and oxygen in the ocean and sediments. It is designed to be used in a wide range of applications, including the study of past, present, and future climates, the carbon cycle, and the biogeochemical cycles of the ocean. cGENIE is a community model, which means that it is developed and maintained by a group of scientists from around the world. It is freely available to anyone who wants to use it. cGENIE is written in Fortran and is designed to be run on a supercomputer. It is a complex model that requires a lot of computational resources to run. cGENIE is a powerful tool for studying the Earth System, but it is also a complex and difficult model to use. It requires a lot of expertise to set up and run, and it can be difficult to understand the results that it produces. This package is designed to make it easier to work with the outputs of cGENIE.


What is this package (cgeniepy)?
--------------------------------

This package provides a handful of tools for dealing with the cGENIE model outputs. It is not a full implementation of the cGENIE model itself, which means one must run the before trying this package.


How was this package created?
------------------------------
I use cGENIE for my PhD research and I kept writing more analysing codes in Python to achieve ideas in my mind. For example, I want to plot the data in a different projection. I want to achieve some clean syntax using method chaining (e.g., data.do_a().do_b().do_c()). I also want to make all my codes sustianble, reproducible, and performs quickly. Finally, I feel it is a good practice to packaging my codes to the cGENIE community.


Why do I need to use this package?
----------------------------------
Actually, you don't need to. You can always write your own codes to achieve the same goal. But, if you mean what are the core features of this package, here are some of them:

- **Python-based**: It is based on Python and many other packages that are open-source. It is easy and free to access compared to other commercial software (e.g., MATLAB).
- **Essentile functionalities**: You can read, compute, and plot the cGENIE model outputs using this package.
- **Model Ensemble**: You can deal with an ensemble of models using this package.
- **Observational data support**: For wet lab scientists, you can use this package to deal with observational data (e.g., core data downloaded from pangea.de).
- **Model-data comparison**: You can compare the cGENIE outputs with observational data using this package, including finding the nearest valid grid point, calculating the skill metric.


