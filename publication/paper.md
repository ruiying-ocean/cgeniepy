---
title: 'cgeniepy: A Python package for analysing cGENIE Earth System Model output'
tags:
  - Python
  - cGENIE
  - Earth System Model
  - Ocean
  - Climate

authors:
  - name: Rui Ying
    orcid: 0000-0001-5811-2388
    corresponding: true
    email: rui.ying@bristol.ac.uk
    affiliation: 1

affiliations:
 - name: School of Earth Sciences, University of Bristol, UK
   index: 1  

date: 17 Feburary 2024
bibliography: paper.bib
---

# Summary

The cGENIE model is an numerical model that simulates the various components of the Earth system (e.g., atmosphere, ocean, land, biosphere, ice sheet) and their interactions [@ridgwell2007]. It has been widely used in studying and reconstructing the past ocean and climate states. Here, I provide a Python package *cgeniepy* for reading, analysing, and visualising the cGENIE model output, and performing the model-data comparison. The package is designed to facilitate the post-simulation analysis for all cGENIE users, as used in my recent studies [@ying2023;@ying2023b]. The package is designed as object-oriented, thus many features are standalone and can be used without cGENIE background.

# Statement of need
Earth System Models (ESM) are the essential tool used to study the mechanism regulating the complex climate and their impacts. cGENIE is such a model with reduced model complexity and enhanced running speed that strengthen its application in paleoceanography studies. For instance, @henehan2019 used the model to study the impact of an extreme event (Cretaceous-Paleogene massive extinction in 66 Million years ago). @pohl2022 used it to study the long-term evolution of ocean oxygen in the past 550 million years. The application of this model has provided invaluable implication to study the evolution of climate and life.

However, despite the power of cGENIE, the analysis of the model output has no consistent and systematic package support yet. Such gap could influence the efficiency and reproducibility of the research, in particular for users who are not familiar with the model output format or need to perform advanced analysis (e.g., model ensemble and model-skill assessment). Python is a open-source modern programming language that has built-in package management system, which makes it easy to install and use. The *cgeniepy* package provides a starting point for the continuously growing cGENIE community to develop and maintain a set of convenient tools, as seen in the other general circulation models (e.g., @romain2023 for the NEMO model and @gael2023 for the MITgcm model).

# Package Design
This package is designed to be modular and object-oriented. The package provide a class `model` to access the cGENIE model output, and two data structure class (`GriddedData` and `ScatterData`) to deal with the model or observational data. These two classes are based on the `xarray.DataArray` and `pandas.DataFrame` respectively, which are the common data structure used in the Python community. The data structure provides corresponding computational methods and visualisation methods. Again, the visualisation is separated from the data structure for decoupling philosophy. Finall, there is a `skill` module to provide the skill score calculation for the model-data comparison and a universal `Interpolator` class to help the users to interpolate the regular or irregular data. Therefore, the package is designed to be general-purpose and not limited to the cGENIE model output.

![A schematic figure showing the structure of `cgeniepy` package and its functionalities. It helps users to access the model output and operate the visualisation and analysis including interpolation, model-data comparison.](fig0.png){width=65%}

# Example

## Access, analyse and visualise the cGENIE model output
Here, I provide an example to show how to use the package to quickly access, analyse and visualise the cGENIE model output.  More examples can be found in the package documentation website (https://cgeniepy.readthedocs.io/en/latest/).

```python
## import the package
from cgeniepy.model import GenieModel

## initialise a model instance bu providing the path to the model output
model = GenieModel("/Users/foo/Downloads/model_experiment_id")

## get the time-slice variable
ocn_sst = model.get_var("ocn_surf_temp")

## plot the last time slice
ocn_sst.isel(time=-1).plot()
```
![The simulated sea surface temperature in the Modern (left) and Paleogene-Eocene Thermal Maximum event (right) in cGENIE and visualised by `cgeniepy`. The data is adapted from @ying2023b and @gutjahr2017 respectively. \label{fig1}](fig1.png)


# Acknowledgements
R.Y. acknowledge the funding from China Scholarship Council (202006380070). R.Y. also thanks Shao Jun for his suggestions on the package.

# References
