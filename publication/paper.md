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

However, despite the power of cGENIE, the analysis of the model output has no consistent and systematic package support yet. Python is a open-source modern programming language that has built-in package management system, which makes it easy to install and use. The *cgeniepy* package provides a starting point for the continuously growing cGENIE community to develop and maintain a set of convenient tools, as seen in the other general circulation models (e.g., @romain2023 for the NEMO model and @gael2023 for the MITgcm model).

# Features
## A general interface to the cGENIE model output
The first-order class for the users is the `GenieModel` class, which can be initialised by the path to the cGENIE model output directory. The class provides a `get_var` method to access the  time-slice variable stored in the netCDF format, and the `get_ts` method to render the time-series data stored in the ascii-based table format.  Once being intialised, the class will automatically index the directory, which allows the users to access the target variable without specifying the exact path to the file. Finally, the class also supports the multiple model runs (i.e., model ensemble) and the multiple variable in a single model run. This could be useful to explore the results with different model parameterisations.

```python
## import the package
from cgeniepy.model import GenieModel

## initialise a model instance bu providing the path to the model output
model = GenieModel("/Users/yingrui/Downloads/muffin.CBE.worlg4.BASESFeTDTL.SPIN")

## get the time-slice variable
ocn_po4 = model.get_var("ocn_PO4")
## get the time-series variable
ocn_temp_ts = model.get_ts("ocn_temp")
```

## Data structure
The accessed data is stored in the `GriddedData` object, a customised data structure based on the `xarray.DataArray`. 

Because the package does not assume any specific model configuration, all the operations are independent to the topographies. Like in the \autoref{fig1}, I read and plot the sea surface temperature (SST) in the modern (0 Ma) and PETM (Paleogene-Eocene Thermal Maximum, 55 Ma) model output, adapted from @ying2023b and @gutjahr2017 respectively. It shows the significant warming in the PETM event, as suggested by the recent data assimilation [@tierney2022], 

![The simulated sea surface temperature in the Modern (left) and Paleogene-Eocene Thermal Maximum event (right) visualised by *cgeniepy* package. \label{fig1}](fig1.png)

Additional calculation methods are provided to enrich the functionalities of the package. For example, the `GriddedData` can be interpolated to finer grid using the `interpolate()` method (\autoref{fig2}). One can mask the ocean basin by using `mask_basin()` method. It also provides a feature to allow the users to find the nearest valid value to a given coordinate using the `search_grid()` method.All the calculation will update the array attribution of the `GriddedData` object, thus the users can easily chain the calculation by using the method chaining (i.e., `data.action_a().action_b()` is supported). Because the `GriddedData` is indenpendent to the cGENIE model, in principal the users can also use the package to analyse the other ocean data (e.g., the observational data from the World Ocean Atlas).


```python
## operation: select the last time slice, mask the basin
atlantic_po4 = ocn_po4.isel(time=-1).mask_basin(base='worjh2',basin='Atlantic', subbasin='')
## operation: calculate the zonal mean
atlantic_zonal_mean_po4 = atlantic_po4.mean(dim='lon')
```

The `GriddedData` is a subclass of `ArrayVis` which provides the `plot` method to visualise the data. The plot will automatically select a method to visualise the data based on the dimension of the data. For example, it will plot a map for data with only latitude and longitude dimension and plot a transect plot for data with only depth and latitude dimension. It will also use the unit string and the variable name (if available) to label the axis and the color bar (\autoref{fig2}).

```python
## no interpolation
atlantic_zonal_mean_po4.plot(contour=True, colorbar=True)

## with interpolation
atlantic_zonal_mean_po4.interpolate().plot(contour=True, colorbar=True)
```
![The zonal mean phosphate concentration in the modern Atlantic (left) before and (right) after the linear interpolation.\label{fig2}](fig2.png){width=80%}

## Model-data comparison
A moduel `cgeniepy.skill` is provided to quantify the model-data comparison. It provides several common skill scores (e.g., the root mean square error, the correlation coefficient, the M-score) to quantify the model performance. In common practice, one could use the `search_grid` method to find the nearest valid value to the observational data, and then intialise a `DataFrameComp` object to calculate the skill scores. Alternatively, one can regrid the observational data (e.g., the World Ocean Atlas data) to the cGENIE model grid and call the `ArrayComp` object to calculate the skill scores. 

# Acknowledgements
R.Y. acknowledge the funding from China Scholarship Council (202006380070). R.Y. also thanks Shao Jun for his suggestions on the package.

# References
