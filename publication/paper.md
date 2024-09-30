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

date: 19 Sep 2024
bibliography: paper.bib
---

# Summary

cGENIE is a numerical model that simulates the Earth system (e.g., atmosphere, ocean, land, biosphere, ice sheet, and their interactions) in different geological ages [@ridgwell2007]. It has been widely used in studying and reconstructing the past ocean and climate states. Here, I provide a Python package *cgeniepy* for reading, analysing, and visualising the cGENIE model output, and performing the model-data comparison. The package is designed to facilitate the post-simulation analysis for all cGENIE users, as used in my recent studies [@ying2023;@ying2023b]. The package is designed as object-oriented, thus many features are standalone and can be used without cGENIE background.

# Statement of need
Earth System Models are the essential tool used to study the mechanisms regulating the complex climate and their impacts. cGENIE is such a model with intermediate model complexity that strengthens its application in paleoceanography studies. For instance, @henehan2019 used cGENIE to study the impact of an extreme climatic event, the Cretaceous-Paleogene massive extinction in 66 Million years ago. @pohl2022 used cGENIE to study the long-term evolution of ocean oxygen over the past 550 million years. The application of this model has promoted our understanding of climate change in the geological past.

Despite the power of cGENIE, the analysis of its model output has relied on a collection of MATLAB scripts developed by the cGENIE maintainer (https://github.com/derpycode/muffinplot). A systematic package has been long missing. Such gap might hamper the efficiency of the research, in particular for users who are not familiar with MATLAB or need to perform customised analysis (e.g., model ensemble based analysis). 

Python is a popular open-source programming language that has a built-in package management system. Therefore, relative to MATLAB, Python packages can be easier to install, use, and demonstrate across platforms based on Jupyter notebooks and Quarto. Many climate and ocean models have their own Python package support (e.g., @romain2023 for the NEMO model and @gael2023 for the MITgcm model). As such, it is useful to develop a similar one for the growing cGENIE community.


# Package Design
This package first provides a class `model` to read the cGENIE model output (\autoref{fig:0}). Then the accessed data will in be stored in corresponding data structure class (`GriddedData` or `ScatterData`). The two data structure classes are based on the `xarray.DataArray` and `pandas.DataFrame` respectively, which are common data structures used in the Python community. The two data structures can also be converted to each other easily using built-in methods.

Once the data classes are initialised (read from GENIE model or not), the users can perform basic operations from `xarray.DataArray` and `pandas.DataFrame`. However, additional features are also provided such as publication-ready visualisations in the `GriddedDataVis` and `ScatterDataVis` classes, which both contain various options to customise plots based on the `matplotlib` and `cartopy` packages.

Another common demand for Earth system model users is model-data comparisons. Thus, I provided a `skill` module to conduct the skill score calculation including the correlation coefficient, root mean square error, and the Taylor diagram (see the Examples section).

For the cGENIE model specifically, its coarse model output can be interpolated using the  `Interpolator` class (\autoref{fig:0}). This is a wrapper of the `scipy.interpolate` subpackage and its purpose is to help increase the grid resolution and create prettier figures. However, a long-term goal is to incorporate more advanced interpolation methods (e.g., the DIVA method) to make cGENIE model output more comparable to the high-resolution model/observational results.

![A schematic figure showing the structure of the `cgeniepy` package and its functionalities. `cgeniepy` helps users to access the model output and operate the visualisation and analysis, including interpolation and model-data comparison. \label{fig:0}](fig0.png){width=65%}

# Examples

In this section, I provide two examples to show the core functionalities of `cgeniepy`. More examples however can be found in the package documentation website (https://cgeniepy.readthedocs.io/en/latest/).

## Access, analyse and visualise the cGENIE model output

The following code example demonstrates using `cgeniepy` in a common use case for cGENIE users. It initialises the cGENIE model instance, reads the sea surface temperature data, and plots the last time slice as a map. The data is adapted from @ying2023b and @gutjahr2017 (\autoref{fig:1}). The users can easily change the variable name to access other model outputs.

```python
## import the package
from cgeniepy.model import GenieModel

## initialise a model instance bu providing the path to the model output
model = GenieModel("/Users/foo/model_experiment_id")

## get the time-slice variable
ocn_sst = model.get_var("ocn_surf_temp")

## plot the last time slice
ocn_sst.isel(time=-1).plot()
```
![The simulated sea surface temperature in the Modern (left) and Paleogene-Eocene Thermal Maximum event (right) in cGENIE and visualised by `cgeniepy`. The data is adapted from @ying2023b and @gutjahr2017 respectively. \label{fig:1}](fig1.png)

## Model-data comparison

In this example, I demonstrate how to conduct a model-data comparison using `cgeniepy`. I compare the ocean carbon isotope in sediment cores (i.e., observation) [@peterson2014] with cGENIE model ouptputs [@ying2023b] in the Last Glacial Maximum (21 ka).  First, both observations and model results are read by the `cgeniepy` package. Then I search the nearest cGENIE model grid boxes for each sediment core and append the matched results to the existing dataframe. Finally, I visualise the model-data comparison by plotting the scatter plot with multiple metrics calculated (\autoref{fig:2}).

```python
from cgeniepy.model import GenieModel
from cgeniepy.table import ScatterData

## The example model and data are archived in
## https://zenodo.org/records/13786013 and
## https://zenodo.org/records/8189647

## initialise a model instance
lgm_model = GenieModel("path/to/the/model/")
## get the variable and select the last time slice of the spin-up model
lgm_d13C = lgm_model.get_var("ocn_DIC_13C").isel(time=-1)

## read in the proxy data and construct ScatterData object
proxy_d13C = ScatterData("path/to/the/data")
proxy_d13C.set_index(["Lat", "Lon", "Depth"])

## find the model data for each core location
model_data = []
for i in proxy_d13C.data.index:
    lat, lon, depth = i
    pos = (depth, lat, lon)
    data = lgm_d13C.search_point(pos, ignore_na=True)
    model_data.append(data)

## add the model data to the dataframe
proxy_d13C.data["GENIE_d13C"] = model_data

## rename the column for y axis label
proxy_d13C.data.rename(columns={"LGM":"Observational d13C"}, inplace=True)

## conduct the model-data comparison and plot the 1:1 lineplot
## by default, model data is in the col, and observational col is in the second
proxy_d13C.compare("GENIE_d13C","LGM").plot()
```

![A `cgeniepy` example of model-data comparison for the Last Glacial Maximum carbon isotope data. The model output is adapted from @ying2023b and the data is from @peterson2014. \label{fig:2}](fig2.png)

# Acknowledgements
R.Y. acknowledges financial support from China Scholarship Council (202006380070). R.Y. also thanks Shao Jun for his suggestions on the package.

# References
