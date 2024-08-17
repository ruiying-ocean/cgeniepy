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

date: 17 August 2024
bibliography: paper.bib
---

# Summary

cGENIE Earth system model is an numerical model that simulates the various components of the Earth system (e.g., atmosphere, ocean, land, biosphere, ice sheet) and their interactions [@ridgwell2007]. It has been widely used in studying and reconstructing the past ocean and climate states. Here, I provide a Python package *cgeniepy* for reading, analysing, and visualising the cGENIE model output, and performing the model-data comparison. The package is designed to facilitate the post-simulation analysis for all cGENIE users, as used in my recent studies [@ying2023;@ying2023b]. The package is designed as object-oriented, thus many features are standalone and can be used without cGENIE background.

# Statement of need
Earth System Models (ESM) are the essential tool used to study the mechanism regulating the complex climate and their impacts. cGENIE is such a model with reduced model complexity and enhanced running speed that strengthen its application in paleoceanography studies. For instance, @henehan2019 used the model to study the impact of an extreme event (Cretaceous-Paleogene massive extinction in 66 Million years ago). @pohl2022 used it to study the long-term evolution of ocean oxygen in the past 550 million years. The application of this model has provided invaluable implication to study the evolution of climate and life.

Despite the power of cGENIE, the analysis of its model output has relied on a collection of Matlab scripts developed by the cGENIE maintainer (https://github.com/derpycode/muffinplot). A systematic package has been long missing. Such gap could influence the efficiency and reproducibility of the research, in particular for users who are not familiar with the Matlab language or need to perform advanced analysis (e.g., model ensemble based analysis). Python is a popular open-source modern programming language that has built-in package management system, which makes it easy to install, use, and advertise. The *cgeniepy* package provides a starting point for the continuously growing cGENIE community to develop and maintain a set of convenient tools, as seen in the other general circulation models (e.g., @romain2023 for the NEMO model and @gael2023 for the MITgcm model).

# Package Design
This package first provides a class `model` to read the cGENIE model output (Figure \autoref{fig:0}). Then the data will in be stored in two data structure class (`GriddedData` and `ScatterData`). The two data structure classes are based on the `xarray.DataArray` and `pandas.DataFrame` respectively, which are common data structure used in the Python community. The two data structures can also be converted to each other easily using built-in methods.

Once the data classes are initialised (read from GENIE model or not), the users can perform basic operations as they do in `xarray.DataArray` and `pandas.DataFrame`. However, additional features are provided such as the publication-ready visualisation. The core of the visualisation is the `GriddedDataVis` and `ScatterDataVis` class, which contains various options to customise the plot. The visualisation is based on the `matplotlib` and `cartopy`.

Another common demand for Earth system model users is the model-data comparison. Thus, I provided a`skill` module to conduct the skill score calculation including the correlation coefficient, root mean square error, and the Taylor diagram (see the Examples).

For cGENIE model specifically, its coarse model output can be interpolated using the  `Interpolator` class (Figure \autoref{fig:0}). This is a wrapper of the `scipy.interpolate` subpackage and its purpose is to help increase the grid resolution and create prettier figure. However, a long-term goal is to incorporate more advanced interpolation methods (e.g., the DIVA method) to make cGENIE model output more comparable to the high-resolution model/observational results.

![A schematic figure showing the structure of `cgeniepy` package and its functionalities. It helps users to access the model output and operate the visualisation and analysis including interpolation, model-data comparison. \label{fig:0}](fig0.png){width=65%}

# Examples

In this section, I provide two examples to show the core functionalities of `cgeniepy`. More examples however can be found in the package documentation website (https://cgeniepy.readthedocs.io/en/latest/).

## Access, analyse and visualise the cGENIE model output

The following codes will be mostly used by the cGENIE users. It initialise the cGENIE model instance, read the sea surface temperature data, and plot the last time slice as map. The data is adapted from @ying2023b and @gutjahr2017 (Figure \autoref{fig:1}). The users can easily change the variable name to access other model output.

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
![The simulated sea surface temperature in the Modern (left) and Paleogene-Eocene Thermal Maximum event (right) in cGENIE and visualised by `cgeniepy`. The data is adapted from @ying2023b and @gutjahr2017 respectively. \label{fig:1}](fig1.png)

## Model-data comparison

In this example, I read the carbon isotope data from the Last Glacial Maximum (21 ka) model, the corresponding data compilation from @peterson2014. Then I find the model data for each core location and add it to the dataframe. Finally, I conduct the model-data comparison and plot the 1:1 lineplot (Figure \autoref{fig:2}). Multiple metrics are calculated and showed in the figure (Figure \autoref{fig:2}).

```python
from cgeniepy.model import GenieModel
from cgeniepy.table import ScatterData

## initialise a model instance
lgm_model = GenieModel("/Users/yingrui/Science/lgm_foram_niche/model/muffin.CBE.GIteiiva.BASESFeTDTL_rb.SPIN", gemflag='biogem')
## get the variable and select the last time slice of the spin-up model
lgm_d13C = lgm_model.get_var("ocn_DIC_13C").isel(time=-1)

## read in the proxy data and construct ScatterData object
proxy_d13C = ScatterData("~/Science/lgm_bcp/data/lgm_d13C.xlsx")
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
