<p align="center">
  <img src="logo.png"/>
</p>

`cgeniepy` is an Python interface to the output of [cGENIE Earth System Model](https://www.seao2.info/mymuffin.html). It is an alternative to the existing [MATLAB libraries](https://github.com/derpycode/muffinplot). It aims to facilitate the post-simulation analysis, including reading and analysing the model output, and observational data.

⚠ `cgeniepy` is still in alpha phase and under active changes.

## Installation

1. Install from [testpypi](https://test.pypi.org/project/cgeniepy/).

```bash
python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy==0.10.1
```

2. Install from GitHub:

```bash
python3 -m pip install git+https://github.com/ruiying-ocean/cgeniepy.git@master
```

## Quickstart
### 0. initialise a model instance
```python
from cgeniepy.model import GenieModel

## single model
model = GenieModel('path_a')

## model ensemble
multi_dirs = ['/path_a/', '/path_b/']
model = GenieModel(path_to_model_output)
```


### 1. Read data
+ netCDF (*.nc)
+ time series (*.res)

```python
## time slice data
model.get_var("XXXXX").array

## timeseries data
model.get_ts("biogem_series_ocn_temp.res")
```

### 2. Data analysis
+ Data subsetting and statistics
+ Model performance
+ Unit-changing operation (e.g., rate to magnitude)

```python
## get zonal average SST of the last model year
zonal_sst = model.get_var("ocn_sur_temp").isel(time=-1).mean(dim='lat')

## North Pacific SST
npac_sst = model.get_var("ocn_sur_temp").select_basin(47).isel(time=-1)
```

### 3. Visualisation
+ 1D line (time series, zonal average)
+ 2D map (including various projections like polar map)
+ 2D cross section
+ 3D (facet)
+ Add a layer of observational data

```python
## simply call `plot` after accessing the data

## map
model.get_var("ocn_sur_temp").isel(time=-1).plot()
```

### 4. Others
+ ECOGEM shortcuts

```python
from cgeniepy.ecology import EcoModel
model = EcoModel(path_to_model_path)
## get all phytoplankton carbon biomass and plot as map
model.get_pft("Phyto", "Biomass", "C").isel(time=-1).plot()
```

## Project Roadmap 🚩

- [ ] Publish the first stable version
- [ ] plot.py 3D facet subplots
- [ ] Show one colorbar in transect plot
- [X] ignore NAs when searching grid 
- [X] use lat/lon/zt from GENIE output
- [X] Documentation Webiste
- [X] figsize influences colorbar length
- [x] create a simple logo
- [X] allow reading an ensemble of models (netcdf & timeseries)
- [X] formatting the ugly unit string
- [ ] observation.py including plot scatter
- [ ] Add global inventory function

## Killer Features
[X] Model-data comparison
[X] Model ensemble support
[X] Search the nearest *valid* grid point

## Citation

```latex
@software{cgeniepy,
  author = {Rui Ying},
  title = {A Python interface to analyse and visualise cGENIE model output},
  url = {https://github.com/ruiying-ocean/cgeniepy/},
  version = {0.10.1},
  date = {2024-02-17},
}
```

## Logo

Logo is designed by me using free **righteous** font.

## Alternative
Alex Phol's [genie_basicdiags](https://github.com/alexpohl/genie_basicdiags/)

## Raise a bug

Please use GitHub's Issues to raise a bug. This makes the issues traceable so that future users having the same problem can find the answer in the public domain.
