`cgeniepy` is an Python interface to the output of [cGENIE Earth System Model](https://www.seao2.info/mymuffin.html). It is an alternative to the exisiting [MATLAB libraries](https://github.com/derpycode/muffinplot). It aims to do three things:

+ Read model output
+ Data analysis
+ Data visualisation
  
## Installation

`cgeniepy` is still under active development. But welcome to try the feature and download cgeniep from [testpypi](https://test.pypi.org/project/cgeniepy/).

```bash
python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy==0.7.1
```

## Usage
### 1. Read data
+ netCDF (*.nc)
+ time series (*.res)

```python
## Example
from cgeniepy.model import GenieModel
model = GenieModel(path_to_model_output)
model.get_var("XXXXX").array
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
zonal_sst.plot()
npac_sst.plot()
```

### 4. Others
+ ECOGEM shortcuts

```python
from cgeniepy.ecology import EcoModel
model = EcoModel(path_to_model_path)
## get all phytoplankton carbon biomass and plot as map
model.get_pft("Phyto", "Biomass", "C").isel(time=-1).plot()
```

## Gallery

### A global biomass map of modelled picophytoplankton (0.6 μm) 

![map](example_map.png)

### A global distribution of basin-level nutrient (PO4) 

![modern_po4](example_transection.png)

## Project Roadmap 🚩

- [ ] Publish the first stable version to `pypi`
- [ ] Examples and Documentation
- [ ] plot.py 3D facet subplots
- [ ] plot.py more dependency of data.dimension
- [ ] plot.py scatter data overlay
- [ ] create a simple logo
- [ ] ignore NAs when searching grid 
- [ ] use lat/lon/zt from GENIE output
- [ ] Move the computation functions in ecology.py and foram.py to array.py
- [X] allow reading an ensemble of models (netcdf & timeseries)
- [X] formatting the ugly unit string

## Citation

```latex
@software{cgeniepy,
  author = {Rui Ying},
  title = {A Python interface to analyse and visualise cGENIE model output},
  url = {https://github.com/ruiying-ocean/cgeniepy/},
  version = {0.7.1},
  date = {2023-10-09},
}
```
