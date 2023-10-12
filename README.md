This package is an interface to the output of [cGENIE Earth System Model](https://www.seao2.info/mymuffin.html). It is an alternative to the exisiting [MATLAB libraries](https://github.com/derpycode/muffinplot). It aims to do three things:

+ Read data
+ Data analysis
+ Visualisation

## Installation 🙂

```bash
python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy==0.7.0
```

## Read data
+ netCDF (*.nc)
+ time series (*.res)

## Data analysis
+ subset basin
+ mean/max/min/median/sum over lat/lon/depth
+ unit conversion (e.g., mmol C/m3 -> Tg C)

## Visualisation
+ plot map/transection/zonal average/polar according to dimension
+ optionally use contour lines
+ add layer of observational data points [TODO]

## Others


## Gallery

### A global biomass map of modelled picophytoplankton (0.6 μm) 

![map](example_map.png)

### A global distribution of basin-level nutrient (PO4) 

![modern_po4](example_transection.png)



### ECOGEM shortcuts

## Project TODO 🚩

- [] Build Documentation
- [] Write Examples


## Citation
@software{cgeniepy,
  author = {Rui Ying},
  title = {A Python interface to analyse and visualise cGENIE model output},
  url = {https://github.com/ruiying-ocean/cgeniepy/},
  version = {0.7.0},
  date = {2023-10-09},
}
