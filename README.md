This is an open-source Python package used in my Ph.D project to analyse and visualise [cGENIE model](https://www.seao2.info/mymuffin.html) output (mostly in netCDF format).


## Main features 🐛

### 1. visualisation based on Matplotlib
- Map for 2D cGENIE array
- Ocean Transection Contour 
- Taylor diagram

### 2. chemistry system
Parse a molecular formula and get its molecular mass

### 3. unit system
unit conversion using `pint`

### 4. ocean basin detection
input latitude/longitude to get where the site is in.


## Installation 🙂

```bash
python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy==0.4.2
```

## Code style

ReST as docstring style, black as code style.

## Project TODO 🚩

- [] Documentation
- [] Examples
- [] Unit Tests
- [X] add decorator for plot object

## Structure

the essential is to pass the model path string until the code finally get data from it.

Class: GenieModel -> EcoModel; GenieVariable -> PlanktonBiomass/PlanktonExport -> ForamBiomass/ForamExport;
