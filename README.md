![logo](logo.png)

This package is an interface to [cGENIE model](https://www.seao2.info/mymuffin.html). Although there are exisiting [MATLAB codes](https://github.com/derpycode/muffinplot) to do this, I feel more comfortable when I have more freedom and control. So I developed this package. I try to use my limited Python knowledge in this package, but it is inevitable that many errors may exist.

## Main features 🐛

### 1. visualisation based on Matplotlib
- Global map
- Polar map (despite low resolution)
- Basin transection
- Taylor diagram (from community)
- Copy some colormaps that was not in matplotlib, e.g., the ODV and MATLAB default color palette (thanks to the `rpal` package in R) 

### 2. chemistry system
Parse a molecular formula and get its molecular mass

### 3. unit system
unit conversion using `pint`

### 4. ocean basin detection
input latitude/longitude to get where the site is in.

### 5. Functional diversity
Calculate community weight mean, functional dispersion, functional specisation rate, over-redundancy.


## Installation 🙂

```bash
python3 -m pip install -i https://test.pypi.org/simple/ cgeniepy==0.6
```

## Code style

ReST as docstring style, black as code style.

## Project TODO 🚩

- [] Documentation built by `sphinx` 
- [] Examples
- [] Unit Tests
- [X] add decorator for plot object

## Structure

the essential is to pass the model path string until the code finally get data from it.

Class: GenieModel -> EcoModel; GenieVariable -> PlanktonBiomass/PlanktonExport -> ForamBiomass/ForamExport;
