<p align="center">
  <img src="logo.png"/>
</p>

[![Documentation Status](https://readthedocs.org/projects/cgeniepy/badge/?version=latest)](https://cgeniepy.readthedocs.io/en/latest/?badge=latest)
![PyPI](https://img.shields.io/pypi/v/PACKAGE?label=pypi%20cgeniepy)
![PyPI - Downloads](https://img.shields.io/pypi/dm/cgeniepy)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`cgeniepy` is an Python package to analyse the output of [cGENIE Earth System Model](https://www.seao2.info/mymuffin.html). It aims to provide a set of convenient tools for the post-simulation analysis, including analysing the model output, plot publication-quality figures, and conducting model-data comparison.

⚠ `cgeniepy` is still in alpha phase and under active changes.



## Installation

1. Install from [PyPI](https://pypi.org/project/cgeniepy/).

```bash
python3 -m pip install cgeniepy
```

2. Install from GitHub:

```bash
python3 -m pip install git+https://github.com/ruiying-ocean/cgeniepy.git@master
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
