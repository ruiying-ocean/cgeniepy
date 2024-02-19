<p align="center">
  <img src="logo.png"/>
</p>

`cgeniepy` is an Python interface to the output of [cGENIE Earth System Model](https://www.seao2.info/mymuffin.html). It is an alternative to the existing [MATLAB libraries](https://github.com/derpycode/muffinplot). It aims to facilitate the post-simulation analysis, including reading and analysing the model output, and observational data.

⚠ `cgeniepy` is still in alpha phase and under active changes.

[![Documentation Status](https://readthedocs.org/projects/cgeniepy/badge/?version=latest)](https://cgeniepy.readthedocs.io/en/latest/?badge=latest)

## Installation

1. Install from [testpypi](https://test.pypi.org/project/cgeniepy/).

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
