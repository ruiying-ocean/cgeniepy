<p align="center">
  <img src="logo.png"/>
</p>

[![Documentation Status](https://readthedocs.org/projects/cgeniepy/badge/?version=latest)](https://cgeniepy.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/cgeniepy.svg)](https://badge.fury.io/py/cgeniepy)
![PyPI - Downloads](https://img.shields.io/pypi/dm/cgeniepy)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![status](https://joss.theoj.org/papers/b08301b8ec79f1da9150cec224da8391/status.svg)](https://joss.theoj.org/papers/b08301b8ec79f1da9150cec224da8391)


`cgeniepy` is a Python package to analyse the output of [cGENIE Earth System Model](https://www.seao2.info/mymuffin.html). It aims to provide a set of convenient tools for the post-simulation work, including analysing the model output, plotting publication-quality figures, and conducting model-data comparison.

⚠ `cgeniepy` is in beta stage with main functionalities being accompolished.



## Installation

1. Install from [PyPI](https://pypi.org/project/cgeniepy/).

```bash
python3 -m pip install cgeniepy
```

2. Install from GitHub:

```bash
python3 -m pip install git+https://github.com/ruiying-ocean/cgeniepy.git@master
```

## Quickstart
I have uploaded a sample model run, which is a preindustrial spinup configuration with marine ecosystem and biogeochemistry enabled.

```python
import cgeniepy

model = cgeniepy.sample_model()
model.get_var('ocn_sur_temp').isel(time=-1).plot(colorbar=True)
```

* If you want to try other cGENIE model runs, you may go to this zenodo record (https://zenodo.org/records/10575295). 
* If you want to try non-cGENIE model, I have also uploaded two example files to here (https://zenodo.org/records/13786014). 


## Documentation

[An online documentation is hosted in readthedoc.](https://cgeniepy.readthedocs.io/en/latest/)



## Citation

Rui Ying. cgeniepy: A Python package for analysing cGENIE Earth System Model output, 20 February 2024, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-3967633/v1]

## Logo

Logo is designed by me using free **righteous** font.

## Alternative
* Prof. Andy Ridgwell's [muffinplot](https://github.com/derpycode/muffinplot)
* Dr. Alex Phol's [genie_basicdiags](https://github.com/alexpohl/genie_basicdiags/)

## Raise a bug

Please use GitHub's Issues to raise a bug. This makes the issues traceable so that future users having the same problem can find the answer in the public domain.

## Contributing

[How to contribute](CONTRIBUTING.md)
