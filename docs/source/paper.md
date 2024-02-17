---
title: 'cgeniepy: A Python package for analysing cGENIE Earth System Model output'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Rui Ying
    orcid: 0000-0001-5811-2388
    affiliation: 1
  
affiliations:
 - name: School of Earth Sciences, University of Bristol, UK
   index: 1
date: 4 Feburary 2024
bibliography: paper.bib

---

# Summary

Earth System Models are powerful tools for understanding the various components of the Earth system (e.g., atmosphere, ocean, land, biosphere, ice sheet) and their interactions. The cGENIE model is one of such with intermediate complexity  and reduced computational cost (ref), which makes it suitable for long-term simulations (e.g., in the whole Phanerozoic, 0-540 million years) (ref). The model is written in Fortran and export simulation output in NetCDF format. The cgeniepy package is a Python package that provides a set of convenient tools to analyse and visualise the model output. The package is designed to be user-friendly and to match the increasing research demands.

# Statement of need

`cgeniepy` is an Python package for reading, analysing, and visualising the cGENIE model output. Whilst the model is widely used in international paleoceanography and paleoclimatology communities, the analysis of the model output has no consistent package support. Python is a open source popular programming languages in the scientific community, 

The `cgeniepy` package aims to provide a set of convenient tools to streamline the analysis process and to make it more user-friendly. The package is designed to be flexible and to match the increasing research demands.

Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

ensemble model
model-data comparison and skill assessment

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References