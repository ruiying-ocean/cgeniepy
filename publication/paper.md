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
    affiliation: 1
  
affiliations:
 - name: School of Earth Sciences, University of Bristol, UK
   index: 1
date: 17 Feburary 2024
bibliography: paper.bib

---

# Summary

The cGENIE model is an Earth system model that is used to study the various components of the Earth system (e.g., atmosphere, ocean, land, biosphere, ice sheet) and their interactions. Here, I provide a Python package `cgeniepy` for reading, analysing, and visualising the cGENIE model output, including conducting the model-data comparison. The package is designed to facilitate the post-simulation analysis for the cGENIE users. However, because it is object oriented, many features are standalone and can be useful for users who lacks cGENIE background.

# Statement of need

The cGENIE model has been widely used in climate and oceanographic science studies, particularly in the paleo-context. For instance, Henehan et al. (2020) used the model to study the impact of an extreme event (Cretaceous-Paleogene massive extinction, 66 Million years ago) (ref). Phol et al. (2021) used it to study the long-term evolution of ocean oxygen in the past 540 million years. However, despite its power, the analysis of the model output has no consistent and systematic package support yet. Python is a open-source modern programming language that has built-in package mangement system, which makes it easy to install and use. The `cgeniepy` package aims to provide a starting point for the continously growing cGENIE community to develop and maintain a set of convenient tools, as seen in the other model community (e.g., [xnemogcm](https://github.com/rcaneill/xnemogcm/) for the NEMO model and [MITgcmTools.jl](https://gaelforget.github.io/MITgcmTools.jl/dev/#MITgcmTools.jl) for the MITgcm model).

# Package Design
## Data format
Two core data structures in the package are: `GriddedData` and `ScatterData`.

## An interface to the cGENIE model
In ``


## Visualisation
## Analysis


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

I acknowledge the funding from China Scholarship Council (202006380070).

# References
