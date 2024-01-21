# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cgeniepy'
copyright = '2023, Rui Ying'
author = 'Rui Ying'
release = '0.7.4'

## source code path
import os
import sys

# Assuming conf.py is inside docs/source/
# Add the path to the directory containing your project to sys.path
sys.path.insert(0, os.path.abspath('../../src'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

## use myst_nb to render html
## use autodoc to generate API documentation (docstrings)
extensions = ["myst_nb", "sphinx.ext.autodoc"]

source_suffix = {
    '.rst': 'restructuredtext',
    '.ipynb': 'myst-nb',
    '.myst': 'myst-nb',
}

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_logo = "logo.png"
html_static_path = ['_static']
