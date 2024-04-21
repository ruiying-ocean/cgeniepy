# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

## -- project path --
import sys, os
import datetime
sys.path.insert(0, os.path.abspath('../src/cgeniepy/'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cgeniepy'
copyright = f"2024-{datetime.datetime.now().year}, Rui Ying"
author = 'Rui Ying'
release = '0.12.0'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx_gallery.gen_gallery',
    "sphinx_design",    
]

source_suffix = {
    '.rst': 'restructuredtext',
}

sphinx_gallery_conf = {
     'examples_dirs': '../../examples',   # path to your example scripts
     'gallery_dirs': 'auto_examples',  # path to where to save gallery generated output
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "pydata_sphinx_theme"
html_static_path = ['_static']

html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"

html_theme_options = {
    "navbar_align": "left",
    "logo": {
        "text": "cgeniepy",
    },
    "navbar_center": ["navbar-nav"],    
}
