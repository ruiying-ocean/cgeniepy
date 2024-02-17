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
release = '0.10.1'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

source_suffix = {
    '.rst': 'restructuredtext',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "furo"
html_static_path = ['_static']
html_logo = "_static/logo.png"
# html_sidebars = {
#     "index": [],
#     "examples/index": [],
#     "**": ["sidebar-nav-bs.html"],
# }

# html_sidebars = {
#     '**': [
#         'globaltoc.html',
#     ]
# }
