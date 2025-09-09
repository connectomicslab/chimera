#!/usr/bin/env python
# Root conf.py for ReadTheDocs - minimal configuration without m2r2

import os
import sys
sys.path.insert(0, os.path.abspath('.'))

import chimera

# Basic configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode', 
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages'
]

templates_path = ['docs/_templates']
source_suffix = '.rst'
master_doc = 'docs/index'

project = 'CHIMERA: An open source framework for combining multiple parcellations'
copyright = "2024, Yasser Aleman Gomez"
author = "Yasser Aleman Gomez"

version = chimera.__version__
release = chimera.__version__

language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'

# HTML output
html_theme = 'sphinx_rtd_theme'
html_static_path = ['docs/_static']

# Extensions configuration
napoleon_google_docstring = True
napoleon_numpy_docstring = True
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
}

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}