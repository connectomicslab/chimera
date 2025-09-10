import os
import sys
sys.path.insert(0, os.path.abspath('..'))

try:
    import chimera
    version = chimera.__version__
    release = chimera.__version__
except ImportError:
    version = '0.2.1'
    release = '0.2.1'

# Project information
project = 'CHIMERA'
copyright = '2024, Yasser Aleman Gomez'
author = 'Yasser Aleman Gomez'

# Extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.githubpages',
    'm2r2',  # Support for Markdown files
]

# m2r2 configuration
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# Templates and static files
templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# HTML theme
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
}

# Language
language = 'en'

# Intersphinx mapping
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
}

# Napoleon settings for Google/NumPy docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# AutoDoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
}