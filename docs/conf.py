import os
import sys
sys.path.insert(0, os.path.abspath('..'))

import chimera

project = 'CHIMERA'
copyright = '2024, Yasser Aleman Gomez'
author = 'Yasser Aleman Gomez'
version = chimera.__version__
release = chimera.__version__

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

language = 'en'