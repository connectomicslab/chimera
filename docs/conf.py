import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# Debug: Print what we're loading
print("DEBUG: Loading chimera module...")
try:
    import chimera
    print(f"DEBUG: Successfully imported chimera, version: {chimera.__version__}")
except Exception as e:
    print(f"DEBUG: Failed to import chimera: {e}")
    # Fallback values
    class FakeChimera:
        __version__ = "0.2.1"
    chimera = FakeChimera()

project = 'CHIMERA'
copyright = '2024, Yasser Aleman Gomez'
author = 'Yasser Aleman Gomez'
version = chimera.__version__
release = chimera.__version__

# Minimal extensions only - NO m2r2!
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
]

print(f"DEBUG: Extensions to load: {extensions}")

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

language = 'en'