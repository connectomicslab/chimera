Installation
============

Install from PyPI (Recommended)
--------------------------------

The easiest way to install CHIMERA is using pip:

.. code-block:: bash

    pip install chimera-brainparcellation

This will automatically install all required dependencies.

Manual Installation
--------------------

Alternatively, you can install from source:

.. code-block:: bash

    git clone https://github.com/connectomicslab/chimera
    cd chimera
    pip install -e .

Requirements
------------

CHIMERA requires the following external tools:

- FreeSurfer (version > 7.2.0)
- FSL
- ANTs

Python Dependencies
-------------------

All Python dependencies are automatically installed with pip:

- numpy
- pandas  
- scipy
- nibabel
- pybids
- templateflow
- clabtoolkit
- rich