**CHIMERA**: An open source framework for combining multiple parcellations
=============================================================================

Creating multi-source parcellations of the human brain is a fundamental task at several steps of the MRI analysis research workflow. **Chimera** facilitates this otherwise difficult operation with an intuitive and flexible interface for humans and machines, thereby assisting in the construction of sophisticated and more reliable processing pipelines.
This repository contains the source code and atlases needed by **Chimera**.

.. image:: Figure1.png
   :alt: CHIMERA Framework Overview
   :align: center

📖 Documentation
=================

Full documentation is available at: **https://chimera-brainparcellation.readthedocs.io**

The documentation includes:

- Complete API reference
- Installation guide  
- Usage examples
- Parcellation methodology details

Parcellations fusion
====================

Chimera defines ten different supra-regions (cortex, subcortical structures, thalamus, amygdala, hippocampus, hypothalamus, cerebellum, brainstem, gyral white matter, and white-matter). Subcortical structures include only the regions that are not labeled as supra-regions. Subdivisions in each supra-region will be populated with the parcellation information of a single source. The available parcellation sources per supra-region, as well as one corresponding parcellation name, and a one-character unique identifier are configured in a JSON (JavaScript Object Notation) file.

**Chimera code**: A sequence of ten one-character identifiers (one per each supra-region) unambiguously denotes a single instance of combined parcellation (Figure 1B). Given the sequence of ten identifier characters, Chimera selects the atlas and/or applies the corresponding methodology to obtain the parcellation for each supra-region. These supra-region-specific parcellations are finally integrated to obtain the combined volumetric parcellation for each input subject, as well as its corresponding tab-separated values table of labels, region names, and rendering colors for visualization.

Chimera uses FreeSurfer to map cortical templates from fsaverage to individual space. It also applies different methods to obtain the hippocampal subfields and brainstem parcellations as well as the thalamic, amygdala and hypothalamic nuclei segmentations. FIRST and ANTs are also used for segmenting subcortical structures and thalamic nuclei respectively.

Requirements
============

Required Python Packages
-------------------------

Standard Library (Built-in, no installation required)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `argparse <https://docs.python.org/3/library/argparse.html>`_ - Command-line argument parsing
- `csv <https://docs.python.org/3/library/csv.html>`_ - CSV file reading and writing
- `datetime <https://docs.python.org/3/library/datetime.html>`_ - Date and time handling
- `json <https://docs.python.org/3/library/json.html>`_ - JSON encoder and decoder
- `operator <https://docs.python.org/3/library/operator.html>`_ - Standard operators as functions
- `os <https://docs.python.org/3/library/os.html>`_ - Operating system interface
- `pathlib <https://docs.python.org/3/library/pathlib.html>`_ - Object-oriented filesystem paths
- `shutil <https://docs.python.org/3/library/shutil.html>`_ - High-level file operations
- `subprocess <https://docs.python.org/3/library/subprocess.html>`_ - Subprocess management
- `sys <https://docs.python.org/3/library/sys.html>`_ - System-specific parameters and functions
- `time <https://docs.python.org/3/library/time.html>`_ - Time access and conversions
- `typing <https://docs.python.org/3/library/typing.html>`_ - Support for type hints

Data Science & Analysis
~~~~~~~~~~~~~~~~~~~~~~~~

- `numpy <https://pypi.org/project/numpy/>`_ - Fundamental package for scientific computing
- `pandas <https://pypi.org/project/pandas/>`_ - Data manipulation and analysis library
- `scipy <https://pypi.org/project/scipy/>`_ - Scientific computing library

Neuroimaging & Medical Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- `nibabel <https://pypi.org/project/nibabel/>`_ - Access to neuroimaging file formats
- `pybids <https://pypi.org/project/pybids/>`_ - BIDS (Brain Imaging Data Structure) toolkit
- `templateflow <https://pypi.org/project/templateflow/>`_ - Neuroimaging template management

CLI & User Interface
~~~~~~~~~~~~~~~~~~~~

- `rich <https://pypi.org/project/rich/>`_ - Rich text and beautiful formatting for terminals

Specialized Tools
~~~~~~~~~~~~~~~~~

- `clabtoolkit <https://pypi.org/project/clabtoolkit/>`_ - Connectomics Lab Toolkit

Installation
============

Install from PyPI (Recommended)
--------------------------------

The easiest way to install CHIMERA is using pip:

.. code-block:: bash

    pip install chimera-brainparcellation

This will automatically install all required dependencies including:

- pandas
- pybids  
- numpy
- nibabel
- rich
- scipy
- templateflow
- clabtoolkit

Manual Installation
-------------------

Alternatively, you can install all required external packages manually:

.. code-block:: bash

    pip install pandas pybids numpy nibabel rich scipy templateflow clabtoolkit

Or using a requirements.txt file:

.. code-block:: bash

    pip install -r requirements.txt

requirements.txt content:
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

    pandas
    pybids
    numpy
    nibabel
    rich
    scipy
    templateflow
    clabtoolkit

Required image processing packages:

- `FreeSurfer (version>7.2.0) <https://surfer.nmr.mgh.harvard.edu/>`_
- `FSL <https://fsl.fmrib.ox.ac.uk/fsl/fslwiki>`_
- `ANTs <http://stnava.github.io/ANTs/>`_

Options
=======

Brief description of input options:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Option
     - Description
   * - ``--regions``, ``-r``
     - List available parcellations for each supra-region.
   * - ``--bidsdir``, ``-b``
     - BIDs dataset folder. Different BIDs directories could be entered separating them by a comma.
   * - ``--derivdir``, ``-d``
     - Derivatives folder. Different directories could be entered separating them by a comma.
   * - ``--parcodes``, ``-p``
     - Sequence of ten one-character identifiers (one per each supra-region).
   * - ``--freesurferdir``, ``-fr``
     - FreeSurfer subjects dir. If the folder does not exist it will be created.
   * - ``--scale``, ``-s``
     - Scale identification. This option should be supplied for multi-resolution cortical parcellations (e.g. Lausanne or Schaeffer).
   * - ``--seg``, ``-e``
     - Segmentation identifier.
   * - ``--nthreads``, ``-n``
     - Number of processes to run in parallel (default= Number of cores - 4).
   * - ``--growwm``, ``-g``
     - Grow of GM labels inside the white matter (mm).
   * - ``--subjids``, ``-ids``
     - Subject IDs. Multiple subject ids can be specified separating them by a comma.
   * - ``--mergectx``, ``-mctx``
     - Join cortical white matter and cortical gray matter regions.
   * - ``--force``, ``-f``
     - Overwrite the results.
   * - ``--verbose``, ``-v``
     - Verbose (**0**, **1** or **2**).
   * - ``--help``, ``-h``
     - Help.

Usage
=====

General command line to use **Chimera**:

.. code-block:: bash

    $ chimera -b <BIDs directory> -d <Derivatives directory> -p <Chimera code>

This command will run Chimera for all the subjects in the BIDs directory.

Simple examples
---------------

1. Running **Chimera** for 3 different parcellation codes (LFMFIIFIF,SFMFIIFIF,CFMFIIFIF). This will obtain the combined parcellations for all the T1-weighted images inside the BIDs dataset.

.. code-block:: bash

    $ chimera -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF,SFMFIIFIF,CFMFIIFI

2. Running **Chimera** for T1-weighted images included in a txt file:

.. code-block:: bash

    $ chimera -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF -ids <t1s.txt>

Example of **t1s.txt** file::

    sub-00001_ses-0001_run-2
    sub-00001_ses-0003_run-1
    sub-00001_ses-post_acq-mprage

3. Cortical volumes will grow 0 and 2 mm respectively inside the white matter for the selected cortical parcellations.

.. code-block:: bash

    $ chimera -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF -g 0,2

Main files in the repository
============================

1. **chimera.py**: Main python library for performing **Chimera** parcellations.
2. **supraregions_dictionary.json**: JSON file specifying the available parcellation sources per supra-region.
3. **annot_atlases** and **gcs_atlases**: Folder containing cortical atlases in *.annot* and *.gcs* file formats.

For detailed information about available parcellations for each supra-region, see the :doc:`parcellations` page.

Results
=======

Chimera parcellations were generated using the following codes: LFMIIIFIF, HFIIIIFIF, BFIIHIFIF (162, 492 and
314 regions respectively). Figure 2A shows the corresponding results of the fused parcellations for a single
subject. By filtering each individual's tractogram with the corresponding Chimera parcellations, we generated
connectivity matrices (Figure 2B).

.. image:: Figure2.png
   :alt: CHIMERA Results
   :align: center

License
=======

.. image:: https://img.shields.io/badge/License-Apache_2.0-blue.svg
   :target: https://opensource.org/licenses/Apache-2.0
   :alt: License