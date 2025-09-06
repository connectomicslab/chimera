===============================
Chimera
===============================

.. image:: https://img.shields.io/pypi/v/chimera.svg
        :target: https://pypi.python.org/pypi/chimera

.. image:: https://img.shields.io/travis/yasseraleman/chimera.svg
        :target: https://travis-ci.com/yasseraleman/chimera

.. image:: https://readthedocs.org/projects/chimera/badge/?version=latest
        :target: https://chimera.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/badge/License-Apache_2.0-blue.svg
        :target: https://opensource.org/licenses/Apache-2.0
        :alt: License

An open source framework for combining multiple parcellations of the human brain.

* Free software: Apache Software License 2.0
* Documentation: https://chimera.readthedocs.io.

Overview
--------

Creating multi-source parcellations of the human brain is a fundamental task at several steps of the MRI analysis research workflow. **Chimera** facilitates this otherwise difficult operation with an intuitive and flexible interface for humans and machines, thereby assisting in the construction of sophisticated and more reliable processing pipelines. This repository contains the source code and atlases needed by **Chimera**.

Parcellations Fusion
--------------------

Chimera defines ten different supra-regions (cortex, basal ganglia, thalamus, amygdala, hippocampus, hypothalamus, cerebellum, brainstem, gyral white matter, and white-matter). Basal ganglia includes only the regions that are not labeled as supra-regions. Subdivisions in each supra-region will be populated with the parcellation information of a single source. The available parcellation sources per supra-region, as well as one corresponding parcellation name, and a one-character unique identifier are configured in a JSON (JavaScript Object Notation) file.

**Chimera code**: A sequence of ten one-character identifiers (one per each supra-region) unambiguously denotes a single instance of combined parcellation. Given the sequence of ten identifier characters, Chimera selects the atlas and/or applies the corresponding methodology to obtain the parcellation for each supra-region. These supra-region-specific parcellations are finally integrated to obtain the combined volumetric parcellation for each input subject, as well as its corresponding tab-separated values table of labels, region names, and rendering colors for visualization.

Chimera uses FreeSurfer to map cortical templates from fsaverage to individual space. It also applies different methods to obtain the hippocampal subfields and brainstem parcellations as well as the thalamic, amygdala and hypothalamic nuclei segmentations. FIRST and ANTs are also used for segmenting subcortical structures and thalamic nuclei respectively.

Installation
------------

Using pip::

    pip install chimera

From source::

    git clone https://github.com/yasseraleman/chimera.git
    cd chimera
    pip install -e .

Requirements
-----------

Required Python packages:

- numpy
- pandas  
- nibabel
- pybids
- scipy
- clabtoolkit
- templateflow
- rich

Required image processing packages:

- FreeSurfer (version>7.2.0)
- FSL
- ANTs

Usage
-----

Brief description of input options:

==================  ============================================
Option              Description
==================  ============================================
``--regions``, ``-r``     List available parcellations for each supra-region
``--bidsdir``, ``-b``     BIDs dataset folder
``--derivdir``, ``-d``    Derivatives folder
``--parcodes``, ``-p``    Sequence of ten one-character identifiers
``--growwm``, ``-g``      Grow of GM labels inside the white matter (mm)
``--t1file``, ``-t``      File containing the basenames of T1w images
``--force``, ``-f``       Overwrite the results
``--verbose``, ``-v``     Verbose (0, 1 or 2)
``--help``, ``-h``        Help
==================  ============================================

General command line usage::

    chimera -b <BIDs directory> -d <Derivatives directory> -p <Chimera code>

Examples
--------

1. Running Chimera for 3 different parcellation codes (LFMFIIFIF,SFMFIIFIF,CFMFIIFIF). This will obtain the combined parcellations for all the T1-weighted images inside the BIDs dataset::

    chimera -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF,SFMFIIFIF,CFMFIIFI

2. Running Chimera for T1-weighted images included in a txt file::

    chimera -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF -t <t1s.txt>

Example of **t1s.txt** file::

    sub-00001_ses-0001_run-2
    sub-00001_ses-0003_run-1
    sub-00001_ses-post_acq-mprage

3. Cortical volumes will grow 0 and 2 mm respectively inside the white matter for the selected cortical parcellations::

    chimera -b <BIDs directory> -d <Derivatives directory> -p LFMFIIFIF -g 0,2

Main Files in the Repository
----------------------------

1. **chimera_parcellation.py**: Main python library for performing **Chimera** parcellations.
2. **parcTypes.json**: JSON file specifying the available parcellation sources per supra-region.
3. **ANNOT_atlases** and **GCS_atlases**: Folder containing cortical atlases in *.annot* and *.gcs* file formats.
4. **mni_icbm152_t1_tal_nlin_asym_09c**: Folder containing the reference atlas used by the MIAL atlas-based thalamic parcellation method. The atlas is referenced in standard MNI (Montreal Neurological Institute) space with a high resolution T1 weighted image (ICBM 2009c Nonlinear Asymmetric).
5. **thalamic_nuclei_MIALatlas**: Folder containing the spatial probabilistic maps of 14 thalamic nuclei.

Parcellations and Methodologies for Each Supra-region
-----------------------------------------------------

1. Cortical parcellation
~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``D``     Desikan et al, 2006            ``X``     Destrieux et al, 2009
``T``     Klein and Tourville, 2012      ``B``     Fan et al, 2016
``R``     Broadmann, 1909                ``C``     Campbell, 1905
``K``     Kleist, 1934                   ``L``     Symmetric version of Cammoun et al, 2012
``H``     Glasser et al, 2016            ``S``     Schaefer et al, 2018
``M``     Smith et al, 1907              ``V``     von Economo and Koskinas, 1925
``Y``     Yeo et al, 2011                ``F``     Flechsig, 1920
========  =============================  ========  ================================

2. Basal ganglia parcellation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``F``     Fischl et al, 2002             ``R``     Patenaude et al, 2011
========  =============================  ========  ================================

3. Thalamus parcellation
~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``F``     Fischl et al, 2002             ``I``     Iglesias et al, 2018
``M``     Najdenovska and Alemán-Gómez   ``R``     Patenaude et al, 2011
          et al, 2018
========  =============================  ========  ================================

4. Amygdala parcellation
~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``F``     Fischl et al, 2002             ``I``     Saygin et al, 2017
``R``     Patenaude et al, 2011
========  =============================  ========  ================================

5. Hippocampus parcellation
~~~~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``F``     Fischl et al, 2002             ``I``     Iglesias et al, 2015
``H``     Iglesias et al, 2015           ``R``     Patenaude et al, 2011
========  =============================  ========  ================================

6. Hypothalamus parcellation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``F``     Based on in-house protocol     ``I``     Billot et al, 2020
========  =============================  ========  ================================

7. Cerebellum parcellation
~~~~~~~~~~~~~~~~~~~~~~~~~

========  =============================
Code      Citation
========  =============================
``F``     Fischl et al, 2002
========  =============================

8. Brainstem parcellation
~~~~~~~~~~~~~~~~~~~~~~~~

========  =============================  ========  ================================
Code      Citation                       Code      Citation
========  =============================  ========  ================================
``F``     Fischl et al, 2002             ``I``     Iglesias et al, 2015
``R``     Patenaude et al, 2011
========  =============================  ========  ================================

9. Gyral white matter parcellation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

========  =============================
Code      Citation
========  =============================
``F``     Cortical (Depends on the cortical parcellation)
========  =============================

Results
-------

Chimera parcellations were generated using the following codes: LFMIIIFIF, HFIIIIFIF, BFIIHIFIF (162, 492 and 314 regions respectively). The corresponding results of the fused parcellations show the integration of multiple atlases for comprehensive brain parcellation. By filtering each individual's tractogram with the corresponding Chimera parcellations, connectivity matrices can be generated for further analysis.

License
-------

This project is licensed under the Apache Software License 2.0.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage