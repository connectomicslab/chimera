Quickstart Guide
================

This guide will help you get started with CHIMERA quickly.

Installation
------------

Install CHIMERA from PyPI:

.. code-block:: bash

    pip install chimera-brainparcellation

Basic Usage
-----------

Create a CHIMERA Parcellation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here's a simple example of how to create a CHIMERA parcellation:

.. code-block:: python

    from chimera import Chimera
    
    # Create a CHIMERA instance with a parcellation code
    # Each character represents a supra-region parcellation choice
    parc_code = "DFMIIIFIF"  # 10-character code for all supra-regions
    
    chimera_obj = Chimera(parc_code=parc_code)
    
    # Prepare templates and build parcellation
    chimera_obj.prepare_templates(fssubj_dir="/path/to/freesurfer/subjects")
    chimera_obj.build_parcellation(
        bids_dir="/path/to/bids",
        deriv_dir="/path/to/derivatives",
        fssubj_dir="/path/to/freesurfer/subjects"
    )

Command Line Usage
~~~~~~~~~~~~~~~~~~

CHIMERA can also be used from the command line:

.. code-block:: bash

    # Basic usage
    chimera -b /path/to/bids -d /path/to/derivatives -p DFMIIIFIF
    
    # Multiple parcellation codes
    chimera -b /path/to/bids -d /path/to/derivatives -p DFMIIIFIF,SFMIIIFIF
    
    # Specific subjects
    chimera -b /path/to/bids -d /path/to/derivatives -p DFMIIIFIF -ids sub-001,sub-002

Understanding Parcellation Codes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CHIMERA uses 10-character codes where each character represents the parcellation choice for a specific supra-region:

1. **Position 1**: Cortical parcellation (A, B, C, D, F, H, K, L, M, R, S, T, V, X, Y)
2. **Position 2**: Subcortical parcellation (F, R)
3. **Position 3**: Thalamus parcellation (F, I, M, R)
4. **Position 4**: Amygdala parcellation (F, I, R)
5. **Position 5**: Hippocampus parcellation (F, I, H, R)
6. **Position 6**: Hypothalamus parcellation (F, I)
7. **Position 7**: Cerebellum parcellation (A, F, S)
8. **Position 8**: Brainstem parcellation (F, I, R)
9. **Position 9**: Gyral White Matter parcellation (F)
10. **Position 10**: White Matter parcellation (F, J)

Example codes:
- ``DFMIIIFIF``: Desikan-Killiany cortical + FreeSurfer subcortical + MIAL thalamus + ...
- ``HFIIIIFIF``: HCP-MMP1 cortical + FreeSurfer subcortical + ...

List Available Parcellations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To see all available parcellations for each supra-region:

.. code-block:: bash

    chimera --regions

Working with Results
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    # Create and export lookup table
    chimera_obj.create_table()
    chimera_obj.export_table(out_basename="my_parcellation", format=["tsv", "json"])
    
    # Build LUT header
    header = chimera_obj.build_lut_header()
    print(header)

Next Steps
----------

- Check the :doc:`api` for detailed function documentation
- See the :doc:`installation` guide for advanced installation options
- Read about parcellation methodologies in the main README