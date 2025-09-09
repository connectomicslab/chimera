Quick Start Guide
=================

This quick start guide will help you get up and running with Chimera in just a few minutes.

Prerequisites
-------------

Before starting, make sure you have:

1. **Python 3.9+** installed
2. **FreeSurfer** (version ≥ 7.2.0) installed and sourced
3. **FSL** installed (optional, for FIRST subcortical segmentation)
4. **ANTs** installed (optional, for advanced registration)
5. A **BIDS-formatted** neuroimaging dataset

Installation
------------

Install Chimera using pip:

.. code-block:: bash

   pip install chimera-brainparcellation

Or install from source:

.. code-block:: bash

   git clone https://github.com/yasseraleman/chimera.git
   cd chimera
   pip install -e .

First Steps
-----------

1. **Check Available Parcellations**
   
   Start by exploring what parcellation methods are available:

   .. code-block:: bash

      chimera --regions

   This will display all available parcellation options for each brain region.

2. **Prepare Your Data**
   
   Ensure your data is in BIDS format:

   .. code-block:: text

      your_dataset/
      ├── sub-001/
      │   └── anat/
      │       └── sub-001_T1w.nii.gz
      ├── sub-002/
      │   └── anat/
      │       └── sub-002_T1w.nii.gz
      └── dataset_description.json

3. **Create Your First Parcellation**
   
   Run Chimera with a simple parcellation code:

   .. code-block:: bash

      chimera -b /path/to/your/bids/dataset \\
              -d /path/to/derivatives \\
              -p LFMFIIFIF

   This command will:
   - Use Lausanne cortical parcellation (L)
   - Use FreeSurfer for basal ganglia (F) 
   - Use Multi-resolution thalamic parcellation (M)
   - Use FreeSurfer for remaining regions (F,I,I,F,I,F)

Understanding the Output
------------------------

After processing, you'll find results in your derivatives directory:

.. code-block:: text

   derivatives/
   └── chimera/
       └── sub-001/
           └── anat/
               ├── sub-001_space-T1w_desc-LFMFIIFIF_dseg.nii.gz
               └── sub-001_space-T1w_desc-LFMFIIFIF_dseg.tsv

The files include:
- ``.nii.gz``: The combined parcellation volume
- ``.tsv``: Label table with region names, indices, and colors

Common Use Cases
----------------

**Process Multiple Subjects:**

.. code-block:: bash

   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -ids sub-001,sub-002,sub-003

**Use Different Scales:**

.. code-block:: bash

   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -s scale60,scale125

**Run in Parallel:**

.. code-block:: bash

   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -n 8

Next Steps
----------

- Read the :doc:`tutorials` for detailed examples
- Check the :doc:`api` for programmatic usage
- See the :doc:`usage` section for all command-line options

Getting Help
------------

If you encounter issues:

1. Check the verbose output: ``chimera ... -v 2``
2. Review the :doc:`tutorials` section
3. Visit the GitHub repository: https://github.com/yasseraleman/chimera
4. Open an issue if you find bugs or have feature requests