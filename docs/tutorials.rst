Tutorials
=========

This section provides comprehensive tutorials for using Chimera to create multi-source brain parcellations.

Getting Started
---------------

Basic Usage
^^^^^^^^^^^

Chimera provides a command-line interface for creating combined brain parcellations from multiple sources. Here's how to get started:

1. **Basic Command Structure**:

   .. code-block:: bash

      chimera -b <BIDs directory> -d <Derivatives directory> -p <Chimera code>

2. **Understanding Chimera Codes**:
   
   A Chimera code is a sequence of ten one-character identifiers, one for each supra-region:
   
   - **Position 1**: Cortex
   - **Position 2**: Basal ganglia  
   - **Position 3**: Thalamus
   - **Position 4**: Amygdala
   - **Position 5**: Hippocampus
   - **Position 6**: Hypothalamus
   - **Position 7**: Cerebellum
   - **Position 8**: Brainstem
   - **Position 9**: Gyral white matter
   - **Position 10**: White matter

Tutorial 1: Simple Parcellation
-------------------------------

**Objective**: Create a basic brain parcellation using default settings.

**Prerequisites**: 
- BIDs-formatted dataset
- FreeSurfer derivatives directory
- FSL and ANTs installed (if using local processing)

**Step 1: Explore Available Parcellations**

First, let's see what parcellations are available for each region:

.. code-block:: bash

   chimera --regions

This command displays all available parcellation methods for each supra-region with their corresponding codes and citations.

**Step 2: Create Your First Parcellation**

Let's create a parcellation using the Lausanne cortical atlas (L), FreeSurfer basal ganglia (F), Multiple-resolution thalamus (M), FreeSurfer amygdala (F), FreeSurfer hippocampus (I), FreeSurfer hypothalamus (I), FreeSurfer cerebellum (F), FreeSurfer brainstem (I), FreeSurfer gyral white matter (F):

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF

**Step 3: Verify Results**

After processing, you'll find the combined parcellation in your derivatives directory under:

.. code-block:: text

   derivatives/chimera/sub-<subject>/anat/
   ├── sub-<subject>_space-<space>_desc-LFMFIIFIF_dseg.nii.gz
   └── sub-<subject>_space-<space>_desc-LFMFIIFIF_dseg.tsv

Tutorial 2: Multi-Resolution Processing
--------------------------------------

**Objective**: Process multiple subjects with different parcellation codes and scales.

**Step 1: Multiple Parcellation Codes**

You can process multiple parcellation combinations in a single command:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF,SFMFIIFIF,HFMFIIFIF

This creates three different parcellations:
- LFMFIIFIF: Lausanne cortical parcellation
- SFMFIIFIF: Schaefer cortical parcellation  
- HFMFIIFIF: HCP Multi-Modal cortical parcellation

**Step 2: Multi-Scale Processing**

For multi-resolution atlases like Lausanne, specify scales:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -s scale33,scale60,scale125

**Step 3: White Matter Growth**

Control how cortical labels extend into white matter:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -g 0,1,2

This creates parcellations with 0mm, 1mm, and 2mm white matter growth.

Tutorial 3: Subject-Specific Processing
---------------------------------------

**Objective**: Process specific subjects or sessions.

**Step 1: Process Specific Subjects**

Process only certain subjects:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -ids sub-001,sub-002,sub-003

**Step 2: Using Subject List File**

Create a text file ``subjects.txt`` with subject identifiers:

.. code-block:: text

   sub-001_ses-baseline_run-1
   sub-001_ses-followup_run-1
   sub-002_ses-baseline_run-2

Then run:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -ids subjects.txt

Tutorial 4: Advanced Configuration
----------------------------------

**Objective**: Customize processing with advanced options.

**Step 1: Parallel Processing**

Control the number of parallel processes:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -n 8

**Step 2: FreeSurfer Integration**

Specify FreeSurfer subjects directory:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -fr /path/to/freesurfer/subjects

**Step 3: Force Reprocessing**

Overwrite existing results:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -f

**Step 4: Merge Cortical Regions**

Join cortical white and gray matter regions:

.. code-block:: bash

   chimera -b /path/to/bids/dataset -d /path/to/derivatives -p LFMFIIFIF -mctx

Tutorial 5: Python API Usage
----------------------------

**Objective**: Use Chimera programmatically from Python.

**Basic Python Usage**:

.. code-block:: python

   from chimera.chimera import Chimera
   from chimera.config_manager import load_parcellations_info

   # Load available parcellations
   parc_dict, supra_dict = load_parcellations_info()
   
   # Create a Chimera object
   chim = Chimera(parc_code='LFMFIIFIF', scale='scale60')
   
   # Access parcellation information
   print(f"Parcellation code: {chim.parc_code}")
   print(f"Available methods: {chim.methods}")

**Working with Configuration**:

.. code-block:: python

   from chimera.config_manager import _pipeline_info, _set_templateflow_home
   
   # Load pipeline configuration
   pipe_config = _pipeline_info()
   
   # Set templateflow directory
   tflow_home = _set_templateflow_home('/path/to/templateflow')

**Creating Custom Parcellations**:

.. code-block:: python

   from chimera.parcellation import create_extra_regions_parc
   
   # Create parcellation with extra regions
   extra_parc = create_extra_regions_parc('/path/to/aparc+aseg.nii.gz')

Troubleshooting
--------------

**Common Issues**:

1. **FreeSurfer Not Found**: Ensure FreeSurfer is properly installed and sourced
2. **Missing Dependencies**: Install FSL and ANTs for full functionality
3. **Permission Issues**: Check write permissions for derivatives directory
4. **Memory Issues**: Reduce parallel processes with ``-n`` option

**Getting Help**:

Use the help option to see all available commands:

.. code-block:: bash

   chimera --help

**Verbose Output**:

Use verbose mode for debugging:

.. code-block:: bash

   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -v 2