=====
Usage
=====

Command Line Interface
----------------------

Chimera provides a powerful command-line interface for creating brain parcellations. The basic syntax is:

.. code-block:: bash

   chimera -b <BIDs directory> -d <Derivatives directory> -p <Chimera code>

Quick Start Example
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   # List available parcellations
   chimera --regions
   
   # Create a basic parcellation
   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF

Common Options
^^^^^^^^^^^^^^

.. code-block:: bash

   # Process specific subjects
   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -ids sub-001,sub-002
   
   # Use multiple scales
   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -s scale60,scale125
   
   # Control white matter growth
   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -g 0,1,2
   
   # Run in parallel with 8 processes
   chimera -b /path/to/bids -d /path/to/derivatives -p LFMFIIFIF -n 8

Python API
----------

You can also use Chimera programmatically in Python:

.. code-block:: python

   import chimera
   from chimera.chimera import Chimera
   from chimera.config_manager import load_parcellations_info
   
   # Load available parcellations
   parc_dict, supra_dict = load_parcellations_info()
   
   # Create a Chimera object
   chim = Chimera(parc_code='LFMFIIFIF', scale='scale60')
   
   # Access parcellation information
   print(f"Parcellation code: {chim.parc_code}")

Understanding Chimera Codes
---------------------------

A Chimera code is a 10-character string where each position corresponds to a specific brain region:

1. **Cortex** (L=Lausanne, S=Schaefer, H=HCP, etc.)
2. **Basal ganglia** (F=FreeSurfer, R=FIRST)
3. **Thalamus** (F=FreeSurfer, M=Multi-resolution, I=Iglesias)
4. **Amygdala** (F=FreeSurfer, I=Iglesias, R=FIRST)
5. **Hippocampus** (F=FreeSurfer, I=Iglesias, H=Hippocampal subfields)
6. **Hypothalamus** (F=FreeSurfer, I=Iglesias)
7. **Cerebellum** (F=FreeSurfer)
8. **Brainstem** (F=FreeSurfer, I=Iglesias, R=FIRST)
9. **Gyral white matter** (F=Based on cortical parcellation)
10. **White matter** (F=FreeSurfer)

For detailed tutorials and examples, see the :doc:`tutorials` section.
