Examples
========

This section provides comprehensive examples of using CHIMERA for various scenarios.

Basic Examples
--------------

Single Subject Processing
~~~~~~~~~~~~~~~~~~~~~~~~~

Process a single subject with a specific parcellation:

.. code-block:: python

    from chimera import Chimera
    import os
    
    # Set up paths
    bids_dir = "/data/my_study"
    derivatives_dir = "/data/my_study/derivatives"
    freesurfer_dir = "/data/my_study/derivatives/freesurfer"
    
    # Create CHIMERA instance
    parcellation_code = "DFMIIIFIF"  # Desikan-Killiany + standard options
    chimera_obj = Chimera(parc_code=parcellation_code)
    
    # Process single subject
    subject_id = "sub-001"
    chimera_obj.build_parcellation(
        bids_dir=bids_dir,
        deriv_dir=derivatives_dir,
        fssubj_dir=freesurfer_dir,
        subj_list=[subject_id]
    )

Batch Processing
~~~~~~~~~~~~~~~

Process multiple subjects in parallel:

.. code-block:: python

    from chimera import chimera_parcellation
    
    # Define multiple parcellation codes
    parc_codes = ["DFMIIIFIF", "HFIIIIFIF", "SFMIIIFIF"]
    
    # Process all subjects in BIDS dataset
    for parc_code in parc_codes:
        chimera_parcellation(
            bids_dir="/data/my_study",
            deriv_dir="/data/my_study/derivatives", 
            parc_codes=[parc_code],
            nthreads=8  # Use 8 parallel threads
        )

Advanced Configuration
----------------------

Custom Parcellation Dictionary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use a custom parcellation configuration:

.. code-block:: python

    from chimera import Chimera
    
    # Path to custom parcellation dictionary
    custom_dict = "/path/to/my_supraregions_dictionary.json"
    
    chimera_obj = Chimera(
        parc_code="DFMIIIFIF",
        parc_dict_file=custom_dict
    )

Multi-scale Parcellations
~~~~~~~~~~~~~~~~~~~~~~~~

Work with multi-scale parcellations like Lausanne or Schaefer:

.. code-block:: python

    from chimera import Chimera
    
    # Lausanne multi-scale (L) with scale specification
    chimera_obj = Chimera(
        parc_code="LFMIIIFIF",  # L = Lausanne cortical parcellation
        scale=["3"]  # Use scale 3 
    )
    
    # Schaefer multi-scale (S) with multiple scales
    chimera_obj = Chimera(
        parc_code="SFMIIIFIF",  # S = Schaefer cortical parcellation
        scale=["400"],
        seg=["7n"]  # 400 parcels, 7 networks
    )

Command Line Examples
---------------------

Basic Command Line Usage
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Process all subjects with Desikan-Killiany cortical parcellation
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF
    
    # Process specific subjects
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -ids sub-001,sub-002,sub-003
    
    # Use multiple parcellation codes
    chimera -b /data/study -d /data/study/derivatives \\
            -p DFMIIIFIF,HFIIIIFIF,SFMIIIFIF

Advanced Command Line Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    # Specify FreeSurfer subjects directory
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -fr /data/study/derivatives/freesurfer
    
    # Control parallel processing
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -n 12  # Use 12 threads
    
    # Grow cortical labels into white matter
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -g 2  # Grow 2mm into white matter
    
    # Merge cortical GM and WM regions
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -mctx
    
    # Force overwrite existing results
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -f
    
    # Verbose output
    chimera -b /data/study -d /data/study/derivatives -p DFMIIIFIF \\
            -v 2  # Maximum verbosity

Integration Examples
--------------------

Integration with neuroimaging workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    """
    Example integration with a typical neuroimaging pipeline
    """
    import os
    from chimera import chimera_parcellation
    import nibabel as nib
    import numpy as np
    
    def process_study_with_chimera(study_path, parcellation_codes):
        """
        Complete processing pipeline with CHIMERA parcellations
        """
        bids_dir = study_path
        deriv_dir = os.path.join(study_path, "derivatives")
        
        # Step 1: Create CHIMERA parcellations
        for parc_code in parcellation_codes:
            print(f"Creating parcellation: {parc_code}")
            chimera_parcellation(
                bids_dir=bids_dir,
                deriv_dir=deriv_dir,
                parc_codes=[parc_code],
                nthreads=8
            )
        
        # Step 2: Load results for further analysis
        results = {}
        for parc_code in parcellation_codes:
            parc_files = []
            subjects_dir = os.path.join(deriv_dir, "chimera", f"parc-{parc_code}")
            
            for subject_dir in os.listdir(subjects_dir):
                if subject_dir.startswith("sub-"):
                    parc_file = os.path.join(
                        subjects_dir, subject_dir, "anat",
                        f"{subject_dir}_space-T1w_desc-{parc_code}_dseg.nii.gz"
                    )
                    if os.path.exists(parc_file):
                        parc_files.append(parc_file)
            
            results[parc_code] = parc_files
        
        return results
    
    # Usage
    study_parcellations = process_study_with_chimera(
        "/data/my_study",
        ["DFMIIIFIF", "HFIIIIFIF"]
    )

Quality Control
~~~~~~~~~~~~~~~

.. code-block:: python

    """
    Quality control checks for CHIMERA parcellations
    """
    import nibabel as nib
    import numpy as np
    from chimera import Chimera
    
    def check_parcellation_quality(parcellation_file, expected_regions=None):
        """
        Perform basic quality checks on a parcellation
        """
        # Load parcellation
        img = nib.load(parcellation_file)
        data = img.get_fdata()
        
        # Basic statistics
        unique_labels = np.unique(data)
        num_regions = len(unique_labels) - 1  # Exclude background (0)
        
        print(f"Parcellation: {parcellation_file}")
        print(f"Number of regions: {num_regions}")
        print(f"Label range: {unique_labels.min()} - {unique_labels.max()}")
        
        if expected_regions and num_regions != expected_regions:
            print(f"WARNING: Expected {expected_regions} regions, found {num_regions}")
        
        # Check for gaps in labeling
        expected_range = set(range(int(unique_labels.max()) + 1))
        actual_labels = set(unique_labels.astype(int))
        missing_labels = expected_range - actual_labels
        
        if missing_labels:
            print(f"Missing labels: {sorted(missing_labels)}")
        
        return {
            'num_regions': num_regions,
            'unique_labels': unique_labels,
            'missing_labels': missing_labels
        }
    
    # Usage
    qc_results = check_parcellation_quality(
        "/data/derivatives/chimera/parc-DFMIIIFIF/sub-001/anat/sub-001_space-T1w_desc-DFMIIIFIF_dseg.nii.gz",
        expected_regions=84  # Expected for Desikan-Killiany
    )