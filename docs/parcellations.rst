Parcellations and Methodologies
===============================

This section describes all available parcellations for each supra-region in CHIMERA.

Overview
--------

CHIMERA defines ten different supra-regions where each can be parcellated using different methodologies. The parcellation code consists of 10 characters, one for each supra-region, allowing you to create custom combinations.

Available Parcellations by Supra-region
---------------------------------------

1. Cortical (Supra-region: Cortical)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The first character of the parcellation code selects the cortical parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 20 40 10 20 40

   * - Code
     - Name
     - Citation
     - Code
     - Name
     - Citation
   * - A
     - AALv2
     - Rolls et al, 2015
     - B
     - Brainnetome
     - Fan et al, 2016
   * - C
     - Campbell
     - Campbell, 1905
     - D
     - Desikan-Killiany
     - Desikan et al, 2006
   * - F
     - Flechsig
     - Flechsig, 1920
     - H
     - HCP-MMP1
     - Glasser et al, 2016
   * - K
     - Kleist
     - Kleist, 1934
     - L
     - Lausanne
     - Symmetric version of Cammoun et al, 2012
   * - M
     - Smith
     - Smith et al, 1907
     - R
     - Broadmann
     - Broadmann, 1909
   * - S
     - Schaefer
     - Schaefer et al, 2018
     - T
     - Desikan-Killiany-Tourville
     - Klein and Tourville, 2012
   * - V
     - vonEconomo
     - von Economo and Koskinas, 1925
     - X
     - Destrieux
     - Destrieux et al, 2009
   * - Y
     - Yeo
     - Yeo et al, 2011
     -
     -
     -

2. Subcortical (Supra-region: Subcortical)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second character selects the subcortical structures parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Fischl et al, 2002
   * - R
     - FIRST
     - Patenaude et al, 2011

3. Thalamus (Supra-region: Thalamus)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The third character selects the thalamic parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 20 40 10 20 40

   * - Code
     - Name
     - Citation
     - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Fischl et al, 2002
     - I
     - FSThalParc
     - Iglesias et al, 2018
   * - M
     - MIALThalParc
     - Najdenovska and Aleman-Gomez et al, 2018
     - R
     - FIRST
     - Patenaude et al, 2011

4. Amygdala (Supra-region: Amygdala)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fourth character selects the amygdala parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Fischl et al, 2002
   * - I
     - FSAmygHippoParc
     - Saygin et al, 2017
   * - R
     - FIRST
     - Patenaude et al, 2011

5. Hippocampus (Supra-region: Hippocampus)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fifth character selects the hippocampal parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 20 40 10 20 40

   * - Code
     - Name
     - Citation
     - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Fischl et al, 2002
     - I
     - FSAmygHippoParc
     - Iglesias et al, 2015
   * - H
     - HBT
     - Iglesias et al, 2015
     - R
     - FIRST
     - Patenaude et al, 2011

6. Hypothalamus (Supra-region: Hypothalamus)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sixth character selects the hypothalamic parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Based on in-house protocol
   * - I
     - FSHypoThalParc
     - Billot et al, 2020

7. Cerebellum (Supra-region: Cerebellum)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The seventh character selects the cerebellar parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - A
     - AALv2
     - Rolls et al, 2015
   * - F
     - Aseg
     - Fischl et al, 2002
   * - S
     - SUIT
     - Diedrichsen, J. et al, 2009

8. Brainstem (Supra-region: Brainstem)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The eighth character selects the brainstem parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Fischl et al, 2002
   * - I
     - FSBrainStemParc
     - Iglesias et al, 2015
   * - R
     - FIRST
     - Patenaude et al, 2011

9. Gyral White Matter (Supra-region: GyralWM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ninth character selects the gyral white matter parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - F
     - Cortical
     - Depends on the cortical parcellation

10. White Matter (Supra-region: WhiteMatter)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tenth character selects the white matter parcellation method.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - Code
     - Name
     - Citation
   * - F
     - Aseg
     - Fischl et al, 2002
   * - J
     - JHU
     - Hua et al, 2008

Example Parcellation Codes
---------------------------

Here are some example codes and their meanings:

**DFMIIIFIF** (Recommended default):
  - D: Desikan-Killiany cortical parcellation
  - F: FreeSurfer subcortical structures
  - M: MIAL thalamic nuclei
  - I: FreeSurfer amygdala nuclei
  - I: FreeSurfer hippocampal subfields
  - I: FreeSurfer hypothalamic regions
  - F: FreeSurfer cerebellum
  - I: FreeSurfer brainstem regions
  - F: Cortical-based gyral white matter
  - F: FreeSurfer white matter

**HFIIIIFIF** (High-resolution cortical):
  - H: HCP-MMP1 high-resolution cortical (360 regions)
  - F: FreeSurfer subcortical structures
  - I: FreeSurfer thalamic nuclei
  - (rest same as above)

**SFMIIIFIF** (Schaefer cortical):
  - S: Schaefer cortical parcellation (customizable scale)
  - F: FreeSurfer subcortical structures
  - M: MIAL thalamic nuclei
  - (rest same as above)

Command Line Usage
------------------

To see all available parcellations:

.. code-block:: bash

    chimera --regions

To use a specific parcellation code:

.. code-block:: bash

    chimera -b /path/to/bids -d /path/to/derivatives -p DFMIIIFIF