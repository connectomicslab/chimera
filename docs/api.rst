API Reference
=============

This section provides detailed documentation for all CHIMERA modules, classes, and functions.

CHIMERA Package Overview
------------------------

.. automodule:: chimera
   :members:
   :undoc-members:
   :show-inheritance:

Core Classes and Functions
===========================

Chimera Class
-------------

The main class for creating and working with CHIMERA parcellations.

.. autoclass:: chimera.chimera.Chimera
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

**Main Methods:**

.. automethod:: chimera.chimera.Chimera.prepare_templates
.. automethod:: chimera.chimera.Chimera.create_table
.. automethod:: chimera.chimera.Chimera.export_table
.. automethod:: chimera.chimera.Chimera.build_lut_header
.. automethod:: chimera.chimera.Chimera.build_parcellation

Core Functions
--------------

.. autofunction:: chimera.chimera.chimera_parcellation
.. autofunction:: chimera.chimera.main

Configuration Management
=========================

Functions for loading and managing parcellation configurations.

.. automodule:: chimera.config_manager
   :members:
   :undoc-members:
   :show-inheritance:

**Key Functions:**

.. autofunction:: chimera.config_manager.load_parcellations_info
.. autofunction:: chimera.config_manager._pipeline_info
.. autofunction:: chimera.config_manager._set_templateflow_home

Parcellation Tools
==================

Tools for creating and manipulating parcellations.

.. automodule:: chimera.parcellation
   :members:
   :undoc-members:
   :show-inheritance:

**Key Functions:**

.. autofunction:: chimera.parcellation.create_extra_regions_parc
.. autofunction:: chimera.parcellation._mix_side_prop
.. autofunction:: chimera.parcellation._print_availab_parcels

Processing Tools
================

External processing tools and utilities.

.. automodule:: chimera.processing
   :members:
   :undoc-members:
   :show-inheritance:

**Key Functions:**

.. autofunction:: chimera.processing.launch_fsl_first

Data Submodules
===============

Configuration Data
------------------

.. automodule:: chimera.config
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: chimera.config.supraregions
   :members:
   :undoc-members:
   :show-inheritance:

Atlas Data
----------

.. automodule:: chimera.data
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: chimera.data.annot_atlases
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: chimera.data.gcs_atlases
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: chimera.data.vol_atlases
   :members:
   :undoc-members:
   :show-inheritance: