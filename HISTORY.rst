=======
History
=======

0.4.0 (2026-07-28)
------------------

* Decomposed the ``chimera.py`` monolith into focused modules: ``core``
  (the ``Chimera`` class), ``template_preparation``, ``region_table``,
  ``parcellation_builder`` and ``cli``.
* Rewrote the 1,780-line ``build_parcellation`` as a ``ParcellationBuilder``
  with named steps, replacing ``"<name>" in locals()`` state tracking with
  explicit attributes and threading ``pipe_dict`` explicitly.
* Moved the console entry point to ``chimera.cli:main`` and removed the
  legacy ``chimera.py``.
* Consolidated the PyPI/TestPyPI publish workflows into one, gated on tags
  (TestPyPI) and GitHub Releases (PyPI); bumped artifact actions to v4.
* Dropped the redundant root ``requirements.txt`` in favour of
  ``pyproject.toml`` and ``docs/requirements.txt``.

0.3.1 (2026-02-20)
------------------

* Corrected typos in the supraregions files.
* Pinned ``clabtoolkit>=0.4.2``.

0.3.0 (2026-02-19)
------------------

* Added opacity support to parcellations and supra-region tables.
* Added atlas entity prefix to parcel names.
* Updated color table export using new clabtoolkit API.
* Replaced deprecated ``codes2mask`` argument with ``mask_codes``.
* Replaced deprecated ``names2look`` argument with ``names2keep``.
* Migrated ``group_by_code`` to use a grouping dictionary interface.
* Renamed ``rearrange_parc`` method to ``rearrange`` throughout.
* Used ``ColorTableLoader`` class for color table manipulation.
* Adapted to updated ``Parcellation`` class that loads color tables as dictionaries.
* Pinned ``clabtoolkit>=0.4.1`` dependency.
* Comprehensive documentation updates.

0.1.0 (2024-01-01)
------------------

* First release on PyPI.
* Initial implementation of brain parcellation fusion framework.
* Support for 10 supra-regions with multiple atlas options.
* Command-line interface for processing BIDS datasets.
* Integration with FreeSurfer, FSL, and ANTs tools.
