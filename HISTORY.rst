=======
History
=======

Unreleased
----------

* Fixed the codes of the cortical white matter regions, which could be wrong
  in the left hemisphere.  ``wm-lh-X`` must have the code of ``ctx-lh-X``
  plus 3000, but its code was computed from the number of regions in the
  right non-cortical table (``nrh_subc``), while the code of ``ctx-lh-X``
  comes from the labels in the image.  When the table listed regions that
  had no voxels in the image, the left white matter codes were shifted (by
  14 in a cohort that had to be corrected by hand).  The codes of
  ``wm-lh-X`` and ``wm-rh-X`` are now set from the codes of their cortical
  regions by name once all the regions are in the image, and the region
  table (``region_table.py``) derives them the same way instead of counting.
* Every parcellation is now checked before it is saved: the color table lists
  have the same length, no name or code is repeated, every code of the image
  is in the color table, and every ``wm-<hemi>-X`` region has a
  ``ctx-<hemi>-X`` region and its code plus 3000.  An inconsistent
  parcellation raises a ``ValueError`` listing the problems and is not
  written.  The checks and the correction are in the new ``validation``
  module (``check_parcellation``, ``validate_parcellation``,
  ``fix_wm_codes``), with the first tests of the package.

0.4.1 (2026-08-11)
------------------

* Fixed ``processing.py`` never importing ``cltmisc``, which made every FSL
  FIRST subcortical parcellation fail with ``NameError``.  The import was
  lost when ``chimera.py`` was decomposed in 0.4.0.
* Fixed ``main()`` performing arithmetic on ``os.cpu_count()`` without
  handling ``None``, which raised ``TypeError`` before any subject was
  processed on platforms where the core count is undeterminable.
* Corrected ``requires-python`` to ``>=3.11``.  0.4.0 advertised ``>=3.9``
  but could not be installed there: ``clabtoolkit`` requires
  ``scipy>=1.16``, which is Python 3.11+ only.  ``scipy>=1.16`` is now
  pinned directly so the constraint is visible.
* Repaired the CI pipeline, which had never passed: dropped the
  unresolvable 3.9/3.10 matrix entries, replaced the ``types-all``
  pre-commit hook (its dependency was removed from PyPI), moved the ruff
  settings under ``[tool.ruff.lint]``, and relaxed mypy to the checks the
  codebase can meet.
* Switched PyPI and TestPyPI uploads to Trusted Publishing.  The workflow
  referenced API token secrets that were never configured, so every
  release since 0.3.1 failed to publish.

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
