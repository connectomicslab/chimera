=======
History
=======

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