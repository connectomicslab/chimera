.. highlight:: shell

============
Installation
============


Stable release
--------------

To install Chimera, run this command in your terminal:

.. code-block:: console

    $ pip install chimera

Alternatively, you can install directly from GitHub for the latest development version:

.. code-block:: console

    $ pip install git+https://github.com/yasseraleman/chimera.git

This is the preferred method to install Chimera, as it will always install the most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From sources
------------

The sources for Chimera can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/yasseraleman/chimera

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/yasseraleman/chimera/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ cd chimera
    $ pip install -e .

For development, install with development dependencies:

.. code-block:: console

    $ pip install -e ".[dev,docs]"

System Requirements
-------------------

Chimera requires several external neuroimaging tools to be installed:

**Required Software:**

- **FreeSurfer** (version ≥ 7.2.0): For cortical surface processing and subcortical segmentation
- **FSL**: For subcortical structure segmentation using FIRST
- **ANTs**: For advanced normalization and registration

**Installation Instructions:**

1. **FreeSurfer**: Download and install from https://surfer.nmr.mgh.harvard.edu/
2. **FSL**: Install from https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation
3. **ANTs**: Install from http://stnava.github.io/ANTs/

Make sure these tools are properly sourced in your environment before using Chimera.

**Container Options:**

Chimera also supports containerized execution using Docker or Singularity, which can simplify dependency management.


.. _Github repo: https://github.com/yasseraleman/chimera
.. _tarball: https://github.com/yasseraleman/chimera/tarball/master
