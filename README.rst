BEAMSpy - Birmingham mEtabolite Annotation for Mass Spectrometry (Python package)
==================================================================================
|Version| |Py versions| |Git| |Bioconda| |Build Status (Travis)| |Build Status (AppVeyor)| |License| |RTD doc| |codecov| |mybinder|

BEAMSpy (Birmingham mEtabolite Annotation for Mass Spectrometry) is a Python package that includes several automated and
seamless computational modules that are applied to putatively annotate metabolites detected in untargeted ultra (high)
performance liquid chromatography-mass spectrometry or untargeted direct infusion mass spectrometry metabolomic assays
in an automated process. All reported metabolites are annotated to level 3 of the Metabolomics Standards
Initiative (MSI) reporting standards (Metabolomics. 2007 Sep; 3(3): 211–221. `doi: 10.1007/s11306-007-0082-2 <https://doi.org/10.1007/s11306-007-0082-2>`_).

The package is highly flexible to suit the diversity of sample types studied and mass spectrometers applied in
untargeted metabolomics studies. The user can use the standard reference files included in the package or can develop
their own reference files.


- `Documentation (Read the Docs) <https://beamspy.readthedocs.io/en/latest/>`_
- `Bug reports <https://github.com/computational-metabolomics/beamspy/issues>`_


.. _quick_installation:

Quick installation
-------------------

Conda_
~~~~~~~

1. `Install Miniconda <https://conda.io/docs/user-guide/install>`_.
2. Run the following commands to install BEAMSpy.

Windows-64, Linux-64 and OSx

::

    $ conda create -n beamspy beamspy -c conda-forge -c bioconda -c computational-metabolomics
    $ activate beamspy

Linux-64 and OSx

::

    $ conda create -n beamspy beamspy -c conda-forge -c bioconda -c computational-metabolomics
    $ source activate beamspy


Usage
------------------------

Command line interface (CLI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    $ beamspy --help

Graphical user interface (GUI)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    $ beamspy start-gui

Bug reports
------------------------

Please report any bugs that you find `here <https://github.com/computational-metabolomics/beamspy/issues>`_.
Or fork the repository on `GitHub <https://github.com/computational-metabolomics/beamspy/>`_
and create a pull request (PR). We welcome all contributions, and we will help you to make the PR if you are new to `git <https://guides.github.com/activities/hello-world/>`_.

License
------------------------

Released under the GNU General Public License v3.0 (see `LICENSE file <https://github.com/computational-metabolomics/beamspy/LICENSE>`_)

.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/beamspy.svg?branch=dev&style=flat&maxAge=3600&label=Travis-CI
   :target: https://travis-ci.com/computational-metabolomics/beamspy

.. |Build Status (AppVeyor)| image:: https://img.shields.io/appveyor/ci/RJMW/beamspy.svg?style=flat&maxAge=3600&label=AppVeyor
   :target: https://ci.appveyor.com/project/RJMW/beamspy

.. |Py versions| image:: https://img.shields.io/pypi/pyversions/beamspy.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/beamspy/

.. |Version| image:: https://img.shields.io/pypi/v/beamspy.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/beamspy/

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/beamspy

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: http://bioconda.github.io/recipes/beamspy/README.html

.. |License| image:: https://img.shields.io/badge/License-GPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/gpl-3.0.html

.. |RTD doc| image:: https://img.shields.io/badge/documentation-RTD-71B360.svg?style=flat&maxAge=3600
   :target: https://beamspy.readthedocs.io/en/latest/

.. |codecov| image:: https://codecov.io/gh/computational-metabolomics/beamspy/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/computational-metabolomics/beamspy

.. |mybinder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/computational-metabolomics/beamspy/master?filepath=notebooks

.. _pip: https://pip.pypa.io/
.. _Conda: http://conda.pydata.org/docs/
