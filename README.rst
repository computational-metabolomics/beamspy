BEAMS - Birmingham mEtabolite Annotation for Mass SpectroMetry
================================================================
|Version| |Py versions| |Git| |Bioconda| |Build Status (Travis)| |Build Status (AppVeyor)| |License| |RTD doc| |codecov| |mybinder|

BEAMS is a Python Package to annotate LC-MS and DIMS data.

- |documentation|
- |source|
- |bug reports|

Installation
------------------------

Conda_
~~~~~~~

1. Install |miniconda|. Follow the steps described |conda_install|.
2. Run the following commands to install BEAMS.

Windows-64, Linux-64 and OSx

::

    $ conda create -n beams beams -c conda-forge -c bioconda -c computational-metabolomics
    $ source activate beams

Linux-64 and OSx

::

    $ conda create -n beams beams -c conda-forge -c bioconda -c computational-metabolomics
    $ source activate beams


Usage
------------------------

Command line
~~~~~~~~~~~~~

::

    $ beams --help

GUI
~~~~~~~~~~~~~

::

    $ beams start-gui

Bugs
------------------------

Please report any bugs that you find `here <https://github.com/computational-metabolomics/beams/issues>`_.
Or fork the repository on `GitHub <https://github.com/computational-metabolomics/beams/>`_
and create a pull request (PR). We welcome all contributions, and we will help you to make the PR if you are new to `git <https://guides.github.com/activities/hello-world/>`_.

License
------------------------

Released under the GNU General Public License v3.0 (see `LICENSE file <https://github.com/computational-metabolomics/beams/LICENSE>`_)

.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/beams.svg?branch=dev&style=flat&maxAge=3600&label=Travis-CI
   :target: https://travis-ci.com/computational-metabolomics/beams

.. |Build Status (AppVeyor)| image:: https://img.shields.io/appveyor/ci/RJMW/beams.svg?style=flat&maxAge=3600&label=AppVeyor
   :target: https://ci.appveyor.com/project/RJMW/beams

.. |Py versions| image:: https://img.shields.io/pypi/pyversions/beams.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/beams/

.. |Version| image:: https://img.shields.io/pypi/v/beams.svg?style=flat&maxAge=3600
   :target: https://pypi.python.org/pypi/beams/

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/beams

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: http://bioconda.github.io/recipes/beams/README.html

.. |License| image:: https://img.shields.io/pypi/l/beams.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html

.. |RTD doc| image:: https://img.shields.io/badge/documentation-RTD-71B360.svg?style=flat&maxAge=3600
   :target: http://beams.readthedocs.io/en/latest/beams/index.html

.. |codecov| image:: https://codecov.io/gh/computational-metabolomics/beams/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/computational-metabolomics/beams

.. |mybinder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/computational-metabolomics/beams/master?filepath=notebooks

.. |documentation| raw:: html

   <a href="http://beams.readthedocs.io/en/latest/beams/index.html" target="_blank">Documentation</a>

.. |source| raw:: html

   <a href="https://github.com/computational-metabolomics/beams/tree/dev/beams" target="_blank">Source</a>

.. |bug reports| raw:: html

   <a href="https://github.com/computational-metabolomics/beams/issues" target="_blank">Bug reports</a>

.. |conda_install| raw:: html

   <a href="https://conda.io/docs/user-guide/install" target="_blank">here</a>

.. |miniconda| raw:: html

   <a href="http://conda.pydata.org/miniconda.html" target="_blank">Miniconda</a>

.. _pip: https://pip.pypa.io/
.. _Conda: http://conda.pydata.org/docs/

