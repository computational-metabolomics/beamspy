beams
========
|Version| |Py versions| |Git| |Bioconda| |Build Status (Travis)| |Build Status (AppVeyor)| |License| |RTD doc|

BEAMS is a Python Package to annotate LC-MS and DIMS data

- **Documentation:** https://computational-metabolomics.github.io/beams
- **Source:** https://github.com/computational-metabolomics/beams
- **Bug reports:** https://github.com/computational-metabolomics/beams/issues

Installation
--------

Conda_
~~~~~~~

1. Install Conda_ (For example: `Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`__).
2. Run the following commands to install BEAMS.

Windows-64

::

    $ conda create -n beams python=2.7 -y --name beams python=2.7 numpy scipy networkx requests pandas -c conda-forge
    $ activate beams
    $ pip install git+https://github.com/computational-metabolomics/beams.git@dev

Linux-64 and OSx

::

    $ conda create -n beams python=2.7 -y --name beams python=2.7 numpy scipy networkx requests pandas pyqt -c conda-forge
    $ source activate beams
    $ pip install git+https://github.com/computational-metabolomics/beams.git@dev

Usage
------

Command line
~~~~~~~~~~~~~

::

    $ beams --help

GUI
~~~~~~~~~~~~~

::

    $ beams start-gui

Bugs
----

Please report any bugs that you find `here <https://github.com/computational-metabolomics/beams/issues>`_.
Or fork the repository on `GitHub <https://github.com/computational-metabolomics/beams/>`_
and create a pull request (PR). We welcome all contributions, and we
will help you to make the PR if you are new to `git` (see `CONTRIBUTING.rst`).

License
-------

Released under the GNU General Public License v3.0 (see `LICENSE` file)::

   Copyright (C) 2017 BEAMS Developers
   Ralf J.M. Weber <r.j.weber@bham.ac.uk>   
   Han Zhang <aynhzhanghan@gmail.com>

.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/beams.svg?style=flat&maxAge=3600&label=Travis-CI
   :target: https://travis-ci.org/computational-metabolomics/beams

.. |Build Status (AppVeyor)| image:: https://img.shields.io/appveyor/ci/computational-metabolomics/mzml2isa.svg?style=flat&maxAge=3600&label=AppVeyor
   :target: https://ci.appveyor.com/project/computational-metabolomics/beams

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

.. _pip: https://pip.pypa.io/
.. _Conda: http://conda.pydata.org/docs/
