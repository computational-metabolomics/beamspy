#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
import sys
import beams
import os
import unittest

def main():

    if sys.version_info[0] != 2 and sys.version_info[1] <= 7:
        sys.exit("Python-2.7.8 is required ")

    if sys.platform == "win32" or sys.platform == "win64":
        install_requires = open('requirements-win.txt').read().splitlines()
    else:
        install_requires = open('requirements.txt').read().splitlines()

    setuptools.setup(name="beams",
        version=beams.__version__,
        description="Python package to annotate LC-MS and DIMS data",
        long_description=open('README.rst').read(),
        author="Han Zhang, Ralf Weber",
        author_email="aynhzhanghan@gmail.com, r.j.weber@bham.ac.uk",
        url="https://github.com/computational-metabolomics/beams",
        license="GPLv3",
        platforms=['Windows, UNIX'],
        keywords=['Metabolomics', 'Mass spectrometry', 'Liquid-Chromatography Mass Spectrometry'],
        packages=setuptools.find_packages(),
        test_suite='tests.suite',
        install_requires=install_requires,
        include_package_data=True,
        classifiers=[
          "Programming Language :: Python :: 2",
          "Programming Language :: Python :: 2.7",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          "Topic :: Scientific/Engineering :: Chemistry",
          "Topic :: Utilities",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
        ],
        entry_points={
         'console_scripts': [
             'beams = beams.__main__:main'
         ]
        }
    )


if __name__ == "__main__":
    main()