#!/usr/bin/python
# -*- coding: utf-8 -*-
import setuptools
import sys
import beams
import os
import unittest

def main():

    install_requires = open('requirements.txt').read().splitlines()

    setuptools.setup(name="beams",
        version=beams.__version__,
        description="Python package to annotate LC-MS and DIMS data",
        long_description=open('README.rst').read(),
        author="Ralf Weber",
        author_email="r.j.weber@bham.ac.uk",
        url="https://github.com/computational-metabolomics/beams",
        license="GPLv3",
        platforms=['Windows, UNIX'],
        keywords=['Metabolomics', 'Mass spectrometry', 'Liquid-Chromatography Mass Spectrometry', 'Metabolite Annotation'],
        packages=setuptools.find_packages(),
        test_suite='tests.suite',
        install_requires=install_requires,
        include_package_data=True,
        classifiers=[
          "Programming Language :: Python :: 2",
          "Programming Language :: Python :: 2.7",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.6",
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