#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools
import beamspy


def main():

    install_requires = open("requirements.txt").read().splitlines()

    setuptools.setup(name="beamspy",
        version=beamspy.__version__,
        description="Putative annotation of metabolites for mass spectrometry-based metabolomics datasets.",
        long_description=open("README.rst").read(),
        long_description_content_type="text/x-rst",
        author="Ralf Weber",
        author_email="r.j.weber@bham.ac.uk",
        url="https://github.com/computational-metabolomics/beamspy",
        license="GPLv3",
        platforms=["Windows, UNIX"],
        keywords=["Metabolomics", "Mass spectrometry", "Liquid-Chromatography Mass Spectrometry", "Metabolite Annotation"],
        packages=setuptools.find_packages(),
        python_requires=">=3.7",
        test_suite="tests.suite",
        install_requires=install_requires,
        include_package_data=True,
        classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Utilities",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
        ],
        entry_points={
            "console_scripts": [
                "beamspy = beamspy.__main__:main"
            ]
        }
    )


if __name__ == "__main__":
    main()
