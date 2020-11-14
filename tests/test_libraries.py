#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import unittest
from beamspy.in_out import *
from beamspy.libraries import *
from beamspy.auxiliary import nist_database_to_pyteomics
from collections import OrderedDict


class LibrariesTestCase(unittest.TestCase):
    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_read_isotopes(self):
        lib_isotopes = read_isotopes(os.path.join(self.path, "beamspy", "data", "isotopes.txt"), "pos")
        self.assertTrue("in library" in lib_isotopes.__str__())

    def test_read_adducts(self):
        lib_adducts = read_adducts(os.path.join(self.path, "beamspy", "data", "adducts.txt"), "pos")
        self.assertTrue("in library" in lib_adducts.__str__())
        lib_adducts.add("test", 100.0, 1)
        self.assertEqual(lib_adducts.lib["test"]["mass"], 100.0)
        self.assertEqual(lib_adducts.lib["test"]["charge"], 1)

        lib_adducts.remove("*")
        self.assertEqual(lib_adducts.lib, OrderedDict())

    # def test_mass_differences(self):
    #     lib_differences = read_mass_differences(os.path.join(self.path, "beamspy", "data", "adducts_differences.txt"), "pos")
    #     self.assertTrue("in library" in lib_differences.__str__())
    #
    #     lib_differences.remove("*", "*")
    #     self.assertEqual(lib_differences.lib, [])

    def test_nist_database_to_pyteomics(self):
        nist_database = nist_database_to_pyteomics(os.path.join(self.path, "beamspy", "data", "nist_database.txt"))
        self.assertEqual(nist_database["C"][0], (12.0, 1.0))
        self.assertEqual(nist_database["H"][0], (1.00782503223, 1.0))
        self.assertEqual(nist_database["N"][0], (14.00307400443, 1.0))
        self.assertEqual(nist_database["O"][0], (15.99491461957, 1.0))
        self.assertEqual(nist_database["P"][0], (30.97376199842, 1.0))
        self.assertEqual(nist_database["S"][0], (31.9720711744, 1.0))


if __name__ == '__main__':
    unittest.main()