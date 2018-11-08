#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import os
import copy
from collections import OrderedDict
import pyteomics
from pyteomics.mass import nist_mass, calculate_mass
from beams.auxiliary import update_and_sort_nist_mass


class AuxiliaryTestCase(unittest.TestCase):

    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_update_and_sort_nist_mass(self):
        exact_mass = calculate_mass(formula="H")
        self.assertEqual(exact_mass, 1.00782503207)

        path_atoms = os.path.join(self.path, "beams", "data", "elements.txt")
        mass_data = update_and_sort_nist_mass(path_atoms, digits=6)
        exact_mass = calculate_mass(formula="H", mass_data=mass_data)
        self.assertEqual(exact_mass, 1.007825)

        exact_mass = calculate_mass(formula="H")
        self.assertEqual(exact_mass, 1.00782503207)


if __name__ == '__main__':
    unittest.main()
