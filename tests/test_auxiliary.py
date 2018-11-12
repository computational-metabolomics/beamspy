#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import unittest
from collections import OrderedDict
from beams.auxiliary import *


class AuxiliaryTestCase(unittest.TestCase):

    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_order_composition_by_hill(self):
        composition = OrderedDict([('C', 6), ('O', 6), ('H', 12)])
        hill = order_composition_by_hill(composition)
        self.assertEqual(list(hill), ['C', 'H', 'O'])

    def test_composition_to_string(self):
        composition = OrderedDict([('C', 6), ('O', 6), ('H', 12)])
        mf = composition_to_string(composition)
        self.assertEqual(mf, "C6H12O6")

    def test_double_bond_equivalents(self):
        composition = OrderedDict([('C', 6), ('H', 12), ('O', 6)])
        dbe = double_bond_equivalents(composition)
        self.assertEqual(dbe, 1)

        composition = OrderedDict([('C', 6), ('H', 24), ('O', 12)])
        dbe = double_bond_equivalents(composition)
        self.assertEqual(dbe, -5.0)

    def test_HC_HNOPS_rules(self):
        molecular_formula = "C6H12O6"
        rules = HC_HNOPS_rules(molecular_formula)
        self.assertEqual(rules, {"HC": 1, "NOPSC": 1})

        molecular_formula = "C6H36O6"
        rules = HC_HNOPS_rules(molecular_formula)
        self.assertEqual(rules, {"HC": 0, "NOPSC": 1})

        molecular_formula = "C6H12O32"
        rules = HC_HNOPS_rules(molecular_formula)
        self.assertEqual(rules, {"HC": 1, "NOPSC": 0})

    def test_lewis_senior_rules(self):
        molecular_formula = "C6H12O6"
        rules = lewis_senior_rules(molecular_formula)
        self.assertEqual(rules, {"lewis": 1, "senior": 1})

        molecular_formula = "C6H24O12"
        rules = lewis_senior_rules(molecular_formula)
        self.assertEqual(rules, {"lewis": 1, "senior": 0})


if __name__ == '__main__':
    unittest.main()
