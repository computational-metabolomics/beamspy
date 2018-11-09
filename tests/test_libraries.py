#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import unittest
from beams.in_out import *
from collections import OrderedDict


class LibrariesTestCase(unittest.TestCase):
    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_read_isotopes(self):
        lib_isotopes = read_isotopes(os.path.join(self.path, "beams", "data", "isotopes.txt"), "pos")

    def test_read_adducts(self):
        lib_adducts = read_adducts(os.path.join(self.path, "beams", "data", "adducts.txt"), "pos")
        self.assertTrue("in library" in lib_adducts.__str__())
        lib_adducts.add("test", 10)
        self.assertEqual(lib_adducts.lib["test"], 10)

        lib_adducts.remove("*")
        self.assertEqual(lib_adducts.lib, OrderedDict())

    def test_multiple_charged_ions(self):
        lib_multiple_charged_ions = read_multiple_charged_ions(os.path.join(self.path, "beams", "data", "multiple_charged_ions.txt"), "pos")
        self.assertTrue("in library" in lib_multiple_charged_ions.__str__())

        lib_multiple_charged_ions.remove("*")
        self.assertEqual(lib_multiple_charged_ions.lib, OrderedDict())


if __name__ == '__main__':
    unittest.main()