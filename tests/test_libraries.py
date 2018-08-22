#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import unittest
from beams.in_out import *


class LibrariesTestCase(unittest.TestCase):
    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_read_isotopes(self):
        self.lib_isotopes = read_isotopes(os.path.join(self.path, "beams", "data", "isotopes.txt"), "pos")

    def test_read_adducts(self):
        self.lib_adducts = read_adducts(os.path.join(self.path, "beams", "data", "adducts.txt"), "pos")

    def test_multiple_charged_ions(self):
        self.lib_multiple_charged_ions = read_multiple_charged_ions(os.path.join(self.path, "beams", "data", "multiple_charged_ions.txt"), "pos")


if __name__ == '__main__':
    unittest.main()