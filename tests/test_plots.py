#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import os
import numpy as np
import pandas as pd
from tests.utils import to_test_data, to_test_results
from beamspy.plots import report


class PlotsTestCase(unittest.TestCase):

    def test_report(self):

        report(to_test_data("results_annotation.sqlite"), to_test_results("test_report_01.pdf"),
               "r_value", "p_value", "ppm_error", "adduct")
        statinfo = os.stat(to_test_results("test_report_01.pdf"))
        # print(statinfo.st_size)
        self.assertTrue(statinfo.st_size > 10000)

        report(to_test_data("results_pearson_all.sqlite"),  to_test_results("test_report_02.pdf"),
               "r_value", "p_value", "ppm_error", "adduct")
        statinfo = os.stat(to_test_results("test_report_02.pdf"))
        # print(statinfo.st_size)
        self.assertTrue(statinfo.st_size > 10000)

if __name__ == '__main__':
    unittest.main()
