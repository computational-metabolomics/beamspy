#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import os
import numpy as np
import pandas as pd
from tests.utils import to_test_results
from beams.plots import report


class PlotsTestCase(unittest.TestCase):

    def test_report(self):

        np.random.seed(0)
        n = 1000
        mu, sigma = 0, 0.1 # mean and standard deviation
        s = np.random.normal(mu, sigma, n)

        lib = ["[M+H]+", "[M+Na]+", "[M+K]+"]
        adducts = [lib[i] for i in np.random.randint(3, size=n)]

        df = pd.DataFrame({'ppm_error': s, "adduct": adducts})

        report(df, "ppm_error", "adduct", to_test_results("test_report_01.pdf"))

        statinfo = os.stat(to_test_results("test_report_01.pdf"))
        # print statinfo.st_size
        self.assertTrue(abs(statinfo.st_size - 18788) < 100)


if __name__ == '__main__':
    unittest.main()
